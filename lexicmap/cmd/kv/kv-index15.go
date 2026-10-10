// Copyright © 2023-2026 Wei Shen <shenwei356@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

package kv

import (
	"bufio"
	"crypto/rand"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"

	mmap "github.com/edsrzf/mmap-go"
	"github.com/pkg/errors"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
)

// KVIndex15FileExt is the file extension of the second-level seed index.
const KVIndex15FileExt = ".idx15"

// IndexMainVersionTagged distinguishes the tagged-offset index from the legacy
// index. Keeping this separate from MainVersion lets the seeds file stay
// byte-for-byte unchanged.
const IndexMainVersionTagged uint8 = 2

const (
	index15HeaderSize      = 48
	index15SuffixBytes     = 5
	index15BlockFixedBytes = 94
	index15OffsetTag       = uint64(1) << 63
	index15RawOffsetMask   = index15OffsetTag - 1
	index15ByteOrderBig    = 1
	index15RelativeFixedBE = 1
	index15FormatVersion   = 1
	index15FormatMinor     = 0
)

// MinIndex15Threshold is the smallest useful seeds-block threshold. The seed
// reader already uses 4 KiB buffers, while each indexed block can fault in an
// idx15 page, so indexing smaller blocks only increases index overhead.
const MinIndex15Threshold uint64 = 4 << 10

var magicIndex15 = [8]byte{'.', 'k', 'v', 'i', 'd', 'x', '1', '5'}

// idx15 binary layout
//
// All multi-byte integers are stored in big-endian byte order. Blocks are
// concatenated without padding and are addressed directly by tagged offsets in
// the primary .idx file.
//
// Header (48 bytes):
//
//	offset  bytes  field
//	0       8      magic (".kvidx15")
//	8       1      format version
//	9       1      format minor version
//	10      1      byte-order code (1 = big endian)
//	11      1      primary prefix length
//	12      1      secondary prefix length (primary + 2)
//	13      1      k-mer length
//	14      1      checkpoint suffix width (5 bytes)
//	15      1      relative-offset encoding (1 = fixed-width big endian)
//	16      8      index of the first mask in the seeds chunk
//	24      8      number of masks in the seeds chunk
//	32      8      seeds file size
//	40      6      reserved, zero
//	46      2      build ID shared with the primary .idx file
//
// The build ID is a nonzero random 16-bit pairing token generated separately
// for each seeds chunk whenever adaptive indexes are created or rebuilt. The same value is written
// to the primary .idx header and this .idx15 header. A reader rejects unequal
// values, which detects a stale or incorrectly paired sidecar, including the
// state left if replacement stops after renaming only one of the two files.
// It is not a format version, database-wide ID, content hash, or cryptographic
// integrity check. The seeds size and chunk metadata provide additional
// consistency checks.
//
// Block (94 + 15*offset_width bytes):
//
//	offset  bytes              field
//	0       1                  offset_width (1..8)
//	1       2                  present_mask, one bit per secondary subprefix
//	3       8                  absolute seeds base offset
//	11      8                  pair-first delta checkpoint at the base offset
//	19      15*5               checkpoint suffixes for subprefixes 1..15
//	94      15*offset_width    offsets relative to the seeds base offset
//
// Subprefix 0 uses the base offset and pair-first checkpoint directly. For
// subprefix i in 1..15, suffix and relative-offset slot i-1 are used. Empty
// subprefixes have a cleared present_mask bit and zero-filled slots. A zero
// relative offset reuses the pair-first checkpoint, which handles two
// subprefixes or primary blocks sharing the same encoded k-mer pair.

// Index15Stats summarizes one adaptively indexed seeds chunk.
type Index15Stats struct {
	TotalBlocks    uint64
	NonEmptyBlocks uint64
	IndexedBlocks  uint64
	WidthCounts    [8]uint64
	Index15Bytes   uint64
}

// Index15ProgressFunc receives the number of completely processed masks and
// the total number of masks in one seeds chunk. The first call reports zero
// processed masks after the seeds header has been validated.
type Index15ProgressFunc func(processedMasks, totalMasks uint64)

// Add merges another per-file summary into s.
func (s *Index15Stats) Add(other Index15Stats) {
	s.TotalBlocks += other.TotalBlocks
	s.NonEmptyBlocks += other.NonEmptyBlocks
	s.IndexedBlocks += other.IndexedBlocks
	s.Index15Bytes += other.Index15Bytes
	for i := range s.WidthCounts {
		s.WidthCounts[i] += other.WidthCounts[i]
	}
}

// IsIndex15Offset reports whether a tagged offset points into idx15.
func IsIndex15Offset(offset uint64) bool {
	return offset&index15OffsetTag != 0
}

// RawIndexOffset removes the idx15 tag. It should only be used for tagged
// offsets; untagged entries retain the legacy pair-position bit.
func RawIndexOffset(offset uint64) uint64 {
	return offset & index15RawOffsetMask
}

// MakeSeedsOffset converts a byte offset to the untagged legacy representation.
func MakeSeedsOffset(offset uint64, isSecond bool) (uint64, error) {
	if offset > index15RawOffsetMask>>1 {
		return 0, fmt.Errorf("seeds offset exceeds 62 bits: %d", offset)
	}
	encoded := offset << 1
	if isSecond {
		encoded |= 1
	}
	return encoded, nil
}

// MakeIndex15Offset validates and tags an idx15 offset.
func MakeIndex15Offset(offset uint64) (uint64, error) {
	if offset&index15OffsetTag != 0 {
		return 0, fmt.Errorf("idx15 offset exceeds 63 bits: %d", offset)
	}
	return offset | index15OffsetTag, nil
}

func validKVIndexVersion(version uint8) bool {
	return version == MainVersion || version == IndexMainVersionTagged
}

func validKVIndexMetadata(version, config uint8) bool {
	return validKVIndexVersion(version) && (version == IndexMainVersionTagged) == (config&MaskHasIndex15 != 0)
}

type index15BuildBlock struct {
	prefix uint64 // Encoded primary prefix identifying this block.

	indexKmer           uint64 // First target k-mer in the primary-prefix block; written to idx.
	pairFirstCheckpoint uint64 // First k-mer of the pair at baseOffset; written to the idx15 block.
	baseOffset          uint64 // Absolute seeds offset of the first pair intersecting this block.
	baseIsSecond        bool   // Whether indexKmer is the second k-mer of the pair at baseOffset.
	endOffset           uint64 // Exclusive seeds offset used to measure the block size.

	presentMask  uint16     // One bit per secondary subprefix that occurs in the block.
	suffixes     [15]uint64 // Pair-first checkpoint suffixes for secondary subprefixes 1 through 15.
	relOffsets   [15]uint64 // Pair offsets relative to baseOffset for subprefixes 1 through 15.
	maxRelOffset uint64     // Largest value in relOffsets; determines the encoded offset width.
}

func (b *index15BuildBlock) add(kmer, pairKmer, pairOffset, pairEnd uint64, k, primaryPrefix, secondaryPrefix uint8) {
	subprefix := uint8(kmer >> ((k - secondaryPrefix) << 1) & 15)
	bit := uint16(1) << subprefix
	if b.presentMask&bit == 0 {
		b.presentMask |= bit
		if subprefix > 0 {
			relative := pairOffset - b.baseOffset
			b.relOffsets[subprefix-1] = relative
			if relative > 0 {
				suffixBits := (k - primaryPrefix) << 1
				b.suffixes[subprefix-1] = pairKmer & (uint64(1)<<suffixBits - 1)
			}
			if relative > b.maxRelOffset {
				b.maxRelOffset = relative
			}
		}
	}
	b.endOffset = pairEnd
}

func putUint40(b []byte, value uint64) {
	_ = b[4]
	b[0] = byte(value >> 32)
	b[1] = byte(value >> 24)
	b[2] = byte(value >> 16)
	b[3] = byte(value >> 8)
	b[4] = byte(value)
}

func uint40(b []byte) uint64 {
	_ = b[4]
	return uint64(b[0])<<32 | uint64(b[1])<<24 | uint64(b[2])<<16 | uint64(b[3])<<8 | uint64(b[4])
}

func putUintWidth(b []byte, value uint64, width uint8) {
	for i := int(width) - 1; i >= 0; i-- {
		b[i] = byte(value)
		value >>= 8
	}
}

func uintWidth(b []byte, width uint8) uint64 {
	var value uint64
	for i := uint8(0); i < width; i++ {
		value = value<<8 | uint64(b[i])
	}
	return value
}

func index15OffsetWidth(value uint64) uint8 {
	return util.ByteLengthUint64(value)
}

func writeIndex15Block(w io.Writer, block *index15BuildBlock, scratch []byte) (uint8, int, error) {
	// The full pair checkpoint is the 8-byte field at [11:19]. It is needed
	// when a pair straddles a primary-prefix boundary because the sparse index
	// must keep the first target k-mer as its anchor-identifying checkpoint.
	width := index15OffsetWidth(block.maxRelOffset)
	size := index15BlockFixedBytes + 15*int(width)
	buf := scratch[:size]
	clear(buf)
	buf[0] = width
	be.PutUint16(buf[1:3], block.presentMask)
	be.PutUint64(buf[3:11], block.baseOffset)
	be.PutUint64(buf[11:19], block.pairFirstCheckpoint)
	for i := 0; i < 15; i++ {
		putUint40(buf[19+i*index15SuffixBytes:], block.suffixes[i])
		putUintWidth(buf[index15BlockFixedBytes+i*int(width):], block.relOffsets[i], width)
	}
	_, err := w.Write(buf)
	return width, size, err
}

func writeIndex15Header(w io.Writer, k, primaryPrefix, secondaryPrefix uint8, chunkIndex, chunkSize, seedsSize uint64, buildID uint16) error {
	buf := make([]byte, index15HeaderSize)
	copy(buf[:8], magicIndex15[:])
	buf[8] = index15FormatVersion
	buf[9] = index15FormatMinor
	buf[10] = index15ByteOrderBig
	buf[11] = primaryPrefix
	buf[12] = secondaryPrefix
	buf[13] = k
	buf[14] = index15SuffixBytes
	buf[15] = index15RelativeFixedBE
	be.PutUint64(buf[16:24], chunkIndex)
	be.PutUint64(buf[24:32], chunkSize)
	be.PutUint64(buf[32:40], seedsSize)
	be.PutUint64(buf[40:48], uint64(buildID))
	_, err := w.Write(buf)
	return err
}

// randomIndexBuildID creates the nonzero pairing token stored in both index
// headers. It is intentionally only an identity check between the two files.
func randomIndexBuildID() (uint16, error) {
	var buf [2]byte
	if _, err := rand.Read(buf[:]); err != nil {
		return 0, err
	}
	id := be.Uint16(buf[:])
	if id == 0 {
		id = 1
	}
	return id, nil
}

func anchorPrefixForPartitions(partitions int) (uint8, error) {
	if partitions < 4 {
		return 0, fmt.Errorf("partitions should be a power of 4 and at least 4: %d", partitions)
	}
	prefix := uint8(0)
	for partitions > 1 {
		if partitions%4 != 0 {
			return 0, fmt.Errorf("partitions should be a power of 4")
		}
		partitions /= 4
		prefix++
	}
	return prefix, nil
}

func index15PrefixLengths(k, maskPrefix, anchorPrefix uint8) (uint8, uint8, error) {
	if maskPrefix > k || anchorPrefix > k-maskPrefix {
		return 0, 0, fmt.Errorf("mask prefix (%d) + anchor prefix (%d) should not exceed k (%d)", maskPrefix, anchorPrefix, k)
	}
	primaryPrefix := maskPrefix + anchorPrefix
	if k-primaryPrefix < 2 {
		return 0, 0, fmt.Errorf("k (%d) should contain at least 2 bases after the primary prefix (%d)", k, primaryPrefix)
	}
	if (k-primaryPrefix)<<1 > index15SuffixBytes*8 {
		return 0, 0, fmt.Errorf("the %d-byte idx15 checkpoint suffix cannot encode k (%d) with primary prefix %d", index15SuffixBytes, k, primaryPrefix)
	}
	return primaryPrefix, primaryPrefix + 2, nil
}

// CheckIndex15Options validates adaptive index settings before seed data is built.
func CheckIndex15Options(k, maskPrefix uint8, partitions int, threshold uint64) error {
	if threshold < MinIndex15Threshold {
		return fmt.Errorf("idx15 threshold should be at least %d bytes", MinIndex15Threshold)
	}
	anchorPrefix, err := anchorPrefixForPartitions(partitions)
	if err != nil {
		return err
	}
	_, _, err = index15PrefixLengths(k, maskPrefix, anchorPrefix)
	return err
}

func discardSeedPositions(r io.Reader, count, bytesPerPosition uint64, offset *uint64) error {
	if count == 0 {
		return nil
	}
	if count > math.MaxUint64/bytesPerPosition {
		return ErrBrokenFile
	}
	n := count * bytesPerPosition
	if n > math.MaxInt64 || *offset > math.MaxUint64-n {
		return ErrBrokenFile
	}
	read, err := io.CopyN(io.Discard, r, int64(n))
	if err != nil {
		return err
	}
	if uint64(read) != n {
		return ErrBrokenFile
	}
	*offset += n
	return nil
}

// CreateKVIndex15 recreates the primary index and the adaptive second-level
// index while leaving the seeds file unchanged.
// Longer secondary prefixes locate seeds within large data blocks in large
// indexes, reducing scanning during prefix matching. Overall gains are most
// noticeable for batch lexicmap search queries when -n/--top-n-genomes limits
// downstream alignments, making seed matching a larger share of runtime.
// Without this limit, sequence alignment often dominates runtime.
func CreateKVIndex15(file string, partitions int, threshold uint64) (Index15Stats, error) {
	return CreateKVIndex15WithProgress(file, partitions, threshold, nil)
}

// CreateKVIndex15WithProgress is CreateKVIndex15 with per-mask progress
// reporting for one seeds chunk.
func CreateKVIndex15WithProgress(file string, partitions int, threshold uint64, progress Index15ProgressFunc) (stats Index15Stats, err error) {
	if threshold < MinIndex15Threshold {
		return stats, fmt.Errorf("idx15 threshold should be at least %d bytes", MinIndex15Threshold)
	}
	anchorPrefix, err := anchorPrefixForPartitions(partitions)
	if err != nil {
		return stats, err
	}

	indexFile := filepath.Clean(file) + KVIndexFileExt
	indexK, indexChunk, indexMasks, maskPrefix, _, err := ReadKVIndexInfo(indexFile)
	if err != nil {
		return stats, errors.Wrap(err, "reading the existing seed index")
	}
	fh, err := os.Open(file)
	if err != nil {
		return stats, errors.Wrap(err, "opening seeds data")
	}
	defer fh.Close()
	fileInfo, err := fh.Stat()
	if err != nil {
		return stats, err
	}
	if fileInfo.Size() < 32 {
		return stats, ErrBrokenFile
	}
	seedsSize := uint64(fileInfo.Size())

	r := bufio.NewReaderSize(fh, 64<<10)
	buf8 := make([]byte, 8)
	buf16 := make([]byte, 16)
	if _, err = io.ReadFull(r, buf8); err != nil {
		return stats, err
	}
	if string(buf8) != string(Magic[:]) {
		return stats, ErrInvalidFileFormat
	}
	if _, err = io.ReadFull(r, buf8); err != nil {
		return stats, err
	}
	if buf8[0] != MainVersion {
		return stats, ErrVersionMismatch
	}
	k := buf8[2]
	config1 := buf8[3]
	primaryPrefix, secondaryPrefix, err := index15PrefixLengths(k, maskPrefix, anchorPrefix)
	if err != nil {
		return stats, err
	}
	bytesPerPosition := uint64(8)
	if config1&MaskUse3BytesForSeedPos != 0 {
		bytesPerPosition = 7
	}
	if _, err = io.ReadFull(r, buf8); err != nil {
		return stats, err
	}
	chunkIndex := be.Uint64(buf8)
	if _, err = io.ReadFull(r, buf8); err != nil {
		return stats, err
	}
	chunkSize := be.Uint64(buf8)
	if indexK != k || indexChunk != int(chunkIndex) || indexMasks != int(chunkSize) {
		return stats, fmt.Errorf("the existing seed index does not match the seeds file")
	}
	stats.TotalBlocks = chunkSize * uint64(partitions)
	if progress != nil {
		progress(0, chunkSize)
	}

	buildID, err := randomIndexBuildID()
	if err != nil {
		return stats, errors.Wrap(err, "creating index build identifier")
	}
	dir := filepath.Dir(indexFile)
	base := filepath.Base(file)
	indexTmp, err := os.CreateTemp(dir, "."+base+".idx-*")
	if err != nil {
		return stats, err
	}
	indexTmpName := indexTmp.Name()
	defer indexTmp.Close()
	keepIndexTmp := false
	defer func() {
		if !keepIndexTmp {
			_ = os.Remove(indexTmpName)
		}
	}()
	index15Tmp, err := os.CreateTemp(dir, "."+base+".idx15-*")
	if err != nil {
		indexTmp.Close()
		return stats, err
	}
	index15TmpName := index15Tmp.Name()
	defer index15Tmp.Close()
	keepIndex15Tmp := false
	defer func() {
		if !keepIndex15Tmp {
			_ = os.Remove(index15TmpName)
		}
	}()

	wi := bufio.NewWriterSize(indexTmp, 64<<10)
	w15 := bufio.NewWriterSize(index15Tmp, 64<<10)
	if _, err = wi.Write(MagicIdx[:]); err != nil {
		return stats, err
	}
	// The final two metadata bytes pair this primary index with the idx15 file.
	idxMeta := [8]byte{IndexMainVersionTagged, 0, k, maskPrefix, anchorPrefix, config1 | MaskHasIndex15, byte(buildID >> 8), byte(buildID)}
	if _, err = wi.Write(idxMeta[:]); err != nil {
		return stats, err
	}
	be.PutUint64(buf16[:8], chunkIndex)
	be.PutUint64(buf16[8:], chunkSize)
	if _, err = wi.Write(buf16); err != nil {
		return stats, err
	}
	if err = writeIndex15Header(w15, k, primaryPrefix, secondaryPrefix, chunkIndex, chunkSize, seedsSize, buildID); err != nil {
		return stats, err
	}
	builder := newIndex15Builder(wi, w15, k, maskPrefix, anchorPrefix, chunkSize, threshold)
	defer func() { stats = builder.stats }()
	var offset uint64 = 32
	var decodedKmers uint64
	for mask := uint64(0); mask < chunkSize; mask++ {
		builder.beginMask()
		if _, err = io.ReadFull(r, buf8); err != nil {
			return stats, err
		}
		nKmers := be.Uint64(buf8)
		offset += 8
		if nKmers == 0 {
			be.PutUint64(buf8, 0)
			if _, err = wi.Write(buf8); err != nil {
				return stats, err
			}
			if progress != nil {
				progress(mask+1, chunkSize)
			}
			continue
		}

		var previousKmer uint64
		decodedKmers = 0

		for {
			pairOffset := offset
			if _, err = io.ReadFull(r, buf8[:1]); err != nil {
				return stats, err
			}
			keyControl := buf8[0]
			lastPair := keyControl&128 != 0
			hasKmer2 := keyControl&64 == 0
			if !hasKmer2 && !lastPair {
				return stats, ErrBrokenFile
			}
			keyControl &= 63
			nKeyBytes := util.CtrlByte2ByteLengthsUint64(keyControl)
			if _, err = io.ReadFull(r, buf16[:nKeyBytes]); err != nil {
				return stats, err
			}
			v1, v2, nDecoded := util.Uint64s(keyControl, buf16[:nKeyBytes])
			if nDecoded == 0 {
				return stats, ErrBrokenFile
			}
			kmer1 := previousKmer + v1
			kmer2 := kmer1 + v2
			previousKmer = kmer2
			offset += uint64(1 + nKeyBytes)

			if _, err = io.ReadFull(r, buf8[:1]); err != nil {
				return stats, err
			}
			valueControl := buf8[0]
			nValueBytes := util.CtrlByte2ByteLengthsUint64(valueControl)
			if _, err = io.ReadFull(r, buf16[:nValueBytes]); err != nil {
				return stats, err
			}
			len1, len2, nDecoded := util.Uint64s(valueControl, buf16[:nValueBytes])
			if nDecoded == 0 {
				return stats, ErrBrokenFile
			}
			offset += uint64(1 + nValueBytes)
			if err = discardSeedPositions(r, len1, bytesPerPosition, &offset); err != nil {
				return stats, err
			}
			if hasKmer2 {
				if err = discardSeedPositions(r, len2, bytesPerPosition, &offset); err != nil {
					return stats, err
				}
			}
			pairEnd := offset

			if err = builder.addPair(kmer1, kmer2, hasKmer2, pairOffset, pairEnd); err != nil {
				return stats, err
			}
			decodedKmers++
			if hasKmer2 {
				decodedKmers++
			}
			if lastPair {
				break
			}
		}
		if decodedKmers != nKmers {
			return stats, fmt.Errorf("number of k-mers mismatch for mask %d: expected %d, got %d", chunkIndex+mask, nKmers, decodedKmers)
		}
		if err = builder.endMask(offset); err != nil {
			return stats, err
		}
		if progress != nil {
			progress(mask+1, chunkSize)
		}
	}
	if offset != seedsSize {
		return stats, fmt.Errorf("unexpected trailing seeds data: decoded %d of %d bytes", offset, seedsSize)
	}
	if err = wi.Flush(); err != nil {
		return stats, err
	}
	if err = w15.Flush(); err != nil {
		return stats, err
	}
	if err = indexTmp.Sync(); err != nil {
		return stats, err
	}
	if err = index15Tmp.Sync(); err != nil {
		return stats, err
	}
	if err = indexTmp.Close(); err != nil {
		return stats, err
	}
	if err = index15Tmp.Close(); err != nil {
		return stats, err
	}
	mode := os.FileMode(0644)
	if oldInfo, statErr := os.Stat(indexFile); statErr == nil {
		mode = oldInfo.Mode().Perm()
	}
	if err = os.Chmod(indexTmpName, mode); err != nil {
		return stats, err
	}
	if err = os.Chmod(index15TmpName, mode); err != nil {
		return stats, err
	}
	index15File := filepath.Clean(file) + KVIndex15FileExt
	if err = os.Rename(index15TmpName, index15File); err != nil {
		return stats, err
	}
	keepIndex15Tmp = true
	if err = os.Rename(indexTmpName, indexFile); err != nil {
		return stats, err
	}
	keepIndexTmp = true
	return stats, nil
}

type index15 struct {
	data            mmap.MMap
	k               uint8
	primaryPrefix   uint8
	secondaryPrefix uint8
	seedsSize       uint64
}

func openIndex15(file string, k, maskPrefix, anchorPrefix uint8, chunkIndex, chunkSize uint64, seedsSize uint64, buildID uint16) (*index15, error) {
	fh, err := os.Open(filepath.Clean(file) + KVIndex15FileExt)
	if err != nil {
		return nil, err
	}
	data, err := mmap.Map(fh, mmap.RDONLY, 0)
	// On Unix the mapping remains valid after closing fh, so idx15 does not
	// consume one of the persistent file descriptors budgeted for seed readers.
	closeErr := fh.Close()
	if err != nil {
		return nil, err
	}
	if closeErr != nil {
		data.Unmap()
		return nil, closeErr
	}
	fail := func(err error) (*index15, error) {
		_ = data.Unmap()
		return nil, err
	}
	if len(data) < index15HeaderSize {
		return fail(ErrBrokenFile)
	}
	if string(data[:8]) != string(magicIndex15[:]) {
		return fail(ErrInvalidFileFormat)
	}
	if data[8] != index15FormatVersion || data[9] != index15FormatMinor || data[10] != index15ByteOrderBig || data[14] != index15SuffixBytes || data[15] != index15RelativeFixedBE {
		return fail(ErrVersionMismatch)
	}
	primaryPrefix, secondaryPrefix, err := index15PrefixLengths(k, maskPrefix, anchorPrefix)
	if err != nil {
		return fail(err)
	}
	if data[11] != primaryPrefix || data[12] != secondaryPrefix || data[13] != k {
		return fail(fmt.Errorf("idx15 parameters do not match the primary index"))
	}
	if be.Uint64(data[16:24]) != chunkIndex || be.Uint64(data[24:32]) != chunkSize {
		return fail(fmt.Errorf("idx15 chunk metadata does not match the primary index"))
	}
	// A matching size and build ID prevent silently combining independently
	// generated primary and secondary indexes for the same seeds chunk.
	if be.Uint64(data[32:40]) != seedsSize || be.Uint64(data[40:48]) != uint64(buildID) {
		return fail(fmt.Errorf("idx15 does not match the seeds and primary index"))
	}
	return &index15{data: data, k: k, primaryPrefix: data[11], secondaryPrefix: data[12], seedsSize: seedsSize}, nil
}

func (idx *index15) close() error {
	if idx == nil || idx.data == nil {
		return nil
	}
	return idx.data.Unmap()
}

// lookup uses a longer prefix to locate the query's sub-block within a large
// primary block in a large index, reducing the seed data scanned during prefix
// matching.
// Shorter queries retain the primary block's starting checkpoint and offset.
func (idx *index15) lookup(taggedOffset, query uint64, prefixLength uint8) (uint64, uint64, bool, error) {
	if !IsIndex15Offset(taggedOffset) {
		return 0, 0, false, ErrInvalidFileFormat
	}
	start := RawIndexOffset(taggedOffset)
	if start > uint64(len(idx.data)) || uint64(len(idx.data))-start < index15BlockFixedBytes {
		return 0, 0, false, ErrBrokenFile
	}
	// Reslicing is metadata-only. The OS faults in mapped pages as the fields
	// below are accessed; this does not eagerly read the rest of idx15.
	block := idx.data[start:]
	width := block[0]
	if width < 1 || width > 8 {
		return 0, 0, false, ErrBrokenFile
	}
	blockSize := uint64(index15BlockFixedBytes + 15*int(width))
	if uint64(len(block)) < blockSize {
		return 0, 0, false, ErrBrokenFile
	}
	baseOffset := be.Uint64(block[3:11])
	baseCheckpoint := be.Uint64(block[11:19])
	if baseOffset >= idx.seedsSize {
		return 0, 0, false, ErrBrokenFile
	}
	if prefixLength < idx.secondaryPrefix {
		return baseCheckpoint, baseOffset, false, nil
	}
	subprefix := uint8(query >> ((idx.k - idx.secondaryPrefix) << 1) & 15)
	if be.Uint16(block[1:3])&(uint16(1)<<subprefix) == 0 {
		return 0, 0, true, nil
	}
	if subprefix == 0 {
		return baseCheckpoint, baseOffset, false, nil
	}
	slot := int(subprefix - 1)
	relative := uintWidth(block[index15BlockFixedBytes+slot*int(width):], width)
	if relative > math.MaxUint64-baseOffset || baseOffset+relative >= idx.seedsSize {
		return 0, 0, false, ErrBrokenFile
	}
	if relative == 0 {
		return baseCheckpoint, baseOffset, false, nil
	}
	suffix := uint40(block[19+slot*index15SuffixBytes:])
	suffixBits := (idx.k - idx.primaryPrefix) << 1
	if suffixBits < 64 && suffix >= uint64(1)<<suffixBits {
		return 0, 0, false, ErrBrokenFile
	}
	checkpoint := query>>suffixBits<<suffixBits | suffix
	if checkpoint>>suffixBits != query>>suffixBits {
		return 0, 0, false, ErrBrokenFile
	}
	return checkpoint, baseOffset + relative, false, nil
}

// readTaggedIndexBuildID reads the 16-bit idx/idx15 pairing token from the
// final two bytes of the tagged primary-index metadata.
func readTaggedIndexBuildID(indexFile string) (uint16, error) {
	fh, err := os.Open(indexFile)
	if err != nil {
		return 0, err
	}
	defer fh.Close()
	var header [16]byte
	if _, err = io.ReadFull(fh, header[:]); err != nil {
		return 0, err
	}
	if string(header[:8]) != string(MagicIdx[:]) || header[8] != IndexMainVersionTagged || header[13]&MaskHasIndex15 == 0 {
		return 0, ErrInvalidFileFormat
	}
	return be.Uint16(header[14:16]), nil
}
