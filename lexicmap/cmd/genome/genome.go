// Copyright © 2023-2024 Wei Shen <shenwei356@gmail.com>
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

package genome

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"strings"
	"sync"
	"time"
)

var be = binary.BigEndian

// Magic number for checking file format
var Magic = [8]byte{'.', 'g', 'e', 'n', 'o', 'm', 'e', 's'}

// Magic number for the index file
var MagicIdx = [8]byte{'.', 'g', 'e', 'n', 'o', 'm', 'e', 'i'}

// KVIndexFileExt is the file extension of k-mer data index file.
var GenomeIndexFileExt = ".idx"

// MainVersion is use for checking compatibility
var MainVersion uint8 = 0

// MinorVersion is less important
var MinorVersion uint8 = 1

// BufferSize is size of reading and writing buffer
var BufferSize = 65536 // os.Getpagesize()

// ErrInvalidFileFormat means invalid file format.
var ErrInvalidFileFormat = errors.New("genome data: invalid binary format")

// ErrEmptySeq means the sequence is empty
var ErrEmptySeq = errors.New("genome data: empty seq")

// ErrInvalidTwoBitData means the length of two bit seq slice does not match the number of bases
var ErrInvalidTwoBitData = errors.New("genome data: invalid two-bit data")

// ErrBrokenFile means the file is not complete.
var ErrBrokenFile = errors.New("genome data: broken file")

// ErrVersionMismatch means version mismatch between files and program
var ErrVersionMismatch = errors.New("genome data: version mismatch")

// Genome represents a reference sequence to insert and a matched subsequence
type Genome struct {
	ID  []byte // genome ID
	Seq []byte // sequence, bases

	GenomeSize int       // bases of all sequences
	Len        int       // length of contatenated sequences
	NumSeqs    int       // number of sequences
	SeqSizes   []int     // sizes of sequences
	SeqIDs     []*[]byte // IDs of all sequences

	// only used in index building
	Kmers     *[]uint64 // lexichash mask result
	Locses    *[][]int  // lexichash mask result
	TwoBit    *[]byte   // bit-packed sequence
	StartTime time.Time

	// seed positions to write to the file
	Locs []uint32
}

func (r Genome) String() string {
	return fmt.Sprintf("%s, genomeSize:%d, len:%d, contigs:%d", r.ID, r.GenomeSize, r.Len, r.NumSeqs)
}

// PoolGenome is the object pool for Genome
var PoolGenome = &sync.Pool{New: func() interface{} {
	return &Genome{
		ID:  make([]byte, 0, 128),
		Seq: make([]byte, 0, 10<<20),

		GenomeSize: 0,
		SeqSizes:   make([]int, 0, 128),
	}
}}

// Reset resets the Genome.
func (r *Genome) Reset() {
	r.ID = r.ID[:0]
	r.Seq = r.Seq[:0]
	r.GenomeSize = 0
	r.Len = 0
	r.NumSeqs = 0
	r.SeqSizes = r.SeqSizes[:0]
	r.SeqIDs = r.SeqIDs[:0]

	if r.Locs != nil {
		r.Locs = r.Locs[:0]
	}
}

// RecycleGenome recycle a Genome
func RecycleGenome(g *Genome) {
	if g.TwoBit != nil {
		RecycleTwoBit(g.TwoBit)
	}
	for _, id := range g.SeqIDs {
		poolID.Put(id)
	}
	PoolGenome.Put(g)
}

var poolID = &sync.Pool{New: func() interface{} {
	tmp := make([]byte, 128)
	return &tmp
}}

// Writer saves a list of DNA sequences into 2bit-encoded format,
// along with its genome information.
type Writer struct {
	batch uint32
	file  string
	fh    *os.File
	w     *bufio.Writer

	bBuf   bytes.Buffer
	buf    []byte // 24 bytes buffer
	offset int

	// offsets
	index [][2]int
}

// NewWriter creates a new Writer.
// Batch is the batch id for this data file.
func NewWriter(file string, batch uint32) (*Writer, error) {
	w := &Writer{
		batch: batch,
		file:  file,
		index: make([][2]int, 0, 1024),
	}
	var err error
	w.fh, err = os.Create(file)
	if err != nil {
		return nil, err
	}
	w.w = bufio.NewWriterSize(w.fh, BufferSize)

	w.buf = make([]byte, 24)

	// 8-byte magic number
	err = binary.Write(w.w, be, Magic)
	if err != nil {
		return nil, err
	}
	w.offset += 8

	// 8-byte meta info
	// actually, only 2 bytes used and the left 6 bytes is preserved.
	err = binary.Write(w.w, be, [8]uint8{MainVersion, MinorVersion})
	if err != nil {
		return nil, err
	}
	w.offset += 8
	return w, nil
}

// Write writes one genome.
// After calling this, you need to call RecycleGenome to recycle the genome.
func (w *Writer) Write(s *Genome) error {
	// collect data for the index file
	w.index = append(w.index, [2]int{w.offset, s.Len})

	// write genome information
	buf := w.buf
	buf0 := w.bBuf
	buf0.Reset()

	// ID length
	be.PutUint16(buf[:2], uint16(len(s.ID)))
	buf0.Write(buf[:2])

	// ID
	buf0.Write(s.ID)

	// meta
	be.PutUint32(buf[:4], uint32(s.GenomeSize))      // genome size
	be.PutUint32(buf[4:8], uint32(s.Len))            // length of contatenated sequences
	be.PutUint32(buf[8:12], uint32(len(s.SeqSizes))) // number of contigs
	buf0.Write(buf[:12])

	var seqid []byte
	for i, size := range s.SeqSizes {
		// seq sizes
		be.PutUint32(buf[:4], uint32(size))
		buf0.Write(buf[:4])

		// seq ids
		seqid = *s.SeqIDs[i]
		be.PutUint16(buf[:2], uint16(len(seqid))) // length of id
		buf0.Write(buf[:2])
		buf0.Write(seqid)
	}

	// write sequence
	b2 := s.TwoBit
	var newTwoBit bool
	if b2 == nil {
		b2 = Seq2TwoBit(s.Seq)
		newTwoBit = true
	}
	bases := len(s.Seq)
	nbytes := len(*b2)
	// possible bases for b2 of n bytes: [n*4-3, n*4]
	if bases < (nbytes<<2)-3 || bases > nbytes<<2 {
		return ErrInvalidTwoBitData
	}

	// the number of bytes and bases
	be.PutUint32(buf[:4], uint32(nbytes))
	be.PutUint32(buf[4:8], uint32(bases))
	buf0.Write(buf[:8])

	// write metadata to file
	_, err := w.w.Write(buf0.Bytes())
	if err != nil {
		return err
	}
	// write 2bit-packed data to file
	_, err = w.w.Write(*b2)
	if err != nil {
		return err
	}
	w.offset += buf0.Len() + nbytes

	if newTwoBit {
		poolTwoBit.Put(b2)
	}
	return err
}

// Close writes the index file and finish writing.
func (w *Writer) Close() error {
	err := w.w.Flush()
	if err != nil {
		return err
	}

	err = w.fh.Close()
	if err != nil {
		return err
	}

	// write the index

	fh, err := os.Create(filepath.Clean(w.file) + GenomeIndexFileExt)
	if err != nil {
		return err
	}
	wtr := bufio.NewWriterSize(fh, BufferSize)

	// magic
	err = binary.Write(wtr, be, MagicIdx)
	if err != nil {
		return err
	}

	// versions
	// actually, only 2 bytes used and the left 6 bytes is preserved.
	err = binary.Write(wtr, be, [8]uint8{MainVersion, MinorVersion})
	if err != nil {
		return err
	}

	buf := w.buf
	buf0 := w.bBuf
	buf0.Reset()

	// batch number
	be.PutUint32(buf[:4], w.batch)
	// the number of records
	be.PutUint32(buf[4:8], uint32(len(w.index)))
	buf0.Write(buf[:8])

	buf = w.buf[:12]
	for _, data := range w.index {
		be.PutUint64(buf[:8], uint64(data[0]))   // offset
		be.PutUint32(buf[8:12], uint32(data[1])) // bases
		buf0.Write(buf)
	}

	_, err = wtr.Write(buf0.Bytes())
	if err != nil {
		return err
	}

	err = wtr.Flush()
	if err != nil {
		return err
	}

	return fh.Close()
}

// Reader is for fast extracting of subsequence of any sequence in the data file.
type Reader struct {
	batch uint32
	nSeqs uint32

	fh     *os.File
	offset int // offset of the first index record

	buf []byte

	fhData *os.File
}

var poolReader = &sync.Pool{New: func() interface{} {
	return &Reader{
		buf: make([]byte, 10<<10), // 10k
	}
}}

// NewReader returns a reader from a genome file.
// The reader is recycled after calling Close().
func NewReader(file string) (*Reader, error) {
	if strings.HasSuffix(file, GenomeIndexFileExt) {
		return nil, fmt.Errorf("genome file, not the index file should be given")
	}

	// ------------ genome index file ----------------

	fileIndex := filepath.Clean(file) + GenomeIndexFileExt
	var err error
	r := poolReader.Get().(*Reader)

	r.fh, err = os.Open(fileIndex)
	if err != nil {
		return nil, err
	}

	buf := r.buf

	// check the magic number
	n, err := io.ReadFull(r.fh, buf[:8])
	if err != nil {
		return nil, err
	}
	if n < 8 {
		return nil, ErrBrokenFile
	}
	same := true
	for i := 0; i < 8; i++ {
		if MagicIdx[i] != buf[i] {
			same = false
			break
		}
	}
	if !same {
		return nil, ErrInvalidFileFormat
	}
	r.offset += 8

	// read metadata
	n, err = io.ReadFull(r.fh, buf[:8])
	if err != nil {
		return nil, err
	}
	if n < 8 {
		return nil, ErrBrokenFile
	}
	r.offset += 8

	// check compatibility
	if MainVersion != buf[0] {
		return nil, ErrVersionMismatch
	}

	// batch number and the number seqs
	n, err = io.ReadFull(r.fh, buf[:8])
	if err != nil {
		return nil, err
	}
	if n < 8 {
		return nil, ErrBrokenFile
	}
	r.offset += 8

	r.batch = be.Uint32(buf[:4])
	r.nSeqs = be.Uint32(buf[4:8])

	// ------------ genome data file ----------------

	r.fhData, err = os.Open(file)
	if err != nil {
		return nil, err
	}

	return r, nil
}

// Close closes and recycles the reader.
func (r *Reader) Close() error {
	err := r.fh.Close()
	if err != nil {
		poolReader.Put(r)
		return err
	}

	err = r.fhData.Close()
	if err != nil {
		poolReader.Put(r)
		return err
	}

	poolReader.Put(r)
	return nil
}

// Seq returns the sequence with index of genome (0-based).
func (r *Reader) Seq(idx int) (*Genome, error) {
	return r.SubSeq(idx, 0, math.MaxUint32)
}

// SubSeq returns the subsequence of a genome (idx is 0-based),
// from start to end (both are 0-based and included).
// Please call RecycleGenome() after using the result.
func (r *Reader) SubSeq(idx int, start int, end int) (*Genome, error) {
	if idx < 0 || idx >= int(r.nSeqs) {
		return nil, fmt.Errorf("sequence index (%d) out of range: [0, %d]", idx, int(r.nSeqs)-1)
	}

	buf := r.buf

	// -----------------------------------------------------------
	// read index information
	// 24 + 12 * idx
	r.fh.Seek(int64(r.offset)+int64(idx)<<3+int64(idx)<<2, 0)

	// offset in the data file and bases
	n, err := io.ReadFull(r.fh, buf[:12])
	if err != nil {
		return nil, err
	}
	if n < 12 {
		return nil, ErrBrokenFile
	}
	offset := int64(be.Uint64(buf[:8]))
	nBases := int(be.Uint32(buf[8:12])) // for check end

	if start < 0 {
		start = 0
	}
	if end >= nBases-1 {
		end = nBases - 1
	}
	if end < start {
		end = start
	}

	// -----------------------------------------------------------
	// get sequence information

	g := PoolGenome.Get().(*Genome)

	r.fhData.Seek(offset, 0)

	// ID length
	n, err = io.ReadFull(r.fhData, buf[:2])
	if err != nil {
		return nil, err
	}
	if n < 2 {
		return nil, ErrBrokenFile
	}
	idLen := be.Uint16(buf[:2])
	offset += 2

	// ID
	n, err = io.ReadFull(r.fhData, buf[:idLen])
	if err != nil {
		return nil, err
	}
	if n < int(idLen) {
		return nil, ErrBrokenFile
	}
	g.ID = g.ID[:0]
	g.ID = append(g.ID, buf[:idLen]...)
	offset += int64(idLen)

	// genome size, Len of concatenated seqs, NumSeqs
	n, err = io.ReadFull(r.fhData, buf[:12])
	if err != nil {
		return nil, err
	}
	if n < 12 {
		return nil, ErrBrokenFile
	}
	g.GenomeSize = int(be.Uint32(buf[:4]))
	g.Len = int(be.Uint32(buf[4:8]))
	g.NumSeqs = int(be.Uint32(buf[8:12]))
	offset += 12

	// SeqSizes and SeqIDs
	g.SeqSizes = g.SeqSizes[:0]
	g.SeqIDs = g.SeqIDs[:0]
	var j, nappend int
	var idLen2 int
	for i := 0; i < g.NumSeqs; i++ {
		n, err = io.ReadFull(r.fhData, buf[:4])
		if err != nil {
			return nil, err
		}
		if n < 4 {
			return nil, ErrBrokenFile
		}
		g.SeqSizes = append(g.SeqSizes, int(be.Uint32(buf[:4])))

		// seq id
		n, err = io.ReadFull(r.fhData, buf[:2])
		if err != nil {
			return nil, err
		}
		if n < 2 {
			return nil, ErrBrokenFile
		}

		idLen2 = int(be.Uint16(buf[:2]))
		id := poolID.Get().(*[]byte)
		if len(*id) >= idLen2 {
			*id = (*id)[:idLen2]
		} else {
			nappend = idLen2 - len(*id)
			for j = 0; j < nappend; j++ {
				*id = append(*id, 0)
			}
		}
		n, err = io.ReadFull(r.fhData, *id)
		if err != nil {
			return nil, err
		}
		if n < idLen2 {
			return nil, ErrBrokenFile
		}
		g.SeqIDs = append(g.SeqIDs, id)

		offset += int64(6 + idLen2)
	}

	// get sequence

	// start of byte, 8 is #bytes+#bases
	offset += 8 + int64(start>>2)
	_, err = r.fhData.Seek(offset, 0)
	if err != nil {
		return nil, err
	}

	nBytes := end>>2 - start>>2 + 1

	// prepair the buf
	if nBytes <= len(r.buf) {
		buf = r.buf[:nBytes]
	} else {
		n := nBytes - len(r.buf)
		for i := 0; i < n; i++ {
			r.buf = append(r.buf, 0)
		}
		buf = r.buf
	}
	n, err = io.ReadFull(r.fhData, buf)
	if err != nil {
		return nil, err
	}

	if n < nBytes {
		return nil, ErrBrokenFile
	}

	l := end - start + 1

	// initialize with l+4 blank values, because if there less than 4 bases
	// to extract, code below would panic.
	s := &g.Seq
	*s = (*s)[:4]

	// -- first byte --
	b := buf[0]
	j = start & 3

	switch j {
	case 0:
		(*s)[3] = bit2base[b&3]
		b >>= 2
		(*s)[2] = bit2base[b&3]
		b >>= 2
		(*s)[1] = bit2base[b&3]
		b >>= 2
		(*s)[0] = bit2base[b&3]
	case 1:
		(*s)[2] = bit2base[b&3]
		b >>= 2
		(*s)[1] = bit2base[b&3]
		b >>= 2
		(*s)[0] = bit2base[b&3]
	case 2:
		(*s)[1] = bit2base[b&3]
		b >>= 2
		(*s)[0] = bit2base[b&3]
	case 3:
		(*s)[0] = bit2base[b&3]
	}
	j = 4 - j
	*s = (*s)[:j]
	if j >= l {
		*s = (*s)[:l]
		return g, nil
	}

	// -- middle byte --
	if nBytes > 2 {
		for _, b = range buf[1 : nBytes-1] {
			*s = append(*s, bit2base[b>>6&3])
			*s = append(*s, bit2base[b>>4&3])
			*s = append(*s, bit2base[b>>2&3])
			*s = append(*s, bit2base[b&3])
		}
	}

	if nBytes > 1 {
		// -- last byte --
		b = buf[nBytes-1]
		j = end & 3
		switch j {
		case 0:
			*s = append(*s, bit2base[b>>6&3])
		case 1:
			*s = append(*s, bit2base[b>>6&3])
			*s = append(*s, bit2base[b>>4&3])
		case 2:
			*s = append(*s, bit2base[b>>6&3])
			*s = append(*s, bit2base[b>>4&3])
			*s = append(*s, bit2base[b>>2&3])
		case 3:
			*s = append(*s, bit2base[b>>6&3])
			*s = append(*s, bit2base[b>>4&3])
			*s = append(*s, bit2base[b>>2&3])
			*s = append(*s, bit2base[b&3])
		}
	}

	*s = (*s)[:l]
	g.Len = len(g.Seq)
	return g, nil
}

// SubSeq2 returns the subsequence of one genome (idx is 0-based),
// from start to end (both are 0-based and included).
// Please call RecycleGenome() after using the result.
func (r *Reader) SubSeq2(idx int, seqid []byte, start int, end int) (*Genome, error) {
	if idx < 0 || idx >= int(r.nSeqs) {
		return nil, fmt.Errorf("genome index (%d) out of range: [0, %d]", idx, int(r.nSeqs)-1)
	}

	buf := r.buf

	// -----------------------------------------------------------
	// read index information
	// 24 + 12 * idx
	r.fh.Seek(int64(r.offset)+int64(idx)<<3+int64(idx)<<2, 0)

	// offset in the data file and bases
	n, err := io.ReadFull(r.fh, buf[:12])
	if err != nil {
		return nil, err
	}
	if n < 12 {
		return nil, ErrBrokenFile
	}
	offset := int64(be.Uint64(buf[:8]))
	nBases := int(be.Uint32(buf[8:12])) // for check end

	if start < 0 {
		start = 0
	}
	if end >= nBases-1 {
		end = nBases - 1
	}
	if end < start {
		end = start
	}

	// -----------------------------------------------------------
	// get sequence information

	g := PoolGenome.Get().(*Genome)

	r.fhData.Seek(offset, 0)

	// ID length
	n, err = io.ReadFull(r.fhData, buf[:2])
	if err != nil {
		return nil, err
	}
	if n < 2 {
		return nil, ErrBrokenFile
	}
	idLen := be.Uint16(buf[:2])
	offset += 2

	// ID
	n, err = io.ReadFull(r.fhData, buf[:idLen])
	if err != nil {
		return nil, err
	}
	if n < int(idLen) {
		return nil, ErrBrokenFile
	}
	g.ID = g.ID[:0]
	g.ID = append(g.ID, buf[:idLen]...)
	offset += int64(idLen)

	// genome size, Len of concatenated seqs, NumSeqs
	n, err = io.ReadFull(r.fhData, buf[:12])
	if err != nil {
		return nil, err
	}
	if n < 12 {
		return nil, ErrBrokenFile
	}
	g.GenomeSize = int(be.Uint32(buf[:4]))
	g.Len = int(be.Uint32(buf[4:8]))
	g.NumSeqs = int(be.Uint32(buf[8:12]))
	offset += 12

	// SeqSizes and SeqIDs
	g.SeqSizes = g.SeqSizes[:0]
	g.SeqIDs = g.SeqIDs[:0]
	var j, nappend int
	var idLen2 int

	// --------------------------------------------------
	var foundSeqID bool
	var seqSize, lenSum int
	var interval int
	if g.NumSeqs > 1 {
		interval = (g.Len - g.GenomeSize) / (g.NumSeqs - 1)
	}
	// --------------------------------------------------

	for i := 0; i < g.NumSeqs; i++ {
		n, err = io.ReadFull(r.fhData, buf[:4])
		if err != nil {
			return nil, err
		}
		if n < 4 {
			return nil, ErrBrokenFile
		}

		seqSize = int(be.Uint32(buf[:4]))
		g.SeqSizes = append(g.SeqSizes, seqSize)

		// seq id
		n, err = io.ReadFull(r.fhData, buf[:2])
		if err != nil {
			return nil, err
		}
		if n < 2 {
			return nil, ErrBrokenFile
		}

		idLen2 = int(be.Uint16(buf[:2]))
		id := poolID.Get().(*[]byte)
		if len(*id) >= idLen2 {
			*id = (*id)[:idLen2]
		} else {
			nappend = idLen2 - len(*id)
			for j = 0; j < nappend; j++ {
				*id = append(*id, 0)
			}
		}
		n, err = io.ReadFull(r.fhData, *id)
		if err != nil {
			return nil, err
		}
		if n < idLen2 {
			return nil, ErrBrokenFile
		}
		g.SeqIDs = append(g.SeqIDs, id)

		offset += int64(6 + idLen2)

		// --------------------------------------------------
		if bytes.Equal(*id, seqid) { // found it!
			foundSeqID = true
			start += lenSum
			end += lenSum
		} else {
			lenSum += seqSize
			lenSum += interval
		}
		// --------------------------------------------------
	}
	// --------------------------------------------------
	if !foundSeqID {
		return nil, fmt.Errorf("seqid not found: %s", seqid)
	}
	// --------------------------------------------------

	// get sequence

	// start of byte, 8 is #bytes+#bases
	offset += 8 + int64(start>>2)
	_, err = r.fhData.Seek(offset, 0)
	if err != nil {
		return nil, err
	}

	nBytes := end>>2 - start>>2 + 1

	// prepair the buf
	if nBytes <= len(r.buf) {
		buf = r.buf[:nBytes]
	} else {
		n := nBytes - len(r.buf)
		for i := 0; i < n; i++ {
			r.buf = append(r.buf, 0)
		}
		buf = r.buf
	}
	n, err = io.ReadFull(r.fhData, buf)
	if err != nil {
		return nil, err
	}

	if n < nBytes {
		return nil, ErrBrokenFile
	}

	l := end - start + 1

	// initialize with l+4 blank values, because if there less than 4 bases
	// to extract, code below would panic.
	s := &g.Seq
	*s = (*s)[:4]

	// -- first byte --
	b := buf[0]
	j = start & 3

	switch j {
	case 0:
		(*s)[3] = bit2base[b&3]
		b >>= 2
		(*s)[2] = bit2base[b&3]
		b >>= 2
		(*s)[1] = bit2base[b&3]
		b >>= 2
		(*s)[0] = bit2base[b&3]
	case 1:
		(*s)[2] = bit2base[b&3]
		b >>= 2
		(*s)[1] = bit2base[b&3]
		b >>= 2
		(*s)[0] = bit2base[b&3]
	case 2:
		(*s)[1] = bit2base[b&3]
		b >>= 2
		(*s)[0] = bit2base[b&3]
	case 3:
		(*s)[0] = bit2base[b&3]
	}
	j = 4 - j
	*s = (*s)[:j]
	if j >= l {
		*s = (*s)[:l]
		return g, nil
	}

	// -- middle byte --
	if nBytes > 2 {
		for _, b = range buf[1 : nBytes-1] {
			*s = append(*s, bit2base[b>>6&3])
			*s = append(*s, bit2base[b>>4&3])
			*s = append(*s, bit2base[b>>2&3])
			*s = append(*s, bit2base[b&3])
		}
	}

	if nBytes > 1 {
		// -- last byte --
		b = buf[nBytes-1]
		j = end & 3
		switch j {
		case 0:
			*s = append(*s, bit2base[b>>6&3])
		case 1:
			*s = append(*s, bit2base[b>>6&3])
			*s = append(*s, bit2base[b>>4&3])
		case 2:
			*s = append(*s, bit2base[b>>6&3])
			*s = append(*s, bit2base[b>>4&3])
			*s = append(*s, bit2base[b>>2&3])
		case 3:
			*s = append(*s, bit2base[b>>6&3])
			*s = append(*s, bit2base[b>>4&3])
			*s = append(*s, bit2base[b>>2&3])
			*s = append(*s, bit2base[b&3])
		}
	}

	*s = (*s)[:l]
	g.Len = len(g.Seq)
	return g, nil
}

var base2bit = [256]uint8{
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 1, 1, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0,
	0, 0, 0, 1, 3, 3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
	0, 0, 1, 1, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0,
	0, 0, 0, 1, 3, 3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
}

var bit2base = [4]byte{'A', 'C', 'G', 'T'}

// RecycleSeq recycles the sequence.
func RecycleTwoBit(b2 *[]byte) {
	poolTwoBit.Put(b2)
}

var poolTwoBit = &sync.Pool{New: func() interface{} {
	tmp := make([]byte, 0, 1<<20)
	return &tmp
}}

// Seq2TwoBit converts a DNA sequence to 2bit-packed sequence.
func Seq2TwoBit(s []byte) *[]byte {
	if s == nil {
		return nil
	}
	if len(s) == 0 {
		return &[]byte{}
	}

	n := len(s) >> 2
	m := len(s) & 3

	// codes := make([]byte, n+1)
	codes := poolTwoBit.Get().(*[]byte)
	*codes = (*codes)[:0]

	var j int
	for i := 0; i < n; i++ {
		j = i << 2

		*codes = append(*codes, base2bit[s[j]]<<6+base2bit[s[j+1]]<<4+base2bit[s[j+2]]<<2+base2bit[s[j+3]])
	}

	if m == 0 {
		tmp := (*codes)[:n]
		return &tmp
	}

	j = n << 2

	switch m {
	case 3:
		*codes = append(*codes, base2bit[s[j]]<<6+base2bit[s[j+1]]<<4+base2bit[s[j+2]]<<2)
	case 2:
		*codes = append(*codes, base2bit[s[j]]<<6+base2bit[s[j+1]]<<4)
	case 1:
		*codes = append(*codes, base2bit[s[j]]<<6)
	}

	return codes
}

// TwoBit2Seq converts a 2bit-packed sequence to DNA.
func TwoBit2Seq(b2 []byte, bases int) ([]byte, error) {
	// possible bases for b2 of n bytes: [n*4-3, n*4]
	if bases < (len(b2)<<2)-3 || bases > len(b2)<<2 {
		return nil, ErrInvalidTwoBitData
	}

	s := make([]byte, bases)
	n := len(s) >> 2
	m := bases & 3
	var b byte
	var j int
	for i := 0; i < n; i++ {
		b = b2[i]
		j = i << 2

		s[j+3] = bit2base[b&3]
		b >>= 2
		s[j+2] = bit2base[b&3]
		b >>= 2
		s[j+1] = bit2base[b&3]
		b >>= 2
		s[j] = bit2base[b&3]
	}
	if m == 0 {
		return s, nil
	}

	b = b2[n]
	j = n << 2
	switch m {
	case 1:
		s[j] = bit2base[b>>6&3]
	case 2:
		b >>= 4
		s[j+1] = bit2base[b&3]
		b >>= 2
		s[j] = bit2base[b&3]
	case 3:
		b >>= 2
		s[j+2] = bit2base[b&3]
		b >>= 2
		s[j+1] = bit2base[b&3]
		b >>= 2
		s[j] = bit2base[b&3]
	}

	return s, nil
}
