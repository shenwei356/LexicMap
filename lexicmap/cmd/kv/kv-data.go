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

package kv

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"sync"

	"github.com/pkg/errors"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
	"github.com/twotwotwo/sorts/sortutil"
)

var be = binary.BigEndian

// Magic number for checking file format
var Magic = [8]byte{'.', 'k', 'v', '-', 'd', 'a', 't', 'a'}

// Magic number for the index file
var MagicIdx = [8]byte{'.', 'k', 'v', 'i', 'n', 'd', 'e', 'x'}

// KVIndexFileExt is the file extension of k-mer data index file.
var KVIndexFileExt = ".idx"

// MainVersion is use for checking compatibility
var MainVersion uint8 = 1

// MinorVersion is less important
var MinorVersion uint8 = 0

// ErrInvalidFileFormat means invalid file format.
var ErrInvalidFileFormat = errors.New("k-mer-value data: invalid binary format")

// ErrBrokenFile means the file is not complete.
var ErrBrokenFile = errors.New("k-mer-value data: broken file")

// ErrKOverflow means K < 1 or K > 32.
var ErrKOverflow = errors.New("k-mer-value data: k-mer size [1, 32] overflow")

// ErrVersionMismatch means version mismatch between files and program
var ErrVersionMismatch = errors.New("k-mer-value data: version mismatch")

// WriteKVData writes k-mer-value data of a chunk of masks to a file.
// At the same time, the index file is also created with the number of
// anchors `nAnchors` (default: sqrt(#kmers)).
//
// Header (32 bytes):
//
//	Magic number, 8 bytes, ".kv-data".
//	Main and minor versions, 2 bytes.
//	K size, 1 byte.
//	Blank, 5 bytes.
//	Mask start index, 8 bytes. The index of the first index.
//	Mask chunk size, 8 bytes. The number of masks in this file.
//
// For each mask:
//
//	Number of k-mers, 8 bytes.
//	k-mer-values pairs, sizes vary.
//		Control byte for 2 k-mers, 1 byte
//		Delta values of the 2 k-mers, 2-16 bytes
//		Control byte for numbers of values, 1 byte
//		Numbers of values of the 2 k-mers, 2-16 bytes, 2 bytes for most cases.
//		Values of the 2 k-mers, 8*n bytes, 16 bytes for most cases.
//
// Index file stores 4^p' k-mers (anchors) and their offsets in
// the kv-data file for fast access, the time complexity would be O(1) instead of previous O(log2N)
//
//	AAAAAAA CCCCC NNNNNNNN
//	------- -----
//	p       p'
//
// Locations of these anchors vary, and some of them might not exist.
//
//	kkvvkkvvkkvvkkvvkkvvkkvvkkvvkkvvkkvvkkvvkkvv
//	AA        AC       AG   AT   CA
//
// Header (40 bytes):
//
//	Magic number, 8 bytes, ".kvindex".
//	Main and minor versions, 2 bytes.
//	K size, 1 byte.
//	Mask prefix length, 1  byte. e.g., 7
//	Anchor prefix length, 1 byte. e.g., 5
//	Blank, 3 bytes.
//	Mask start index, 8 bytes. The index of the first index.
//	Mask chunk size, 8 bytes. The number of masks in this file.
//
// For each mask:
//
//	Number of anchors, 8 bytes.
//	kmer-offset data:
//
//		k-mer: 8 bytes
//		offset: 8 bytes
func WriteKVData(k uint8, MaskOffset int, data []*map[uint64]*[]uint64, file string, maskPrefix uint8, anchorPrefix uint8) (int, error) {
	if len(data) == 0 {
		return 0, errors.New("k-mer-value data: no data given")
	}

	// some masks might do not capture any k-mers from a short genome

	// if len(*data[0]) == 0 {
	// 	return 0, errors.New("k-mer-value data: no data given")
	// }

	wtr, err := NewWriter(k, MaskOffset, len(data), file, maskPrefix, anchorPrefix)
	if err != nil {
		return 0, err
	}

	for _, m := range data {
		err = wtr.WriteDataOfAMask(*m)
		if err != nil {
			return 0, err
		}
	}
	err = wtr.Close()
	if err != nil {
		return 0, err
	}

	return wtr.N, nil
}

// Writer is used for k-mer-value data for multiple mask
type Writer struct {
	K          uint8 // kmer size
	ChunkIndex int   // index of the first mask in this chunk
	ChunkSize  int   // the number of masks in this chunk

	// bufers
	bufVar []byte // needs at most 8+8=16
	buf    []byte // needs at most 1+16+1+16=34
	buf8   []byte // for writing uint8

	// for kv data
	N  int // the number of bytes.
	fh *os.File
	w  *bufio.Writer

	// for index file
	fhi          *os.File
	wi           *bufio.Writer
	maskPrefix   uint8
	anchorPrefix uint8
	poolP2O      *sync.Pool
	getAnchor    func(uint64) uint64
}

// Close is very important
func (wtr *Writer) Close() (err error) {
	err = wtr.w.Flush()
	if err != nil {
		return err
	}
	err = wtr.fh.Close()
	if err != nil {
		return err
	}
	err = wtr.wi.Flush()
	if err != nil {
		return err
	}
	err = wtr.fhi.Close()
	if err != nil {
		return err
	}
	return nil
}

// NewWriter returns a new writer
func NewWriter(k uint8, MaskOffset int, chunkSize int, file string, maskPrefix uint8, anchorPrefix uint8) (*Writer, error) {
	if maskPrefix+anchorPrefix > k {
		return nil, fmt.Errorf("maskPrefix + anchorPrefix should be <= k")
	}
	if anchorPrefix == 0 {
		return nil, fmt.Errorf("anchorPrefix could not be 0")
	}

	// file handlers
	fh, err := os.Create(file)
	if err != nil {
		return nil, err
	}
	w := bufio.NewWriter(fh)
	fhi, err := os.Create(filepath.Clean(file) + KVIndexFileExt)
	if err != nil {
		return nil, err
	}
	wi := bufio.NewWriter(fhi)

	wtr := &Writer{
		K:          k,
		ChunkIndex: MaskOffset,
		ChunkSize:  chunkSize,
		fh:         fh,
		w:          w,
		fhi:        fhi,
		wi:         wi,

		maskPrefix:   maskPrefix,
		anchorPrefix: anchorPrefix,
		poolP2O: &sync.Pool{New: func() interface{} {
			// 2 is for recording the offset of the first k-mer whose prefix might not be AA
			tmp := make([]uint64, 2+int(math.Pow(4, float64(anchorPrefix)))<<1)
			return &tmp
		}},
		getAnchor: AnchorExtracter(k, maskPrefix, anchorPrefix),

		bufVar: make([]byte, 16),
		buf:    make([]byte, 36),
		buf8:   make([]byte, 8),
	}

	// ---------------------------------------------------------------------------

	var N int // the number of bytes.

	// 8-byte magic number
	err = binary.Write(w, be, Magic)
	if err != nil {
		return nil, err
	}
	N += 8

	// 8-byte meta info
	err = binary.Write(w, be, [8]uint8{MainVersion, MinorVersion, k})
	if err != nil {
		return nil, err
	}
	N += 8

	// 16-byte the MaskOffset and the chunk size
	err = binary.Write(w, be, [2]uint64{uint64(MaskOffset), uint64(chunkSize)})
	if err != nil {
		return nil, err
	}
	N += 16

	// --------------------- Index file

	// 8-byte magic number
	err = binary.Write(wi, be, MagicIdx)
	if err != nil {
		return nil, err
	}

	// 8-byte meta info
	err = binary.Write(wi, be, [8]uint8{MainVersion, MinorVersion, k, maskPrefix, anchorPrefix})
	if err != nil {
		return nil, err
	}

	// 16-byte the MaskOffset and the chunk size
	err = binary.Write(wi, be, [2]uint64{uint64(MaskOffset), uint64(chunkSize)})
	if err != nil {
		return nil, err
	}

	// ---------------------------------------------------------------------------

	wtr.N = N
	return wtr, nil
}

// AnchorExtracter returns the function for extracting anchors, i.e., CCCCC below
//
//	maskPrefix
//	-------
//	AAAAAAA CCCCC NNNNNNNN
//	        -----
//	        anchorPrefix
func AnchorExtracter(k uint8, maskPrefix uint8, anchorPrefix uint8) func(uint64) uint64 {
	shift := uint64(k-maskPrefix-anchorPrefix) << 1
	mask := uint64(1<<(anchorPrefix<<1)) - 1
	return func(kmer uint64) uint64 {
		return kmer >> shift & mask
	}
}

// WriteDataOfAMask writes data of one mask.
func (wtr *Writer) WriteDataOfAMask(m map[uint64]*[]uint64) (err error) {
	var hasPrev bool
	var preKey, key, _v uint64
	var preVal, v *[]uint64
	var offset uint64
	var ctrlByteKey, ctrlByteVal byte
	var nBytesKey, nBytesVal, n int
	bufVar := wtr.bufVar // needs at most 8+8=16
	buf := wtr.buf       // needs at most 1+16+1+16=34
	buf8 := wtr.buf8     // for writing uint8
	bufVals := poolBytesBuffer.Get().(*bytes.Buffer)
	defer poolBytesBuffer.Put(bufVals)
	var even bool
	var i, nm1 int
	var j int

	nKmers := len(m)

	w := wtr.w
	wi := wtr.wi

	hasPrev = false
	offset = 0

	// 8-byte the number of k-mers
	err = binary.Write(w, be, uint64(nKmers))
	if err != nil {
		return err
	}
	wtr.N += 8

	if nKmers == 0 { // this hapens when no captured k-mer for a mask
		// 8-byte the number of anchors
		err = binary.Write(wi, be, uint64(0))
		if err != nil {
			return err
		}

		return nil
	}

	p2o := wtr.poolP2O.Get().(*[]uint64)
	clear(*p2o)
	(*p2o)[1] = uint64(wtr.N) << 1 // offset of the first k-mer

	keys := poolUint64s.Get().(*[]uint64)
	// sort keys
	*keys = (*keys)[:0]
	for key = range m {
		*keys = append(*keys, key)
	}
	sortutil.Uint64s(*keys)

	// for decide should we set flag for the last control byte of the last k-mer
	even = len(*keys)&1 == 0 // the number of kmers is even
	nm1 = len(*keys) - 1     // idx of the last element

	j = 0

	getAnchor := wtr.getAnchor
	var prefix, prefixPre uint64
	first := true

	for i, key = range *keys {
		v = m[key]

		if !hasPrev { // write it later
			preKey = key
			preVal = v
			hasPrev = true

			continue
		}

		// ------------------------------------------------------------------------
		// index anchor

		// key 1
		prefix = getAnchor(preKey)
		// if preKey == 2233842599699997050 {
		// 	fmt.Printf("i:%d, key1:%s, prefix:%s, v:%d, offset:%d\n", i, lexichash.MustDecode(preKey, wtr.K), lexichash.MustDecode(prefix, wtr.anchorPrefix), preVal, wtr.N)
		// }
		if first || prefix != prefixPre { // the first new prefix
			first = false

			j = int(prefix<<1) + 2
			(*p2o)[j], (*p2o)[j+1] = preKey, uint64(wtr.N)<<1
			// fmt.Printf("  %d, record %s, %d\n", j, lexichash.MustDecode(preKey, wtr.K), wtr.N)

			prefixPre = prefix
		}

		// key 2
		prefix = getAnchor(key)
		// if key == 2233842599699997050 {
		// 	fmt.Printf("i:%d, key2:%s, prefix:%s, v:%d, offset:%d\n", i, lexichash.MustDecode(key, wtr.K), lexichash.MustDecode(prefix, wtr.anchorPrefix), v, wtr.N)
		// }
		if prefix != prefixPre { // the first new prefix
			j = int(prefix<<1) + 2
			(*p2o)[j], (*p2o)[j+1] = key, uint64(wtr.N)<<1|1 // add a flag to mark it's the second k-mer
			// fmt.Printf("  %d, record %s, %d\n", j, lexichash.MustDecode(key, wtr.K), wtr.N)

			prefixPre = prefix
		}

		// ------------------------------------------------------------------------

		// 2 k-mers and numbers of values

		// only save key2 - key1, which is small so it could be saved in few bytes
		ctrlByteKey, nBytesKey = util.PutUint64s(bufVar, preKey-offset, key-preKey)
		if even && i == nm1 {
			// fmt.Printf("write last two kmers: %s, %s\n",
			// 	lexichash.MustDecode(preKey, k), lexichash.MustDecode(key, k))
			ctrlByteKey |= 1 << 7 // it means this is the last record(s) for this mask
		}
		buf[0] = ctrlByteKey
		copy(buf[1:nBytesKey+1], bufVar[:nBytesKey])
		n = nBytesKey + 1

		// save lengths of values
		ctrlByteVal, nBytesVal = util.PutUint64s(bufVar, uint64(len(*preVal)), uint64(len(*v)))
		buf[n] = ctrlByteVal
		copy(buf[n+1:n+nBytesVal+1], bufVar[:nBytesVal])
		n += nBytesVal + 1

		_, err = w.Write(buf[:n])
		if err != nil {
			return err
		}
		wtr.N += n

		// values

		bufVals.Reset()

		for _, _v = range *preVal {
			be.PutUint64(buf8, _v)
			bufVals.Write(buf8)
		}
		for _, _v = range *v {
			be.PutUint64(buf8, _v)
			bufVals.Write(buf8)
		}

		_, err = w.Write(bufVals.Bytes())
		if err != nil {
			return err
		}
		wtr.N += bufVals.Len()

		// update

		offset = key
		hasPrev = false
	}

	if hasPrev { // the last single one
		// ------------------------------------------------------------------------
		// index anchor

		// key 1
		prefix = getAnchor(preKey)
		if first || prefix != prefixPre { // the first new prefix
			first = false

			j = int(prefix<<1) + 2
			(*p2o)[j], (*p2o)[j+1] = preKey, uint64(wtr.N)<<1

			prefixPre = prefix
		}

		// ------------------------------------------------------------------------

		// fmt.Printf("write the last two kmer: %s\n",
		// 	lexichash.MustDecode(preKey, k))

		// 2 k-mers and numbers of values

		// only save key2 - key1, which is small so it could be saved in few bytes
		ctrlByteKey, nBytesKey = util.PutUint64s(bufVar, preKey-offset, 0)
		ctrlByteKey |= 1 << 7 // it means this is the last record(s) for this mask.
		ctrlByteKey |= 1 << 6 // it means this is the last single record
		buf[0] = ctrlByteKey
		copy(buf[1:nBytesKey+1], bufVar[:nBytesKey])
		n = nBytesKey + 1

		// save lengths of values
		ctrlByteVal, nBytesVal = util.PutUint64s(bufVar, uint64(len(*preVal)), 0)
		buf[n] = ctrlByteVal
		copy(buf[n+1:n+nBytesVal+1], bufVar[:nBytesVal])
		n += nBytesVal + 1

		_, err = w.Write(buf[:n])
		if err != nil {
			return err
		}
		wtr.N += n

		// values

		bufVals.Reset()
		for _, _v = range *preVal {
			be.PutUint64(buf8, _v)
			bufVals.Write(buf8)
		}

		_, err = w.Write(bufVals.Bytes())
		if err != nil {
			return err
		}
		wtr.N += bufVals.Len()
	}

	poolUint64s.Put(keys)

	// -----------------------------------------
	// save index

	var kmer uint64
	var nAnchors uint64
	e := len(*p2o) >> 1
	for i := 0; i < e; i++ {
		offset = (*p2o)[i<<1+1]
		if offset > 0 {
			nAnchors++
		}
	}
	// 8-byte the number of anchors
	err = binary.Write(wi, be, nAnchors) // including the extra 1
	if err != nil {
		return err
	}

	(*p2o)[0] = nAnchors // might be useful

	// k-mer and offset
	for i := 0; i < e; i++ {
		j = i << 1
		kmer, offset = (*p2o)[j], (*p2o)[j+1]
		if offset > 0 {
			be.PutUint64(buf[:8], kmer)     // k-mer
			be.PutUint64(buf[8:16], offset) // offset
			_, err = wi.Write(buf[:16])
			if err != nil {
				return err
			}
		}
	}

	return nil
}

// ReadKVIndex parses the k-mer-value index file.
//
// Returned:
//
//	k-mer size
//	Index (0-based) of the first Mask in current chunk.
//	Index data of masks saved in a list, the list size equals to the number of masks.
//	error
//
// A list of k-mer and offset pairs are intermittently saved in a []uint64.
// e.g., [k1, o1, k2, o2].
func ReadKVIndex(file string) (uint8, int, [][]uint64, uint8, uint8, error) {
	fh, err := os.Open(file)
	if err != nil {
		return 0, -1, nil, 0, 0, err
	}
	// r := bufio.NewReader(fh)
	r := poolBufReader.Get().(*bufio.Reader)
	r.Reset(fh)
	defer func() {
		poolBufReader.Put(r)
		fh.Close()
	}()

	// ---------------------------------------------

	buf := make([]byte, 8)
	buf16 := make([]byte, 16)

	var n int

	// check the magic number
	n, err = io.ReadFull(r, buf)
	if err != nil {
		return 0, -1, nil, 0, 0, err
	}
	if n < 8 {
		return 0, -1, nil, 0, 0, ErrBrokenFile
	}
	same := true
	for i := 0; i < 8; i++ {
		if MagicIdx[i] != buf[i] {
			same = false
			break
		}
	}
	if !same {
		return 0, -1, nil, 0, 0, ErrInvalidFileFormat
	}
	// read version information
	n, err = io.ReadFull(r, buf)
	if err != nil {
		return 0, -1, nil, 0, 0, err
	}
	if n < 8 {
		return 0, -1, nil, 0, 0, ErrBrokenFile
	}
	// check compatibility
	if MainVersion != buf[0] {
		return 0, -1, nil, 0, 0, ErrVersionMismatch
	}
	k := buf[2] // k-mer size
	maskPrefix := buf[3]
	anchorPrefix := buf[4]

	// index of the first mask in current chunk.
	var iFirstMask int
	_, err = io.ReadFull(r, buf)
	if err != nil {
		return 0, -1, nil, 0, 0, err
	}
	iFirstMask = int(be.Uint64(buf))

	// mask chunk size
	var nMasks int
	_, err = io.ReadFull(r, buf)
	if err != nil {
		return 0, -1, nil, 0, 0, err
	}
	nMasks = int(be.Uint64(buf))

	// the number of anchors
	var nAnchors int

	// ---------------------------------------------

	indexSize := 2 + int(math.Pow(4, float64(anchorPrefix)))<<1
	getAnchor := AnchorExtracter(k, maskPrefix, anchorPrefix)

	data := make([][]uint64, nMasks)

	var kmer, offset uint64
	var prefix uint64
	var j int
	var _j uint64
	for i := 0; i < nMasks; i++ {
		_, err = io.ReadFull(r, buf)
		if err != nil {
			return 0, -1, nil, 0, 0, err
		}
		nAnchors = int(be.Uint64(buf))

		if nAnchors == 0 { // this hapens when no captured k-mer for a mask
			data[i] = make([]uint64, 0)
			continue
		}

		// index := make([]uint64, 0, nAnchors<<1)
		index := make([]uint64, indexSize)
		for j = 0; j < nAnchors; j++ {
			_, err = io.ReadFull(r, buf16)
			if err != nil {
				return 0, -1, nil, 0, 0, err
			}
			kmer = be.Uint64(buf16[:8])

			// _, err = io.ReadFull(r, buf)
			// if err != nil {
			// 	return 0, -1, nil, 0, 0, err
			// }
			offset = be.Uint64(buf16[8:16])

			// index = append(index, kmer)
			// index = append(index, offset)
			if j == 0 {
				_j = 0
			} else {
				prefix = getAnchor(kmer)
				_j = prefix<<1 + 2
			}
			index[_j] = kmer
			index[_j+1] = offset
		}
		data[i] = index
	}

	return k, iFirstMask, data, maskPrefix, anchorPrefix, nil
}

// ReadKVIndexInfo read the information.
func ReadKVIndexInfo(file string) (uint8, int, int, uint8, uint8, error) {
	fh, err := os.Open(file)
	if err != nil {
		return 0, -1, 0, 0, 0, err
	}
	r := bufio.NewReader(fh)
	defer fh.Close()

	// ---------------------------------------------

	buf := make([]byte, 8)

	var n int

	// check the magic number
	n, err = io.ReadFull(r, buf)
	if err != nil {
		return 0, -1, 0, 0, 0, err
	}
	if n < 8 {
		return 0, -1, 0, 0, 0, ErrBrokenFile
	}
	same := true
	for i := 0; i < 8; i++ {
		if MagicIdx[i] != buf[i] {
			same = false
			break
		}
	}
	if !same {
		return 0, -1, 0, 0, 0, ErrInvalidFileFormat
	}
	// read version information
	n, err = io.ReadFull(r, buf)
	if err != nil {
		return 0, -1, 0, 0, 0, err
	}
	if n < 8 {
		return 0, -1, 0, 0, 0, ErrBrokenFile
	}
	// check compatibility
	if MainVersion != buf[0] {
		return 0, -1, 0, 0, 0, ErrVersionMismatch
	}
	k := buf[2] // k-mer size
	maskPrefix := buf[3]
	anchorPrefix := buf[4]

	// index of the first mask in current chunk.
	var iFirstMask int
	_, err = io.ReadFull(r, buf)
	if err != nil {
		return 0, -1, 0, 0, 0, err
	}
	iFirstMask = int(be.Uint64(buf))

	// mask chunk size
	var nMasks int
	_, err = io.ReadFull(r, buf)
	if err != nil {
		return 0, -1, 0, 0, 0, err
	}
	nMasks = int(be.Uint64(buf))

	return k, iFirstMask, nMasks, maskPrefix, anchorPrefix, nil
}

var poolBytesBuffer = &sync.Pool{New: func() interface{} {
	return &bytes.Buffer{}
}}

var poolUint64s = &sync.Pool{New: func() interface{} {
	tmp := make([]uint64, 0, 1<<20)
	return &tmp
}}

// CreateKVIndex recreates kv index file for the kv-data file.
func CreateKVIndex(file string, nAnchors int) error {
	fh, err := os.Open(file)
	if err != nil {
		return errors.Wrapf(err, "reading kv-data file")
	}

	r := bufio.NewReader(fh)

	// --------------------------------------------------------------------
	// header information of kv-data file

	var offset int

	buf8 := make([]byte, 8)
	buf := make([]byte, 64)

	var n int

	// check the magic number
	n, err = io.ReadFull(r, buf8)
	if err != nil {
		return err
	}
	if n < 8 {
		return ErrBrokenFile
	}
	same := true
	for i := 0; i < 8; i++ {
		if Magic[i] != buf8[i] {
			same = false
			break
		}
	}
	if !same {
		return ErrInvalidFileFormat
	}
	offset += 8

	// read version information
	n, err = io.ReadFull(r, buf8)
	if err != nil {
		return err
	}
	if n < 8 {
		return ErrBrokenFile
	}
	// check compatibility
	if MainVersion != buf8[0] {
		return ErrVersionMismatch
	}

	K := buf8[2] // k-mer size
	offset += 8

	// index of the first mask in current chunk.
	_, err = io.ReadFull(r, buf8)
	if err != nil {
		return err
	}
	ChunkIndex := be.Uint64(buf8)
	offset += 8

	// mask chunk size
	_, err = io.ReadFull(r, buf8)
	if err != nil {
		return err
	}
	ChunkSize := be.Uint64(buf8)
	offset += 8

	// --------------------------------------------------------------------
	// writer of kv-index file

	indexFile := filepath.Clean(file) + KVIndexFileExt
	_, _, _, maskPrefix, _, err := ReadKVIndexInfo(indexFile)
	if err != nil {
		return err
	}

	anchorPrefix := 0
	partitions := nAnchors
	for partitions > 0 {
		partitions >>= 2
		anchorPrefix++
	}
	anchorPrefix--
	if anchorPrefix < 1 {
		anchorPrefix = 1
	}

	getAnchor := AnchorExtracter(K, maskPrefix, uint8(anchorPrefix))

	fhi, err := os.Create(indexFile)
	if err != nil {
		return err
	}
	wi := bufio.NewWriter(fhi)

	// 8-byte magic number
	err = binary.Write(wi, be, MagicIdx)
	if err != nil {
		return err
	}

	// 8-byte meta info
	err = binary.Write(wi, be, [8]uint8{MainVersion, MinorVersion, K, maskPrefix, uint8(anchorPrefix)})
	if err != nil {
		return err
	}

	// 16-byte the MaskOffset and the chunk size
	err = binary.Write(wi, be, [2]uint64{ChunkIndex, ChunkSize})
	if err != nil {
		return err
	}

	// --------------------------------------------------------------------
	// read kv data

	var nKmers int
	var ctrlByte byte
	var lastPair bool  // check if this is the last pair
	var hasKmer2 bool  // check if there's a kmer2
	var _offset uint64 // offset of kmer
	var nBytes int
	var nReaded, nDecoded int
	var v1, v2 uint64
	var kmer1, kmer2 uint64
	var lenVal1, lenVal2 uint64
	var i, j uint64

	var _j int
	var prefix, prefixPre uint64
	var first bool

	p2o := make([]uint64, 2+int(math.Pow(4, float64(anchorPrefix)))<<1)
	var offset1 uint64

	for i = 0; i < ChunkSize; i++ { // for chunkSize masks
		// fmt.Printf("chunk: %d/%d\n", i, ChunkSize)

		_offset = 0

		// 8-byte the number of k-mers
		nReaded, err = io.ReadFull(r, buf8)
		if err != nil {
			return err
		}
		if nReaded < 8 {
			return ErrBrokenFile
		}
		nKmers = int(be.Uint64(buf8))
		offset += 8

		if nKmers == 0 { // this hapens when no captured k-mer for a mask
			// 8-byte the number of anchors
			err = binary.Write(wi, be, uint64(0))
			if err != nil {
				return err
			}

			continue
		}

		first = true
		clear(p2o)
		p2o[1] = uint64(offset) << 1 // offset of the first k-mer

		// fmt.Printf("nKmers: %d, nAnchors: %d, offset: %d\n", nKmers, _nAnchors, offset)
		for {
			// ------------------------------------------------------------------------

			// read the control byte
			_, err = io.ReadFull(r, buf[:1])
			if err != nil {
				return err
			}
			ctrlByte = buf[0]

			lastPair = ctrlByte&128 > 0 // 1<<7
			hasKmer2 = ctrlByte&64 == 0 // 1<<6

			ctrlByte &= 63

			// parse the control byte
			nBytes = util.CtrlByte2ByteLengthsUint64(ctrlByte)

			// read encoded bytes
			nReaded, err = io.ReadFull(r, buf[:nBytes])
			if err != nil {
				return err
			}
			if nReaded < nBytes {
				return ErrBrokenFile
			}

			v1, v2, nDecoded = util.Uint64s(ctrlByte, buf[:nBytes])
			if nDecoded == 0 {
				return ErrBrokenFile
			}

			kmer1 = v1 + _offset
			kmer2 = kmer1 + v2
			_offset = kmer2

			// ------------------------------------------------------------------------
			// index anchor

			// key 1
			prefix = getAnchor(kmer1)
			if first || prefix != prefixPre { // the first new prefix
				first = false

				_j = int(prefix<<1) + 2
				p2o[_j], p2o[_j+1] = kmer1, uint64(offset)<<1

				prefixPre = prefix
			}

			offset1 = uint64(offset) << 1 // for kmer2, just for keeping compatibility

			// ------------------ lengths of values -------------------

			offset += 1 + nBytes

			// read the control byte
			_, err = io.ReadFull(r, buf[:1])
			if err != nil {
				return err
			}
			ctrlByte = buf[0]

			// parse the control byte
			nBytes = util.CtrlByte2ByteLengthsUint64(ctrlByte)

			// read encoded bytes
			nReaded, err = io.ReadFull(r, buf[:nBytes])
			if err != nil {
				return err
			}
			if nReaded < nBytes {
				return ErrBrokenFile
			}

			offset += 1 + nBytes

			lenVal1, lenVal2, nDecoded = util.Uint64s(ctrlByte, buf[:nBytes])
			if nDecoded == 0 {
				return ErrBrokenFile
			}

			// ------------------ values -------------------

			for j = 0; j < lenVal1; j++ {
				nReaded, err = io.ReadFull(r, buf8)
				if err != nil {
					return err
				}
				if nReaded < 8 {
					return ErrBrokenFile
				}

			}

			offset += int(lenVal1) << 3

			if lastPair && !hasKmer2 {
				break
			}

			// key 2
			prefix = getAnchor(kmer2)
			if prefix != prefixPre { // the first new prefix
				_j = int(prefix<<1) + 2
				p2o[_j], p2o[_j+1] = kmer2, offset1|1

				prefixPre = prefix
			}

			for j = 0; j < lenVal2; j++ {
				nReaded, err = io.ReadFull(r, buf8)
				if err != nil {
					return err
				}
				if nReaded < 8 {
					return ErrBrokenFile
				}

			}

			offset += int(lenVal2) << 3

			if lastPair {
				break
			}
		}

		// -----------------------------------------
		// save index

		var kmer uint64
		var nAnchors uint64
		var offset2 uint64
		e := len(p2o) >> 1
		for i := 0; i < e; i++ {
			offset2 = p2o[i<<1+1]
			if offset2 > 0 {
				nAnchors++
			}
		}
		// 8-byte the number of anchors
		err = binary.Write(wi, be, nAnchors) // including the extra 1
		if err != nil {
			return err
		}

		p2o[0] = nAnchors // might be useful

		// k-mer and offset
		for i := 0; i < e; i++ {
			_j = i << 1
			kmer, offset2 = p2o[_j], p2o[_j+1]
			if offset2 > 0 {
				be.PutUint64(buf[:8], kmer)      // k-mer
				be.PutUint64(buf[8:16], offset2) // offset
				_, err = wi.Write(buf[:16])
				if err != nil {
					return err
				}
			}
		}
	}

	// --------------------------------------------------------------------
	// close reader and writer

	err = wi.Flush()
	if err != nil {
		return err
	}
	err = fhi.Close()
	if err != nil {
		return err
	}

	return nil
}
