// Copyright © 2023-2025 Wei Shen <shenwei356@gmail.com>
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
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"sync"

	"github.com/pkg/errors"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
)

// Reader provides methods for reading kv data of a mask, used in kv-data merging.
type Reader struct {
	K          uint8 // kmer size
	ChunkIndex int   // index of the first mask in this chunk
	ChunkSize  int   // the number of masks in this chunk

	file string
	fh   *os.File // file handler of the kv-data file
	r    *bufio.Reader

	buf     []byte
	buf8    []uint8
	buf2048 []uint8 // for parsing seed data

	// index of seed data
	readIndexInfo bool
	maskPrefix    uint8
	anchorPrefix  uint8

	Use3BytesForSeedPos bool
}

// NewReader creates a reader.
func NewReader(file string) (*Reader, error) {
	fh, err := os.Open(file)
	if err != nil {
		return nil, errors.Wrapf(err, "reading kv-data file")
	}

	r := bufio.NewReader(fh)

	rdr := &Reader{
		file: file,
		fh:   fh,
		r:    r,
		buf:  make([]byte, 64),
		buf8: make([]uint8, 8),

		buf2048: make([]uint8, 2048),
	}

	// ---------------------------------------------

	buf := rdr.buf8

	var n int

	// check the magic number
	n, err = io.ReadFull(r, buf)
	if err != nil {
		return nil, err
	}
	if n < 8 {
		return nil, ErrBrokenFile
	}
	same := true
	for i := 0; i < 8; i++ {
		if Magic[i] != buf[i] {
			same = false
			break
		}
	}
	if !same {
		return nil, ErrInvalidFileFormat
	}
	// read version information
	n, err = io.ReadFull(r, buf)
	if err != nil {
		return nil, err
	}
	if n < 8 {
		return nil, ErrBrokenFile
	}
	// check compatibility
	if MainVersion != buf[0] {
		return nil, ErrVersionMismatch
	}
	rdr.K = buf[2] // k-mer size

	rdr.Use3BytesForSeedPos = buf[3]&MaskUse3BytesForSeedPos > 0

	// index of the first mask in current chunk.
	_, err = io.ReadFull(r, buf)
	if err != nil {
		return nil, err
	}
	rdr.ChunkIndex = int(be.Uint64(buf))

	// mask chunk size
	_, err = io.ReadFull(r, buf)
	if err != nil {
		return nil, err
	}
	rdr.ChunkSize = int(be.Uint64(buf))

	return rdr, nil
}

var PoolKmerData = &sync.Pool{New: func() interface{} {
	m := make(map[uint64]*[]uint64, 1024) // can't be too big
	return &m
}}

// RecycleKmerData recycles a k-mer data object.
func RecycleKmerData(m *map[uint64]*[]uint64) {
	clear(*m) // clear the map before recyle it
	PoolKmerData.Put(m)
}

// Close closes the reader
func (rdr *Reader) Close() error {
	return rdr.fh.Close()
}

// ReadDataOfAMaskAsMap reads data of a mask.
// Please remember to recycle the result.
func (rdr *Reader) ReadDataOfAMaskAsMap() (*map[uint64]*[]uint64, error) {
	buf := rdr.buf
	buf8 := rdr.buf8
	buf2048 := rdr.buf2048
	r := rdr.r

	var ctrlByte byte
	var lastPair bool  // check if this is the last pair
	var hasKmer2 bool  // check if there's a kmer2
	var _offset uint64 // offset of kmer
	var nBytes int
	var nReaded int
	var v1, v2 uint64
	var kmer1, kmer2 uint64
	var lenVal, lenVal1, lenVal2 uint64
	var j uint64
	var values *[]uint64
	var v uint64
	var ok bool
	var n uint64

	m := PoolKmerData.Get().(*map[uint64]*[]uint64)
	// clear(*m)  // unnecessary as it's cleared
	var err error

	// 8-byte the number of k-mers
	nReaded, _ = io.ReadFull(r, buf8)
	if nReaded < 8 {
		return nil, ErrBrokenFile
	}
	nKmers := int(be.Uint64(buf8))

	if nKmers == 0 {
		return m, nil
	}

	var nSeedPosBytes uint64 = 8
	fUint64 := be.Uint64
	if rdr.Use3BytesForSeedPos {
		nSeedPosBytes = 7
		fUint64 = Uint64ThreeBytes
	}
	batchSizeBytes := uint64(seedPosBatchSize) * nSeedPosBytes

	for {
		// read the control byte
		_, err = io.ReadFull(r, buf[:1])
		if err != nil {
			return nil, err
		}
		ctrlByte = buf[0]

		lastPair = ctrlByte&128 > 0 // 1<<7
		hasKmer2 = ctrlByte&64 == 0 // 1<<6

		ctrlByte &= 63

		// parse the control byte
		nBytes = util.CtrlByte2ByteLengthsUint64(ctrlByte)

		// read encoded bytes
		nReaded, _ = io.ReadFull(r, buf[:nBytes])
		if nReaded < nBytes {
			return nil, ErrBrokenFile
		}

		v1, v2, _ = util.Uint64s(ctrlByte, buf[:nBytes])

		kmer1 = v1 + _offset
		kmer2 = kmer1 + v2
		_offset = kmer2

		// fmt.Printf("%s, %s\n", lexichash.MustDecode(kmer1, rdr.K), lexichash.MustDecode(kmer2, rdr.K))

		// ------------------ lengths of values -------------------

		// read the control byte
		_, err = io.ReadFull(r, buf[:1])
		if err != nil {
			return nil, err
		}
		ctrlByte = buf[0]

		// parse the control byte
		nBytes = util.CtrlByte2ByteLengthsUint64(ctrlByte)

		// read encoded bytes
		nReaded, _ = io.ReadFull(r, buf[:nBytes])
		if nReaded < nBytes {
			return nil, ErrBrokenFile
		}

		lenVal1, lenVal2, _ = util.Uint64s(ctrlByte, buf[:nBytes])

		// ------------------ values -------------------

		if values, ok = (*m)[kmer1]; !ok {
			values = &[]uint64{}
			(*m)[kmer1] = values
		}

		// for j = 0; j < lenVal1; j++ {
		// 	nReaded, err = io.ReadFull(r, buf8)
		// 	if err != nil {
		// 		return nil, err
		// 	}
		// 	if nReaded < 8 {
		// 		return nil, ErrBrokenFile
		// 	}

		// 	v = be.Uint64(buf8)

		// 	*values = append(*values, v)
		// }

		lenVal = lenVal1
		for lenVal >= seedPosBatchSize {
			nReaded, _ = io.ReadFull(r, buf2048[:batchSizeBytes])
			if nReaded < int(batchSizeBytes) {
				return nil, ErrBrokenFile
			}

			for j = 0; j < batchSizeBytes; j += nSeedPosBytes {
				v = fUint64(buf2048[j : j+nSeedPosBytes])

				*values = append(*values, v)
			}

			lenVal -= seedPosBatchSize
		}
		if lenVal > 0 {
			n = lenVal * uint64(nSeedPosBytes)

			nReaded, _ = io.ReadFull(r, buf2048[:n])
			if nReaded < int(n) {
				return nil, ErrBrokenFile
			}

			for j = 0; j < n; j += nSeedPosBytes {
				v = fUint64(buf2048[j : j+nSeedPosBytes])

				*values = append(*values, v)
			}
		}

		if lastPair && !hasKmer2 {
			break
		}

		if values, ok = (*m)[kmer2]; !ok {
			values = &[]uint64{}
			(*m)[kmer2] = values
		}

		// for j = 0; j < lenVal2; j++ {
		// 	nReaded, err = io.ReadFull(r, buf8)
		// 	if err != nil {
		// 		return nil, err
		// 	}
		// 	if nReaded < 8 {
		// 		return nil, ErrBrokenFile
		// 	}

		// 	v = be.Uint64(buf8)

		// 	*values = append(*values, v)
		// }

		lenVal = lenVal2
		for lenVal >= seedPosBatchSize {
			nReaded, _ = io.ReadFull(r, buf2048[:batchSizeBytes])
			if nReaded < int(batchSizeBytes) {
				return nil, ErrBrokenFile
			}

			for j = 0; j < batchSizeBytes; j += nSeedPosBytes {
				v = fUint64(buf2048[j : j+nSeedPosBytes])

				*values = append(*values, v)
			}

			lenVal -= seedPosBatchSize
		}
		if lenVal > 0 {
			n = lenVal * uint64(nSeedPosBytes)

			nReaded, _ = io.ReadFull(r, buf2048[:n])
			if nReaded < int(n) {
				return nil, ErrBrokenFile
			}

			for j = 0; j < n; j += nSeedPosBytes {
				v = fUint64(buf2048[j : j+nSeedPosBytes])

				*values = append(*values, v)
			}
		}

		if lastPair {
			break
		}
	}

	if len(*m) != nKmers {
		return m, fmt.Errorf("number of k-mers mismatch. expected: %d, got: %d", nKmers, len(*m))
	}

	return m, nil
}

// ReadDataOfAMaskAsMap reads data of a mask.
// Please remember to recycle the result.
func (rdr *Reader) ReadDataOfAMaskAndAppendToMap(m *map[uint64]*[]uint64) error {
	buf := rdr.buf
	buf8 := rdr.buf8
	buf2048 := rdr.buf2048
	r := rdr.r

	var ctrlByte byte
	var lastPair bool  // check if this is the last pair
	var hasKmer2 bool  // check if there's a kmer2
	var _offset uint64 // offset of kmer
	var nBytes int
	var nReaded int
	var v1, v2 uint64
	var kmer1, kmer2 uint64
	var lenVal, lenVal1, lenVal2 uint64
	var j uint64
	var values *[]uint64
	var v uint64
	var ok bool
	var n uint64

	var err error

	// 8-byte the number of k-mers
	nReaded, _ = io.ReadFull(r, buf8)
	if nReaded < 8 {
		return ErrBrokenFile
	}
	nKmers := int(be.Uint64(buf8))

	if nKmers == 0 {
		return nil
	}
	var N int

	var nSeedPosBytes uint64 = 8
	fUint64 := be.Uint64
	if rdr.Use3BytesForSeedPos {
		nSeedPosBytes = 7
		fUint64 = Uint64ThreeBytes
	}
	batchSizeBytes := uint64(seedPosBatchSize) * nSeedPosBytes

	for {
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
		nReaded, _ = io.ReadFull(r, buf[:nBytes])
		if nReaded < nBytes {
			return ErrBrokenFile
		}

		v1, v2, _ = util.Uint64s(ctrlByte, buf[:nBytes])

		kmer1 = v1 + _offset
		kmer2 = kmer1 + v2
		_offset = kmer2

		// fmt.Printf("%s, %s\n", lexichash.MustDecode(kmer1, rdr.K), lexichash.MustDecode(kmer2, rdr.K))

		// ------------------ lengths of values -------------------

		// read the control byte
		_, err = io.ReadFull(r, buf[:1])
		if err != nil {
			return err
		}
		ctrlByte = buf[0]

		// parse the control byte
		nBytes = util.CtrlByte2ByteLengthsUint64(ctrlByte)

		// read encoded bytes
		nReaded, _ = io.ReadFull(r, buf[:nBytes])
		if nReaded < nBytes {
			return ErrBrokenFile
		}

		lenVal1, lenVal2, _ = util.Uint64s(ctrlByte, buf[:nBytes])

		// ------------------ values -------------------

		N++
		if values, ok = (*m)[kmer1]; !ok {
			values = &[]uint64{}
			(*m)[kmer1] = values
		}

		lenVal = lenVal1
		for lenVal >= seedPosBatchSize {
			nReaded, _ = io.ReadFull(r, buf2048[:batchSizeBytes])
			if nReaded < int(batchSizeBytes) {
				return ErrBrokenFile
			}

			for j = 0; j < batchSizeBytes; j += nSeedPosBytes {
				v = fUint64(buf2048[j : j+nSeedPosBytes])

				*values = append(*values, v)
			}

			lenVal -= seedPosBatchSize
		}
		if lenVal > 0 {
			n = lenVal * uint64(nSeedPosBytes)

			nReaded, _ = io.ReadFull(r, buf2048[:n])
			if nReaded < int(n) {
				return ErrBrokenFile
			}

			for j = 0; j < n; j += nSeedPosBytes {
				v = fUint64(buf2048[j : j+nSeedPosBytes])

				*values = append(*values, v)
			}
		}

		if lastPair && !hasKmer2 {
			break
		}

		N++
		if values, ok = (*m)[kmer2]; !ok {
			values = &[]uint64{}
			(*m)[kmer2] = values
		}

		lenVal = lenVal2
		for lenVal >= seedPosBatchSize {
			nReaded, _ = io.ReadFull(r, buf2048[:batchSizeBytes])
			if nReaded < int(batchSizeBytes) {
				return ErrBrokenFile
			}

			for j = 0; j < batchSizeBytes; j += nSeedPosBytes {
				v = fUint64(buf2048[j : j+nSeedPosBytes])

				*values = append(*values, v)
			}

			lenVal -= seedPosBatchSize
		}
		if lenVal > 0 {
			n = lenVal * uint64(nSeedPosBytes)

			nReaded, _ = io.ReadFull(r, buf2048[:n])
			if nReaded < int(n) {
				return ErrBrokenFile
			}

			for j = 0; j < n; j += nSeedPosBytes {
				v = fUint64(buf2048[j : j+nSeedPosBytes])

				*values = append(*values, v)
			}
		}

		if lastPair {
			break
		}
	}

	if N != nKmers {
		return fmt.Errorf("number of k-mers mismatch. expected: %d, got: %d", nKmers, N)
	}

	return nil
}

// ReadDataOfAMaskAsList reads data of a mask
// Returned: a list of k-mer and value pairs are intermittently saved in a []uint64.
func (rdr *Reader) ReadDataOfAMaskAsList() ([]uint64, error) {
	buf := rdr.buf
	buf8 := rdr.buf8
	buf2048 := rdr.buf2048
	r := rdr.r

	var ctrlByte byte
	var lastPair bool  // check if this is the last pair
	var hasKmer2 bool  // check if there's a kmer2
	var _offset uint64 // offset of kmer
	var nBytes int
	var nReaded int
	var v1, v2 uint64
	var kmer, kmer1, kmer2 uint64
	var lenVal, lenVal1, lenVal2 uint64
	var j uint64
	var v uint64
	var n uint64

	var err error

	// 8-byte the number of k-mers
	nReaded, _ = io.ReadFull(r, buf8)
	if nReaded < 8 {
		return nil, ErrBrokenFile
	}
	nKmers := int(be.Uint64(buf8))

	if nKmers == 0 { // this hapens when no captured k-mer for a mask
		return make([]uint64, 0), nil
	}

	// A list of k-mer and value pairs are intermittently saved in a []uint64
	m := make([]uint64, 0, nKmers<<1)
	// multiping 2.2 is because that some k-mers would have more than one locations,
	// it help to reduce slice growing, but it's slightly slower in batch querying, interesting.
	// m := make([]uint64, 0, int(float64(nKmers)*2.2))

	var nSeedPosBytes uint64 = 8
	fUint64 := be.Uint64
	if rdr.Use3BytesForSeedPos {
		nSeedPosBytes = 7
		fUint64 = Uint64ThreeBytes
	}
	batchSizeBytes := uint64(seedPosBatchSize) * nSeedPosBytes

	for {
		// read the control byte
		_, err = io.ReadFull(r, buf[:1])
		if err != nil {
			return nil, err
		}
		ctrlByte = buf[0]

		lastPair = ctrlByte&128 > 0 // 1<<7
		hasKmer2 = ctrlByte&64 == 0 // 1<<6

		ctrlByte &= 63

		// parse the control byte
		nBytes = util.CtrlByte2ByteLengthsUint64(ctrlByte)

		// read encoded bytes
		nReaded, _ = io.ReadFull(r, buf[:nBytes])
		if nReaded < nBytes {
			return nil, ErrBrokenFile
		}

		v1, v2, _ = util.Uint64s(ctrlByte, buf[:nBytes])

		kmer1 = v1 + _offset
		kmer2 = kmer1 + v2
		_offset = kmer2

		// fmt.Printf("%s, %s\n", lexichash.MustDecode(kmer1, rdr.K), lexichash.MustDecode(kmer2, rdr.K))

		// ------------------ lengths of values -------------------

		// read the control byte
		_, err = io.ReadFull(r, buf[:1])
		if err != nil {
			return nil, err
		}
		ctrlByte = buf[0]

		// parse the control byte
		nBytes = util.CtrlByte2ByteLengthsUint64(ctrlByte)

		// read encoded bytes
		nReaded, _ = io.ReadFull(r, buf[:nBytes])
		if nReaded < nBytes {
			return nil, ErrBrokenFile
		}

		lenVal1, lenVal2, _ = util.Uint64s(ctrlByte, buf[:nBytes])

		// ------------------ values -------------------

		// for j = 0; j < lenVal1; j++ {
		// 	nReaded, err = io.ReadFull(r, buf8)
		// 	if err != nil {
		// 		return nil, err
		// 	}
		// 	if nReaded < 8 {
		// 		return nil, ErrBrokenFile
		// 	}

		// 	v = be.Uint64(buf8)

		// 	m = append(m, kmer1)
		// 	m = append(m, v)
		// }

		lenVal = lenVal1
		kmer = kmer1
		for lenVal >= seedPosBatchSize {
			nReaded, _ = io.ReadFull(r, buf2048[:batchSizeBytes])
			if nReaded < int(batchSizeBytes) {
				return nil, ErrBrokenFile
			}

			for j = 0; j < batchSizeBytes; j += nSeedPosBytes {
				v = fUint64(buf2048[j : j+nSeedPosBytes])

				m = append(m, kmer)
				m = append(m, v)
			}

			lenVal -= seedPosBatchSize
		}
		if lenVal > 0 {
			n = lenVal * uint64(nSeedPosBytes)

			nReaded, _ = io.ReadFull(r, buf2048[:n])
			if nReaded < int(n) {
				return nil, ErrBrokenFile
			}

			for j = 0; j < n; j += nSeedPosBytes {
				v = fUint64(buf2048[j : j+nSeedPosBytes])

				m = append(m, kmer)
				m = append(m, v)
			}
		}

		if lastPair && !hasKmer2 {
			break
		}

		// for j = 0; j < lenVal2; j++ {
		// 	nReaded, err = io.ReadFull(r, buf8)
		// 	if err != nil {
		// 		return nil, err
		// 	}
		// 	if nReaded < 8 {
		// 		return nil, ErrBrokenFile
		// 	}

		// 	v = be.Uint64(buf8)

		// 	m = append(m, kmer2)
		// 	m = append(m, v)
		// }

		lenVal = lenVal2
		kmer = kmer2
		for lenVal >= seedPosBatchSize {
			nReaded, _ = io.ReadFull(r, buf2048[:batchSizeBytes])
			if nReaded < int(batchSizeBytes) {
				return nil, ErrBrokenFile
			}

			for j = 0; j < batchSizeBytes; j += nSeedPosBytes {
				v = fUint64(buf2048[j : j+nSeedPosBytes])

				m = append(m, kmer)
				m = append(m, v)
			}

			lenVal -= seedPosBatchSize
		}
		if lenVal > 0 {
			n = lenVal * uint64(nSeedPosBytes)

			nReaded, _ = io.ReadFull(r, buf2048[:n])
			if nReaded < int(n) {
				return nil, ErrBrokenFile
			}

			for j = 0; j < n; j += nSeedPosBytes {
				v = fUint64(buf2048[j : j+nSeedPosBytes])

				m = append(m, kmer)
				m = append(m, v)
			}
		}

		if lastPair {
			break
		}
	}

	if len(m)>>1 < nKmers {
		return m, fmt.Errorf("number of k-mers mismatch. expected: >=%d, got: %d", nKmers, len(m)>>1)
	}

	return m, nil
}

// ReadDataOfAMaskAsListAndCreateIndex reads data of a mask,
// and create a new index with n anchors.
// Returned: a list of k-mer and value pairs are intermittently saved in a []uint64.
func (rdr *Reader) ReadDataOfAMaskAsListAndCreateIndex() ([]uint64, []int, uint8, uint8, error) {
	if !rdr.readIndexInfo {
		var err error
		_, _, _, rdr.maskPrefix, rdr.anchorPrefix, err = ReadKVIndexInfo(filepath.Clean(rdr.file) + KVIndexFileExt)
		if err != nil {
			return nil, nil, 0, 0, errors.Wrapf(err, "reading kv-data index file")
		}
		rdr.readIndexInfo = true
	}

	buf := rdr.buf
	buf8 := rdr.buf8
	buf2048 := rdr.buf2048
	r := rdr.r

	var ctrlByte byte
	var lastPair bool  // check if this is the last pair
	var hasKmer2 bool  // check if there's a kmer2
	var _offset uint64 // offset of kmer
	var nBytes int
	var nReaded int
	var v1, v2 uint64
	var kmer, kmer1, kmer2 uint64
	var lenVal, lenVal1, lenVal2 uint64
	var j uint64
	var v uint64
	var n uint64

	var err error

	// 8-byte the number of k-mers
	nReaded, _ = io.ReadFull(r, buf8)
	if nReaded < 8 {
		return nil, nil, 0, 0, ErrBrokenFile
	}
	nKmers := int(be.Uint64(buf8))

	if nKmers == 0 { // this hapens when no captured k-mer for a mask
		return make([]uint64, 0), nil, 0, 0, nil
	}

	// A list of k-mer and value pairs are intermittently saved in a []uint64
	m := make([]uint64, 0, nKmers<<1)
	// multiping 2.2 is because that some k-mers would have more than one locations,
	// it help to reduce slice growing, but it's slightly slower in batch querying, interesting.
	// m := make([]uint64, 0, int(float64(nKmers)*2.2))

	var iOffset uint64 // offset of kmer

	index := make([]int, int(math.Pow(4, float64(rdr.anchorPrefix))))
	for i := range index {
		index[i] = -1
	}
	getAnchor := AnchorExtracter(rdr.K, rdr.maskPrefix, rdr.anchorPrefix)
	var prefix, prefixPre uint64
	first := true

	var nSeedPosBytes uint64 = 8
	fUint64 := be.Uint64
	if rdr.Use3BytesForSeedPos {
		nSeedPosBytes = 7
		fUint64 = Uint64ThreeBytes
	}
	batchSizeBytes := uint64(seedPosBatchSize) * nSeedPosBytes

	for {
		// read the control byte
		_, err = io.ReadFull(r, buf[:1])
		if err != nil {
			return nil, nil, 0, 0, err
		}
		ctrlByte = buf[0]

		lastPair = ctrlByte&128 > 0 // 1<<7
		hasKmer2 = ctrlByte&64 == 0 // 1<<6

		ctrlByte &= 63

		// parse the control byte
		nBytes = util.CtrlByte2ByteLengthsUint64(ctrlByte)

		// read encoded bytes
		nReaded, _ = io.ReadFull(r, buf[:nBytes])
		if nReaded < nBytes {
			return nil, nil, 0, 0, ErrBrokenFile
		}

		v1, v2, _ = util.Uint64s(ctrlByte, buf[:nBytes])

		kmer1 = v1 + _offset
		kmer2 = kmer1 + v2
		_offset = kmer2

		// ------------------------------------------------------------------------
		// index anchor

		// key 1
		prefix = getAnchor(kmer1)
		if first || prefix != prefixPre { // the first new prefix
			first = false

			index[prefix] = int(iOffset)

			prefixPre = prefix
		}

		// fmt.Printf("%s, %s\n", lexichash.MustDecode(kmer1, rdr.K), lexichash.MustDecode(kmer2, rdr.K))

		// ------------------ lengths of values -------------------

		// read the control byte
		_, err = io.ReadFull(r, buf[:1])
		if err != nil {
			return nil, nil, 0, 0, err
		}
		ctrlByte = buf[0]

		// parse the control byte
		nBytes = util.CtrlByte2ByteLengthsUint64(ctrlByte)

		// read encoded bytes
		nReaded, _ = io.ReadFull(r, buf[:nBytes])
		if nReaded < nBytes {
			return nil, nil, 0, 0, ErrBrokenFile
		}

		lenVal1, lenVal2, _ = util.Uint64s(ctrlByte, buf[:nBytes])

		// ------------------ values -------------------

		// for j = 0; j < lenVal1; j++ {
		// 	nReaded, err = io.ReadFull(r, buf8)
		// 	if err != nil {
		// 		return nil, nil, 0, 0, err
		// 	}
		// 	if nReaded < 8 {
		// 		return nil, nil, 0, 0, ErrBrokenFile
		// 	}

		// 	v = be.Uint64(buf8)

		// 	m = append(m, kmer1)
		// 	m = append(m, v)
		// 	iOffset += 2
		// }

		lenVal = lenVal1
		kmer = kmer1
		for lenVal >= seedPosBatchSize {
			nReaded, _ = io.ReadFull(r, buf2048[:batchSizeBytes])
			if nReaded < int(batchSizeBytes) {
				return nil, nil, 0, 0, ErrBrokenFile
			}

			for j = 0; j < batchSizeBytes; j += nSeedPosBytes {
				v = fUint64(buf2048[j : j+nSeedPosBytes])

				m = append(m, kmer)
				m = append(m, v)
			}

			iOffset += seedPosBatchSize << 1
			lenVal -= seedPosBatchSize
		}
		if lenVal > 0 {
			n = lenVal * uint64(nSeedPosBytes)

			nReaded, _ = io.ReadFull(r, buf2048[:n])
			if nReaded < int(n) {
				return nil, nil, 0, 0, ErrBrokenFile
			}

			for j = 0; j < n; j += nSeedPosBytes {
				v = fUint64(buf2048[j : j+nSeedPosBytes])

				m = append(m, kmer)
				m = append(m, v)
			}

			iOffset += lenVal << 1
		}

		if lastPair && !hasKmer2 {
			break
		}

		// key 2
		prefix = getAnchor(kmer2)
		if prefix != prefixPre { // the first new prefix
			index[prefix] = int(iOffset)

			prefixPre = prefix
		}

		// for j = 0; j < lenVal2; j++ {
		// 	nReaded, err = io.ReadFull(r, buf8)
		// 	if err != nil {
		// 		return nil, nil, 0, 0, err
		// 	}
		// 	if nReaded < 8 {
		// 		return nil, nil, 0, 0, ErrBrokenFile
		// 	}

		// 	v = be.Uint64(buf8)

		// 	m = append(m, kmer2)
		// 	m = append(m, v)
		// 	iOffset += 2
		// }

		lenVal = lenVal2
		kmer = kmer2
		for lenVal >= seedPosBatchSize {
			nReaded, _ = io.ReadFull(r, buf2048[:batchSizeBytes])
			if nReaded < int(batchSizeBytes) {
				return nil, nil, 0, 0, ErrBrokenFile
			}

			for j = 0; j < batchSizeBytes; j += nSeedPosBytes {
				v = fUint64(buf2048[j : j+nSeedPosBytes])

				m = append(m, kmer)
				m = append(m, v)
			}

			iOffset += seedPosBatchSize << 1
			lenVal -= seedPosBatchSize
		}
		if lenVal > 0 {
			n = lenVal * uint64(nSeedPosBytes)

			nReaded, _ = io.ReadFull(r, buf2048[:n])
			if nReaded < int(n) {
				return nil, nil, 0, 0, ErrBrokenFile
			}

			for j = 0; j < n; j += nSeedPosBytes {
				v = fUint64(buf2048[j : j+nSeedPosBytes])

				m = append(m, kmer)
				m = append(m, v)
			}

			iOffset += lenVal << 1
		}

		if lastPair {
			break
		}
	}

	if len(m)>>1 < nKmers {
		return m, nil, 0, 0, fmt.Errorf("number of k-mers mismatch. expected: >=%d, got: %d", nKmers, len(m)>>1)
	}
	// if int(index[len(index)-1]) >= len(m) {
	// 	fmt.Println(_nAnchors, len(m), len(index), index)
	// }

	return m, index, rdr.maskPrefix, rdr.anchorPrefix, nil
}

// --------------------------------------------------------------------

// IndexReader provides methods for reading kv-data index data.
type IndexReader struct {
	K          uint8 // kmer size
	ChunkIndex int   // index of the first mask in this chunk
	ChunkSize  int   // the number of masks in this chunk
	NAnchors   int

	fh *os.File // file handler of the kv-data file
	r  *bufio.Reader

	buf  []byte
	buf8 []uint8

	Use3BytesForSeedPos bool
}

// NewIndexReader creates a index reader
func NewIndexReader(file string) (*IndexReader, error) {
	fh, err := os.Open(file)
	if err != nil {
		return nil, errors.Wrapf(err, "reading kv-data file")
	}

	r := bufio.NewReader(fh)

	rdr := &IndexReader{
		fh:   fh,
		r:    r,
		buf:  make([]byte, 64),
		buf8: make([]uint8, 8),
	}

	// ---------------------------------------------

	buf := rdr.buf8

	var n int

	// check the magic number
	n, err = io.ReadFull(r, buf)
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
	// read version information
	n, err = io.ReadFull(r, buf)
	if err != nil {
		return nil, err
	}
	if n < 8 {
		return nil, ErrBrokenFile
	}
	// check compatibility
	if MainVersion != buf[0] {
		return nil, ErrVersionMismatch
	}
	rdr.K = buf[2] // k-mer size
	rdr.Use3BytesForSeedPos = buf[5]&MaskUse3BytesForSeedPos > 0

	// index of the first mask in current chunk.
	_, err = io.ReadFull(r, buf)
	if err != nil {
		return nil, err
	}
	rdr.ChunkIndex = int(be.Uint64(buf))

	// mask chunk size
	_, err = io.ReadFull(r, buf)
	if err != nil {
		return nil, err
	}
	rdr.ChunkSize = int(be.Uint64(buf))

	// the number of anchors
	_, err = io.ReadFull(r, buf)
	if err != nil {
		return nil, err
	}
	rdr.NAnchors = int(be.Uint64(buf))

	return rdr, nil
}

// Close closes the reader
func (rdr *IndexReader) Close() error {
	return rdr.fh.Close()
}
