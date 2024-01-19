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
	"errors"
	"io"
	"math"
	"os"
	"path/filepath"
	"sync"

	"github.com/shenwei356/LexicMap/lexicmap/util"
	"github.com/twotwotwo/sorts/sortutil"
)

var be = binary.BigEndian

// Magic number for checking file format
var Magic = [8]byte{'.', 'k', 'v', '-', 'd', 'a', 't', 'a'}

// Magic number for the index file
var MagicIdx = [8]byte{'.', 'k', 'v', 'i', 'n', 'd', 'e', 'x'}

// MainVersion is use for checking compatibility
var MainVersion uint8 = 0

// MinorVersion is less important
var MinorVersion uint8 = 1

// ErrInvalidFileFormat means invalid file format.
var ErrInvalidFileFormat = errors.New("k-mer-value data: invalid binary format")

// ErrBrokenFile means the file is not complete.
var ErrBrokenFile = errors.New("k-mer-value data: broken file")

// ErrKOverflow means K < 1 or K > 32.
var ErrKOverflow = errors.New("k-mer-value data: k-mer size [1, 32] overflow")

// ErrVersionMismatch means version mismatch between files and program
var ErrVersionMismatch = errors.New("k-mer-value data: version mismatch")

// KVIndexFileExt is the file extension of k-mer data index file.
var KVIndexFileExt = ".idx"

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
// Index file stores a fix number of k-mers (anchors) and their offsets in
// the kv-data file for fast access.
//
// Locations of these anchors, e.g., 5 anchors.
//
//	kkvvkkvvkkvvkkvvkkvvkkvvkkvvkkvvkkvvkkvvkkvv
//	A       A       A       A       A
//
// Header (40 bytes):
//
//	Magic number, 8 bytes, ".kvindex".
//	Main and minor versions, 2 bytes.
//	K size, 1 byte.
//	Blank, 5 bytes.
//	Mask start index, 8 bytes. The index of the first index.
//	Mask chunk size, 8 bytes. The number of masks in this file.
//	Number of anchors, 8 bytes, default: $(squre root of ref genomes).
//
// For each mask:
//
//	kmer-offset data:
//
//		k-mer: 8 bytes
//		offset: 8 bytes
func WriteKVData(k uint8, MaskOffset int, data []map[uint64]*[]uint64, file string, nAnchors int) (int, error) {
	if len(data) == 0 {
		return 0, errors.New("k-mer-value data: no data given")
	}

	var N int // the number of bytes.

	// file handlers
	fh, err := os.Create(file)
	if err != nil {
		return N, err
	}
	w := bufio.NewWriter(fh)
	fhi, err := os.Create(filepath.Clean(file) + KVIndexFileExt)
	if err != nil {
		return N, err
	}
	wi := bufio.NewWriter(fhi)
	defer func() {
		w.Flush()
		fh.Close()
		wi.Flush()
		fhi.Close()
	}()

	// check nAnchors
	nKmers := len(data[0])
	if nAnchors <= 0 {
		nAnchors = int(math.Sqrt(float64(nKmers)))
	} else if nAnchors > nKmers {
		nAnchors = nKmers >> 1
	}
	idxChunkSize := (nKmers / nAnchors) >> 1

	// ---------------------------------------------------------------------------

	// 8-byte magic number
	err = binary.Write(w, be, Magic)
	if err != nil {
		return N, err
	}
	N += 8

	// 8-byte meta info
	err = binary.Write(w, be, [8]uint8{MainVersion, MinorVersion, k})
	if err != nil {
		return N, err
	}
	N += 8

	// 16-byte the MaskOffset and the chunk size
	err = binary.Write(w, be, [2]uint64{uint64(MaskOffset), uint64(len(data))})
	if err != nil {
		return N, err
	}
	N += 16

	// --------------------- Index file

	// 8-byte magic number
	err = binary.Write(wi, be, MagicIdx)
	if err != nil {
		return N, err
	}

	// 8-byte meta info
	err = binary.Write(wi, be, [8]uint8{MainVersion, MinorVersion, k})
	if err != nil {
		return N, err
	}

	// 16-byte the MaskOffset and the chunk size
	err = binary.Write(wi, be, [2]uint64{uint64(MaskOffset), uint64(len(data))})
	if err != nil {
		return N, err
	}

	// 8-byte the number of anchors
	err = binary.Write(wi, be, uint64(nAnchors))
	if err != nil {
		return N, err
	}

	// ---------------------------------------------------------------------------

	var hasPrev bool
	var preKey, key, _v uint64
	var preVal, v *[]uint64
	var offset uint64
	var ctrlByteKey, ctrlByteVal byte
	var nBytesKey, nBytesVal, n int
	bufVar := make([]byte, 16) // needs at most 8+8=16
	buf := make([]byte, 36)    // needs at most 1+16+1+16=34
	buf8 := make([]byte, 8)    // for writing uint8
	bufVals := poolBytesBuffer.Get().(*bytes.Buffer)
	defer poolBytesBuffer.Put(bufVals)
	var even bool
	var i, nm1 int
	var j int
	var recordedAnchors int

	for _, m := range data {
		hasPrev = false
		offset = 0

		// 8-byte the number of k-mers
		err = binary.Write(w, be, uint64(len(m)))
		if err != nil {
			return N, err
		}
		N += 8

		// sort keys
		keys := poolUint64s.Get().(*[]uint64)
		*keys = (*keys)[:0]
		for key = range m {
			*keys = append(*keys, key)
		}
		sortutil.Uint64s(*keys)

		// for decide should we set flag for the last control byte of the last k-mer
		even = len(*keys)&1 == 0 // the number of kmers is even
		nm1 = len(*keys) - 1     // idx of the last element

		j = 0
		recordedAnchors = 0

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
			if j%idxChunkSize == 0 {
				if recordedAnchors < nAnchors {
					recordedAnchors++

					// fmt.Printf("[%d] %d, %d, %d\n", recordedAnchors, i, preKey, N)
					be.PutUint64(buf[:8], preKey)      // k-mer
					be.PutUint64(buf[8:16], uint64(N)) // offset
					_, err = wi.Write(buf[:16])
					if err != nil {
						return N, err
					}
				}

				j = 0
			}
			j++

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
				return N, err
			}
			N += n

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
				return N, err
			}
			N += bufVals.Len()

			// update

			offset = key
			hasPrev = false
		}

		if hasPrev { // the last single one
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
				return N, err
			}
			N += n

			// values

			bufVals.Reset()
			for _, _v = range *preVal {
				be.PutUint64(buf8, _v)
				bufVals.Write(buf8)
			}

			_, err = w.Write(bufVals.Bytes())
			if err != nil {
				return N, err
			}
			N += bufVals.Len()
		}

	}

	return N, nil
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
func ReadKVIndex(file string) (uint8, int, [][]uint64, error) {
	fh, err := os.Open(file)
	if err != nil {
		return 0, -1, nil, err
	}
	r := bufio.NewReader(fh)
	defer fh.Close()

	// ---------------------------------------------

	buf := make([]byte, 8)

	var n int

	// check the magic number
	n, err = io.ReadFull(r, buf)
	if err != nil {
		return 0, -1, nil, err
	}
	if n < 8 {
		return 0, -1, nil, ErrBrokenFile
	}
	same := true
	for i := 0; i < 8; i++ {
		if MagicIdx[i] != buf[i] {
			same = false
			break
		}
	}
	if !same {
		return 0, -1, nil, ErrInvalidFileFormat
	}

	// read version information
	n, err = io.ReadFull(r, buf)
	if err != nil {
		return 0, -1, nil, err
	}
	if n < 8 {
		return 0, -1, nil, ErrBrokenFile
	}
	// check compatibility
	if MainVersion != buf[0] {
		return 0, -1, nil, ErrVersionMismatch
	}
	k := buf[2] // k-mer size

	// index of the first mask in current chunk.
	var iFirstMask int
	_, err = io.ReadFull(r, buf)
	if err != nil {
		return 0, -1, nil, err
	}
	iFirstMask = int(be.Uint64(buf))

	// mask chunk size
	var nMasks int
	_, err = io.ReadFull(r, buf)
	if err != nil {
		return 0, -1, nil, err
	}
	nMasks = int(be.Uint64(buf))

	// the number of anchors
	var nAnchors int
	_, err = io.ReadFull(r, buf)
	if err != nil {
		return 0, -1, nil, err
	}
	nAnchors = int(be.Uint64(buf))

	// ---------------------------------------------

	data := make([][]uint64, nMasks)

	var kmer, offset uint64
	var j int
	for i := 0; i < nMasks; i++ {
		index := make([]uint64, 0, nAnchors<<1)
		for j = 0; j < nAnchors; j++ {
			_, err = io.ReadFull(r, buf)
			if err != nil {
				return 0, -1, nil, err
			}
			kmer = be.Uint64(buf)

			_, err = io.ReadFull(r, buf)
			if err != nil {
				return 0, -1, nil, err
			}
			offset = be.Uint64(buf)

			index = append(index, kmer)
			index = append(index, offset)
		}
		data[i] = index
	}

	return k, iFirstMask, data, nil
}

var poolBytesBuffer = &sync.Pool{New: func() interface{} {
	return &bytes.Buffer{}
}}

var poolUint64s = &sync.Pool{New: func() interface{} {
	tmp := make([]uint64, 0, 1<<20)
	return &tmp
}}
