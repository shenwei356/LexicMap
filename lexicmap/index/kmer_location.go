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
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OFTestSerializationTestSerialization ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

package index

import (
	"bytes"
	"fmt"
	"io"
	"path/filepath"
	"sync"

	"github.com/shenwei356/LexicMap/lexicmap/util"
	"github.com/shenwei356/xopen"
	"github.com/twotwotwo/sorts/sortutil"
)

// ExtractKmerLocations interate all mask trees to summarize k-mer locations of all ref seq.
func (idx *Index) ExtractKmerLocations() {
	var refpos uint64
	var idIdx int
	// var pos int
	// var rc uint8
	var m uint64 = (1 << 38) - 1
	locs := make([][]uint64, len(idx.IDs))
	for i := range locs {
		locs[i] = make([]uint64, 0, 1024)
	}

	// visit all trees
	for _, t := range idx.Trees {
		t.Walk(func(key uint64, v []uint64) bool {
			for _, refpos = range v {
				idIdx = int(refpos >> 38)
				// pos = int(refpos << 26 >> 28)
				// rc = uint8(refpos & 1)
				// store the pos+rc
				locs[idIdx] = append(locs[idIdx], refpos&m)
			}
			return false
		})
	}

	// sort
	var wg sync.WaitGroup
	tokens := make(chan int, Threads)
	for i := range locs {
		wg.Add(1)
		tokens <- 1
		go func(i int) {
			sortutil.Uint64s(locs[i])
			util.UniqUint64s(&locs[i])
			wg.Done()
			<-tokens
		}(i)
	}
	wg.Wait()

	idx.KmerLocations = locs
}

// MaskLocationFile is the name of the mask file.
const MaskLocationFile = "kmer_locations.bin"

// WriteKmerLocations writes kmer locations to the binary file.
func (idx *Index) WriteKmerLocations() error {
	return idx.writeKmerLocations(filepath.Join(idx.path, MaskLocationFile))
}

func (idx *Index) writeKmerLocations(file string) error {
	outfh, err := xopen.Wopen(file)
	if err != nil {
		return err
	}
	defer outfh.Close()

	var buf bytes.Buffer
	var hasPrev bool
	var v, preVal uint64
	var offset uint64
	var ctrlByte byte
	var nBytes, n int
	bufVar := make([]byte, 16) // needs at most 8+8=16
	buf24 := make([]byte, 24)  // needs at most 1+16 = 17
	buf8 := make([]byte, 8)    // for writing uint8

	var _n int

	// the number of refs
	be.PutUint64(buf8, uint64(len(idx.IDs)))
	_n, err = outfh.Write(buf8)
	if err != nil {
		return fmt.Errorf("write kmers locations error: %s", err)
	}
	if _n != 8 {
		return fmt.Errorf("write kmers locations error: unexpected written bytes")
	}

	for _, locs := range idx.KmerLocations {
		buf.Reset()

		// reset variables
		offset = 0
		hasPrev = false

		// the number of records
		be.PutUint64(buf8, uint64(len(locs)))
		buf.Write(buf8)

		for _, v = range locs {
			if !hasPrev { // write it later
				preVal = v
				hasPrev = true
				continue
			}

			ctrlByte, nBytes = util.PutUint64s(bufVar, preVal-offset, v-preVal)
			buf24[0] = ctrlByte
			n = nBytes + 1
			copy(buf24[1:n], bufVar[:nBytes])

			buf.Write(buf24[:n])

			offset = v
			hasPrev = false
		}
		if hasPrev { // the last single one
			ctrlByte, nBytes = util.PutUint64s(bufVar, preVal-offset, 0)
			buf24[0] = ctrlByte
			n = nBytes + 1
			copy(buf24[1:n], bufVar[:nBytes])

			buf.Write(buf24[:n])
		}

		_n, err = outfh.Write(buf.Bytes())
		if err != nil {
			return fmt.Errorf("write kmer locations error: %s", err)
		}
		if _n != buf.Len() {
			return fmt.Errorf("write kmer locations error: unexpected written bytes")
		}
	}

	return nil
}

// ReadKmerLocations reads kmer locations from the default binary file.
func (idx *Index) ReadKmerLocations() error {
	var err error
	idx.KmerLocations, err = ReadKmerLocationsFromFile(filepath.Join(idx.path, MaskLocationFile))
	return err
}

// ReadKmerLocationsFromFile reads kmer locations from the binary file.
func ReadKmerLocationsFromFile(file string) ([][]uint64, error) {
	fh, err := xopen.Ropen(file)
	if err != nil {
		return nil, err
	}
	defer fh.Close()

	var nVals int
	var nReadRounds int
	var hasSingle bool
	var offset uint64
	var ctrlByte byte
	var bytes [2]uint8
	var nBytes int
	var nReaded, nDecoded int
	var decodedVals [2]uint64
	var i int
	var v uint64

	buf := make([]byte, 24) // needs at most 1+16 = 17
	buf8 := make([]byte, 8)

	// the number of refs
	_, err = io.ReadFull(fh, buf8)
	if err != nil {
		return nil, err
	}
	nRefs := int(be.Uint64(buf8))

	kmerLocations := make([][]uint64, nRefs)

	for j := 0; j < nRefs; j++ {
		// reset variables
		offset = 0

		// the number of values
		_, err = io.ReadFull(fh, buf8)
		if err != nil {
			return nil, err
		}
		nVals = int(be.Uint64(buf8))
		nReadRounds = nVals / 2
		hasSingle = nVals%2 > 0

		vals := make([]uint64, 0, nVals)

		for i = 0; i < int(nReadRounds); i++ {
			// read the control byte
			_, err = io.ReadFull(fh, buf[:1])
			if err != nil {
				return nil, err
			}
			ctrlByte = buf[0]

			// parse the control byte
			bytes = util.CtrlByte2ByteLengths[ctrlByte]
			nBytes = int(bytes[0] + bytes[1])

			// read encoded bytes
			nReaded, err = io.ReadFull(fh, buf[:nBytes])
			if err != nil {
				return nil, err
			}
			if nReaded < nBytes {
				return nil, ErrBrokenFile
			}

			// decoding
			decodedVals, nDecoded = util.Uint64s(ctrlByte, buf[:nBytes])
			if nDecoded == 0 {
				return nil, ErrBrokenFile
			}

			v = offset + decodedVals[0]
			vals = append(vals, v)

			v = v + decodedVals[1]
			vals = append(vals, v)

			offset = v
		}

		if hasSingle {
			// read the control byte
			_, err = io.ReadFull(fh, buf[:1])
			if err != nil {
				return nil, err
			}
			ctrlByte = buf[0]

			// parse the control byte
			bytes = util.CtrlByte2ByteLengths[ctrlByte]
			nBytes = int(bytes[0] + bytes[1])

			// read encoded bytes
			nReaded, err = io.ReadFull(fh, buf[:nBytes])
			if err != nil {
				return nil, err
			}
			if nReaded < nBytes {
				return nil, ErrBrokenFile
			}

			// decoding
			decodedVals, nDecoded = util.Uint64s(ctrlByte, buf[:nBytes])
			if nDecoded == 0 {
				return nil, ErrBrokenFile
			}

			v = offset + decodedVals[0]
			vals = append(vals, v)
		}

		kmerLocations[j] = vals
	}

	return kmerLocations, nil
}
