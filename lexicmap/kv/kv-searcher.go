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
	"fmt"
	"io"
	"math"
	"math/bits"
	"os"
	"path/filepath"
	"sync"

	"github.com/pkg/errors"
	"github.com/shenwei356/LexicMap/lexicmap/util"
)

// Searcher provides searching service of querying k-mer values in a k-mer-value file.
type Searcher struct {
	K          uint8 // kmer size
	ChunkIndex int   // index of the first mask in this chunk
	ChunkSize  int   // the number of masks in this chunk

	fh *os.File // file handler of the kv-data file

	Indexes [][]uint64 // indexes of the ChunkSize masks

	maxKmer uint64
	buf     []byte
	buf8    []uint8
}

// NewSearcher creates a new Searcher for the given kv-data file.
func NewSearcher(file string) (*Searcher, error) {
	k, chunkIndex, indexes, err := ReadKVIndex(filepath.Clean(file) + KVIndexFileExt)
	if err != nil {
		return nil, errors.Wrapf(err, "reading kv-data file")
	}

	fh, err := os.Open(file)
	if err != nil {
		return nil, errors.Wrapf(err, "reading kv-data index file")
	}

	scr := &Searcher{
		K:          k,
		ChunkIndex: chunkIndex,
		ChunkSize:  len(indexes),
		Indexes:    indexes,
		fh:         fh,

		maxKmer: 1<<(k<<1) - 1,
		buf:     make([]byte, 64),
		buf8:    make([]uint8, 8),
	}
	return scr, nil
}

// SearchResult represents a search result.
type SearchResult struct {
	Kmer      uint64   // searched kmer
	LenPrefix uint8    // length of common prefix between the query and this k-mer
	Mismatch  uint8    // number of mismatch, it has meanning only when checking mismatch!
	Values    []uint64 // value of this key
}

// Reset just resets the stats of a SearchResult
func (r *SearchResult) Reset() {
	r.Values = r.Values[:0]
}

var poolSearchResults = &sync.Pool{New: func() interface{} {
	tmp := make([]*SearchResult, 0, 128)
	return &tmp
}}

var poolSearchResult = &sync.Pool{New: func() interface{} {
	return &SearchResult{Values: make([]uint64, 0, 1)}
}}

// RecycleSearchResults recycles search results objects.
func RecycleSearchResults(sr *[]*SearchResult) {
	for _, r := range *sr {
		poolSearchResult.Put(r)
	}
	poolSearchResults.Put(sr)
}

// Search queries a k-mer and returns k-mers with a minimum prefix of p,
// and maximum m mismatches.
// For m <0 or m >= k-p, mismatch will not be checked.
//
// Please remember to recycle the results object with RecycleSearchResults().
func (scr *Searcher) Search(kmer uint64, p uint8, m int) (*[]*SearchResult, error) {
	if kmer > scr.maxKmer {
		return nil, fmt.Errorf("invalid kmer for k=%d: %d", scr.K, kmer)
	}
	k := scr.K
	if p < 1 || p > k {
		p = k
	}

	checkMismatch := m >= 0 && m < int(k-p)
	m8 := uint8(m)

	// ----------------------------------------------------------
	// scope to search
	// e.g., For a query ACGAC and m=3,
	// kmers shared >=3 prefix are: ACGAA ... ACGTT.

	var suffix2 uint8
	var leftBound, rightBound uint64
	var mask uint64
	if k > p {
		suffix2 = (k - p) << 1
		mask = (1 << suffix2) - 1                  // 1111
		leftBound = kmer & (math.MaxUint64 - mask) // kmer & 1111110000
		rightBound = kmer>>suffix2<<suffix2 + mask // kmer with last 4bits being 1
	} else {
		leftBound = kmer
		rightBound = kmer
	}
	// fmt.Printf("k:%d, m:%d\n", k, m)
	// fmt.Printf("%s\n", lexichash.MustDecode(kmer, k))
	// fmt.Printf("%s\n", lexichash.MustDecode(leftBound, k))
	// fmt.Printf("%s\n", lexichash.MustDecode(rightBound, k))

	// ----------------------------------------------------------
	var last, begin, middle, end int
	var i int
	var offset uint64 // offset in kv-data file

	var first bool    // the first kmer has a different way to comput the value
	var lastPair bool // check if this is the last pair
	var hasKmer2 bool // check if there's a kmer2

	var _offset uint64 // offset of kmer
	var ctrlByte byte
	var bytes [2]uint8
	var nBytes int
	var nReaded, nDecoded int
	var decodedVals [2]uint64
	var kmer1, kmer2 uint64
	var lenVal1, lenVal2 uint64
	var j uint64
	buf8 := scr.buf8
	buf := scr.buf

	var err error

	results := poolSearchResults.Get().(*[]*SearchResult)
	*results = (*results)[:0]
	var found, saveKmer bool
	var mismatch uint8
	var v1, v2 *SearchResult

	for _, index := range scr.Indexes {
		// -----------------------------------------------------
		// find the nearest anchor

		last = len(index) - 2
		// fmt.Printf("len: %d, last: %d\n", len(index), last)
		begin, end = 0, last
		for {
			middle = begin + (end-begin)>>1
			if middle&1 > 0 {
				middle--
			}
			// fmt.Printf("[%d, %d] %d: %d %s\n", begin, end, middle,
			// 	index[middle], lexichash.MustDecode(index[middle], k))
			if leftBound < index[middle] {
				// fmt.Printf(" left\n")
				end = middle // new end
			} else {
				// fmt.Printf(" right\n")
				begin = middle // new start
			}
			if begin+2 == end { // next to eacher
				i = begin
				break
			}
		}
		offset = index[i+1]

		// fmt.Printf("i: %d, kmer:%s, offset: %d\n", i, lexichash.MustDecode(index[i], k), offset)

		// -----------------------------------------------------
		// check one by one

		r := scr.fh

		r.Seek(int64(offset), 0)

		first = true
		found = false

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
			bytes = util.CtrlByte2ByteLengths[ctrlByte]
			nBytes = int(bytes[0] + bytes[1])

			// read encoded bytes
			nReaded, err = io.ReadFull(r, buf[:nBytes])
			if err != nil {
				return nil, err
			}
			if nReaded < nBytes {
				return nil, ErrBrokenFile
			}

			decodedVals, nDecoded = util.Uint64s(ctrlByte, buf[:nBytes])
			if nDecoded == 0 {
				return nil, ErrBrokenFile
			}

			if first {
				kmer1 = index[i] // from the index
				first = false
			} else {
				kmer1 = decodedVals[0] + _offset
			}
			kmer2 = kmer1 + decodedVals[1]
			_offset = kmer2

			if kmer1 > rightBound { // finished
				// fmt.Printf("  kmer1 out of scope: %s\n", lexichash.MustDecode(kmer1, k))
				break
			}

			if kmer1 >= leftBound || kmer2 >= leftBound {
				// fmt.Printf("  found: %v, %v\n", kmer1 >= leftBound, kmer2 >= leftBound)
				found = true
			}

			// ------------------ lengths of values -------------------

			// read the control byte
			_, err = io.ReadFull(r, buf[:1])
			if err != nil {
				return nil, err
			}
			ctrlByte = buf[0]

			// parse the control byte
			bytes = util.CtrlByte2ByteLengths[ctrlByte]
			nBytes = int(bytes[0] + bytes[1])

			// read encoded bytes
			nReaded, err = io.ReadFull(r, buf[:nBytes])
			if err != nil {
				return nil, err
			}
			if nReaded < nBytes {
				return nil, ErrBrokenFile
			}

			decodedVals, nDecoded = util.Uint64s(ctrlByte, buf[:nBytes])
			if nDecoded == 0 {
				return nil, ErrBrokenFile
			}

			lenVal1 = decodedVals[0]
			lenVal2 = decodedVals[1]

			// ------------------ values -------------------

			saveKmer = false
			if found && kmer1 >= leftBound {
				if checkMismatch {
					mismatch = util.MustSharingPrefixKmersMismatch(kmer, kmer1, k, p)
					if mismatch <= m8 {
						saveKmer = true
					}
				} else {
					saveKmer = true
				}
			}
			if saveKmer {
				v1 = poolSearchResult.Get().(*SearchResult)
				v1.Kmer = kmer1
				v1.LenPrefix = uint8(bits.LeadingZeros64(kmer^kmer1)>>1) + k - 32
				v1.Mismatch = mismatch
				v1.Values = v1.Values[:0]

				for j = 0; j < lenVal1; j++ {
					nReaded, err = io.ReadFull(r, buf8)
					if err != nil {
						return nil, err
					}
					if nReaded < 8 {
						return nil, ErrBrokenFile
					}

					v1.Values = append(v1.Values, be.Uint64(buf8))
				}
				*results = append(*results, v1)
			} else {
				for j = 0; j < lenVal1; j++ {
					nReaded, err = io.ReadFull(r, buf8)
					if err != nil {
						return nil, err
					}
					if nReaded < 8 {
						return nil, ErrBrokenFile
					}
				}
			}

			if kmer2 > rightBound { // only record kmer1
				// fmt.Printf("  kmer2 out of scope: %s\n", lexichash.MustDecode(kmer2, k))
				break
			}

			if lastPair && !hasKmer2 {
				// fmt.Printf("  last pair without kmer2: %s\n",
				// 	lexichash.MustDecode(kmer1, k))
				break
			}

			saveKmer = false
			if found {
				if checkMismatch {
					mismatch = util.MustSharingPrefixKmersMismatch(kmer, kmer2, k, p)
					if mismatch <= m8 {
						saveKmer = true
					}
				} else {
					saveKmer = true
				}
			}

			if saveKmer {
				v2 = poolSearchResult.Get().(*SearchResult)
				v2.Kmer = kmer2
				v2.LenPrefix = uint8(bits.LeadingZeros64(kmer^kmer2)>>1) + k - 32
				v2.Mismatch = mismatch
				v2.Values = v2.Values[:0]

				for j = 0; j < lenVal2; j++ {
					nReaded, err = io.ReadFull(r, buf8)
					if err != nil {
						return nil, err
					}
					if nReaded < 8 {
						return nil, ErrBrokenFile
					}

					v2.Values = append(v2.Values, be.Uint64(buf8))
				}

				*results = append(*results, v2)
			} else {
				for j = 0; j < lenVal2; j++ {
					nReaded, err = io.ReadFull(r, buf8)
					if err != nil {
						return nil, err
					}
					if nReaded < 8 {
						return nil, ErrBrokenFile
					}
				}
			}

			if lastPair {
				// fmt.Printf("  last pair: %s, %s\n",
				// 	lexichash.MustDecode(kmer1, k), lexichash.MustDecode(kmer2, k))
				break
			}

		}
	}

	return results, nil
}
