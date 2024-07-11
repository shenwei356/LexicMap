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
	"math"
	"math/bits"

	"github.com/pkg/errors"
)

// Searcher provides searching service of querying k-mer values in a k-mer-value file.
type InMemorySearcher2 struct {
	K          uint8 // kmer size
	ChunkIndex int   // index of the first mask in this chunk
	ChunkSize  int   // the number of masks in this chunk

	rdr *Reader // reader of the kv-data file

	// kv data of the ChunkSize masks.
	// A list of k-mer and value pairs are intermittently saved in a []uint64
	KVdata [][]uint64 // indexes of the ChunkSize masks

	Indexes [][]uint64

	maxKmer uint64
}

// NewSearcher creates a new Searcher for the given kv-data file.
func NewInMemomrySearcher2(file string) (*InMemorySearcher2, error) {
	rdr, err := NewReader(file)
	if err != nil {
		return nil, errors.Wrapf(err, "reading kv-data file")
	}

	kvdata := make([][]uint64, rdr.ChunkSize)
	indexes := make([][]uint64, rdr.ChunkSize)
	for i := 0; i < rdr.ChunkSize; i++ {
		m, index, err := rdr.ReadDataOfAMaskAsListAndCreateIndex(512)
		if err != nil {
			return nil, errors.Wrapf(err, "reading kv-data")
		}
		kvdata[i] = m
		indexes[i] = index
	}

	scr := &InMemorySearcher2{
		K:          rdr.K,
		ChunkIndex: rdr.ChunkIndex,
		ChunkSize:  rdr.ChunkSize,
		rdr:        rdr,
		KVdata:     kvdata,
		Indexes:    indexes,

		maxKmer: 1<<(rdr.K<<1) - 1,
	}
	return scr, nil
}

// Search queries a k-mer and returns k-mers with a minimum prefix of p,
// and maximum m mismatches.
// For m <0 or m >= k-p, mismatch will not be checked.
//
// Please remember to recycle the results object with RecycleSearchResults().
func (scr *InMemorySearcher2) Search(kmers []uint64, p uint8) (*[]*SearchResult, error) {
	// func (scr *InMemorySearcher2) Search(kmers []uint64, p uint8, m int) (*[]*SearchResult, error) {
	if len(kmers) != scr.ChunkSize {
		return nil, fmt.Errorf("number of query kmers (%d) != number of masks (%d)", len(kmers), len(scr.KVdata))
	}
	// if kmer > scr.maxKmer {
	// 	return nil, fmt.Errorf("invalid kmer for k=%d: %d", scr.K, kmer)
	// }
	k := scr.K
	if p < 1 || p > k {
		p = k
	}

	// checkMismatch := m >= 0 && m < int(k-p)
	// m8 := uint8(m)

	// ----------------------------------------------------------

	var suffix2 uint8
	var leftBound, rightBound uint64
	var mask uint64

	var last, begin, middle, end int
	var i int

	var kmer0 uint64 // previous one
	var kmer1 uint64 // current one

	results := poolSearchResults.Get().(*[]*SearchResult)
	*results = (*results)[:0]

	var found, saveKmer bool
	// var mismatch uint8
	var sr1 *SearchResult

	var kmer uint64
	prefixSearch := p < k
	chunkIndex := scr.ChunkIndex
	var first bool
	var index []uint64
	var lenData int
	ttt := (uint64(1) << (k << 1)) - 1

	for iQ, data := range scr.KVdata {
		// scope to search
		// e.g., For a query ACGAC and p=3,
		// kmers shared >=3 prefix are: ACGAA ... ACGTT.
		kmer = kmers[iQ]

		if kmer == 0 || kmer == ttt { // skip AAAAAAAAAA and TTTTTTTTT
			continue
		}

		if prefixSearch {
			suffix2 = (k - p) << 1
			mask = (1 << suffix2) - 1                  // 1111
			leftBound = kmer & (math.MaxUint64 - mask) // kmer & 1111110000
			rightBound = kmer>>suffix2<<suffix2 + mask // kmer with last 4bits being 1
		} else {
			leftBound = kmer
			rightBound = kmer
		}

		// fmt.Printf("k:%d, m:%d\n", k, m)
		// fmt.Printf("kmer: %s\n", lexichash.MustDecode(kmer, k))
		// fmt.Printf("left bound: %s\n", lexichash.MustDecode(leftBound, k))
		// fmt.Printf("right bound: %s\n", lexichash.MustDecode(rightBound, k))

		// -----------------------------------------------------
		// find the nearest anchor

		index = scr.Indexes[iQ]
		if len(index) == 2 {
			i = 0
			last = 0
			// fmt.Printf("len: %d\n", len(data))
		} else {
			last = len(index) - 2
			// fmt.Printf("len: %d, last: %d\n", len(data), last)
			begin, end = 0, last
			for {
				middle = begin + (end-begin)>>1
				if middle&1 > 0 {
					middle--
				}
				if middle == begin { // when there are only two indexes, middle = 1 and then middle = 0
					i = begin
					break
				}
				// fmt.Printf("[%d, %d] %d: %d %s\n", begin, end, middle,
				// 	data[middle], lexichash.MustDecode(data[middle], k))
				// we still need to check even if they are equal,
				// because the k-mer in the middle might be duplicated
				if leftBound <= index[middle] {
					// fmt.Printf(" left\n")
					end = middle // new end
				} else {
					// fmt.Printf(" right\n")
					begin = middle // new start
				}
				if begin+2 == end { // next to each other
					i = begin
					break
				}
			}
			// i is the target
		}
		i = int(index[i+1]) // the offset in data

		// fmt.Printf("i: %d, kmer:%s\n", i, lexichash.MustDecode(data[i], k))

		// -----------------------------------------------------
		// check one by one

		found = false
		first = true
		sr1 = nil
		lenData = len(data)
		for {
			kmer1 = data[i]

			if kmer1 > rightBound { // finished
				// fmt.Printf("  kmer1 out of scope: %s\n", lexichash.MustDecode(kmer1, k))
				break
			}

			if kmer1 >= leftBound {
				// fmt.Printf("  found: %s\n", lexichash.MustDecode(kmer1, k))
				found = true
			}

			saveKmer = false
			if found && kmer1 >= leftBound {
				// if checkMismatch {
				// 	mismatch = util.MustSharingPrefixKmersMismatch(kmer, kmer1, k, p)
				// 	if mismatch <= m8 {
				// 		saveKmer = true
				// 	}
				// } else {
				saveKmer = true
				// }
			}
			if saveKmer {
				// fmt.Printf("  save: %s\n", lexichash.MustDecode(kmer1, k))
				if kmer1 != kmer0 || first { // new kmer
					if sr1 != nil {
						// fmt.Printf("  record new result: %p\n", sr1)
						*results = append(*results, sr1) // previous one
					}

					sr1 = poolSearchResult.Get().(*SearchResult)
					sr1.IQuery = iQ + chunkIndex // do not forget to add mask offset
					// sr1.Kmer = kmer1
					sr1.LenPrefix = uint8(bits.LeadingZeros64(kmer^kmer1)>>1) + k - 32
					// sr1.Mismatch = mismatch
					sr1.Values = sr1.Values[:0]

					//	fmt.Printf("  create new result: %p\n", sr1)

					first = false
				}

				sr1.Values = append(sr1.Values, data[i+1])

				kmer0 = kmer1
			} else {
				sr1 = nil
			}

			if i == last {
				break
			}

			i += 2

			if i >= lenData {
				break
			}
		}
		if sr1 != nil {
			*results = append(*results, sr1)
		}
	}

	return results, nil
}

// Close closes the searcher.
func (scr *InMemorySearcher2) Close() error {
	return scr.rdr.Close()
}
