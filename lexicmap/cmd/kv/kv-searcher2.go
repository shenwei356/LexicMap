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
type InMemorySearcher struct {
	K          uint8 // kmer size
	ChunkIndex int   // index of the first mask in this chunk
	ChunkSize  int   // the number of masks in this chunk

	rdr *Reader // reader of the kv-data file

	// kv data of the ChunkSize masks.
	// A list of k-mer and value pairs are intermittently saved in a []uint64
	KVdata    [][]uint64
	Indexes   [][]int
	getAnchor func(uint64) uint64

	maxKmer uint64
}

// NewSearcher creates a new Searcher for the given kv-data file.
func NewInMemomrySearcher(file string) (*InMemorySearcher, error) {
	rdr, err := NewReader(file)
	if err != nil {
		return nil, errors.Wrapf(err, "reading kv-data file")
	}

	kvdata := make([][]uint64, rdr.ChunkSize)
	indexes := make([][]int, rdr.ChunkSize)
	var getAnchor func(uint64) uint64
	once := true
	for i := 0; i < rdr.ChunkSize; i++ {
		m, index, maskPrefix, anchorPrefix, err := rdr.ReadDataOfAMaskAsListAndCreateIndex()
		if err != nil {
			return nil, errors.Wrapf(err, "reading kv-data")
		}

		kvdata[i] = m
		indexes[i] = index

		if once {
			once = false
			getAnchor = AnchorExtracter(rdr.K, maskPrefix, anchorPrefix)
		}
	}

	scr := &InMemorySearcher{
		K:          rdr.K,
		ChunkIndex: rdr.ChunkIndex,
		ChunkSize:  rdr.ChunkSize,
		rdr:        rdr,
		KVdata:     kvdata,
		Indexes:    indexes,
		getAnchor:  getAnchor,

		maxKmer: 1<<(rdr.K<<1) - 1,
	}
	return scr, nil
}

// Search queries a k-mer and returns k-mers with a minimum prefix of p,
// and maximum m mismatches.
// For m <0 or m >= k-p, mismatch will not be checked.
//
// Please remember to recycle the results object with RecycleSearchResults().
func (scr *InMemorySearcher) Search(kmers []uint64, p uint8, checkFlag bool, reversedKmer bool) (*[]*SearchResult, error) {
	// func (scr *InMemorySearcher) Search(kmers []uint64, p uint8, m int) (*[]*SearchResult, error) {
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

	var rvflag uint64
	if reversedKmer {
		rvflag = MASK_REVERSE
	}

	// ----------------------------------------------------------

	var suffix2 uint8
	var leftBound, rightBound uint64
	var mask uint64

	var last, begin, middle, end int
	var i, j int

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
	ttt := (uint64(1) << (k << 1)) - 1

	var index []int
	getAnchor := scr.getAnchor
	var anchor, anchorNext uint64
	var lastNext uint64

	for iQ, data := range scr.KVdata {
		if len(data) == 0 { // this hapens when no captured k-mer for a mask
			continue
		}

		last = len(data) - 2
		index = scr.Indexes[iQ]

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

		anchor = getAnchor(leftBound)
		i = index[anchor]
		if i < 0 {
			continue
		}
		// fmt.Printf("\nquery: i: %d, kmer:%s\n", i, lexichash.MustDecode(leftBound, k))
		// fmt.Printf("before: i: %d, kmer:%s\n", i, lexichash.MustDecode(data[i], k))

		lastNext = uint64(len(index))
		anchorNext = anchor + 1
		for j = index[anchorNext]; anchorNext < lastNext && j < 0; {
			anchorNext++
		}
		if j > 0 && // use binary search between i and j
			j-i > 2 { // more than one seeds for the anchor

			// fmt.Printf("%s: %d, %s: %d\n", lexichash.MustDecode(anchor, 5), i, lexichash.MustDecode(anchorNext, 5), j)

			begin, end = i, j
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
				if leftBound <= data[middle] { // when they are equal, we still need to check
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
			// fmt.Printf("after: i: %d, kmer:%s\n", i, lexichash.MustDecode(data[i], k))
		} // else: no next anchor available, reasons: 1) only 1 reference, 2) the anchor is TTTTT

		// fmt.Printf("i: %d, kmer:%s\n", i, lexichash.MustDecode(data[i], k))

		// -----------------------------------------------------
		// check one by one

		found = false
		first = true
		sr1 = nil
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
					sr1.Len = uint8(bits.LeadingZeros64(kmer^kmer1)>>1) + k - 32
					sr1.IsSuffix = reversedKmer
					// sr1.Mismatch = mismatch
					sr1.Values = sr1.Values[:0]

					//	fmt.Printf("  create new result: %p\n", sr1)

					first = false
				}

				if !checkFlag || data[i+1]&MASK_REVERSE == rvflag {
					sr1.Values = append(sr1.Values, data[i+1])
				}

				kmer0 = kmer1
			} else {
				sr1 = nil
			}

			if i == last {
				break
			}

			i += 2
		}
		if sr1 != nil {
			*results = append(*results, sr1)
		}
	}

	return results, nil
}

// Search2 is very similar to Search, only the data structure of input kmers is different.
func (scr *InMemorySearcher) Search2(kmers []*[]uint64, p uint8, checkFlag bool, reversedKmer bool) (*[]*SearchResult, error) {
	// func (scr *InMemorySearcher) Search(kmers []uint64, p uint8, m int) (*[]*SearchResult, error) {
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

	var rvflag uint64
	if reversedKmer {
		rvflag = MASK_REVERSE
	}

	// ----------------------------------------------------------

	var suffix2 uint8
	var leftBound, rightBound uint64
	var mask uint64

	var last, begin, middle, end int
	var i, j int

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
	ttt := (uint64(1) << (k << 1)) - 1

	var iKmer int

	var index []int
	getAnchor := scr.getAnchor
	var anchor, anchorNext uint64
	var lastNext uint64

	for iQ, data := range scr.KVdata {
		if len(data) == 0 { // this hapens when no captured k-mer for a mask
			continue
		}

		last = len(data) - 2
		index = scr.Indexes[iQ]

		// scope to search
		// e.g., For a query ACGAC and p=3,
		// kmers shared >=3 prefix are: ACGAA ... ACGTT.
		// kmer = kmers[iQ]
		for iKmer, kmer = range *kmers[iQ] {

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

			// fmt.Printf("i: %d, k:%d, m:%d\n", iK, k, iQ+chunkIndex+1)
			// fmt.Printf("kmer: %s\n", lexichash.MustDecode(kmer, k))
			// fmt.Printf("left bound: %s\n", lexichash.MustDecode(leftBound, k))
			// fmt.Printf("right bound: %s\n", lexichash.MustDecode(rightBound, k))

			// -----------------------------------------------------
			// find the nearest anchor

			anchor = getAnchor(leftBound)
			i = index[anchor]
			if i < 0 {
				continue
			}
			// fmt.Printf("\nquery: i: %d, kmer:%s\n", i, lexichash.MustDecode(leftBound, k))
			// fmt.Printf("before: i: %d, kmer:%s\n", i, lexichash.MustDecode(data[i], k))

			lastNext = uint64(len(index))
			anchorNext = anchor + 1
			for j = index[anchorNext]; anchorNext < lastNext && j < 0; {
				anchorNext++
			}
			if j > 0 && // use binary search between i and j
				j-i > 2 { // more than one seeds for the anchor

				// fmt.Printf("%s: %d, %s: %d\n", lexichash.MustDecode(anchor, 5), i, lexichash.MustDecode(anchorNext, 5), j)

				begin, end = i, j
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
					if leftBound <= data[middle] { // when they are equal, we still need to check
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
				// fmt.Printf("after: i: %d, kmer:%s\n", i, lexichash.MustDecode(data[i], k))
			} // else: no next anchor available, reasons: 1) only 1 reference, 2) the anchor is TTTTT

			// fmt.Printf("i: %d, kmer:%s\n", i, lexichash.MustDecode(data[i], k))

			// -----------------------------------------------------
			// check one by one

			found = false
			first = true
			sr1 = nil
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
						sr1.Len = uint8(bits.LeadingZeros64(kmer^kmer1)>>1) + k - 32
						sr1.IsSuffix = reversedKmer
						sr1.IQuery2 = iKmer
						// sr1.Mismatch = mismatch
						sr1.Values = sr1.Values[:0]
						// fmt.Printf("  create new result: %p\n", sr1)

						first = false
					}

					if !checkFlag || data[i+1]&MASK_REVERSE == rvflag {
						// fmt.Printf("  save: %s, %d\n", lexichash.MustDecode(kmer1, k), data[i+1])
						sr1.Values = append(sr1.Values, data[i+1])
					}

					kmer0 = kmer1
				} else {
					sr1 = nil
				}

				if i == last {
					break
				}

				i += 2
			}
			if sr1 != nil {
				*results = append(*results, sr1)
			}
		}
	}

	return results, nil
}

// Close closes the searcher.
func (scr *InMemorySearcher) Close() error {
	return scr.rdr.Close()
}
