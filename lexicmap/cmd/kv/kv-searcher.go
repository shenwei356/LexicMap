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
	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
)

// MASK_REVERSE is the mask of reversed flag
const MASK_REVERSE = 1

// Searcher provides searching service of querying k-mer values in a k-mer-value file.
type Searcher struct {
	K          uint8 // kmer size
	ChunkIndex int   // index of the first mask in this chunk
	ChunkSize  int   // the number of masks in this chunk

	fh *os.File // file handler of the kv-data file

	// indexes of the ChunkSize masks.
	// A list of k-mer and offset pairs are intermittently saved in a []uint64
	Indexes [][]uint64

	maxKmer uint64
	buf     []byte
	buf8    []uint8
}

// NewSearcher creates a new Searcher for the given kv-data file.
func NewSearcher(file string) (*Searcher, error) {
	k, chunkIndex, indexes, err := ReadKVIndex(filepath.Clean(file) + KVIndexFileExt)
	if err != nil {
		return nil, errors.Wrapf(err, "reading kv-data index file")
	}

	fh, err := os.Open(file)
	if err != nil {
		return nil, errors.Wrapf(err, "reading kv-data file")
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
	IQuery int // index of the query kmer, i.e., index of mask
	// Kmer      uint64 // matched kmer
	IQuery2  int   // index of the query of the mask, cause a mask would have multiple k-mers when matchinged by suffx
	Len      uint8 // length of common prefix/suffix between the query and this k-mer
	IsSuffix bool  // if matched by suffix
	// Mismatch  uint8    // number of mismatch, it has meanning only when checking mismatch!
	Values []uint64 // value of this key
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
func (scr *Searcher) Search(kmers []uint64, p uint8, checkFlag bool, reversedKmer bool) (*[]*SearchResult, error) {
	// func (scr *Searcher) Search(kmers []uint64, p uint8, m int) (*[]*SearchResult, error) {
	if len(kmers) != len(scr.Indexes) {
		return nil, fmt.Errorf("number of query kmers (%d) != number of masks (%d)", len(kmers), len(scr.Indexes))
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
	var i int
	var offset uint64 // offset in kv-data file

	var first bool    // the first kmer has a different way to comput the value
	var lastPair bool // check if this is the last pair
	var hasKmer2 bool // check if there's a kmer2

	var _offset uint64 // offset of kmer
	var ctrlByte byte
	var nBytes int
	var nReaded, nDecoded int
	var v1, v2 uint64
	var kmer1, kmer2 uint64
	var lenVal1, lenVal2 uint64
	var j uint64
	buf8 := scr.buf8
	buf := scr.buf

	var err error

	results := poolSearchResults.Get().(*[]*SearchResult)
	*results = (*results)[:0]
	var found, saveKmer bool
	// var mismatch uint8
	var sr1, sr2 *SearchResult

	var kmer uint64
	prefixSearch := p < k
	chunkIndex := scr.ChunkIndex
	ttt := (uint64(1) << (k << 1)) - 1

	for iQ, index := range scr.Indexes {
		if len(index) == 0 { // this hapens when no captured k-mer for a mask
			continue
		}

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
		// fmt.Printf("%s\n", lexichash.MustDecode(kmer, k))
		// fmt.Printf("%s\n", lexichash.MustDecode(leftBound, k))
		// fmt.Printf("%s\n", lexichash.MustDecode(rightBound, k))

		// -----------------------------------------------------
		// find the nearest anchor

		if len(index) == 2 {
			i = 0
			offset = index[1]
		} else {
			last = len(index) - 2
			// fmt.Printf("len: %d, last: %d\n", len(index), last)
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
				// 	index[middle], lexichash.MustDecode(index[middle], k))
				if leftBound <= index[middle] { // when they are equal, we still need to check
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
			offset = index[i+1]
		}

		// fmt.Printf("i: %d, kmer:%s, offset: %d\n", i, lexichash.MustDecode(index[i], k), offset)

		// -----------------------------------------------------
		// check one by one

		r := scr.fh

		r.Seek(int64(offset), 0)

		first = true
		found = false
		var v uint64

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
			nReaded, err = io.ReadFull(r, buf[:nBytes])
			if err != nil {
				return nil, err
			}
			if nReaded < nBytes {
				return nil, ErrBrokenFile
			}

			v1, v2, nDecoded = util.Uint64s(ctrlByte, buf[:nBytes])
			if nDecoded == 0 {
				return nil, ErrBrokenFile
			}

			if first {
				kmer1 = index[i] // from the index
				first = false
			} else {
				kmer1 = v1 + _offset
			}
			kmer2 = kmer1 + v2
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
			nBytes = util.CtrlByte2ByteLengthsUint64(ctrlByte)

			// read encoded bytes
			nReaded, err = io.ReadFull(r, buf[:nBytes])
			if err != nil {
				return nil, err
			}
			if nReaded < nBytes {
				return nil, ErrBrokenFile
			}

			lenVal1, lenVal2, nDecoded = util.Uint64s(ctrlByte, buf[:nBytes])
			if nDecoded == 0 {
				return nil, ErrBrokenFile
			}

			// ------------------ values -------------------

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
				sr1 = poolSearchResult.Get().(*SearchResult)
				sr1.IQuery = iQ + chunkIndex // do not forget to add mask offset
				// sr1.Kmer = kmer1
				sr1.Len = uint8(bits.LeadingZeros64(kmer^kmer1)>>1) + k - 32
				sr1.IsSuffix = reversedKmer
				// sr1.Mismatch = mismatch
				sr1.Values = sr1.Values[:0]

				for j = 0; j < lenVal1; j++ {
					nReaded, err = io.ReadFull(r, buf8)
					if err != nil {
						return nil, err
					}
					if nReaded < 8 {
						return nil, ErrBrokenFile
					}

					v = be.Uint64(buf8)
					if !checkFlag || v&MASK_REVERSE == rvflag {
						sr1.Values = append(sr1.Values, v)
					}
				}

				*results = append(*results, sr1)
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
				// if checkMismatch {
				// 	mismatch = util.MustSharingPrefixKmersMismatch(kmer, kmer2, k, p)
				// 	if mismatch <= m8 {
				// 		saveKmer = true
				// 	}
				// } else {
				saveKmer = true
				// }
			}

			if saveKmer {
				sr2 = poolSearchResult.Get().(*SearchResult)
				sr2.IQuery = iQ + chunkIndex // do not forget to add mask offset
				// sr2.Kmer = kmer2
				sr2.Len = uint8(bits.LeadingZeros64(kmer^kmer2)>>1) + k - 32
				sr2.IsSuffix = reversedKmer
				// sr2.Mismatch = mismatch
				sr2.Values = sr2.Values[:0]

				for j = 0; j < lenVal2; j++ {
					nReaded, err = io.ReadFull(r, buf8)
					if err != nil {
						return nil, err
					}
					if nReaded < 8 {
						return nil, ErrBrokenFile
					}

					v = be.Uint64(buf8)
					if !checkFlag || v&MASK_REVERSE == rvflag {
						sr2.Values = append(sr2.Values, v)
					}
				}

				*results = append(*results, sr2)
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

// Search2 is very similar to Search, only the data structure of input kmers is different.
func (scr *Searcher) Search2(kmers []*[]uint64, p uint8, checkFlag bool, reversedKmer bool) (*[]*SearchResult, error) {
	// func (scr *Searcher) Search(kmers []uint64, p uint8, m int) (*[]*SearchResult, error) {
	if len(kmers) != len(scr.Indexes) {
		return nil, fmt.Errorf("number of query kmers (%d) != number of masks (%d)", len(kmers), len(scr.Indexes))
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
	var i int
	var offset uint64 // offset in kv-data file

	var first bool    // the first kmer has a different way to comput the value
	var lastPair bool // check if this is the last pair
	var hasKmer2 bool // check if there's a kmer2

	var _offset uint64 // offset of kmer
	var ctrlByte byte
	var nBytes int
	var nReaded, nDecoded int
	var v1, v2 uint64
	var kmer1, kmer2 uint64
	var lenVal1, lenVal2 uint64
	var j uint64
	buf8 := scr.buf8
	buf := scr.buf

	var err error

	results := poolSearchResults.Get().(*[]*SearchResult)
	*results = (*results)[:0]
	var found, saveKmer bool
	// var mismatch uint8
	var sr1, sr2 *SearchResult

	var kmer uint64
	prefixSearch := p < k
	chunkIndex := scr.ChunkIndex
	ttt := (uint64(1) << (k << 1)) - 1

	var iKmer int

	for iQ, index := range scr.Indexes {
		if len(index) == 0 { // this hapens when no captured k-mer for a mask
			continue
		}

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

			// fmt.Printf("k:%d, m:%d\n", k, m)
			// fmt.Printf("%s\n", lexichash.MustDecode(kmer, k))
			// fmt.Printf("%s\n", lexichash.MustDecode(leftBound, k))
			// fmt.Printf("%s\n", lexichash.MustDecode(rightBound, k))

			// -----------------------------------------------------
			// find the nearest anchor

			if len(index) == 2 {
				i = 0
				offset = index[1]
			} else {
				last = len(index) - 2
				// fmt.Printf("len: %d, last: %d\n", len(index), last)
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
					// 	index[middle], lexichash.MustDecode(index[middle], k))
					if leftBound <= index[middle] { // when they are equal, we still need to check
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
				offset = index[i+1]
			}

			// fmt.Printf("i: %d, kmer:%s, offset: %d\n", i, lexichash.MustDecode(index[i], k), offset)

			// -----------------------------------------------------
			// check one by one

			r := scr.fh

			r.Seek(int64(offset), 0)

			first = true
			found = false
			var v uint64

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
				nReaded, err = io.ReadFull(r, buf[:nBytes])
				if err != nil {
					return nil, err
				}
				if nReaded < nBytes {
					return nil, ErrBrokenFile
				}

				v1, v2, nDecoded = util.Uint64s(ctrlByte, buf[:nBytes])
				if nDecoded == 0 {
					return nil, ErrBrokenFile
				}

				if first {
					kmer1 = index[i] // from the index
					first = false
				} else {
					kmer1 = v1 + _offset
				}
				kmer2 = kmer1 + v2
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
				nBytes = util.CtrlByte2ByteLengthsUint64(ctrlByte)

				// read encoded bytes
				nReaded, err = io.ReadFull(r, buf[:nBytes])
				if err != nil {
					return nil, err
				}
				if nReaded < nBytes {
					return nil, ErrBrokenFile
				}

				lenVal1, lenVal2, nDecoded = util.Uint64s(ctrlByte, buf[:nBytes])
				if nDecoded == 0 {
					return nil, ErrBrokenFile
				}

				// ------------------ values -------------------

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
					sr1 = poolSearchResult.Get().(*SearchResult)
					sr1.IQuery = iQ + chunkIndex // do not forget to add mask offset
					// sr1.Kmer = kmer1
					sr1.Len = uint8(bits.LeadingZeros64(kmer^kmer1)>>1) + k - 32
					sr1.IsSuffix = reversedKmer
					sr1.IQuery2 = iKmer
					// sr1.Mismatch = mismatch
					sr1.Values = sr1.Values[:0]

					for j = 0; j < lenVal1; j++ {
						nReaded, err = io.ReadFull(r, buf8)
						if err != nil {
							return nil, err
						}
						if nReaded < 8 {
							return nil, ErrBrokenFile
						}

						v = be.Uint64(buf8)
						if !checkFlag || v&MASK_REVERSE == rvflag {
							sr1.Values = append(sr1.Values, v)
						}
					}

					*results = append(*results, sr1)
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
					// if checkMismatch {
					// 	mismatch = util.MustSharingPrefixKmersMismatch(kmer, kmer2, k, p)
					// 	if mismatch <= m8 {
					// 		saveKmer = true
					// 	}
					// } else {
					saveKmer = true
					// }
				}

				if saveKmer {
					sr2 = poolSearchResult.Get().(*SearchResult)
					sr2.IQuery = iQ + chunkIndex // do not forget to add mask offset
					// sr2.Kmer = kmer2
					sr2.Len = uint8(bits.LeadingZeros64(kmer^kmer2)>>1) + k - 32
					sr2.IsSuffix = reversedKmer
					sr2.IQuery2 = iKmer
					// sr2.Mismatch = mismatch
					sr2.Values = sr2.Values[:0]

					for j = 0; j < lenVal2; j++ {
						nReaded, err = io.ReadFull(r, buf8)
						if err != nil {
							return nil, err
						}
						if nReaded < 8 {
							return nil, ErrBrokenFile
						}

						v = be.Uint64(buf8)
						if !checkFlag || v&MASK_REVERSE == rvflag {
							sr2.Values = append(sr2.Values, v)
						}
					}

					*results = append(*results, sr2)
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
	}

	return results, nil
}

// Close closes the searcher.
func (scr *Searcher) Close() error {
	return scr.fh.Close()
}

// func kmerValueString(v uint64) string {
// 	return fmt.Sprintf("batchIdx: %d, genomeIdx: %d, pos: %d, rc: %v",
// 		int(v>>47), int(v<<17>>47), int(v<<34>>35), v&1 > 0)
// }
