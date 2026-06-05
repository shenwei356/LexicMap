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

package cmd

import (
	"slices"
	"sync"

	rtree "github.com/shenwei356/LexicMap/lexicmap/cmd/tree"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
	"github.com/shenwei356/lexichash/iterator"
	"github.com/twotwotwo/sorts"
	"github.com/twotwotwo/sorts/sortutil"
)

// FragmentComparatorOptions defines the options for comparing two sets
// of genome fragments by counting shared k-mers. See FragmentComparator for details.
type FragmentComparatorOptions struct {
	K uint8
	// MinPrefix uint8

	MinSharedKmers uint16

	TopNFragments int // if > 0, for a, only return the top N subject b with the most shared k-mers

	// FracMinHash scaling factor: only k-mers with Hash64(canonical) % Scaled == 0
	// are kept. 0 or 1 disables filtering. Larger values reduce memory and CPU
	// roughly proportionally; MinSharedKmers should be retuned accordingly
	// (e.g. Scaled=10 -> divide MinSharedKmers by ~10).
	Scaled uint32
}

// FragmentComparator compares two sets of genome fragments by counting shared
// k-mers. It is designed for the case where one set of fragments (fragsA) is
// fixed and compared against many sets of fragments (fragsB). The first step
// is to pre-index fragsA into a sorted slice of canonical k-mer entries, which
// can be reused across many comparisons. Each comparison then collects and
// sorts entries for fragsB and merges the two sorted streams to count shared
// k-mers per fragment pair.
type FragmentComparator struct {
	// options
	options *FragmentComparatorOptions

	ccc, ggg, ttt uint64

	// reusable buffer for collecting (kmer, fragID|strand|genome) entries.
	// after collection, sorted by Key and scanned linearly for equal-key runs.
	entries []rtree.BatchEntry

	poolChainers *sync.Pool
}

var poolMapUint64Uint16 = &sync.Pool{New: func() interface{} {
	tmp := make(map[uint64]uint16, 1<<20)
	return &tmp
}}

// NewFragmentComparator creates a new FragmentComparator with
// the given options and pool of chainers.
func NewFragmentComparator(options *FragmentComparatorOptions, poolChainers *sync.Pool) *FragmentComparator {
	cpr := &FragmentComparator{
		options:      options,
		poolChainers: poolChainers,

		ccc: util.Ns(0b01, options.K),
		ggg: util.Ns(0b10, options.K),
		ttt: util.Ns(0b11, options.K),

		entries: make([]rtree.BatchEntry, 0, 1<<20),
	}

	return cpr
}

// Compare compares fragsA and fragsB by counting shared k-mers.
// Do not forget to call RecycleFragmentCompareResult on the result to recycle the internal buffer.
func (cpr *FragmentComparator) Compare(fragsA, fragsB *[][]byte) (*[]uint64, error) {
	entriesA, err := cpr.IndexA(fragsA)
	if err != nil {
		return nil, err
	}
	pairs, err := cpr.CompareWithIndexedA(entriesA, fragsB)
	RecycleResultOfIndexA(entriesA)
	return pairs, err
}

var poolBatchEntries = &sync.Pool{New: func() interface{} {
	return make([]rtree.BatchEntry, 0, 10<<20)
}}

// IndexA pre-computes canonical k-mer entries for fragsA, sorted by Key, so
// that the result can be reused across many CompareWithIndexedA calls (one
// query vs. many subject genomes). The returned slice is owned by the caller
// and read concurrently by workers; it must not be mutated after this call.
func (cpr *FragmentComparator) IndexA(fragsA *[][]byte) ([]rtree.BatchEntry, error) {
	entries := poolBatchEntries.Get().([]rtree.BatchEntry)
	entries, err := cpr.collectEntries(fragsA, 0, entries[:0])
	if err != nil {
		RecycleResultOfIndexA(entries)
		return nil, err
	}
	sorts.ByUint64(rtree.BatchEntries(entries))
	return entries, nil
}

// RecycleResultOfIndexA recycles the result of IndexA.
func RecycleResultOfIndexA(entries []rtree.BatchEntry) {
	if entries != nil {
		entries = entries[:0]
		poolBatchEntries.Put(entries)
	}
}

// CompareWithIndexedA compares pre-indexed entries of genome A (sorted by Key,
// produced by IndexA) against fragsB of genome B. entriesA is read-only.
// Avoids copying entriesA by collecting fragsB into a separate buffer
// (cpr.entries), sorting only that buffer, and merging the two sorted streams.
func (cpr *FragmentComparator) CompareWithIndexedA(entriesA []rtree.BatchEntry, fragsB *[][]byte) (*[]uint64, error) {
	entriesB, err := cpr.collectEntries(fragsB, 1, cpr.entries[:0])
	if err != nil {
		return nil, err
	}
	cpr.entries = entriesB

	sorts.ByUint64(rtree.BatchEntries(entriesB))

	return cpr.scanPairsMerged(entriesA, entriesB), nil
}

// collectEntries iterates k-mers in frags, applies the same filters as Compare,
// and appends one canonical entry per surviving position to dst. genomeBit
// (0 or 1) is OR'ed into Val's lowest bit. If dst is nil a fresh buffer is
// allocated.
func (cpr *FragmentComparator) collectEntries(frags *[][]byte, genomeBit uint32, dst []rtree.BatchEntry) ([]rtree.BatchEntry, error) {
	k := cpr.options.K
	scaled := uint64(cpr.options.Scaled)
	useSketching := scaled > 1

	ccc := cpr.ccc
	ggg := cpr.ggg
	ttt := cpr.ttt

	if dst == nil {
		dst = make([]rtree.BatchEntry, 0, 1<<20)
	}

	var kmer, kmerRC, canon uint64
	var ok bool
	for i, s := range *frags {
		iter, err := iterator.NewKmerIterator(s, int(k))
		if err != nil {
			return nil, err
		}

		for {
			kmer, kmerRC, ok, _ = iter.NextKmer()
			if !ok {
				break
			}

			if kmer == 0 || kmer == ccc || kmer == ggg || kmer == ttt ||
				util.IsLowComplexityDust(kmer, k) {
				continue
			}

			canon = min(kmer, kmerRC)

			if useSketching && util.Hash64(canon)%scaled != 0 {
				continue
			}

			dst = append(dst, rtree.BatchEntry{Key: canon, Val: uint32(i)<<1 | genomeBit})
		}
	}
	return dst, nil
}

// scanPairsMerged walks the two sorted streams entriesA (genome A) and entriesB
// (genome B) in lockstep. For every shared Key, it enumerates the cross product
// of the A-run and B-run and accumulates per-pair counts in a pooled map.
func (cpr *FragmentComparator) scanPairsMerged(entriesA, entriesB []rtree.BatchEntry) *[]uint64 {
	counter := poolMapUint64Uint16.Get().(*map[uint64]uint16)

	nA, nB := len(entriesA), len(entriesB)
	var ia, ib int
	for ia < nA && ib < nB {
		ka := entriesA[ia].Key
		kb := entriesB[ib].Key
		if ka < kb {
			ia++
			continue
		}
		if kb < ka {
			ib++
			continue
		}
		// equal: find end of run on each side
		ja := ia + 1
		for ja < nA && entriesA[ja].Key == ka {
			ja++
		}
		jb := ib + 1
		for jb < nB && entriesB[jb].Key == ka {
			jb++
		}
		// enumerate cross pairs
		for p := ia; p < ja; p++ {
			fa := uint64(entriesA[p].Val >> 1)
			rowKey := fa << 32
			for q := ib; q < jb; q++ {
				fb := uint64(entriesB[q].Val >> 1)
				(*counter)[rowKey|fb]++
			}
		}
		ia = ja
		ib = jb
	}

	pairs := poolUint64s.Get().(*[]uint64)
	*pairs = (*pairs)[:0]
	threshold := cpr.options.MinSharedKmers
	for key, v := range *counter {
		if v >= threshold {
			*pairs = append(*pairs, key)
		}
	}

	topNFragments := cpr.options.TopNFragments
	if topNFragments > 0 {
		sortutil.Uint64s(*pairs)

		ma := poolFragPairMap.Get().(*map[uint64]*[]uint64)

		var ia, ib uint64
		var ls *[]uint64
		var ok bool
		for _, p := range *pairs {
			ia, ib = p>>32, p&4294967295

			if ls, ok = (*ma)[ia]; !ok {
				ls = poolKmerAndLocs.Get().(*[]uint64)
				*ls = (*ls)[:0]
				(*ma)[ia] = ls
			}
			*ls = append(*ls, ib)
		}

		*pairs = (*pairs)[:0]
		for ia, ls := range *ma {
			if len(*ls) > topNFragments {
				slices.SortFunc(*ls, func(a, b uint64) int {
					return int((*counter)[ia<<32|b]) - int((*counter)[ia<<32|a])
				})

				*ls = (*ls)[:topNFragments]
			}
			for _, ib := range *ls {
				*pairs = append(*pairs, ia<<32|ib)
			}
		}

		for _, ls := range *ma {
			*ls = (*ls)[:0]
			poolKmerAndLocs.Put(ls)
		}
		clear(*ma)
		poolFragPairMap.Put(ma)
	}

	clear(*counter)
	poolMapUint64Uint16.Put(counter)

	return pairs
}

func RecycleFragmentCompareResult(pairs *[]uint64) {
	if pairs != nil {
		*pairs = (*pairs)[:0]
		poolUint64s.Put(pairs)
	}
}

var poolFragPairMap = &sync.Pool{New: func() interface{} {
	tmp := make(map[uint64]*[]uint64, 4096)
	return &tmp
}}
