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
	"fmt"

	rtree "github.com/shenwei356/LexicMap/lexicmap/cmd/tree"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
	"github.com/shenwei356/lexichash/iterator"
)

type FragmentComparator struct {
	k             uint8
	ccc, ggg, ttt uint64

	flat []uint32
}

func NewFragmentComparator(k uint8) *FragmentComparator {
	cpr := &FragmentComparator{
		k: k,

		ccc: util.Ns(0b01, k),
		ggg: util.Ns(0b10, k),
		ttt: util.Ns(0b11, k),

		flat: make([]uint32, 0, 128),
	}

	return cpr
}

func (cpr *FragmentComparator) Compare(fragsA, fragsB *[][]byte) error {
	// ----------------------------------------------------
	// 1. Find similar fragment pairs

	// 1.1 build a prefix tree

	k := cpr.k

	// a reusable Radix tree for searching k-mers sharing at least n-base prefixes.
	t := rtree.NewTree(k)

	var kmer, kmerRC uint64
	var ok bool

	ccc := cpr.ccc
	ggg := cpr.ggg
	ttt := cpr.ttt

	//  The last two bits of the uint32 tree node value are used to mark genome and kmer info
	//  .... 0 1
	//         | --- 0 for genome A,        1 for genome B
	//       |------ 0 for positive strand, 1 for negative strand
	//  ____-------- index of fragment
	//
	for i, s := range *fragsA {
		iter, err := iterator.NewKmerIterator(s, int(k))
		if err != nil {
			return err
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

			t.Insert(kmer, uint32(i<<2))
			t.Insert(kmerRC, uint32(i<<2|2))
		}
	}

	for i, s := range *fragsB {
		iter, err := iterator.NewKmerIterator(s, int(k))
		if err != nil {
			return err
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

			t.Insert(kmer, uint32(i<<2|1))
			t.Insert(kmerRC, uint32(i<<2|3))

		}
	}

	// 1.2 find pairs

	// key: uint32(frag index in genome A) | uint32 (frag index in genome B)
	m := make(map[uint64]uint16)

	// Group all k-mers sharing a prefix of at least X bases (singletons too).
	// Flatten each group's vals and enumerate all unordered value pairs:
	// this covers both same-k-mer pairs (originally a separate Walk pass)
	// and shared-prefix pairs in one tree traversal.
	flat := cpr.flat
	t.WalkGroups(15, func(keys []uint64, vals [][]uint32, lenPrefix uint8) bool {
		flat = flat[:0]
		for _, v := range vals {
			flat = append(flat, v...)
		}
		n := len(flat)
		if n < 2 {
			return false
		}

		var i, j int
		var a, b uint32
		for i = 0; i < n; i++ {
			for j = i + 1; j < n; j++ {
				a, b = flat[i], flat[j]
				if a&1 == b&1 { // belong to the same genome
					continue
				}
				if a&1 == 1 { // a belongs to genome B
					a, b = b, a
				}
				m[uint64(a>>2)<<32|uint64(b>>2)]++
			}
		}
		return false
	})

	// 1.3 filter pairs

	for key, v := range m {
		fmt.Printf("%d\t%d\t%d\n", key>>32, (key & 4294967295), v)
	}

	// ----------------------------------------------------

	rtree.RecycleTree(t)

	return nil
}
