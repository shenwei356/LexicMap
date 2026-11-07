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
	"strings"
	"sync"
	"testing"

	rtree "github.com/shenwei356/LexicMap/lexicmap/cmd/tree"
)

// Exercise grouped multiplicity, saturation, sparse fallback, non-adjacent runs,
// and the same count-only Top-N ties as the positional reference.
func TestGroupedFragmentPairsMatchReference(t *testing.T) {
	for _, sparse := range []bool{false, true} {
		var a, b []rtree.BatchEntry
		for repeat := range 3 {
			for fragment := range 5 {
				id := uint32(fragment)
				if sparse {
					id += 1 << 30
				}
				count := 4 + fragment + repeat
				if fragment == 0 {
					count = 100
				}
				for range count {
					a = append(a, rtree.BatchEntry{Key: 17, Val: id << 1})
					b = append(b, rtree.BatchEntry{Key: 17, Val: uint32(fragment)<<1 | 1})
				}
			}
		}
		originalA, originalB := slices.Clone(a), slices.Clone(b)
		for _, threshold := range []uint16{1, 324, 325, 6300, 65535} {
			for _, topN := range []int{0, 1, 3, 5} {
				opts := &FragmentComparatorOptions{MinSharedKmers: threshold, TopNFragments: topN}
				cpr := &FragmentComparator{options: opts}
				for range 2 {
					got := cpr.scanPairsMerged(a, b)
					slices.Sort(*got)
					want := referenceFragmentPairs(a, b, opts)
					if !slices.Equal(*got, want) {
						t.Fatalf("sparse=%v threshold=%d topN=%d: different pairs", sparse, threshold, topN)
					}
					RecycleFragmentCompareResult(got)
				}
			}
		}
		if !slices.Equal(a, originalA) || !slices.Equal(b, originalB) {
			t.Fatal("shared input changed")
		}
	}
}

// Run counts larger than uint16 must saturate before multiplication.
func TestFragmentOccurrenceSaturation(t *testing.T) {
	entries := make([]rtree.BatchEntry, 70000)
	runs := appendFragmentRuns(nil, make(map[uint32]int), entries)
	if len(runs) != 1 || runs[0].count != 65535 {
		t.Fatal("run count did not saturate")
	}
	if got := addFragmentOccurrences(7, runs[0].count, runs[0].count); got != 65535 {
		t.Fatalf("wrapped count: %d", got)
	}
	if got := addFragmentOccurrences(7, 9, 11); got != 106 {
		t.Fatalf("lost multiplicity: %d", got)
	}
}

// Canonical k-mers must select the same fragment pairs after reversing both
// genomes; reverse fragment order explicitly to compare original coordinates.
func TestGroupedFragmentPairsReverseComplement(t *testing.T) {
	cpr := &FragmentComparator{options: &FragmentComparatorOptions{K: 11, Scaled: 1, MinSharedKmers: 30, TopNFragments: 0}}
	motif := "ACGTTGCATGTCGCATGATGCATGAGAGCT"
	fragments := make([][]byte, 16)
	for i := range fragments {
		if i == 8 {
			motif = "GATCCGATAGCTTCGACCGTATGACCTG"
		}
		fragments[i] = []byte(strings.Repeat(motif, 36))
	}
	collect := func(fragments [][]byte) []rtree.BatchEntry {
		entries, err := cpr.IndexA(&fragments)
		if err != nil {
			t.Fatal(err)
		}
		t.Cleanup(func() { RecycleResultOfIndexA(entries) })
		return entries
	}
	a := collect(fragments)
	before := cpr.scanPairsMerged(a, a)
	defer RecycleFragmentCompareResult(before)
	reversed := make([][]byte, len(fragments))
	for i, fragment := range fragments {
		rc := make([]byte, len(fragment))
		for j, base := range fragment {
			switch base {
			case 'A':
				rc[len(rc)-1-j] = 'T'
			case 'T':
				rc[len(rc)-1-j] = 'A'
			case 'C':
				rc[len(rc)-1-j] = 'G'
			case 'G':
				rc[len(rc)-1-j] = 'C'
			}
		}
		reversed[len(fragments)-1-i] = rc
	}
	b := collect(reversed)
	after := cpr.scanPairsMerged(b, b)
	defer RecycleFragmentCompareResult(after)
	for i, pair := range *after {
		(*after)[i] = (uint64(len(fragments)-1)-(pair>>32))<<32 | (uint64(len(fragments)-1) - (pair & 0xffffffff))
	}
	slices.Sort(*before)
	slices.Sort(*after)
	if !slices.Equal(*before, *after) || len(*before) != 128 {
		t.Fatal("reverse complementation changed fragment pairs")
	}
	if len(cpr.runsA) == 0 {
		t.Fatal("fixture did not activate grouped counting")
	}
}

// Each comparator may group the same prepared entries concurrently.
func TestGroupedFragmentPairsConcurrent(t *testing.T) {
	entries := make([]rtree.BatchEntry, 512)
	for i := range entries {
		entries[i] = rtree.BatchEntry{Key: 1, Val: uint32(i%8) << 1}
	}
	original := slices.Clone(entries)
	opts := &FragmentComparatorOptions{MinSharedKmers: 1, TopNFragments: 3}
	want := referenceFragmentPairs(entries, entries, opts)
	var wg sync.WaitGroup
	for range 8 {
		wg.Add(1)
		go func() {
			defer wg.Done()
			cpr := &FragmentComparator{options: opts}
			for range 20 {
				got := cpr.scanPairsMerged(entries, entries)
				slices.Sort(*got)
				if !slices.Equal(*got, want) {
					t.Error("concurrent grouping changed selected pairs")
				}
				RecycleFragmentCompareResult(got)
			}
		}()
	}
	wg.Wait()
	if !slices.Equal(entries, original) {
		t.Fatal("shared entries changed")
	}
}
