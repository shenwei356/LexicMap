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
	"os"
	"path/filepath"
	"reflect"
	"sync"
	"testing"

	"github.com/shenwei356/bio/taxdump"
)

func TestKeepGenomeByTaxIDDoesNotDropFirstAllowedHit(t *testing.T) {
	idx := &Index{
		opt: &IndexSearchingOptions{
			NegativeTaxIds: []uint32{3},
		},
		filterByNegativeTaxId: true,
		genomeIdx2TaxId: map[uint64]uint32{
			10: 2,
			11: 4,
		},
		Taxonomy: &taxdump.Taxonomy{
			Nodes: map[uint32]uint32{
				1: 1,
				2: 1,
				3: 1,
				4: 3,
			},
			MergeNodes: make(map[uint32]uint32),
		},
	}
	cache := make(map[uint64]bool)

	if !idx.keepGenomeByTaxID(&cache, 10) {
		t.Fatal("the first hit outside the negative TaxId subtree should be kept")
	}
	if !idx.keepGenomeByTaxID(&cache, 10) {
		t.Fatal("the cached allowed hit should be kept")
	}
	if idx.keepGenomeByTaxID(&cache, 11) {
		t.Fatal("a hit inside the negative TaxId subtree should be rejected")
	}
}

func TestCompactMaskSelection(t *testing.T) {
	selection := []bool{false, true, false, true, true, false}
	maskIndexes := map[int]struct{}{1: {}, 4: {}, 5: {}}

	gotSlots, gotCount := compactMaskSelection(len(selection), selection, maskIndexes)
	wantSlots := []int32{-1, 0, -1, -1, 1, -1}
	if !reflect.DeepEqual(gotSlots, wantSlots) {
		t.Fatalf("mask slots: got %v, want %v", gotSlots, wantSlots)
	}
	if gotCount != 2 {
		t.Fatalf("selected mask count: got %d, want 2", gotCount)
	}

	gotSlots, gotCount = compactMaskSelection(6, nil, nil)
	if gotSlots != nil || gotCount != 6 {
		t.Fatalf("identity mask selection: got slots=%v count=%d, want slots=<nil> count=6", gotSlots, gotCount)
	}
}

func TestMergeGSearchScreenResultDetail(t *testing.T) {
	dst := &GSearchScreenResultDetail{
		BatchGenomeIndex: []uint64{1},
		SumPrefix:        45,
		LongestMatches:   []uint8{20, 0, 25},
	}
	src := &GSearchScreenResultDetail{
		BatchGenomeIndex: []uint64{2},
		SumPrefix:        50,
		LongestMatches:   []uint8{21, 29, 0},
	}

	mergeGSearchScreenResultDetail(dst, src)

	if got, want := dst.BatchGenomeIndex, []uint64{1, 2}; !reflect.DeepEqual(got, want) {
		t.Fatalf("batch genome indexes: got %v, want %v", got, want)
	}
	if got, want := dst.LongestMatches, []uint8{21, 29, 25}; !reflect.DeepEqual(got, want) {
		t.Fatalf("longest matches: got %v, want %v", got, want)
	}
	if got, want := dst.SumPrefix, uint64(75); got != want {
		t.Fatalf("sum prefix: got %d, want %d", got, want)
	}
}

func TestResetGSearchScreenResultDetailRetainsClearedLongestMatches(t *testing.T) {
	r := &GSearchScreenResultDetail{
		BatchGenomeIndex: []uint64{1},
		SumPrefix:        31,
		LongestMatches:   make([]uint8, 1024),
	}
	for i := range r.LongestMatches {
		r.LongestMatches[i] = 31
	}

	resetGSearchScreenResultDetail(r)

	if len(r.BatchGenomeIndex) != 0 || r.SumPrefix != 0 {
		t.Fatalf("result was not reset: %+v", r)
	}
	if len(r.LongestMatches) != 1024 || cap(r.LongestMatches) != 1024 {
		t.Fatalf("longest matches not retained: len=%d cap=%d", len(r.LongestMatches), cap(r.LongestMatches))
	}
	for i, v := range r.LongestMatches {
		if v != 0 {
			t.Fatalf("longest match %d not cleared: %d", i, v)
		}
	}
}

func TestGSearchWindowsRejectsTooManyWindows(t *testing.T) {
	if _, err := gsearchWindows(2_000, 10_000_000, 31); err == nil {
		t.Fatal("expected an error when screening windows are shorter than k")
	}
}

func TestWindowSkipRegionsClipsAndTranslatesCoordinates(t *testing.T) {
	regions := []int{5, 12, 25, 35, 50, 60}
	got := windowSkipRegions(regions, 10, 30)
	want := []int{0, 2, 15, 19}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("window skip regions: got %v, want %v", got, want)
	}
}

// TestTrimGSearchScreenResultsKeepsCutoffTies covers counting genomes, retaining
// the complete cutoff score group, stable ordering, and recycling discarded candidates.
func TestTrimGSearchScreenResultsKeepsCutoffTies(t *testing.T) {
	tests := []struct {
		name    string
		scores  []uint64
		topN    int
		wantIDs []uint64
	}{
		{"cutoff in top group", []uint64{10, 10, 10, 10, 9, 9, 9, 8, 7, 6}, 3, []uint64{1, 2, 3, 4}},
		{"cutoff in second group", []uint64{10, 10, 10, 10, 9, 9, 9, 8, 7, 6}, 5, []uint64{1, 2, 3, 4, 5, 6, 7}},
		{"cutoff after earlier ties", []uint64{10, 10, 9, 9, 8}, 3, []uint64{1, 2, 3, 4}},
		{"exact group boundary", []uint64{10, 10, 9, 9, 8}, 2, []uint64{1, 2}},
		{"top group only", []uint64{10, 10, 9, 8}, 1, []uint64{1, 2}},
		{"all tied", []uint64{10, 10, 10}, 2, []uint64{1, 2, 3}},
		{"multiple groups with cutoff ties", []uint64{10, 10, 9, 9, 9}, 3, []uint64{1, 2, 3, 4, 5}},
		{"distinct scores", []uint64{6, 8, 10, 9, 7}, 3, []uint64{3, 4, 2}},
		{"unlimited", []uint64{8, 10, 9, 9}, 0, []uint64{2, 3, 4, 1}},
		{"larger than candidates", []uint64{9, 10}, 10, []uint64{2, 1}},
		{"equal to candidates", []uint64{9, 10, 10}, 3, []uint64{2, 3, 1}},
		{"empty", nil, 3, []uint64{}},
	}
	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			rs := make([]*GSearchScreenResultDetail, len(test.scores))
			for i, score := range test.scores {
				// Reverse input to test sorting by score and internal genome index.
				rs[len(rs)-1-i] = &GSearchScreenResultDetail{
					// Multiple chunks of one genome still count as one candidate.
					BatchGenomeIndex: []uint64{uint64(i + 1), uint64(i + 101)},
					SumPrefix:        score,
					LongestMatches:   []uint8{uint8(score)},
				}
			}
			original := append([]*GSearchScreenResultDetail(nil), rs...)
			backing := rs // retain the full slice to check cleared discarded slots
			idx := &Index{poolGSearchDetailResult: &sync.Pool{}}
			trimGSearchScreenResults(&rs, test.topN, idx)
			got := make([]uint64, len(rs))
			kept := make(map[*GSearchScreenResultDetail]bool, len(rs))
			for i, r := range rs {
				got[i] = r.BatchGenomeIndex[0]
				kept[r] = true
				if r.SumPrefix != test.scores[got[i]-1] {
					t.Fatalf("retained candidate %d was reset", got[i])
				}
			}
			if !reflect.DeepEqual(got, test.wantIDs) {
				t.Fatalf("ordered candidates: got %v, want %v", got, test.wantIDs)
			}
			for _, r := range original {
				if !kept[r] && (r.SumPrefix != 0 || len(r.BatchGenomeIndex) != 0 || r.LongestMatches[0] != 0) {
					t.Fatalf("discarded candidate was not recycled: %+v", r)
				}
			}
			for _, r := range backing[len(rs):] {
				if r != nil {
					t.Fatal("discarded slot still retains a candidate")
				}
			}
		})
	}
}

// TestTrimGSearchResultsKeepsDistinctANI checks printed ANI groups,
// rounded-ANI sorting, and clearing and recycling discarded results.
func TestTrimGSearchResultsKeepsDistinctANI(t *testing.T) {
	tests := []struct {
		name    string
		anis    []float64
		topN    int
		wantIDs []uint64
	}{
		{"three ANI groups", []float64{1, 1, 1, 1, .99, .99, .99, .98, .97, .96}, 3, []uint64{1, 2, 3, 4, 5, 6, 7, 8}},
		{"ties above cutoff", []float64{1, 1, .99, .99, .98}, 2, []uint64{1, 2, 3, 4}},
		{"top group only", []float64{1, 1, .99}, 1, []uint64{1, 2}},
		{"all tied", []float64{1, 1, 1}, 2, []uint64{1, 2, 3}},
		{"fewer ANI groups", []float64{1, 1, .99, .99, .99}, 3, []uint64{1, 2, 3, 4, 5}},
		{"distinct ANI", []float64{.96, .98, 1, .99, .97}, 3, []uint64{3, 4, 2}},
		{"equal after output rounding", []float64{.9999900, .9999901, .9999900}, 1, []uint64{1, 2, 3}},
		{"different integer widths", []float64{.09, .1, 1, .99}, 3, []uint64{3, 4, 2}},
		{"distinct at four decimals", []float64{.999960, .999964}, 1, []uint64{2}},
		{"rounding boundary", []float64{.99996449, .99996451, .99996448}, 1, []uint64{2}},
		{"lower rounded group", []float64{.99996449, .99996451, .99996448}, 2, []uint64{2, 1, 3}},
		{"zero ANI", []float64{0, 0, 0}, 1, []uint64{1, 2, 3}},
		{"unlimited", []float64{.98, 1, .99, .99}, 0, []uint64{2, 3, 4, 1}},
		{"larger than candidates", []float64{.99, 1}, 10, []uint64{2, 1}},
		{"empty", nil, 3, []uint64{}},
	}
	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			rs := make([]*GSearchResult, len(test.anis))
			for i, ani := range test.anis {
				// Reverse input to check stable ordering within an ANI group.
				rs[len(rs)-1-i] = &GSearchResult{BatchGenomeIndex: uint64(i + 1), ANI: ani, AFq: .8, AFs: .7}
			}
			original := append([]*GSearchResult(nil), rs...)
			backing := rs // full slice for checking discarded pointer slots
			trimGSearchResults(&rs, test.topN)
			got := make([]uint64, len(rs))
			kept := make(map[*GSearchResult]bool, len(rs))
			printedGroups := make(map[string]struct{}) // group keys from the actual output format
			for i, r := range rs {
				got[i] = r.BatchGenomeIndex
				kept[r] = true
				printedGroups[fmt.Sprintf("%.4f", r.ANI*100)] = struct{}{}
				if r.ANI != test.anis[got[i]-1] || r.AFq != .8 || r.AFs != .7 {
					t.Fatalf("retained result was reset: %+v", r)
				}
			}
			if !reflect.DeepEqual(got, test.wantIDs) {
				t.Fatalf("ordered results: got %v, want %v", got, test.wantIDs)
			}
			if test.topN > 0 {
				allGroups := make(map[string]struct{})
				for _, ani := range test.anis {
					allGroups[fmt.Sprintf("%.4f", ani*100)] = struct{}{}
				}
				if want := min(test.topN, len(allGroups)); len(printedGroups) != want {
					t.Fatalf("printed ANI groups: got %d, want %d", len(printedGroups), want)
				}
			}
			for _, r := range original {
				if !kept[r] && (r.ANI != 0 || r.AFq != 0 || r.AFs != 0 || r.BatchGenomeIndex != ^uint64(0)) {
					t.Fatalf("discarded result was not reset: %+v", r)
				}
			}
			for _, r := range backing[len(rs):] {
				if r != nil {
					t.Fatal("discarded slot still retains a result")
				}
			}
		})
	}
}

// TestTrimGSearchResultsSortsANIQueryAFSubjectAF checks every sorting key and
// compares AF within a printed ANI group despite differences in unrounded ANI.
func TestTrimGSearchResultsSortsANIQueryAFSubjectAF(t *testing.T) {
	rs := []*GSearchResult{
		{BatchGenomeIndex: 6, ANI: .99, AFq: 1, AFs: 1},
		{BatchGenomeIndex: 5, ANI: .9999604, AFq: .8, AFs: .9},
		{BatchGenomeIndex: 4, ANI: .9999603, AFq: .9, AFs: .7},
		{BatchGenomeIndex: 3, ANI: .9999602, AFq: .9, AFs: .8},
		{BatchGenomeIndex: 2, ANI: .9999601, AFq: .9, AFs: .8},
		{BatchGenomeIndex: 1, ANI: .9999600, AFq: 1, AFs: .6},
	}
	trimGSearchResults(&rs, 1)
	got := make([]uint64, len(rs))
	for i, r := range rs {
		got[i] = r.BatchGenomeIndex
	}
	if want := []uint64{1, 2, 3, 4, 5}; !reflect.DeepEqual(got, want) {
		t.Fatalf("ordered results: got %v, want %v", got, want)
	}
}

func TestRecycleGSearchScreenDetailResultsAcceptsNil(t *testing.T) {
	(&Index{}).RecycleGSearchScreenDetailResults(nil)
}

func TestTaxidFileExclusions(t *testing.T) {
	for _, tc := range []struct {
		name, text         string
		args               []string
		positive, negative []uint32
	}{
		{"negative-only", "-2\n-2\n-3\n", nil, nil, []uint32{2, 3}},
		{"mixed", "2\n-2\n0\n3\n-4\n", nil, []uint32{2, 3}, []uint32{2, 4}},
		{"combined", "2\n-3\n4\n-5\n", []string{"2", "-3"}, []uint32{2, 4}, []uint32{3, 5}},
	} {
		t.Run(tc.name, func(t *testing.T) {
			file := filepath.Join(t.TempDir(), "taxids.txt")
			if err := os.WriteFile(file, []byte(tc.text), 0600); err != nil {
				t.Fatal(err)
			}
			positive, negative := parseTaxids("taxdump", "mapping", tc.args, file)
			if !reflect.DeepEqual(positive, tc.positive) || !reflect.DeepEqual(negative, tc.negative) {
				t.Fatalf("got %v / %v, want %v / %v", positive, negative, tc.positive, tc.negative)
			}
		})
	}
}
