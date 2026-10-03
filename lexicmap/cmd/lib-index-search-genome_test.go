package cmd

import (
	"reflect"
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

func TestTrimGSearchScreenResultsKeepsCutoffTies(t *testing.T) {
	rs := []*GSearchScreenResultDetail{
		{BatchGenomeIndex: []uint64{8}, SumPrefix: 80},
		{BatchGenomeIndex: []uint64{4}, SumPrefix: 90},
		{BatchGenomeIndex: []uint64{3}, SumPrefix: 90},
		{BatchGenomeIndex: []uint64{1}, SumPrefix: 100},
	}

	trimGSearchScreenResults(&rs, 2, nil)
	if got, want := len(rs), 3; got != want {
		t.Fatalf("result count: got %d, want %d", got, want)
	}
	got := []uint64{rs[0].BatchGenomeIndex[0], rs[1].BatchGenomeIndex[0], rs[2].BatchGenomeIndex[0]}
	want := []uint64{1, 3, 4}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("ordered candidates: got %v, want %v", got, want)
	}
}

func TestTrimGSearchResultsKeepsStrictTopN(t *testing.T) {
	rs := []*GSearchResult{
		{BatchGenomeIndex: 8, ANI: 0.90, AFq: 0.8, AFs: 0.7},
		{BatchGenomeIndex: 4, ANI: 0.95, AFq: 0.8, AFs: 0.7},
		{BatchGenomeIndex: 3, ANI: 0.95, AFq: 0.8, AFs: 0.7},
	}

	trimGSearchResults(&rs, 2)
	got := []uint64{rs[0].BatchGenomeIndex, rs[1].BatchGenomeIndex}
	want := []uint64{3, 4}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("ordered ANI results: got %v, want %v", got, want)
	}
}

func TestRecycleGSearchScreenDetailResultsAcceptsNil(t *testing.T) {
	(&Index{}).RecycleGSearchScreenDetailResults(nil)
}
