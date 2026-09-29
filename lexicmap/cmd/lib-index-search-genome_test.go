package cmd

import (
	"reflect"
	"testing"
)

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
