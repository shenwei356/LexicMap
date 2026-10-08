package tree

import (
	"math/rand"
	"reflect"
	"slices"
	"testing"
)

func TestInsertSortedBatchPreservesLookupAndEntries(t *testing.T) {
	rng := rand.New(rand.NewSource(23))
	for _, k := range []uint8{6, 31, 32} {
		mask := ^uint64(0) >> (64 - 2*k)
		entries := []BatchEntry{{Key: 5, Val: 3}, {Key: 10, Val: 2}, {Key: 5, Val: 1}}
		for i := range 2000 {
			entries = append(entries, BatchEntry{Key: rng.Uint64() & mask, Val: uint32(i)})
		}
		ref := NewTree(k)
		ref.InsertBatch(entries)
		before := slices.Clone(entries)
		got := NewTree(k)
		got.InsertSortedBatch(entries)
		if !reflect.DeepEqual(entries, before) {
			t.Fatal("sorted bulk loading mutated shared entries")
		}
		queries := []uint64{0, mask, 256}
		for _, entry := range entries {
			queries = append(queries, entry.Key, entry.Key^(1<<(2*k-2)))
		}
		for _, p := range []uint8{0, 1, 4, 5, 11, k, 255} {
			for _, key := range queries {
				var want, actual []SearchResult
				wf := ref.SearchEach(key, p, func(key uint64, length uint8, values []uint32) {
					want = append(want, SearchResult{Kmer: key, LenPrefix: length, Values: slices.Clone(values)})
				})
				gf := got.SearchEach(key, p, func(key uint64, length uint8, values []uint32) {
					actual = append(actual, SearchResult{Kmer: key, LenPrefix: length, Values: slices.Clone(values)})
				})
				if wf != gf || !reflect.DeepEqual(want, actual) {
					t.Fatalf("k=%d p=%d key=%x: sorted bulk lookup differs", k, p, key)
				}
			}
		}
		RecycleTree(ref)
		RecycleTree(got)
	}

	// Preserve the current lookup behavior at short compressed edges as well.
	// This optimization must not silently change legacy prefix matching.
	entries := []BatchEntry{{Key: 69, Val: 3}, {Key: 74, Val: 2}}
	ref, got := NewTree(6), NewTree(6)
	defer RecycleTree(ref)
	defer RecycleTree(got)
	ref.InsertBatch(entries)
	got.InsertSortedBatch(entries)
	var want, actual []SearchResult
	ref.SearchEach(0, 5, func(key uint64, length uint8, values []uint32) {
		want = append(want, SearchResult{Kmer: key, LenPrefix: length, Values: slices.Clone(values)})
	})
	got.SearchEach(0, 5, func(key uint64, length uint8, values []uint32) {
		actual = append(actual, SearchResult{Kmer: key, LenPrefix: length, Values: slices.Clone(values)})
	})
	if len(want) != 2 || !reflect.DeepEqual(want, actual) {
		t.Fatal("short compressed-edge lookup changed")
	}
}
