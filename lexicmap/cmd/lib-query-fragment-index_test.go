package cmd

import (
	"math/rand"
	"reflect"
	"strings"
	"sync"
	"testing"
	"unsafe"

	rtree "github.com/shenwei356/LexicMap/lexicmap/cmd/tree"
)

// fragmentIndexTestSequence includes repeats and ambiguous bases in a short
// deterministic query, exercising both canonical strands and duplicate keys.
func fragmentIndexTestSequence() []byte {
	rng := rand.New(rand.NewSource(21))
	s := make([]byte, 1020)
	for i := range s {
		s[i] = "ACGT"[rng.Intn(4)]
	}
	copy(s[600:700], s[100:200])
	copy(s[450:460], "NNNNNNNNNN")
	return s
}

func TestCachedSeqComparisonMatchesUncached(t *testing.T) {
	idx := compareTestIndex(t)
	sequence := fragmentIndexTestSequence()
	subject := append([]byte(nil), sequence...)
	for i := 20; i < len(subject); i += 71 {
		subject[i] = "ACGT"[(i/71)%4]
	}
	cache := newQueryFragmentIndexes(idx, 1, 4)
	defer cache.close()
	ref := idx.poolSeqComparator.Get().(*SeqComparator)
	defer idx.poolSeqComparator.Put(ref)
	if err := ref.Index(sequence); err != nil {
		t.Fatal(err)
	}
	defer ref.RecycleIndex()
	positiveMatches := 0
	for _, reverse := range []bool{false, true} {
		if reverse {
			RC(subject)
		}
		for _, bounds := range [][2]uint32{{0, uint32(len(sequence))}, {80, 810}} {
			want, err := ref.Compare(bounds[0], bounds[1], subject, len(sequence))
			if err != nil {
				t.Fatalf("reference comparison reverse=%v bounds=%v: result=%v err=%v", reverse, bounds, want, err)
			}
			if !reverse && want != nil {
				positiveMatches++
			}
			cpr := idx.poolSeqComparator.Get().(*SeqComparator)
			if err := cache.index(cpr, 0, sequence); err != nil {
				t.Fatal(err)
			}
			got, err := cpr.Compare(bounds[0], bounds[1], subject, len(sequence))
			if err != nil || !reflect.DeepEqual(got, want) {
				t.Fatalf("reverse=%v bounds=%v: cached result differs: got=%+v want=%+v err=%v", reverse, bounds, got, want, err)
			}
			if want != nil {
				RecycleSeqComparatorResult(want)
				RecycleSeqComparatorResult(got)
			}
			cpr.RecycleIndex()
			idx.poolSeqComparator.Put(cpr)
		}
	}
	if positiveMatches != 2 {
		t.Fatal("full and partial forward comparisons must exercise nonempty alignments")
	}
}

func TestQueryFragmentIndexesConcurrentOwnership(t *testing.T) {
	idx := compareTestIndex(t)
	sequence := fragmentIndexTestSequence()
	cache := newQueryFragmentIndexes(idx, 1, 8)
	var wg sync.WaitGroup
	for range 8 {
		wg.Add(1)
		go func() {
			defer wg.Done()
			cpr := idx.poolSeqComparator.Get().(*SeqComparator)
			defer idx.poolSeqComparator.Put(cpr)
			if err := cache.index(cpr, 0, sequence); err != nil {
				t.Error(err)
				return
			}
			borrowed := cache.fragments[0].entries
			if borrowed == nil || (len(cpr.entries) > 0 && &cpr.entries[0] == &borrowed[0]) {
				t.Error("shared entries remained mutable comparator scratch")
			}
			result, err := cpr.Compare(0, uint32(len(sequence)), sequence, len(sequence))
			if err != nil || result == nil {
				t.Errorf("shared comparison: result=%v err=%v", result, err)
			} else {
				RecycleSeqComparatorResult(result)
			}
			cpr.RecycleIndex()
			if err := cpr.Index(sequence[:100]); err != nil {
				t.Error(err)
			}
			if &cpr.entries[0] == &borrowed[0] {
				t.Error("private indexing overwrote a borrowed index")
			}
			cpr.RecycleIndex()
		}()
	}
	wg.Wait()
	if want := int64(len(cache.fragments[0].entries)) * int64(unsafe.Sizeof(rtree.BatchEntry{})); idx.fragmentIndexBytes.Load() != want || cache.reserved.Load() != want {
		t.Fatal("concurrent users charged the same fragment more than once")
	}
	cache.close()
	if idx.fragmentIndexBytes.Load() != 0 || cache.fragments != nil {
		t.Fatal("closed query retained shared indexes or budget")
	}
}

func TestQueryFragmentIndexesShareBudgetAcrossQueries(t *testing.T) {
	idx := compareTestIndex(t)
	sequence := fragmentIndexTestSequence()
	first := newQueryFragmentIndexes(idx, 1, 2)
	cpr := idx.poolSeqComparator.Get().(*SeqComparator)
	if err := first.index(cpr, 0, sequence); err != nil {
		t.Fatal(err)
	}
	charge := first.reserved.Load()
	cpr.RecycleIndex()
	idx.poolSeqComparator.Put(cpr)
	first.close()

	// Leave space for exactly one fragment among eight concurrent queries.
	base := maxGenomeSearchFragmentIndexMemory - charge
	idx.fragmentIndexBytes.Store(base)
	caches := make([]*queryFragmentIndexes, 8)
	var wg sync.WaitGroup
	for i := range caches {
		caches[i] = newQueryFragmentIndexes(idx, 1, 2)
		wg.Add(1)
		go func(c *queryFragmentIndexes) {
			defer wg.Done()
			cpr := idx.poolSeqComparator.Get().(*SeqComparator)
			defer idx.poolSeqComparator.Put(cpr)
			if err := c.index(cpr, 0, sequence); err != nil {
				t.Error(err)
			}
			cpr.RecycleIndex()
		}(caches[i])
	}
	wg.Wait()
	if idx.fragmentIndexBytes.Load() != maxGenomeSearchFragmentIndexMemory {
		t.Fatal("concurrent queries did not respect their shared budget")
	}
	retained := 0
	for _, c := range caches {
		if c.fragments[0].entries != nil {
			retained++
		}
		c.close()
	}
	if retained != 1 || idx.fragmentIndexBytes.Load() != base {
		t.Fatal("cache admission or release charged another query's memory")
	}
}

func BenchmarkQueryFragmentIndex(b *testing.B) {
	sequence := fragmentIndexTestSequence()
	for _, variant := range []string{"original", "cached", "full"} {
		b.Run(variant, func(b *testing.B) {
			idx := compareTestIndex(b)
			cpr := idx.poolSeqComparator.Get().(*SeqComparator)
			defer idx.poolSeqComparator.Put(cpr)
			var cache *queryFragmentIndexes
			if variant != "original" {
				cache = newQueryFragmentIndexes(idx, 1, 2)
				defer cache.close()
				if variant == "full" {
					idx.fragmentIndexBytes.Store(maxGenomeSearchFragmentIndexMemory)
				}
			}
			if err := cache.index(cpr, 0, sequence); err != nil {
				b.Fatal(err)
			}
			cpr.RecycleIndex()
			b.ReportAllocs()
			for b.Loop() {
				if err := cache.index(cpr, 0, sequence); err != nil {
					b.Fatal(err)
				}
				cpr.RecycleIndex()
			}
		})
	}
}

func TestQueryFragmentIndexesFullBudgetAndSingleton(t *testing.T) {
	idx := compareTestIndex(t)
	sequence := fragmentIndexTestSequence()
	idx.fragmentIndexBytes.Store(maxGenomeSearchFragmentIndexMemory)
	for _, subjects := range []int{1, 2} {
		cache := newQueryFragmentIndexes(idx, 1, subjects)
		for range 2 {
			cpr := idx.poolSeqComparator.Get().(*SeqComparator)
			if err := cache.index(cpr, 0, sequence); err != nil {
				t.Fatal(err)
			}
			if cpr.tree == nil || cache.reserved.Load() != 0 {
				t.Fatal("full budget or singleton did not build a private tree")
			}
			cpr.RecycleIndex()
			idx.poolSeqComparator.Put(cpr)
		}
		cache.close()
		if idx.fragmentIndexBytes.Load() != maxGenomeSearchFragmentIndexMemory {
			t.Fatal("fallback changed another query's budget")
		}
	}
}

func TestQueryFragmentIndexesPublishError(t *testing.T) {
	idx := compareTestIndex(t)
	cache := newQueryFragmentIndexes(idx, 1, 2)
	defer cache.close()
	for range 2 {
		cpr := idx.poolSeqComparator.Get().(*SeqComparator)
		if err := cache.index(cpr, 0, []byte("ACGT")); err == nil {
			t.Fatal("expected the short fragment's indexing error")
		}
		cpr.RecycleIndex()
		idx.poolSeqComparator.Put(cpr)
	}
	if idx.fragmentIndexBytes.Load() != 0 {
		t.Fatal("failed index consumed cache memory")
	}
}

func TestQueryFragmentIndexesEmptyFilteredEntries(t *testing.T) {
	idx := compareTestIndex(t)
	cache := newQueryFragmentIndexes(idx, 1, 2)
	sequence := []byte(strings.Repeat("A", 1020))
	for range 2 {
		cpr := idx.poolSeqComparator.Get().(*SeqComparator)
		if err := cache.index(cpr, 0, sequence); err != nil {
			t.Fatal(err)
		}
		result, err := cpr.Compare(0, uint32(len(sequence)), sequence, len(sequence))
		if err != nil || result != nil {
			t.Fatal("filtered homopolymer produced a match")
		}
		cpr.RecycleIndex()
		idx.poolSeqComparator.Put(cpr)
	}
	if cache.fragments[0].entries == nil || len(cache.fragments[0].entries) != 0 || idx.fragmentIndexBytes.Load() != 0 {
		t.Fatal("empty prepared entries were not shared without charging memory")
	}
	cache.close()
	cache.close()
	if idx.fragmentIndexBytes.Load() != 0 {
		t.Fatal("closing an empty cache changed the memory budget")
	}
}

func TestQueryFragmentIndexesBoundLargeQueryMetadata(t *testing.T) {
	idx := compareTestIndex(t)
	cache := newQueryFragmentIndexes(idx, 1<<20, 2)
	defer cache.close()
	if len(cache.fragments) != maxGenomeSearchFragmentIndexSlots {
		t.Fatal("long query allocated unbounded fragment metadata")
	}
	cpr := idx.poolSeqComparator.Get().(*SeqComparator)
	defer idx.poolSeqComparator.Put(cpr)
	if err := cache.index(cpr, 1<<20-1, fragmentIndexTestSequence()); err != nil {
		t.Fatal(err)
	}
	defer cpr.RecycleIndex()
	if cpr.tree == nil || idx.fragmentIndexBytes.Load() != 0 {
		t.Fatal("fragment beyond the cache slots did not use uncached indexing")
	}
}

// A 10 Mb query must also reuse its final fragment, which exceeded the old
// 8192-slot limit even when shared entry memory was available.
func TestQueryFragmentIndexesCacheFinalFragmentOf10MbQuery(t *testing.T) {
	const genomeLength = 10_000_000
	const fragmentLength = 1020
	fragments := (genomeLength + fragmentLength - 1) / fragmentLength
	idx := compareTestIndex(t)
	cache := newQueryFragmentIndexes(idx, fragments, 2)
	defer cache.close()
	if len(cache.fragments) != fragments {
		t.Fatal("10 Mb query's final fragment has no cache slot")
	}
	sequence := fragmentIndexTestSequence()[:genomeLength%fragmentLength]
	cpr := idx.poolSeqComparator.Get().(*SeqComparator)
	defer idx.poolSeqComparator.Put(cpr)
	if err := cpr.Index(sequence); err != nil {
		t.Fatal(err)
	}
	want := append([]rtree.BatchEntry(nil), cpr.entries...)
	cpr.RecycleIndex()
	charge := int64(len(want)) * int64(unsafe.Sizeof(rtree.BatchEntry{}))
	for range 2 {
		if err := cache.index(cpr, fragments-1, sequence); err != nil {
			t.Fatal(err)
		}
		if !reflect.DeepEqual(cache.fragments[fragments-1].entries, want) {
			t.Fatal("10 Mb query's final fragment was not cached with the original entries")
		}
		if cpr.tree == nil || cache.reserved.Load() != charge || idx.fragmentIndexBytes.Load() != charge {
			t.Fatal("reusing the final fragment changed its cache charge or failed to build a tree")
		}
		cpr.RecycleIndex()
	}
	cache.close()
	if idx.fragmentIndexBytes.Load() != 0 {
		t.Fatal("10 Mb query's cache did not release its shared budget")
	}
}
