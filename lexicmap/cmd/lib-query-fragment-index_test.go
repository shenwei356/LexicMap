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
	"math/rand"
	"reflect"
	"strings"
	"sync"
	"testing"

	rtree "github.com/shenwei356/LexicMap/lexicmap/cmd/tree"
	"github.com/twotwotwo/sorts"
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
	metadataCharge := cache.reserved.Load()
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
			borrowed := cache.fragments[0].tree
			if borrowed == nil || cpr.compact != borrowed || cpr.tree != nil {
				t.Error("comparator did not borrow the complete cached tree")
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
			if cpr.compact != nil {
				t.Error("private indexing reused a borrowed tree")
			}
			cpr.RecycleIndex()
		}()
	}
	wg.Wait()
	if want := metadataCharge + cache.fragments[0].tree.RetainedBytes(); idx.fragmentIndexBytes.Load() != want || cache.reserved.Load() != want {
		t.Fatal("concurrent users charged the same fragment more than once")
	}
	cache.close()
	if idx.fragmentIndexBytes.Load() != 0 || cache.fragments != nil {
		t.Fatal("closed query retained shared indexes or budget")
	}
}

// Reservations from concurrent queries share one limit, including slot metadata.
func TestQueryFragmentIndexesShareBudgetAcrossQueries(t *testing.T) {
	idx := compareTestIndex(t)
	sequence := fragmentIndexTestSequence()
	probe := newQueryFragmentIndexes(idx, 1, 2)
	metadata := probe.reserved.Load()
	cpr := idx.poolSeqComparator.Get().(*SeqComparator)
	if err := probe.index(cpr, 0, sequence); err != nil {
		t.Fatal(err)
	}
	charge := probe.reserved.Load() - metadata
	cpr.RecycleIndex()
	idx.poolSeqComparator.Put(cpr)
	probe.close()

	// Leave room for one real tree after reserving all eight queries' slots.
	caches := make([]*queryFragmentIndexes, 8)
	var metadataCharge int64
	for i := range caches {
		caches[i] = newQueryFragmentIndexes(idx, 1, 2)
		metadataCharge += caches[i].reserved.Load()
	}
	base := idx.opt.MaxGenomeCacheMemory - metadataCharge - charge
	idx.fragmentIndexBytes.Add(base)
	var wg sync.WaitGroup
	for _, c := range caches {
		wg.Add(1)
		go func(c *queryFragmentIndexes) {
			defer wg.Done()
			cpr := idx.poolSeqComparator.Get().(*SeqComparator)
			defer idx.poolSeqComparator.Put(cpr)
			if err := c.index(cpr, 0, sequence); err != nil {
				t.Error(err)
				return
			}
			defer cpr.RecycleIndex()
			result, err := cpr.Compare(0, uint32(len(sequence)), sequence, len(sequence))
			if err != nil || result == nil {
				t.Errorf("cached or fallback comparison: result=%v err=%v", result, err)
			}
			if result != nil {
				RecycleSeqComparatorResult(result)
			}
		}(c)
	}
	wg.Wait()
	if idx.fragmentIndexBytes.Load() != idx.opt.MaxGenomeCacheMemory {
		t.Fatal("concurrent queries did not respect the shared budget")
	}
	retained := 0
	for _, c := range caches {
		if c.fragments[0].tree != nil {
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
				if variant == "full" {
					idx.fragmentIndexBytes.Store(idx.opt.MaxGenomeCacheMemory)
				}
				cache = newQueryFragmentIndexes(idx, 1, 2)
				defer cache.close()
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
	idx.fragmentIndexBytes.Store(idx.opt.MaxGenomeCacheMemory)
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
		if idx.fragmentIndexBytes.Load() != idx.opt.MaxGenomeCacheMemory {
			t.Fatal("fallback changed another query's budget")
		}
	}
}

func TestQueryFragmentIndexesPublishError(t *testing.T) {
	idx := compareTestIndex(t)
	cache := newQueryFragmentIndexes(idx, 1, 2)
	metadataCharge := cache.reserved.Load()
	defer cache.close()
	for range 2 {
		cpr := idx.poolSeqComparator.Get().(*SeqComparator)
		if err := cache.index(cpr, 0, []byte("ACGT")); err == nil {
			t.Fatal("expected the short fragment's indexing error")
		}
		cpr.RecycleIndex()
		idx.poolSeqComparator.Put(cpr)
	}
	if idx.fragmentIndexBytes.Load() != metadataCharge || cache.fragments[0].tree != nil {
		t.Fatal("failed index consumed tree cache memory")
	}
}

func TestQueryFragmentIndexesEmptyFilteredEntries(t *testing.T) {
	idx := compareTestIndex(t)
	cache := newQueryFragmentIndexes(idx, 1, 2)
	metadataCharge := cache.reserved.Load()
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
	if cache.fragments[0].tree == nil || idx.fragmentIndexBytes.Load() != metadataCharge+cache.fragments[0].tree.RetainedBytes() {
		t.Fatal("empty tree and its metadata were not fully charged")
	}
	cache.close()
	cache.close()
	if idx.fragmentIndexBytes.Load() != 0 {
		t.Fatal("closing an empty cache changed the memory budget")
	}
}

func TestQueryFragmentIndexesBoundLargeQueryMetadata(t *testing.T) {
	idx := compareTestIndex(t)
	idx.opt.MaxGenomeCacheMemory = 1024
	cache := newQueryFragmentIndexes(idx, 1<<20, 2)
	defer cache.close()
	if len(cache.fragments) != 0 || cache.reserved.Load() != 0 {
		t.Fatal("long query allocated unbounded fragment metadata")
	}
	cpr := idx.poolSeqComparator.Get().(*SeqComparator)
	defer idx.poolSeqComparator.Put(cpr)
	if err := cache.index(cpr, 1<<20-1, fragmentIndexTestSequence()); err != nil {
		t.Fatal(err)
	}
	defer cpr.RecycleIndex()
	if cpr.tree == nil || (cpr.compact != nil) || idx.fragmentIndexBytes.Load() != cache.reserved.Load() {
		t.Fatal("fragment beyond the cache slots did not use uncached indexing")
	}
}

// A 20 Mb query must also reuse its final fragment, beyond the former 16384
// slots, when the configured shared budget can accommodate it.
func TestQueryFragmentIndexesCacheFinalFragmentOf20MbQuery(t *testing.T) {
	const genomeLength = 20_000_000
	const fragmentLength = 1020
	fragments := (genomeLength + fragmentLength - 1) / fragmentLength
	idx := compareTestIndex(t)
	cache := newQueryFragmentIndexes(idx, fragments, 2)
	defer cache.close()
	if len(cache.fragments) != fragments {
		t.Fatal("20 Mb query's final fragment has no cache slot")
	}
	sequence := fragmentIndexTestSequence()[:genomeLength%fragmentLength]
	cpr := idx.poolSeqComparator.Get().(*SeqComparator)
	defer idx.poolSeqComparator.Put(cpr)
	if err := cpr.Index(sequence); err != nil {
		t.Fatal(err)
	}
	want := append([]rtree.BatchEntry(nil), cpr.entries...)
	cpr.RecycleIndex()
	metadataCharge := cache.reserved.Load()
	for range 2 {
		if err := cache.index(cpr, fragments-1, sequence); err != nil {
			t.Fatal(err)
		}
		var got []rtree.BatchEntry
		cache.fragments[fragments-1].tree.Walk(func(key uint64, values []uint32) bool {
			for _, value := range values {
				got = append(got, rtree.BatchEntry{Key: key, Val: value})
			}
			return false
		})
		if !reflect.DeepEqual(got, want) {
			t.Fatal("20 Mb query's final fragment was not cached with the original entries")
		}
		charge := metadataCharge + cache.fragments[fragments-1].tree.RetainedBytes()
		if cpr.compact != cache.fragments[fragments-1].tree || cpr.tree != nil || cache.reserved.Load() != charge || idx.fragmentIndexBytes.Load() != charge {
			t.Fatal("reusing the final fragment changed its cache charge or failed to build a tree")
		}
		cpr.RecycleIndex()
	}
	cache.close()
	if idx.fragmentIndexBytes.Load() != 0 {
		t.Fatal("20 Mb query's cache did not release its shared budget")
	}
}

// A tree that cannot fit after slot admission must remain private for both the
// first builder and later callers. Closing releases only this cache's metadata.
func TestQueryFragmentIndexesTreeBudgetFallback(t *testing.T) {
	idx := compareTestIndex(t)
	cache := newQueryFragmentIndexes(idx, 1, 2)
	metadataCharge := cache.reserved.Load()
	if metadataCharge == 0 || len(cache.fragments) != 1 {
		t.Fatal("slot metadata was not charged")
	}
	base := idx.opt.MaxGenomeCacheMemory - metadataCharge - 1
	idx.fragmentIndexBytes.Add(base)
	cpr := idx.poolSeqComparator.Get().(*SeqComparator)
	defer idx.poolSeqComparator.Put(cpr)
	sequence := fragmentIndexTestSequence()
	for range 2 {
		if err := cache.index(cpr, 0, sequence); err != nil {
			t.Fatal(err)
		}
		if cpr.tree == nil || (cpr.compact != nil) || cache.fragments[0].tree != nil || cache.reserved.Load() != metadataCharge {
			t.Fatal("unfunded tree was retained")
		}
		result, err := cpr.Compare(0, uint32(len(sequence)), sequence, len(sequence))
		if err != nil || result == nil {
			t.Fatalf("private fallback failed: %v", err)
		}
		RecycleSeqComparatorResult(result)
		cpr.RecycleIndex()
	}
	cache.close()
	if idx.fragmentIndexBytes.Load() != base {
		t.Fatal("fallback released another owner's charge")
	}
}

// Configured budgets must control real tree admission, including zero and exact
// boundaries. A larger budget admits bytes beyond the former fixed 1 GiB limit.
func TestQueryFragmentIndexesConfiguredBudget(t *testing.T) {
	sequence := fragmentIndexTestSequence()
	probe := compareTestIndex(t)
	cpr := probe.poolSeqComparator.Get().(*SeqComparator)
	entries, err := cpr.collectIndexEntries(sequence)
	if err != nil {
		t.Fatal(err)
	}
	sorts.ByUint64(rtree.BatchEntries(entries))
	charge := queryFragmentSlotBytes(1) + rtree.CompactTreeBytes(cpr.options.K, entries)
	probe.poolSeqComparator.Put(cpr)
	for _, tc := range []struct {
		name         string
		budget, used int64
		cached       bool
	}{
		{"disabled", 0, 0, false},
		{"one-byte-short", charge - 1, 0, false},
		{"exact", charge, 0, true},
		{"above-old-limit", 2 << 30, 1 << 30, true},
		{"default-budget", DefaultIndexSearchingOptions.MaxGenomeCacheMemory, 1 << 30, true},
	} {
		t.Run(tc.name, func(t *testing.T) {
			idx := compareTestIndex(t)
			idx.opt.MaxGenomeCacheMemory = tc.budget
			idx.fragmentIndexBytes.Store(tc.used)
			cache := newQueryFragmentIndexes(idx, 1, 2)
			cpr := idx.poolSeqComparator.Get().(*SeqComparator)
			defer idx.poolSeqComparator.Put(cpr)
			for range 2 {
				if err := cache.index(cpr, 0, sequence); err != nil {
					t.Fatal(err)
				}
				if (cpr.compact != nil) != tc.cached || (cpr.tree != nil) == tc.cached {
					t.Fatal("wrong cached/private index selection")
				}
				cpr.RecycleIndex()
			}
			if idx.fragmentIndexBytes.Load() > tc.budget {
				t.Fatal("configured budget exceeded")
			}
			if tc.cached && cache.reserved.Load() != charge {
				t.Fatal("incomplete cache charge")
			}
			cache.close()
			if idx.fragmentIndexBytes.Load() != tc.used {
				t.Fatal("released another query's charge")
			}
		})
	}
}

func TestQueryFragmentIndexesRejectNegativeBudget(t *testing.T) {
	opt := DefaultIndexSearchingOptions
	opt.MaxGenomeCacheMemory = -1
	if err := CheckIndexSearchingOptions(&opt); err == nil {
		t.Fatal("negative budget passed option validation")
	}
	for _, noIndex := range []bool{false, true} {
		opt.NoIndex = noIndex
		if idx, err := NewIndexSearcher("", &opt); idx != nil || err == nil || !strings.Contains(err.Error(), "MaxGenomeCacheMemory") {
			t.Fatal("negative budget passed constructor validation")
		}
	}
}
