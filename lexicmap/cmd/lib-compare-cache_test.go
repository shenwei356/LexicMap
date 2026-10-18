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
	"errors"
	"math"
	"math/rand"
	"reflect"
	"sync"
	"sync/atomic"
	"testing"
	"unsafe"

	rtree "github.com/shenwei356/LexicMap/lexicmap/cmd/tree"
)

func TestCompareCacheEvictionAndPinnedFallback(t *testing.T) {
	c := newCompareGenomeCache(200, func(key string) (*preparedCompareGenome, error) {
		return &preparedCompareGenome{bytes: 120}, nil
	})
	a := c.acquire("A")
	b := c.acquire("B") // A is pinned; B must run uncached without waiting.
	if c.used != 120 || b.retained {
		t.Fatalf("pinned fallback: used=%d retained=%v", c.used, b.retained)
	}
	c.release(b)
	c.release(a)
	b = c.acquire("B") // Evicts idle A.
	if c.entries["A"] != nil || c.evictions != 1 || c.used > c.budget {
		t.Fatalf("eviction: entries=%v evictions=%d used=%d", c.entries, c.evictions, c.used)
	}
	c.release(b)
	b2 := c.acquire("B")
	if b2 != b || c.loads != 3 || c.hits != 1 {
		t.Fatalf("reuse: loads=%d hits=%d", c.loads, c.hits)
	}
	c.release(b2)
	c.close()
	if c.used != 0 || len(c.entries) != 0 || c.idle.Len() != 0 {
		t.Fatal("close retained cache entries")
	}
}

func TestCompareCacheConcurrentLoadAndError(t *testing.T) {
	for _, loadErr := range []error{nil, errors.New("read failed")} {
		var loads atomic.Int32
		started, finish := make(chan struct{}), make(chan struct{})
		c := newCompareGenomeCache(100, func(string) (*preparedCompareGenome, error) {
			loads.Add(1)
			close(started)
			<-finish
			return &preparedCompareGenome{bytes: 80}, loadErr
		})
		const workers = 8
		got := make(chan *compareGenomeCacheEntry, workers)
		var wg sync.WaitGroup
		for range workers {
			wg.Add(1)
			go func() { defer wg.Done(); got <- c.acquire("same") }()
		}
		<-started
		close(finish)
		wg.Wait()
		close(got)
		var first *compareGenomeCacheEntry
		for e := range got {
			if first == nil {
				first = e
			}
			if e != first || !errors.Is(e.err, loadErr) {
				t.Fatalf("load was not shared: entry=%p first=%p err=%v", e, first, e.err)
			}
			c.release(e)
		}
		if loads.Load() != 1 {
			t.Fatalf("loads=%d", loads.Load())
		}
		if loadErr != nil && len(c.entries) != 0 {
			t.Fatal("failed entry remained cached")
		}
		c.close()
	}
}

func TestCompareCacheOversizedAndEmpty(t *testing.T) {
	for _, value := range []*preparedCompareGenome{nil, {bytes: 101}} {
		c := newCompareGenomeCache(100, func(string) (*preparedCompareGenome, error) { return value, nil })
		a, b := c.acquire("same"), c.acquire("same")
		if a != b || a.retained || c.used != 0 {
			t.Fatal("oversized/empty admission")
		}
		c.release(a)
		c.release(b)
		if len(c.entries) != 0 {
			t.Fatal("uncached entry retained after final release")
		}
	}
}

func TestRecycleCachedPairKeepsSharedSequences(t *testing.T) {
	seq := []byte("ACGTACGT")
	g := &GQuery{seqs: []*[]byte{&seq}, genomeSize: len(seq)}
	c := newCompareGenomeCache(100, func(string) (*preparedCompareGenome, error) {
		return &preparedCompareGenome{genome: g, bytes: 80}, nil
	})
	q := poolGPair.Get().(*GPair)
	q.cache = c
	q.entries = [2]*compareGenomeCacheEntry{c.acquire("self"), c.acquire("self")}
	q.views = [2]GQuery{*g, *g}
	q.g1, q.g2 = &q.views[0], &q.views[1]
	RecycleGPair(q)
	if len(seq) != 8 || len(g.seqs) != 1 || g.result != nil {
		t.Fatal("pair recycling modified shared genome")
	}
	e := c.acquire("self")
	if e.value.genome != g || c.loads != 1 {
		t.Fatal("self-pair did not preserve its cached genome")
	}
	c.release(e)
	c.close()
}

func TestCompareCacheLRUOrder(t *testing.T) {
	c := newCompareGenomeCache(160, func(string) (*preparedCompareGenome, error) {
		return &preparedCompareGenome{bytes: 80}, nil
	})
	for _, key := range []string{"A", "B", "A", "C"} {
		c.release(c.acquire(key))
	}
	if c.entries["A"] == nil || c.entries["B"] != nil || c.entries["C"] == nil {
		t.Fatal("cache did not evict the least recently used idle entry")
	}
	c.close()
}

// compareTestIndex uses the CLI's alignment settings without reading an index.
func compareTestIndex(t testing.TB) *Index {
	t.Helper()
	idx, err := NewGenomeComparator("", &IndexSearchingOptions{NoIndex: true, ExtendLength2: 50, MaxEvalue: 1e-15})
	if err != nil {
		t.Fatal(err)
	}
	idx.SetSeqCompareOptions(&SeqComparatorOptions{K: 31, MinPrefix: 11,
		Chaining2Options:   Chaining2Options{MaxGap: 100, MinScore: 21, MinAlignLen: 30, MinIdentity: 70, BandBase: 100, BandCount: 50},
		MinAlignedFraction: 30, MinIdentity: 70})
	idx.SetFragmentCompareOptions(&FragmentComparatorOptions{K: 11, MinSharedKmers: 3, Scaled: 4, TopNFragments: 5})
	return idx
}

func TestPreparedGenomeComparisonMatchesUncached(t *testing.T) {
	idx := compareTestIndex(t)
	rng := rand.New(rand.NewSource(7))
	seq := make([]byte, 4400)
	for i := range seq {
		seq[i] = "ACGT"[rng.Intn(4)]
	}
	other := append([]byte(nil), seq...)
	for i := 15; i < len(other); i += 101 {
		other[i] = "ACGT"[rng.Intn(4)]
	}
	copy(seq[1150:1160], "NNNNNNNNNN")
	copy(other[1150:1160], "NNNNNNNNNN")
	a1, a2 := seq[:2200], seq[2200:]
	b1, b2 := other[:2200], other[2200:]
	q := &GQuery{seqs: []*[]byte{&a1, &a2}, genomeSize: len(seq)}
	s := &GQuery{seqs: []*[]byte{&b1, &b2}, genomeSize: len(other)}
	for _, ortho := range []bool{false, true} {
		minSeedLen := gsa3SampledK
		if ortho {
			minSeedLen = orthoANISeedK
		}
		qp, err := idx.prepareCompareGenome(q, 1020, minSeedLen, ortho)
		if err != nil {
			t.Fatal(err)
		}
		sp, err := idx.prepareCompareGenome(s, 1020, minSeedLen, ortho)
		if err != nil {
			t.Fatal(err)
		}
		cache := newCompareGenomeCache(10<<20, func(key string) (*preparedCompareGenome, error) {
			if key == "query" {
				return qp, nil
			}
			return sp, nil
		})
		cache.cacheFragmentIndexes = true
		defer cache.close()
		qe, se := cache.acquire("query"), cache.acquire("subject")
		defer cache.release(qe)
		defer cache.release(se)
		for _, reverse := range []bool{false, true} {
			if reverse {
				qp, sp = sp, qp
			}
			oldQ, oldS := *qp.genome, *sp.genome
			newQ, newS := oldQ, oldS
			if ortho {
				err = idx.CompareTwoGenomesOrthoANI(&oldQ, &oldS, 1020, minSeedLen, 0, .7)
				if err == nil {
					err = idx.compareTwoGenomesOrthoANIPrepared(&newQ, &newS, qp, sp, 1020, minSeedLen, 0, .7)
				}
			} else {
				err = idx.CompareTwoGenomes(&oldQ, &oldS, 1020, minSeedLen, 0, .7)
				if err == nil {
					err = idx.compareTwoGenomesPrepared(&newQ, &newS, qp, sp, 1020, minSeedLen, 0, .7)
				}
			}
			if err != nil {
				t.Fatal(err)
			}
			if oldQ.result == nil || newQ.result == nil || len(*oldQ.result) != 1 || len(*newQ.result) != 1 {
				t.Fatal("expected aligned genomes")
			}
			want, got := *(*oldQ.result)[0], *(*newQ.result)[0]
			// RBH summation order can vary with Go map iteration; compare floats
			// within rounding noise and require all discrete alignment fields equal.
			wv, gv := reflect.ValueOf(want), reflect.ValueOf(got)
			for i := 0; i < wv.NumField(); i++ {
				x, y := wv.Field(i), gv.Field(i)
				if x.Kind() == reflect.Float64 {
					if math.Abs(x.Float()-y.Float()) > 1e-9 {
						t.Fatalf("ortho=%v reverse=%v field=%s: %v != %v", ortho, reverse, wv.Type().Field(i).Name, x, y)
					}
				} else if !reflect.DeepEqual(x.Interface(), y.Interface()) {
					t.Fatalf("field %s differs", wv.Type().Field(i).Name)
				}
			}
			RecycleGSearchResults(oldQ.result)
			RecycleGSearchResults(newQ.result)
		}
		if q.result != nil || s.result != nil {
			t.Fatal("shared genome acquired pair-local results")
		}
		// Independent pair views concurrently borrow the same seeds/maps/entries.
		var wg sync.WaitGroup
		errs := make(chan error, 4)
		for range 4 {
			wg.Add(1)
			go func() {
				defer wg.Done()
				cq, cs := *qp.genome, *sp.genome
				var err error
				if ortho {
					err = idx.compareTwoGenomesOrthoANIPrepared(&cq, &cs, qp, sp, 1020, minSeedLen, 0, .7)
				} else {
					err = idx.compareTwoGenomesPrepared(&cq, &cs, qp, sp, 1020, minSeedLen, 0, .7)
				}
				if err == nil && (cq.result == nil || len(*cq.result) != 1) {
					err = errors.New("missing concurrent alignment")
				}
				if cq.result != nil {
					RecycleGSearchResults(cq.result)
				}
				errs <- err
			}()
		}
		wg.Wait()
		close(errs)
		for err := range errs {
			if err != nil {
				t.Fatal(err)
			}
		}
	}
}

func TestCompareFineIndexCacheBudgetAndOwnership(t *testing.T) {
	idx := compareTestIndex(t)
	sequence := fragmentIndexTestSequence()
	fragments := [][]byte{sequence}
	const budget = 1 << 20
	c := newCompareGenomeCache(budget, func(key string) (*preparedCompareGenome, error) {
		if key == "large" {
			return &preparedCompareGenome{bytes: budget}, nil
		}
		return &preparedCompareGenome{bytes: 512, fragments: &fragments}, nil
	})
	c.cacheFragmentIndexes = true
	e := c.acquire("A")
	indexes := e.value.indexes
	if indexes == nil {
		t.Fatal("retained genome has no fine-index slots")
	}
	metadataCharge := c.used
	var wg sync.WaitGroup
	for range 8 {
		wg.Add(1)
		go func() {
			defer wg.Done()
			cpr := idx.poolSeqComparator.Get().(*SeqComparator)
			defer idx.poolSeqComparator.Put(cpr)
			if err := indexes.index(cpr, 0, sequence); err != nil {
				t.Error(err)
				return
			}
			result, err := cpr.Compare(0, uint32(len(sequence)), sequence, len(sequence))
			if err != nil || result == nil {
				t.Errorf("cached fine alignment: result=%v err=%v", result, err)
			}
			if result != nil {
				RecycleSeqComparatorResult(result)
			}
			cpr.RecycleIndex()
		}()
	}
	wg.Wait()
	charge := int64(len(indexes.fragments[0].entries)) * int64(unsafe.Sizeof(rtree.BatchEntry{}))
	if charge == 0 || c.used != metadataCharge+charge || e.value.bytes != c.used || c.peak > budget {
		t.Fatalf("lazy charge: used=%d metadata=%d entries=%d", c.used, metadataCharge, charge)
	}
	if idx.fragmentIndexBytes.Load() != 0 {
		t.Fatal("comparison borrowed the separate genome-search cache budget")
	}
	c.release(e)
	second := c.acquire("A")
	if second != e || second.value.indexes != indexes || c.used != metadataCharge+charge {
		t.Fatal("repeated pair did not reuse the same charged entries")
	}
	c.release(second)
	large := c.acquire("large")
	if c.evictions != 1 || e.value != nil || c.used != budget {
		t.Fatal("eviction did not release the base and lazy entry charges")
	}
	c.release(large)
	c.close()
	if c.used != 0 {
		t.Fatal("close retained fine-index charges")
	}
}

func TestCompareFineIndexCacheFullFallsBack(t *testing.T) {
	idx := compareTestIndex(t)
	sequence := fragmentIndexTestSequence()
	fragments := [][]byte{sequence}
	c := newCompareGenomeCache(2048, func(string) (*preparedCompareGenome, error) {
		return &preparedCompareGenome{bytes: 512, fragments: &fragments}, nil
	})
	c.cacheFragmentIndexes = true
	defer c.close()
	e := c.acquire("A")
	defer c.release(e)
	indexes := e.value.indexes
	if indexes == nil {
		t.Fatal("expected enough budget for slot metadata")
	}
	before := c.used
	for range 2 { // the first builder and later cache-full callers both work
		cpr := idx.poolSeqComparator.Get().(*SeqComparator)
		if err := indexes.index(cpr, 0, sequence); err != nil {
			t.Fatal(err)
		}
		result, err := cpr.Compare(0, uint32(len(sequence)), sequence, len(sequence))
		if err != nil || result == nil {
			t.Fatalf("fallback alignment: result=%v err=%v", result, err)
		}
		RecycleSeqComparatorResult(result)
		cpr.RecycleIndex()
		idx.poolSeqComparator.Put(cpr)
	}
	if c.used != before || indexes.fragments[0].entries != nil {
		t.Fatal("full cache retained uncharged entries")
	}
}

func TestCompareFineIndexCacheDisabled(t *testing.T) {
	fragments := [][]byte{fragmentIndexTestSequence()}
	c := newCompareGenomeCache(1<<20, func(string) (*preparedCompareGenome, error) {
		return &preparedCompareGenome{bytes: 512, fragments: &fragments}, nil
	})
	e := c.acquire("A")
	if e.value.indexes != nil || c.used != 512 {
		t.Fatal("single-pair cache retained fine-index metadata")
	}
	c.release(e)
	c.close()
}
