// Copyright © 2026 Wei Shen <shenwei356@gmail.com>

package cmd

import (
	"errors"
	"math"
	"math/rand"
	"reflect"
	"sync"
	"sync/atomic"
	"testing"
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
func compareTestIndex(t *testing.T) *Index {
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
		qp, err := idx.prepareCompareGenome(q, 1020, 100, ortho)
		if err != nil {
			t.Fatal(err)
		}
		sp, err := idx.prepareCompareGenome(s, 1020, 100, ortho)
		if err != nil {
			t.Fatal(err)
		}
		for _, reverse := range []bool{false, true} {
			if reverse {
				qp, sp = sp, qp
			}
			oldQ, oldS := *qp.genome, *sp.genome
			newQ, newS := oldQ, oldS
			if ortho {
				err = idx.CompareTwoGenomesOrthoANI(&oldQ, &oldS, 1020, 100, 0, .7)
				if err == nil {
					err = idx.compareTwoGenomesOrthoANIPrepared(&newQ, &newS, qp, sp, 1020, 100, 0, .7)
				}
			} else {
				err = idx.CompareTwoGenomes(&oldQ, &oldS, 1020, 100, 0, .7)
				if err == nil {
					err = idx.compareTwoGenomesPrepared(&newQ, &newS, qp, sp, 1020, 100, 0, .7)
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
					err = idx.compareTwoGenomesOrthoANIPrepared(&cq, &cs, qp, sp, 1020, 100, 0, .7)
				} else {
					err = idx.compareTwoGenomesPrepared(&cq, &cs, qp, sp, 1020, 100, 0, .7)
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
