package cmd

import (
	"fmt"
	"reflect"
	"slices"
	"sync"
	"testing"
)

func TestTrimSeedSearchResultsKeepsCutoffTies(t *testing.T) {
	for _, tt := range []struct {
		name   string
		scores []float32
		topN   int
		want   []uint64
	}{
		{name: "empty", topN: 1, want: []uint64{}},
		{name: "unlimited", scores: []float32{80, 100, 90}, want: []uint64{0, 1, 2}},
		{name: "below limit", scores: []float32{80, 100}, topN: 3, want: []uint64{0, 1}},
		{name: "at limit", scores: []float32{80, 100}, topN: 2, want: []uint64{0, 1}},
		{name: "unique cutoff", scores: []float32{80, 100, 90, 70}, topN: 2, want: []uint64{1, 2}},
		{name: "cutoff ties", scores: []float32{80, 90, 100, 90, 70, 90}, topN: 2, want: []uint64{1, 2, 3, 5}},
		{name: "top one ties", scores: []float32{90, 100, 80, 100}, topN: 1, want: []uint64{1, 3}},
		{name: "all tied", scores: []float32{90, 90, 90}, topN: 1, want: []uint64{0, 1, 2}},
		{name: "ties reach end", scores: []float32{90, 100, 90}, topN: 2, want: []uint64{0, 1, 2}},
		{name: "near scores are distinct", scores: []float32{90, 90.00001, 80}, topN: 1, want: []uint64{1}},
	} {
		t.Run(tt.name, func(t *testing.T) {
			// Different arrival orders must retain the same candidate set.
			for shift := range max(1, len(tt.scores)) {
				rs := make([]*SearchResult, len(tt.scores))
				for i := range rs {
					id := (i + shift) % len(rs)
					rs[i] = &SearchResult{BatchGenomeIndex: uint64(id), Score: tt.scores[id]}
				}
				backing := rs
				trimSeedSearchResults(&rs, tt.topN)
				got := make([]uint64, len(rs))
				for i, r := range rs {
					got[i] = r.BatchGenomeIndex
					if r.Score != tt.scores[r.BatchGenomeIndex] {
						t.Fatal("candidate score changed")
					}
					if tt.topN > 0 && len(backing) > tt.topN && i > 0 && rs[i-1].Score < r.Score {
						t.Fatal("candidates are not sorted by descending chaining score")
					}
				}
				slices.Sort(got)
				if !reflect.DeepEqual(got, tt.want) {
					t.Fatalf("arrival shift %d: got candidates %v, want %v", shift, got, tt.want)
				}
				for _, r := range backing[len(rs):] {
					if r != nil {
						t.Fatal("discarded candidate remains referenced by the backing slice")
					}
				}
			}
		})
	}
}

func seedChainingTestIndex(workers int) *Index {
	idx := &Index{k: 31, opt: &IndexSearchingOptions{
		NumCPUs: workers, MinSinglePrefix: 20, MaxGap: 50, MaxDistance: 1000,
	}}
	idx.initChainingResources()
	return idx
}

func seedChainingTestResults(n int) [][]*SearchResult {
	results := make([][]*SearchResult, 3)
	for i := range n {
		var anchors []SubstrPair
		switch i % 5 {
		case 0: // Below the minimum single-anchor score.
			anchors = []SubstrPair{{Len: 15}}
		case 1: // A single reverse-strand anchor.
			anchors = []SubstrPair{{QBegin: 10, TBegin: 100, Len: 31, TRC: true}}
		case 2:
			anchors = []SubstrPair{{TBegin: 1000, Len: 31}, {QBegin: 60, TBegin: 1061, Len: 31}}
		case 3:
			anchors = []SubstrPair{{TBegin: 2000, Len: 31, TRC: true}, {QBegin: 60, TBegin: 1940, Len: 31, TRC: true}}
		case 4: // Unsorted, duplicate and nested anchors.
			anchors = []SubstrPair{{QBegin: 60, TBegin: 3065, Len: 31},
				{TBegin: 3000, Len: 31}, {QBegin: 1, TBegin: 3001, Len: 20}, {TBegin: 3000, Len: 31}}
		}
		subs := anchors
		r := &SearchResult{BatchGenomeIndex: uint64(i), Subs: &subs}
		results[i%len(results)] = append(results[i%len(results)], r)
	}
	return results
}

func TestChainSeedResultsMatchesSerial(t *testing.T) {
	for _, n := range []int{0, 1, 17, 500, 2100} {
		for _, workers := range []int{1, 4, 32} {
			t.Run(fmt.Sprintf("targets=%d/workers=%d", n, workers), func(t *testing.T) {
				idx := seedChainingTestIndex(workers)
				serial := seedChainingTestResults(n)
				want := make(map[uint64]*SearchResult)
				for _, targets := range serial {
					for _, r := range targets {
						ClearSubstrPairs(r.Subs, idx.k)
						chainer := idx.poolChainers.Get().(*Chainer)
						r.resetChainRegions()
						r.chainRegions, r.Score = chainer.Chain(r.Subs, r.chainRegions)
						RecycleChainer(idx.poolChainers, chainer)
						if r.Score < idx.chainingOptions.MinScore {
							idx.RecycleSearchResult(r)
							continue
						}
						RecycleSubstrPairs(poolSubs, r.Subs)
						r.Subs = nil
						want[r.BatchGenomeIndex] = r
					}
				}
				got := idx.chainSeedResults(seedChainingTestResults(n), n)
				if len(*got) != len(want) {
					t.Fatalf("got %d targets, want %d", len(*got), len(want))
				}
				seen := make(map[uint64]bool)
				for _, r := range *got {
					w := want[r.BatchGenomeIndex]
					if w == nil || seen[r.BatchGenomeIndex] {
						t.Fatalf("unexpected or duplicate target %d", r.BatchGenomeIndex)
					}
					seen[r.BatchGenomeIndex] = true
					if r.Score != w.Score || !reflect.DeepEqual(r.chainRegions, w.chainRegions) {
						t.Fatalf("target %d: score or chain bounds differ", r.BatchGenomeIndex)
					}
					if r.Subs != nil {
						t.Fatalf("target %d retains anchors or chain paths", r.BatchGenomeIndex)
					}
				}
				idx.RecycleSearchResults(got)
				for _, r := range want {
					idx.RecycleSearchResult(r)
				}
			})
		}
	}
}

func TestChainSeedResultsDropsOversizedChainer(t *testing.T) {
	idx := seedChainingTestIndex(1)
	nCreated := 0 // The sole worker finishes before this counter is read.
	idx.poolChainers = &sync.Pool{New: func() any {
		nCreated++
		return NewChainer(idx.chainingOptions)
	}}
	subs := make([]SubstrPair, chainerInitSize+1)
	for i := range subs {
		subs[i] = SubstrPair{QBegin: int32(i * 64), TBegin: int32(i * 64), Len: 31}
	}
	single := []SubstrPair{{Len: 31}}
	got := idx.chainSeedResults([][]*SearchResult{{{Subs: &subs}, {Subs: &single}}}, 2)
	if nCreated < 2 {
		t.Fatalf("oversized chainer was retained: created %d chainers", nCreated)
	}
	if len(*got) != 2 {
		t.Fatalf("got %d targets, want 2", len(*got))
	}
	idx.RecycleSearchResults(got)
}
