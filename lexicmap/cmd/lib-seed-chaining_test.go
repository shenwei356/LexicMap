package cmd

import (
	"fmt"
	"reflect"
	"sync"
	"testing"
)

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
		subs := make([]*SubstrPair, len(anchors))
		for j := range anchors {
			subs[j] = &anchors[j]
		}
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
						ClearSubstrPairs(poolSub, r.Subs, idx.k)
						chainer := idx.poolChainers.Get().(*Chainer)
						r.Chains, r.Score = chainer.Chain(r.Subs)
						RecycleChainer(idx.poolChainers, chainer)
						if r.Score < idx.chainingOptions.MinScore {
							idx.RecycleSearchResult(r)
							continue
						}
						r.prepareAlignmentRegions()
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
					if r.Subs != nil || r.Chains != nil {
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
	subs := make([]*SubstrPair, chainerInitSize+1)
	for i := range subs {
		subs[i] = &SubstrPair{QBegin: int32(i * 64), TBegin: int32(i * 64), Len: 31}
	}
	single := []*SubstrPair{{Len: 31}}
	got := idx.chainSeedResults([][]*SearchResult{{{Subs: &subs}, {Subs: &single}}}, 2)
	if nCreated < 2 {
		t.Fatalf("oversized chainer was retained: created %d chainers", nCreated)
	}
	if len(*got) != 2 {
		t.Fatalf("got %d targets, want 2", len(*got))
	}
	idx.RecycleSearchResults(got)
}
