package kv

import "testing"

func TestRecycleSearchResultsBoundsSeedPositions(t *testing.T) {
	for _, capacity := range []int{maxPooledSeedPositions, maxPooledSeedPositions + 1} {
		results := []SearchResult{{IQuery: 3, Values: make([]uint64, 1, capacity)}, {IQuery: 4, Values: make([]uint64, 1, capacity)}}
		backing := results
		results = results[:1]
		RecycleSearchResults(&results)
		for _, r := range backing {
			if r.IQuery != 0 || len(r.Values) != 0 {
				t.Fatal("result was not reset")
			}
			if capacity > maxPooledSeedPositions {
				if r.Values != nil {
					t.Fatal("oversized seed data array was retained")
				}
			} else if cap(r.Values) != capacity {
				t.Fatal("small seed data array was not retained")
			}
		}
		if len(results) != 0 {
			t.Fatal("result list was not reset")
		}
	}
}

func TestRecycleSearchResultsDropsShortenedOversizedSlice(t *testing.T) {
	results := make([]SearchResult, 1, maxPooledSearchResults+1)
	results[0].Values = make([]uint64, 1, maxPooledSeedPositions+1)
	backing := results[:cap(results)]
	backing[len(backing)-1].Values = []uint64{9}
	RecycleSearchResults(&results)
	if results != nil || backing[0].Values != nil || backing[len(backing)-1].Values != nil {
		t.Fatal("oversized arrays were retained")
	}
}

func TestSearchResultBuffersSurviveGrowthAndReuse(t *testing.T) {
	results := make([]SearchResult, 0, 1)
	for i := range 1000 {
		r := appendSearchResult(&results)
		r.IQuery = i
		r.Values = append(r.Values, uint64(i))
	}
	for i, r := range results {
		if r.IQuery != i || len(r.Values) != 1 || r.Values[0] != uint64(i) {
			t.Fatal("growth changed a previous result")
		}
	}
	bufferStart := &results[0].Values[0]
	results = results[:0]
	r := appendSearchResult(&results)
	r.Values = append(r.Values, 1234)
	if &r.Values[0] != bufferStart {
		t.Fatal("seed data buffer was not reused")
	}
	RecycleSearchResults(&results)
}

func TestTransferredSearchResultsOwnSeedDataBuffers(t *testing.T) {
	for _, capacity := range []int{1, 4} {
		source := []SearchResult{{IQuery: 7, Values: []uint64{11, 12}}}
		destination := make([]SearchResult, 1, capacity)
		destination[0] = SearchResult{IQuery: 3, Values: []uint64{4}}
		var unusedSeedData *uint64
		if capacity > 1 {
			backing := destination[:capacity]
			backing[1].Values = []uint64{21, 22}
			unusedSeedData = &backing[1].Values[0]
		}
		AppendSearchResults(&destination, &source)
		if len(source) != 0 || len(destination) != 2 || destination[0].Values[0] != 4 {
			t.Fatal("transfer did not append results and reset the source")
		}
		r := appendSearchResult(&source)
		r.Values = append(r.Values, 99)
		if unusedSeedData != nil && &r.Values[0] != unusedSeedData {
			t.Fatal("transfer discarded the destination's unused seed data buffer")
		}
		if destination[1].IQuery != 7 || destination[1].Values[0] != 11 || destination[1].Values[1] != 12 {
			t.Fatal("source reuse overwrote transferred seed data")
		}
		RecycleSearchResults(&source)
		RecycleSearchResults(&destination)
	}
}
