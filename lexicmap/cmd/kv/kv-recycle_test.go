package kv

import "testing"

func TestRecycleSearchResultsBoundsSeedPositionsAndClearsReferences(t *testing.T) {
	for _, capacity := range []int{maxPooledSeedPositions, maxPooledSeedPositions + 1} {
		r := &SearchResult{Values: make([]uint64, 1, capacity)}
		results := []*SearchResult{r, nil, r}
		backing := results
		results = results[:1]
		RecycleSearchResults(&results)
		for _, stale := range backing {
			if stale != nil {
				t.Fatal("recycled slice retains a seed result")
			}
		}
		if len(r.Values) != 0 {
			t.Fatal("seed positions were not reset")
		}
		if capacity > maxPooledSeedPositions {
			if r.Values != nil {
				t.Fatal("oversized seed-position array was retained")
			}
		} else if cap(r.Values) != capacity {
			t.Fatal("small seed-position array was not kept for reuse")
		}
	}
}

func TestRecycleSearchResultsDropsShortenedOversizedSlice(t *testing.T) {
	results := make([]*SearchResult, 1, maxPooledSearchResults+1)
	results[0] = &SearchResult{Values: make([]uint64, 1, maxPooledSeedPositions+1)}
	r := results[0]
	RecycleSearchResults(&results)
	if results != nil || r.Values != nil {
		t.Fatal("oversized seed arrays were retained after recycling")
	}
}
