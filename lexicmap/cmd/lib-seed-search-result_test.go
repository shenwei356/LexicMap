package cmd

import (
	"slices"
	"testing"
)

func TestSeedSearchResultPagesKeepAddressesStable(t *testing.T) {
	var arena seedSearchResultArena
	defer arena.recycle()
	first := arena.add(7)
	first.Subs = append(first.Subs, SubstrPair{QBegin: 11, TBegin: 17, Len: 31})
	if &first.Subs[0] != &first.inlineSub[0] {
		t.Fatal("single anchor did not use inline storage")
	}
	for i := range seedSearchResultPageSize*2 + 3 {
		r := arena.add(uint64(i + 100))
		r.Subs = append(r.Subs, SubstrPair{QBegin: int32(i), Len: 31})
	}
	if first != &arena.pages[0][0] || first.BatchGenomeIndex != 7 ||
		first.Subs[0] != (SubstrPair{QBegin: 11, TBegin: 17, Len: 31}) {
		t.Fatal("growing the arena moved or overwrote an active candidate")
	}
}

func TestSeedSearchResultPromotionOwnsRegions(t *testing.T) {
	for _, n := range []int{1, 3} {
		var arena seedSearchResultArena
		candidate := arena.add(uint64(2)<<BITS_GENOME_IDX | 17)
		candidate.Score = 123.5
		for i := range n {
			candidate.chainRegions = append(candidate.chainRegions, seedChainRegion{qBegin: int32(i), tBegin: 100, tEnd: 131, rc: i%2 == 0})
		}
		want := slices.Clone(candidate.chainRegions)
		result := candidate.alignmentResult()
		if result.BatchGenomeIndex != candidate.BatchGenomeIndex || result.GenomeBatch != 2 ||
			result.GenomeIndex != 17 || result.Score != candidate.Score || !slices.Equal(result.chainRegions, want) {
			t.Fatal("promotion changed candidate metadata or chain bounds")
		}
		if &result.chainRegions[0] == &candidate.chainRegions[0] {
			t.Fatal("alignment result still borrows a candidate region buffer")
		}
		clear(candidate.chainRegions)
		arena.recycle()
		if !slices.Equal(result.chainRegions, want) {
			t.Fatal("recycling candidate pages changed a pending alignment")
		}
		(&Index{}).RecycleSearchResult(result)
	}
}

func TestSeedSearchResultPagesDropOversizedUnusedBuffers(t *testing.T) {
	var arena seedSearchResultArena
	r := arena.add(7)
	r.Subs = make([]SubstrPair, 1, thresholdNSubs+1)
	r.chainRegions = make([]seedChainRegion, 1, thresholdNSubs+1)
	unused := &arena.pages[0][seedSearchResultPageSize-1]
	unused.Subs = make([]SubstrPair, 1, thresholdNSubs+1)
	unused.chainRegions = make([]seedChainRegion, 1, thresholdNSubs+1)
	arena.recycle()
	for _, result := range []*seedSearchResult{r, unused} {
		if len(result.Subs) != 0 || cap(result.Subs) != 1 || len(result.chainRegions) != 0 || cap(result.chainRegions) != 1 {
			t.Fatal("recycling a page retained oversized buffers")
		}
	}
}
