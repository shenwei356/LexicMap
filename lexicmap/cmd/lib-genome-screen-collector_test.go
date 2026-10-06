package cmd

import (
	"sync"
	"testing"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
)

func newGSearchScreenCollectorTestIndex(genomeBatches int) *Index {
	idx := &Index{info: &IndexInfo{GenomeBatches: genomeBatches}}
	idx.poolGSearchDetailResult = &sync.Pool{New: func() any {
		return &GSearchScreenResultDetail{}
	}}
	idx.poolGSearchDetailResultsMap = &sync.Pool{New: func() any {
		m := make(map[uint64]*GSearchScreenResultDetail)
		return &m
	}}
	return idx
}

func TestGSearchScreenCollectorCollectsEachGenomeOnce(t *testing.T) {
	idx := newGSearchScreenCollectorTestIndex(3)
	collector := newGSearchScreenCollector(idx, 4, 2)

	genomes := []uint64{
		7,
		uint64(1)<<BITS_GENOME_IDX | 11,
		uint64(1)<<BITS_GENOME_IDX | 2049,
		uint64(2)<<BITS_GENOME_IDX | 13,
	}
	const nProducers = 4
	var wg sync.WaitGroup
	for producer := range nProducers {
		wg.Add(1)
		go func(producer int) {
			defer wg.Done()
			buffers := collector.newBuffers()
			for _, genome := range genomes {
				for mask := range 4 {
					collector.add(buffers, gsearchScreenHit{
						batchGenomeIndex: genome,
						mask:             int32(mask),
						length:           uint8(20 + producer),
					})
				}
			}
			collector.flush(buffers)
		}(producer)
	}
	wg.Wait()
	collector.finish()

	if got, want := collector.matchCount(), uint64(nProducers*len(genomes)*4); got != want {
		t.Fatalf("screen matches: got %d, want %d", got, want)
	}
	seen := make(map[uint64]bool, len(genomes))
	for _, results := range collector.results {
		for _, r := range results {
			genome := r.BatchGenomeIndex[0]
			if seen[genome] {
				t.Fatalf("duplicate genome entry: %d", genome)
			}
			seen[genome] = true
			if got, want := r.SumPrefix, uint64(4*23); got != want {
				t.Fatalf("sum prefix for genome %d: got %d, want %d", genome, got, want)
			}
			for mask, got := range r.LongestMatches {
				if want := uint8(23); got != want {
					t.Fatalf("longest match for genome %d, mask %d: got %d, want %d", genome, mask, got, want)
				}
			}
			idx.RecycleGSearchScreenDetailResult(r)
		}
	}
	if got, want := len(seen), len(genomes); got != want {
		t.Fatalf("genome entries: got %d, want %d", got, want)
	}
	clearGSearchScreenResults(collector.results)
}

func TestCollectGSearchScreenResultsSerial(t *testing.T) {
	idx := newGSearchScreenCollectorTestIndex(2)
	genomes := []uint64{7, uint64(1)<<BITS_GENOME_IDX | 2049}
	srs := []*kv.SearchResult{
		{IQuery: 0, Len: 21, Values: []uint64{genomes[0] << BITS_NONE_IDX, genomes[1] << BITS_NONE_IDX}},
		{IQuery: 1, Len: 25, Values: []uint64{genomes[0] << BITS_NONE_IDX, genomes[1] << BITS_NONE_IDX}},
	}
	ch := make(chan *[]*kv.SearchResult, 1)
	ch <- &srs
	close(ch)

	results, nMatches := collectGSearchScreenResultsSerial(idx, ch, nil, 2)
	if got, want := nMatches, uint64(4); got != want {
		t.Fatalf("screen matches: got %d, want %d", got, want)
	}
	if got, want := len(results), len(genomes); got != want {
		t.Fatalf("genome entries: got %d, want %d", got, want)
	}
	for _, r := range results {
		if got, want := r.SumPrefix, uint64(46); got != want {
			t.Fatalf("sum prefix for genome %d: got %d, want %d", r.BatchGenomeIndex[0], got, want)
		}
		if got, want := r.LongestMatches, []uint8{21, 25}; got[0] != want[0] || got[1] != want[1] {
			t.Fatalf("longest matches for genome %d: got %v, want %v", r.BatchGenomeIndex[0], got, want)
		}
		idx.RecycleGSearchScreenDetailResult(r)
	}
}
