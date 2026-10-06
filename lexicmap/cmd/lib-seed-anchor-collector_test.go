package cmd

import (
	"sync"
	"testing"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
)

func TestSeedAnchorCollectorCollectsEachGenomeOnce(t *testing.T) {
	idx := &Index{
		info: &IndexInfo{GenomeBatches: 3},
		opt:  &IndexSearchingOptions{NumCPUs: 4},
	}
	collector := newSeedAnchorCollector(idx, 2)

	const nProducers = 4
	const nPerProducer = 100
	genomes := []uint64{
		7,
		uint64(1)<<BITS_GENOME_IDX | 11,
		uint64(1)<<BITS_GENOME_IDX | 2049,
		uint64(2)<<BITS_GENOME_IDX | 13,
	}
	var wg sync.WaitGroup
	for producer := range nProducers {
		wg.Add(1)
		go func(producer int) {
			defer wg.Done()
			buffers := collector.newBuffers()
			for i := range nPerProducer {
				for _, genome := range genomes {
					collector.add(buffers, seedAnchor{
						batchGenomeIndex: genome,
						qBegin:           int32(producer*nPerProducer + i),
						tBegin:           int32(i),
						length:           21,
						qrc:              producer&1 != 0,
						trc:              i&1 != 0,
					})
				}
			}
			collector.flush(buffers)
		}(producer)
	}
	wg.Wait()
	collector.finish()

	if got, want := collector.resultCount(), len(genomes); got != want {
		t.Fatalf("genome entries: got %d, want %d", got, want)
	}
	if got, want := collector.anchorCount(), uint64(nProducers*nPerProducer*len(genomes)); got != want {
		t.Fatalf("anchors: got %d, want %d", got, want)
	}

	seen := make(map[uint64]bool, len(genomes))
	for _, results := range collector.results {
		for _, r := range results {
			if seen[r.BatchGenomeIndex] {
				t.Fatalf("duplicate genome entry: %d", r.BatchGenomeIndex)
			}
			seen[r.BatchGenomeIndex] = true
			if got, want := len(*r.Subs), nProducers*nPerProducer; got != want {
				t.Fatalf("anchors for genome %d: got %d, want %d", r.BatchGenomeIndex, got, want)
			}
			idx.RecycleSearchResult(r)
		}
	}
	collector.clearResults()
}

func TestCollectSeedAnchorsSerial(t *testing.T) {
	idx := &Index{k: 31}
	locses := [][]int{{5 << BITS_STRAND}}
	reverseLocses := make([]*[]int, 1)
	genome := uint64(7)
	srs := []*kv.SearchResult{{
		IQuery: 0,
		Len:    21,
		Values: []uint64{
			genome<<BITS_NONE_IDX | uint64(13)<<BITS_FLAGS,
			genome<<BITS_NONE_IDX | uint64(29)<<BITS_FLAGS,
		},
	}}
	ch := make(chan *[]*kv.SearchResult, 1)
	ch <- &srs
	close(ch)

	results := collectSeedAnchorsSerial(idx, ch, &locses, &reverseLocses, nil)
	if got, want := len(results), 1; got != want {
		t.Fatalf("genome entries: got %d, want %d", got, want)
	}
	if got, want := len(*results[0].Subs), 2; got != want {
		t.Fatalf("anchors: got %d, want %d", got, want)
	}
	if got, want := (*results[0].Subs)[0].QBegin, int32(5); got != want {
		t.Fatalf("query begin: got %d, want %d", got, want)
	}
	if got, want := (*results[0].Subs)[1].TBegin, int32(29); got != want {
		t.Fatalf("target begin: got %d, want %d", got, want)
	}

	idx.RecycleSearchResult(results[0])
}
