package cmd

import (
	"reflect"
	"runtime"
	"sync"
	"testing"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
)

// Keep anchors live until the heap measurement so this covers accumulation,
// before deduplication or chaining can release them.
func BenchmarkSeedAnchorCollectorAccumulation(b *testing.B) {
	const nAnchors = 1 << 20
	const nGenomes = 64
	idx := &Index{
		info: &IndexInfo{GenomeBatches: 1},
		opt:  &IndexSearchingOptions{NumCPUs: 1},
	}
	b.ReportAllocs()
	var liveBytes uint64
	for range b.N {
		b.StopTimer()
		runtime.GC()
		var before, after runtime.MemStats
		runtime.ReadMemStats(&before)
		b.StartTimer()
		collector := newSeedAnchorCollector(idx, 1)
		buffers := collector.newBuffers()
		for i := range nAnchors {
			collector.add(buffers, seedAnchor{
				batchGenomeIndex: uint64(i % nGenomes),
				qBegin:           int32(i / nGenomes),
				tBegin:           int32(i),
				length:           31,
				qrc:              i&1 != 0,
				trc:              i&2 != 0,
			})
		}
		collector.flush(buffers)
		collector.finish()
		b.StopTimer()
		if collector.anchorCount() != nAnchors || collector.resultCount() != nGenomes {
			b.Fatal("collector lost anchors or genomes")
		}
		runtime.GC()
		runtime.ReadMemStats(&after)
		if after.HeapAlloc > before.HeapAlloc {
			liveBytes += after.HeapAlloc - before.HeapAlloc
		}
		collector.clearResults()
		for i := range collector.arenas {
			collector.arenas[i].recycle()
		}
		b.StartTimer()
	}
	b.ReportMetric(float64(liveBytes)/float64(b.N*nAnchors), "live-B/anchor")
}

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
			if got, want := len(r.Subs), nProducers*nPerProducer; got != want {
				t.Fatalf("anchors for genome %d: got %d, want %d", r.BatchGenomeIndex, got, want)
			}
			seenAnchors := make(map[int32]bool, nProducers*nPerProducer)
			for _, sub := range r.Subs {
				producer, i := int(sub.QBegin)/nPerProducer, int(sub.QBegin)%nPerProducer
				want := SubstrPair{
					QBegin: sub.QBegin, TBegin: int32(i), Len: 21,
					QRC: producer&1 != 0, TRC: i&1 != 0,
				}
				if producer < 0 || producer >= nProducers || sub != want || seenAnchors[sub.QBegin] {
					t.Fatalf("genome %d: corrupted or duplicate anchor: %+v", r.BatchGenomeIndex, sub)
				}
				seenAnchors[sub.QBegin] = true
			}

		}
	}
	collector.clearResults()
	for i := range collector.arenas {
		collector.arenas[i].recycle()
	}
}

func TestCollectSeedAnchorsSerial(t *testing.T) {
	for _, tt := range []struct {
		name              string
		qrc, trc, reverse bool
		qBegin, tBegin    int32
	}{
		{"forward/forward/prefix", false, false, false, 5, 13},
		{"forward/forward/suffix", false, false, true, 15, 23},
		{"forward/reverse/prefix", false, true, false, 5, 23},
		{"forward/reverse/suffix", false, true, true, 15, 13},
		{"reverse/forward/prefix", true, false, false, 15, 13},
		{"reverse/forward/suffix", true, false, true, 5, 23},
		{"reverse/reverse/prefix", true, true, false, 15, 23},
		{"reverse/reverse/suffix", true, true, true, 5, 13},
	} {
		t.Run(tt.name, func(t *testing.T) {
			idx := &Index{k: 31}
			qpos := 5 << BITS_STRAND
			var flags uint64
			if tt.qrc {
				qpos |= 1
			}
			if tt.trc {
				flags |= 1 << BITS_REVERSE
			}
			if tt.reverse {
				flags |= 1
			}
			locses := [][]int{{qpos}}
			reverseQueryIndices := []int{0}
			reverseLocses := [][]int{reverseQueryIndices}
			genome := uint64(7)
			srs := []kv.SearchResult{{
				IQuery: 0, IsSuffix: tt.reverse, IQuery2: 0, Len: 21,
				Values: []uint64{
					genome<<BITS_NONE_IDX | uint64(13)<<BITS_FLAGS | flags,
					genome<<BITS_NONE_IDX | uint64(29)<<BITS_FLAGS | flags,
				},
			}}
			ch := make(chan *[]kv.SearchResult, 1)
			ch <- &srs
			close(ch)

			var arena seedSearchResultArena
			results := collectSeedAnchorsSerial(idx, ch, &locses, &reverseLocses, nil, &arena)
			if len(results) != 1 || results[0].BatchGenomeIndex != genome {
				t.Fatalf("genome entries: got %+v, want genome %d", results, genome)
			}
			want := []SubstrPair{
				{QBegin: tt.qBegin, TBegin: tt.tBegin, Len: 21, QRC: tt.qrc, TRC: tt.trc},
				{QBegin: tt.qBegin, TBegin: tt.tBegin + 16, Len: 21, QRC: tt.qrc, TRC: tt.trc},
			}
			if !reflect.DeepEqual(results[0].Subs, want) {
				t.Fatalf("anchors: got %+v, want %+v", results[0].Subs, want)
			}
			arena.recycle()
		})
	}
}
