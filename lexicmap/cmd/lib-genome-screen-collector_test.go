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
	"sync"
	"testing"
)

func newGSearchScreenCollectorTestIndex(genomeBatches int) *Index {
	idx := &Index{info: &IndexInfo{GenomeBatches: genomeBatches}}
	idx.poolGSearchDetailResult = &sync.Pool{New: func() any {
		return &GSearchScreenResultDetail{}
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
