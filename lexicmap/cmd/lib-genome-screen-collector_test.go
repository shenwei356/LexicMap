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
	"fmt"
	"reflect"
	"sync"
	"testing"
)

// GSearchScreen can choose a different mask subset on each call to one index.
func TestGSearchScreenCollectorReusesDifferentMaskCounts(t *testing.T) {
	for _, counts := range [][]int{{8, 1}, {1, 8}, {8, 3, 8}} {
		t.Run(fmt.Sprint(counts), func(t *testing.T) {
			idx := newGSearchScreenCollectorTestIndex(1)
			reused := &GSearchScreenResultDetail{}
			for _, count := range counts {
				// Guarantee reuse even if sync.Pool would discard its contents.
				// Each completed collector owns exactly one result at a time.
				idx.poolGSearchDetailResult = &sync.Pool{New: func() any { return reused }}
				collector := newGSearchScreenCollector(idx, count, 1, nil)
				buffers := collector.newBuffers()
				var wantSum uint64
				for mask := range count {
					length := uint8(20 + mask)
					collector.add(buffers, gsearchScreenHit{batchGenomeIndex: 7, mask: int32(mask), length: length})
					wantSum += uint64(length)
				}
				collector.flush(buffers)
				collector.finish()
				if len(collector.results[0]) != 1 {
					t.Fatalf("got %d results, want one", len(collector.results[0]))
				}
				r := collector.results[0][0]
				if len(r.LongestMatches) != count || r.SumPrefix != wantSum || r.HitMasks != count {
					t.Fatalf("%d masks: got %d match slots, sum %d; want %d, %d", count, len(r.LongestMatches), r.SumPrefix, count, wantSum)
				}
				for mask, length := range r.LongestMatches {
					if length != uint8(20+mask) {
						t.Fatalf("mask %d: got prefix length %d", mask, length)
					}
				}
				idx.RecycleGSearchScreenDetailResult(r)
				clearGSearchScreenResults(collector.results)
			}
		})
	}
}

// Ordered readers may interleave at any mask. Repeated hits must still count
// only the longest prefix, and split genomes must retain their merge details.
func TestGSearchScreenCollectorRollingScores(t *testing.T) {
	const nMasks, nSources = 64, 4
	sources := make([]int32, nMasks)
	for mask := range sources {
		sources[mask] = int32(mask / (nMasks / nSources))
	}
	genomes := []uint64{0, 7, 2049, 1<<BITS_GENOME_IDX | 3, 2<<BITS_GENOME_IDX | 13}
	streams := make([][]gsearchScreenHit, nSources)
	want := make(map[uint64][]uint8)
	for _, genome := range genomes {
		want[genome] = make([]uint8, nMasks)
	}
	for mask, source := range sources {
		for repeat := range 130 { // cross collector batch boundaries within a mask
			for i, genome := range genomes {
				if (mask+i)%7 == 0 { // a genome can skip masks in any reader
					continue
				}
				length := uint8(15 + (mask+repeat+i)%18) // includes k=32
				streams[source] = append(streams[source], gsearchScreenHit{genome, int32(mask), length})
				want[genome][mask] = max(want[genome][mask], length)
			}
		}
	}
	for _, workers := range []int{1, 3} {
		t.Run(fmt.Sprint(workers), func(t *testing.T) {
			idx := newGSearchScreenCollectorTestIndex(3)
			idx.genomeChunksIdx = map[uint64]genomeChunkLocation{7: {}, genomes[3]: {}}
			for _, slots := range [][]int32{nil, sources, sources, nil} {
				collector := newGSearchScreenCollector(idx, nMasks, workers, slots)
				var wg sync.WaitGroup
				var count uint64
				for _, hits := range streams {
					count += uint64(len(hits))
					wg.Add(1)
					go func() {
						defer wg.Done()
						buffers := collector.newBuffers()
						for _, hit := range hits {
							collector.add(buffers, hit)
						}
						collector.flush(buffers)
					}()
				}
				wg.Wait()
				collector.finish()
				if collector.matchCount() != count {
					t.Fatalf("match count: got %d, want %d", collector.matchCount(), count)
				}
				seen := make(map[uint64]bool)
				for _, batch := range collector.results {
					for _, r := range batch {
						genome := r.BatchGenomeIndex[0]
						if seen[genome] {
							t.Fatalf("duplicate genome %d", genome)
						}
						seen[genome] = true
						var sum uint64
						var hitMasks int
						for _, length := range want[genome] {
							sum += uint64(length)
							if length > 0 {
								hitMasks++
							}
						}
						if r.SumPrefix != sum {
							t.Fatalf("genome %d score: got %d, want %d", genome, r.SumPrefix, sum)
						}
						if r.HitMasks != hitMasks {
							t.Fatalf("genome %d matched masks: got %d, want %d", genome, r.HitMasks, hitMasks)
						}
						_, split := idx.genomeChunksIdx[genome]
						if slots == nil || split {
							if !reflect.DeepEqual(r.LongestMatches, want[genome]) || len(r.lastMatches) != 0 {
								t.Fatalf("lost mask details for genome %d", genome)
							}
						} else if len(r.LongestMatches) != 0 || len(r.lastMatches) != nSources {
							t.Fatalf("genome %d did not use rolling maxima", genome)
						}
						idx.RecycleGSearchScreenDetailResult(r)
					}
				}
				if len(seen) != len(want) {
					t.Fatalf("got %d genomes, want %d", len(seen), len(want))
				}
				clearGSearchScreenResults(collector.results)
			}
		})
	}
}

func TestGSearchScreenCollectorReusesRollingStorage(t *testing.T) {
	idx := newGSearchScreenCollectorTestIndex(1)
	reused := &GSearchScreenResultDetail{}
	for _, nSources := range []int{0, 4, 1, 8, 0, 2} {
		idx.poolGSearchDetailResult = &sync.Pool{New: func() any { return reused }}
		var sources []int32
		if nSources > 0 {
			sources = make([]int32, 64)
			for mask := range sources {
				sources[mask] = int32(mask / (64 / nSources))
			}
		}
		collector := newGSearchScreenCollector(idx, 64, 1, sources)
		buffers := collector.newBuffers()
		collector.add(buffers, gsearchScreenHit{7, 0, 20})
		collector.add(buffers, gsearchScreenHit{7, 63, 31})
		collector.flush(buffers)
		collector.finish()
		r := collector.results[0][0]
		if r.SumPrefix != 51 {
			t.Fatalf("%d readers: got score %d, want 51", nSources, r.SumPrefix)
		}
		if r.HitMasks != 2 {
			t.Fatalf("%d readers: got %d matched masks, want 2", nSources, r.HitMasks)
		}
		idx.RecycleGSearchScreenDetailResult(r)
		if r.HitMasks != 0 {
			t.Fatal("recycled matched-mask count was not cleared")
		}
		for _, v := range r.lastMatches {
			if v != 0 {
				t.Fatal("recycled rolling maximum was not cleared")
			}
		}
		clearGSearchScreenResults(collector.results)
	}
}

func newGSearchScreenCollectorTestIndex(genomeBatches int) *Index {
	idx := &Index{info: &IndexInfo{GenomeBatches: genomeBatches}}
	idx.poolGSearchDetailResult = &sync.Pool{New: func() any {
		return &GSearchScreenResultDetail{}
	}}
	return idx
}

func TestGSearchScreenCollectorCollectsEachGenomeOnce(t *testing.T) {
	idx := newGSearchScreenCollectorTestIndex(3)
	collector := newGSearchScreenCollector(idx, 4, 2, nil)

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
