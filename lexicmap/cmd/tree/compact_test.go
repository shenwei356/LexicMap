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

package tree

import (
	"fmt"
	"math/rand"
	"reflect"
	"slices"
	"testing"
)

func TestCompactTreeMatchesMutable(t *testing.T) {
	for _, k := range []uint8{1, 2, 6, 15, 31, 32} {
		for _, n := range []int{0, 1, 96, 256, 1024, 4096} {
			t.Run(fmt.Sprintf("k=%d/n=%d", k, n), func(t *testing.T) {
				source := NewTree(k)
				rng := rand.New(rand.NewSource(62))
				mask := ^uint64(0) >> (64 - 2*k)
				queries := []uint64{0, 1, mask, mask >> 1}
				entries := make([]BatchEntry, 0, n*2)
				for i := 0; i < n; i++ {
					key := rng.Uint64() & mask
					entries = append(entries, BatchEntry{key, uint32(i * 2)}, BatchEntry{key, uint32(i*2 + 1)})
					if i < 100 {
						queries = append(queries, key)
					}
				}
				source.InsertBatch(entries)
				for i := 0; i < 100; i++ {
					queries = append(queries, rng.Uint64()&mask)
				}
				for _, directory := range []bool{false, true} {
					if directory {
						source.BuildPrefixDirectory()
					}
					frozen := NewCompactTree(k, entries)
					if got, want := frozen.RetainedBytes(), CompactTreeBytes(k, entries); got != want {
						t.Fatalf("charge=%d want=%d", got, want)
					}
					for p := uint8(0); p <= k+1; p++ {
						for _, query := range queries {
							var want, got []SearchResult
							found := source.SearchEach(query, p, func(key uint64, p uint8, v []uint32) { want = append(want, SearchResult{key, p, slices.Clone(v)}) })
							found2 := frozen.SearchEach(query, p, func(key uint64, p uint8, v []uint32) { got = append(got, SearchResult{key, p, slices.Clone(v)}) })
							if found != found2 || !reflect.DeepEqual(got, want) {
								t.Fatalf("directory=%v query=%x p=%d got=%v want=%v", directory, query, p, got, want)
							}
						}
					}
				}
				// Copy ownership must survive reused entry scratch and pooled source arenas.
				frozen := NewCompactTree(k, entries)
				var want, got []BatchEntry
				frozen.Walk(func(key uint64, v []uint32) bool {
					for _, value := range v {
						want = append(want, BatchEntry{key, value})
					}
					return false
				})
				clear(entries)
				source.Insert(mask, 99999)
				RecycleTree(source)
				replacement := NewTree(k)
				replacement.Insert(0, 777)
				frozen.Walk(func(key uint64, v []uint32) bool {
					for _, value := range v {
						got = append(got, BatchEntry{key, value})
					}
					return false
				})
				RecycleTree(replacement)
				if !reflect.DeepEqual(got, want) {
					t.Fatal("mutable source changed compact copy")
				}
			})
		}
	}
}

// Preserve the original short-edge behavior even when the directory misses.
func TestCompactTreeLegacyShortEdge(t *testing.T) {
	for _, n := range []int{2, 256} {
		source := NewTree(6)
		entries := make([]BatchEntry, n)
		for i := range entries {
			entries[i] = BatchEntry{256 + uint64(i), uint32(i)}
		} // ACxxxx
		source.InsertSortedBatch(entries)
		source.BuildPrefixDirectory()
		compact := NewCompactTree(6, entries)
		var want, got []SearchResult
		source.SearchEach(0, 6, func(key uint64, p uint8, v []uint32) { want = append(want, SearchResult{key, p, slices.Clone(v)}) })
		compact.SearchEach(0, 6, func(key uint64, p uint8, v []uint32) { got = append(got, SearchResult{key, p, slices.Clone(v)}) })
		RecycleTree(source)
		if len(want) != n || !reflect.DeepEqual(got, want) {
			t.Fatalf("n=%d: got=%v want=%v", n, got, want)
		}
	}
}

func BenchmarkCompactLookup(b *testing.B) {
	source := NewTree(31)
	defer RecycleTree(source)
	rng := rand.New(rand.NewSource(71))
	queries := make([]uint64, 1980)
	entries := make([]BatchEntry, 0, len(queries))
	for i := range queries {
		queries[i] = rng.Uint64() >> 2
		entries = append(entries, BatchEntry{queries[i], uint32(i)})
	}
	source.InsertBatch(entries)
	source.BuildPrefixDirectory()
	frozen := NewCompactTree(31, entries)
	missing := make([]uint64, 0, len(queries))
	for len(missing) < len(queries) {
		key := rng.Uint64() >> 2
		if !source.SearchEach(key, 11, func(uint64, uint8, []uint32) {}) {
			missing = append(missing, key)
		}
	}
	for _, sample := range []struct {
		name string
		keys []uint64
	}{{"hits", queries}, {"misses", missing}} {
		b.Run(sample.name, func(b *testing.B) {
			expected := 0
			for _, key := range sample.keys {
				source.SearchEach(key, 11, func(_ uint64, _ uint8, v []uint32) { expected += len(v) })
			}
			for _, variant := range []string{"mutable", "compact"} {
				b.Run(variant, func(b *testing.B) {
					count := 0
					visit := func(_ uint64, _ uint8, v []uint32) { count += len(v) }
					b.ReportAllocs()
					for b.Loop() {
						for _, query := range sample.keys {
							if variant == "mutable" {
								source.SearchEach(query, 11, visit)
							} else {
								frozen.SearchEach(query, 11, visit)
							}
						}
					}
					if count != b.N*expected {
						b.Fatalf("got %d values; want %d", count, b.N*expected)
					}
				})
			}
		})
	}
}
