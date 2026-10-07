package tree

import (
	"fmt"
	"math/rand"
	"reflect"
	"slices"
	"testing"
)

func TestSearchEachMatchesLeaves(t *testing.T) {
	for _, k := range []uint8{2, 6, 31, 32} {
		for _, n := range []int{0, 1, 96} {
			t.Run(fmt.Sprintf("k=%d/leaves=%d", k, n), func(t *testing.T) {
				tree := NewTree(k)
				defer RecycleTree(tree)
				rng := rand.New(rand.NewSource(42))
				mask := ^uint64(0) >> (64 - 2*k)
				leaves := make(map[uint64][]uint32)
				queries := []uint64{0, 1, mask, mask >> 1}
				for i := range n {
					key := rng.Uint64() & mask
					// Preserve the order of multiple values on the same leaf.
					for _, v := range []uint32{uint32(i*2 + 1), uint32(i*2 + 2)} {
						tree.Insert(key, v)
						leaves[key] = append(leaves[key], v)
					}
					queries = append(queries, key)
				}
				for range 16 {
					queries = append(queries, rng.Uint64()&mask)
				}
				keys := make([]uint64, 0, len(leaves))
				for key := range leaves {
					keys = append(keys, key)
				}
				slices.Sort(keys)
				for _, p := range []uint8{0, 1, 2, k / 2, k - 1, k, k + 1, 255} {
					for _, query := range queries {
						var want []SearchResult
						for _, key := range keys {
							var prefix uint8
							for i := k; i > 0; i-- {
								shift := (i - 1) * 2
								if (query>>shift)&3 != (key>>shift)&3 {
									break
								}
								prefix++
							}
							if prefix >= min(max(p, 1), k) {
								want = append(want, SearchResult{Kmer: key, LenPrefix: prefix, Values: leaves[key]})
							}
						}
						var got []SearchResult
						found := tree.SearchEach(query, p, func(key uint64, prefix uint8, values []uint32) {
							got = append(got, SearchResult{Kmer: key, LenPrefix: prefix, Values: slices.Clone(values)})
						})
						if found != (len(want) > 0) || !reflect.DeepEqual(got, want) {
							t.Fatalf("query=%x p=%d: found=%v got=%v want=%v", query, p, found, got, want)
						}
						// The allocating API must retain the same order and borrowed data.
						results, found := tree.Search(query, p)
						got = nil
						if found {
							for _, r := range *results {
								got = append(got, *r)
							}
							if !reflect.DeepEqual(got, want) {
								t.Fatalf("Search query=%x p=%d: got=%v want=%v", query, p, got, want)
							}
							tree.RecycleSearchResult(results)
						} else if len(want) != 0 {
							t.Fatal("Search missed matching leaves")
						}
					}
				}
			})
		}
	}
}

func BenchmarkSearchEach(b *testing.B) {
	for _, n := range []int{1, 8, 64} {
		b.Run(fmt.Sprintf("leaves=%d", n), func(b *testing.B) {
			tree := NewTree(31)
			defer RecycleTree(tree)
			for i := range n {
				tree.Insert(uint64(i), uint32(i))
			}
			b.Run("results", func(b *testing.B) {
				results, _ := tree.Search(0, 15)
				tree.RecycleSearchResult(results)
				count := 0
				b.ReportAllocs()
				b.ResetTimer()
				for b.Loop() {
					results, _ := tree.Search(0, 15)
					for _, r := range *results {
						count += len(r.Values)
					}
					tree.RecycleSearchResult(results)
				}
				if count != b.N*n {
					b.Fatal("incorrect number of leaf values")
				}
			})
			b.Run("visitor", func(b *testing.B) {
				count := 0
				visit := func(_ uint64, _ uint8, values []uint32) { count += len(values) }
				b.ReportAllocs()
				b.ResetTimer()
				for b.Loop() {
					tree.SearchEach(0, 15, visit)
				}
				if count != b.N*n {
					b.Fatal("incorrect number of leaf values")
				}
			})
		})
	}
}
