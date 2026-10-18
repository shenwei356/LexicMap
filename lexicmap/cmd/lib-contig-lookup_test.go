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
	"math/rand"
	"testing"
)

// referenceContigLookup preserves the original inlined scan as an oracle.
func referenceContigLookup(sizes []int, interval, begin, end, k int) (int, int, int) {
	start, last := 0, 0
	for j, length := range sizes {
		last += length - 1
		if begin+k >= start && end-k <= last {
			return j, start, last
		} else if end < start {
			break
		}
		last += interval + 1
		start = last
	}
	return -1, 0, 0
}

func TestContigLookupMatchesOriginal(t *testing.T) {
	rng := rand.New(rand.NewSource(47))
	for _, n := range []int{0, 1, 2, 31, 32, 33, 256, 1024} {
		for _, interval := range []int{0, 1, 31, 1000} {
			sizes := make([]int, n)
			total := 0
			for i := range sizes {
				sizes[i] = 1 + rng.Intn(400)
				total += sizes[i] + interval
			}
			for _, k := range []int{0, 1, 17, 31} {
				var lookup contigLookup
				check := func(begin, end int) {
					t.Helper()
					want, wb, we := referenceContigLookup(sizes, interval, begin, end, k)
					got, gb, ge := lookup.locate(sizes, interval, begin, end, k)
					if got != want || got >= 0 && (gb != wb || ge != we) {
						t.Fatalf("n=%d interval=%d k=%d query=%d..%d: got %d [%d,%d], want %d [%d,%d]", n, interval, k, begin, end, got, gb, ge, want, wb, we)
					}
				}
				start := 0
				for _, size := range sizes {
					for _, delta := range []int{-k - 1, -k, -1, 0, 1, k, k + 1} {
						check(start+delta, start+delta)
						check(start+delta, start+size-1-delta)
						check(start+size-1+delta, start+size-1+delta)
					}
					start += size + interval
				}
				for range 2000 {
					begin := rng.Intn(total+200) - 100
					check(begin, begin+rng.Intn(800))
				}
			}
		}
	}
}

func TestContigLookupFallbackAndFirstTie(t *testing.T) {
	for _, interval := range []int{-1, 0, 20} {
		for _, size := range []int{0, -1, 1, int(^uint(0) >> 1)} {
			sizes := make([]int, 40)
			for i := range sizes {
				sizes[i] = size
			}
			var lookup contigLookup
			for _, bounds := range [][2]int{{0, 0}, {10, 10}, {70, 100}, {-10, -5}} {
				want, wb, we := referenceContigLookup(sizes, interval, bounds[0], bounds[1], 31)
				got, gb, ge := lookup.locate(sizes, interval, bounds[0], bounds[1], 31)
				if got != want || got >= 0 && (gb != wb || ge != we) {
					t.Fatalf("size=%d interval=%d: got %d [%d,%d], want %d [%d,%d]", size, interval, got, gb, ge, want, wb, we)
				}
			}
		}
	}
	sizes := make([]int, 64)
	for i := range sizes {
		sizes[i] = 1
	}
	var lookup contigLookup
	if got, _, _ := lookup.locate(sizes, 0, 40, 40, 31); got != 9 {
		t.Fatalf("K tolerance must retain the first qualifying contig: got %d, want 9", got)
	}
}

func BenchmarkContigLookup(b *testing.B) {
	for _, n := range []int{4, 32, 256, 4096} {
		sizes := make([]int, n)
		for i := range sizes {
			sizes[i] = 1000
		}
		for _, indexed := range []bool{false, true} {
			b.Run(fmt.Sprintf("contigs=%d/indexed=%t", n, indexed), func(b *testing.B) {
				var lookup contigLookup
				lookup.locate(sizes, 1000, 0, 100, 31)
				lookup.locate(sizes, 1000, 0, 100, 31)
				b.ReportAllocs()
				b.ResetTimer()
				for i := range b.N {
					pos := (i%n)*2000 + 100
					var got int
					if indexed {
						got, _, _ = lookup.locate(sizes, 1000, pos, pos+500, 31)
					} else {
						got, _, _ = referenceContigLookup(sizes, 1000, pos, pos+500, 31)
					}
					if got != i%n {
						b.Fatalf("unexpected contig %d", got)
					}
				}
			})
		}
	}
}
