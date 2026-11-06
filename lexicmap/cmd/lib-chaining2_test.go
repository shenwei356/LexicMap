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
	"math"
	"math/rand"
	"reflect"
	"slices"
	"testing"
)

// An anchor whose length cancels its gap penalty still extends the same chain.
// Exercise the initial maximum and recursive regions on both sides of it.
func TestChainer2ZeroScoreExtension(t *testing.T) {
	for _, recursive := range []bool{false, true} {
		for _, reverse := range []bool{false, true} {
			for _, terminalAnchors := range []int{1, 2} {
				t.Run(fmt.Sprintf("recursive=%t/reverse=%t/ends=%d", recursive, reverse, terminalAnchors), func(t *testing.T) {
					subs := []SubstrPair{{QBegin: 200, TBegin: 200, Len: 31}, {QBegin: 300, TBegin: 325, Len: 25}}
					want := [5]int{2, 200, 324, 200, 349}
					if terminalAnchors == 2 {
						subs = append(subs, SubstrPair{QBegin: 400, TBegin: 450, Len: 25})
						want = [5]int{3, 200, 424, 200, 474}
					}
					if recursive {
						// A higher independent score forces the tested chain into a
						// recursive region whose packed predecessor indexes are global.
						subs = append([]SubstrPair{{TBegin: 2000, Len: 41}}, subs...)
					}
					if reverse {
						for i := range subs {
							subs[i].QBegin = 3000 - subs[i].QBegin - int32(subs[i].Len)
							subs[i].TBegin = 3000 - subs[i].TBegin - int32(subs[i].Len)
						}
						slices.Reverse(subs)
						want = [5]int{want[0], 2999 - want[2], 2999 - want[1], 2999 - want[4], 2999 - want[3]}
					}
					opt := Chaining2Options{MaxGap: 100, MinScore: 22, BandBase: 100, BandCount: 50}
					chains, _, _, _, _, _, _, _ := NewChainer2(&opt).Chain(&subs)
					if chains == nil {
						t.Fatal("missing chain")
					}
					defer RecycleChaining2Result(chains)
					wantCount, position := 1, 0
					if recursive {
						wantCount, position = 2, 1
					}
					if len(*chains) != wantCount {
						t.Fatalf("got %d chains, want %d: %+v", len(*chains), wantCount, denseChain2Summary(chains))
					}
					c := (*chains)[position]
					if got := [5]int{c.NAnchors, c.QBegin, c.QEnd, c.TBegin, c.TEnd}; got != want {
						t.Fatalf("chain = %v, want %v", got, want)
					}
				})
			}
		}
	}
}

// An equal score cannot justify shortening either side of the selected chain.
func TestChainer2ZeroScoreExtensionPreservesBounds(t *testing.T) {
	for _, swap := range []bool{false, true} {
		subs := []SubstrPair{{Len: 31}, {QBegin: 10, TBegin: 23, Len: 13}}
		if swap {
			subs[1].QBegin, subs[1].TBegin = subs[1].TBegin, subs[1].QBegin
		}
		opt := Chaining2Options{MaxGap: 100, MinScore: 22, BandBase: 100, BandCount: 50}
		chains, _, _, _, _, _, _, _ := NewChainer2(&opt).Chain(&subs)
		if chains == nil {
			t.Fatal("missing chain")
		}
		got := denseChain2Summary(chains)
		RecycleChaining2Result(chains)
		want := []denseChainSummary{{anchors: 1, matched: 31, alignedQ: 31, alignedT: 31, qe: 30, te: 30, identity: 100}}
		if !reflect.DeepEqual(got, want) {
			t.Fatalf("swap=%t: got %+v, want %+v", swap, got, want)
		}
	}
}

// Each linked anchor contributes its own length. Reversing a chain must not
// replace its shorter terminal anchor with another copy of the longer one.
func TestChainer2VariableLengthScore(t *testing.T) {
	for _, gap := range []int32{0, 3} {
		for _, reverse := range []bool{false, true} {
			for _, extra := range []int{0, 1} {
				t.Run(fmt.Sprintf("gap=%d/reverse=%t/threshold+%d", gap, reverse, extra), func(t *testing.T) {
					opt := DefaultChaining2Options
					wantScore := 13 + 31 - int(gap)
					opt.MinScore, opt.MinAlignLen = wantScore+extra, 1
					subs := []SubstrPair{{QBegin: 20, TBegin: 20, Len: 13}, {QBegin: 60, TBegin: 60 + gap, Len: 31}}
					if reverse {
						for i := range subs {
							subs[i].QBegin = 200 - subs[i].QBegin - int32(subs[i].Len)
							subs[i].TBegin = 200 - subs[i].TBegin - int32(subs[i].Len)
						}
						slices.Reverse(subs)
					}
					chainer := NewChainer2(&opt)
					chains, _, _, _, _, _, _, _ := chainer.Chain(&subs)
					if chains != nil {
						defer RecycleChaining2Result(chains)
					}
					if got := chainer.maxscoresIdxs[1]; got != uint64(wantScore)<<32 {
						t.Errorf("score/predecessor=%d/%d, want %d/0", got>>32, uint32(got), wantScore)
					}
					if (chains != nil) != (extra == 0) {
						t.Fatalf("retained=%t at minimum score %d", chains != nil, opt.MinScore)
					}
					if chains != nil && (len(*chains) != 1 || (*chains)[0].NAnchors != 2) {
						t.Fatalf("want one two-anchor chain, got %+v", denseChain2Summary(chains))
					}
				})
			}
		}
	}
}

// Reuse one chainer across changing group sizes, orders and band limits. Compare
// every packed score/predecessor and chain field, including reverse complements.
func TestChainer2CursorReuseMatchesReference(t *testing.T) {
	rng := rand.New(rand.NewSource(281))
	opt := DefaultChaining2Options
	opt.MinScore, opt.MinAlignLen = 1, 1
	got, ref := NewChainer2(&opt), NewChainer2(&opt)
	for trial := 0; trial < 48; trial++ {
		opt.BandBase = []int{0, 30, 10000}[trial%3]
		opt.BandCount = []int{0, 1, 10, 64}[trial%4]
		opt.MaxGap = []int{0, 5, 200}[trial%3]
		subs := make([]SubstrPair, 0, 2048)
		for q := 0; q < 24; q++ {
			for copy := 0; copy < 9+rng.Intn(80); copy++ {
				length := uint8(13)
				if trial%2 == 1 {
					length += uint8(rng.Intn(19))
				}
				a := SubstrPair{QBegin: int32(q * 25), TBegin: int32(rng.Intn(80) * 25), Len: length}
				if trial%4 >= 2 {
					a.QBegin = 3000 - a.QBegin - int32(a.Len)
					a.TBegin = 3000 - a.TBegin - int32(a.Len)
				}
				subs = append(subs, a)
			}
		}
		slices.SortFunc(subs, func(a, b SubstrPair) int {
			if a.QBegin != b.QBegin {
				return int(a.QBegin) - int(b.QBegin)
			}
			return int(a.TBegin) - int(b.TBegin)
		})
		if trial%8 == 7 {
			rng.Shuffle(len(subs), func(i, j int) { subs[i], subs[j] = subs[j], subs[i] })
		}
		g, gm, gq, gt, gqb, gqe, gtb, gte := got.Chain(&subs)
		r, rm, rq, rt, rqb, rqe, rtb, rte := referenceChaining2(ref, &subs)
		if !slices.Equal(got.maxscoresIdxs, ref.maxscoresIdxs) ||
			[7]int{gm, gq, gt, gqb, gqe, gtb, gte} != [7]int{rm, rq, rt, rqb, rqe, rtb, rte} ||
			!reflect.DeepEqual(denseChain2Summary(g), denseChain2Summary(r)) {
			t.Fatalf("trial %d differs from reference", trial)
		}
		if g != nil {
			RecycleChaining2Result(g)
		}
		if r != nil {
			RecycleChaining2Result(r)
		}
	}
}

// Disjoint chains force repeated maximum searches. A zero-score terminal anchor
// in each group also exercises continuation ties in indexed subregions.
func disjointChain2Anchors(groups int) []SubstrPair {
	subs := make([]SubstrPair, 0, groups*3)
	for i := range groups {
		q, t := int32(i*500), int32(i*1000)
		subs = append(subs, SubstrPair{QBegin: q, TBegin: t, Len: 31},
			SubstrPair{QBegin: q + 31, TBegin: t + 31, Len: 31},
			SubstrPair{QBegin: q + 100, TBegin: t + 125, Len: 25})
	}
	return subs
}

func TestChainer2DisjointRegionsMatchReference(t *testing.T) {
	for _, reverse := range []bool{false, true} {
		subs := disjointChain2Anchors(2000)
		if reverse {
			for i := range subs {
				subs[i].QBegin = 2000000 - subs[i].QBegin - int32(subs[i].Len)
				subs[i].TBegin = 3000000 - subs[i].TBegin - int32(subs[i].Len)
			}
			slices.Reverse(subs)
		}
		opt := DefaultChaining2Options
		opt.MinScore = 38
		got, ref := NewChainer2(&opt), NewChainer2(&opt)
		g, gm, gq, gt, gqb, gqe, gtb, gte := got.Chain(&subs)
		r, rm, rq, rt, rqb, rqe, rtb, rte := referenceChaining2(ref, &subs)
		if g == nil || r == nil {
			t.Fatal("missing disjoint chains")
		}
		if len(*g) != 2000 || !reflect.DeepEqual(denseChain2Summary(g), denseChain2Summary(r)) ||
			[7]int{gm, gq, gt, gqb, gqe, gtb, gte} != [7]int{rm, rq, rt, rqb, rqe, rtb, rte} {
			t.Fatalf("reverse=%t: independent chains changed", reverse)
		}
		RecycleChaining2Result(g)
		RecycleChaining2Result(r)
	}
}

// Partial blocks and same-score branches must select exactly the endpoint of
// the original left-to-right scan, including nonzero global predecessor offsets.
func TestChain2PeakIndexMatchesScan(t *testing.T) {
	rng := rand.New(rand.NewSource(1906))
	for _, n := range []int{63, 64, 65, 1024, 2049} {
		subs := make([]SubstrPair, n)
		scores := make([]uint64, n)
		const offset = 17
		for i := range subs {
			subs[i] = SubstrPair{QBegin: int32(i * 3), TBegin: int32(i * 5), Len: uint8(1 + rng.Intn(31))}
			previous := i
			if i > 0 && rng.Intn(3) != 0 {
				previous = rng.Intn(i)
			}
			scores[i] = uint64(rng.Intn(5))<<32 | uint64(offset+previous)
		}
		index := chain2PeakIndex{subs: subs, scores: scores, offset: offset}
		index.best(0, n)
		if index.tree != nil {
			t.Fatal("a single scan allocated the optional index")
		}
		if n < chain2PeakBlockSize {
			index.build() // Exercise partial blocks smaller than the activation limit.
		} else {
			for range 5 {
				index.best(0, n)
			}
			if index.tree == nil || index.scanned > 5*n {
				t.Fatal("repeated scans did not activate the index")
			}
		}
		for range 2000 {
			begin := rng.Intn(n)
			end := begin + 1 + rng.Intn(n-begin)
			want, score := begin, uint64(0)
			for i := begin; i < end; i++ {
				m := scores[i] >> 32
				if m > score || m == score && int(uint32(scores[i])) == offset+want && extendsChain2End(subs[want], subs[i]) {
					want, score = i, m
				}
			}
			got, s := index.best(begin, end)
			if got != want || s != float64(score) {
				t.Fatalf("n=%d range=[%d,%d): got %d/%g, want %d/%d", n, begin, end, got, s, want, score)
			}
		}
	}
}

// Compare the complete traversal with the recursive oracle, including rejected
// chains, partial predecessor paths, bounds and depth-first output order.
func TestChainARegionMatchesReference(t *testing.T) {
	rng := rand.New(rand.NewSource(8494))
	for trial := range 100 {
		n := 64 + rng.Intn(512)
		subs := make([]SubstrPair, n)
		scores := make([]uint64, n)
		offset := rng.Intn(100)
		for i := range subs {
			subs[i] = SubstrPair{QBegin: int32(i * 10), TBegin: int32(rng.Intn(n * 10)), Len: uint8(1 + rng.Intn(31))}
			previous := i
			if i > 0 && rng.Intn(3) == 0 {
				previous = rng.Intn(i)
			}
			scores[i] = uint64(rng.Intn(100))<<32 | uint64(offset+previous)
		}
		var expected []denseChainSummary
		var want [7]int
		var wantScore float64
		minimum, align, identity := float64(rng.Intn(110)), rng.Intn(50), float64(rng.Intn(80))
		for _, reference := range []bool{true, false} {
			paths := poolChains2.Get().(*[]*Chain2Result)
			var matched, q, tgt int
			fn := chainARegion
			if reference {
				fn = referenceChainARegion
			}
			s, qb, qe, tb, te := fn(&subs, &scores, offset, minimum, align, 0, paths, &matched, &q, &tgt, -1, nil, identity)
			got := denseChain2Summary(paths)
			stats := [7]int{matched, q, tgt, qb, qe, tb, te}
			if reference {
				expected, want, wantScore = got, stats, s
			} else if s != wantScore || stats != want || !reflect.DeepEqual(got, expected) {
				t.Fatalf("trial=%d: iterative traversal differs from reference", trial)
			}
			RecycleChaining2Result(paths)
		}
	}
}

func BenchmarkChainARegionDisjoint(b *testing.B) {
	for _, groups := range []int{1000, 4000, 16000} {
		subs := disjointChain2Anchors(groups)
		scores := make([]uint64, len(subs))
		for i := range groups {
			scores[3*i] = 31<<32 | uint64(3*i)
			scores[3*i+1] = 62<<32 | uint64(3*i)
			scores[3*i+2] = 62<<32 | uint64(3*i+1)
		}
		for _, reference := range []bool{true, false} {
			b.Run(fmt.Sprintf("groups=%d/reference=%t", groups, reference), func(b *testing.B) {
				fn := chainARegion
				if reference {
					fn = referenceChainARegion
				}
				// Warm the shared result pools equally; large legacy runs may
				// otherwise calibrate to b.N=1 while the new path runs many times.
				warm := poolChains2.Get().(*[]*Chain2Result)
				var wm, wq, wt int
				fn(&subs, &scores, 0, 38, 50, 0, warm, &wm, &wq, &wt, 2, nil, 15)
				RecycleChaining2Result(warm)
				b.ReportAllocs()
				b.ResetTimer()
				for range b.N {
					paths := poolChains2.Get().(*[]*Chain2Result)
					var m, q, tgt int
					fn(&subs, &scores, 0, 38, 50, 0, paths, &m, &q, &tgt, 2, nil, 15)
					if len(*paths) != groups {
						b.Fatal("lost chains")
					}
					RecycleChaining2Result(paths)
				}
			})
		}
	}
}

// Custom zero-length/zero-threshold inputs retain the legacy complete-path
// comparison, including NaN identity, while normal callers use positive lengths.
func TestChainARegionZeroLength(t *testing.T) {
	subs := []SubstrPair{{}, {QBegin: 1, TBegin: 1}}
	scores := []uint64{0, 1}
	for _, reference := range []bool{true, false} {
		fn := chainARegion
		if reference {
			fn = referenceChainARegion
		}
		paths := poolChains2.Get().(*[]*Chain2Result)
		var m, q, tgt int
		fn(&subs, &scores, 0, 0, 0, 0, paths, &m, &q, &tgt, 0, nil, 0)
		if len(*paths) != 2 || !math.IsNaN((*paths)[0].PIdent) || !math.IsNaN((*paths)[1].PIdent) {
			t.Fatalf("reference=%t: zero-length behavior changed", reference)
		}
		RecycleChaining2Result(paths)
	}
}
