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
