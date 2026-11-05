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

type denseChainSummary struct {
	anchors, matched, alignedQ, alignedT int
	qb, qe, tb, te                       int
	identity                             float64
}

func denseChain2Summary(chains *[]*Chain2Result) []denseChainSummary {
	if chains == nil {
		return nil
	}
	out := make([]denseChainSummary, len(*chains))
	for i, c := range *chains {
		out[i] = denseChainSummary{c.NAnchors, c.MatchedBases, c.AlignedBasesQ, c.AlignedBasesT,
			c.QBegin, c.QEnd, c.TBegin, c.TEnd, c.PIdent}
	}
	return out
}

// TestChainer2FirstAnchorSelection checks independent optima and tie order,
// including reversed coordinates and reuse after a rejected chain.
func TestChainer2FirstAnchorSelection(t *testing.T) {
	opt := DefaultChaining2Options
	chainer := NewChainer2(&opt)
	for _, tc := range []struct {
		name                string
		first, second       uint8
		minScore, minLength int
	}{
		{"first maximum", 31, 11, 22, 20},
		{"first tied maximum", 31, 31, 22, 20},
		{"below score", 31, 11, 32, 20},
		{"later maximum", 11, 31, 22, 20},
		{"below length", 31, 11, 22, 32},
	} {
		for _, reverse := range []bool{false, true} {
			t.Run(fmt.Sprintf("%s/reverse=%t", tc.name, reverse), func(t *testing.T) {
				opt.MinScore, opt.MinAlignLen = tc.minScore, tc.minLength
				// The 250-base diagonal difference prevents these anchors from linking.
				subs := []SubstrPair{{QBegin: 20, TBegin: 20, Len: tc.first}, {QBegin: 200, TBegin: 450, Len: tc.second}}
				if reverse {
					for i := range subs {
						subs[i].QBegin = 600 - subs[i].QBegin - int32(subs[i].Len)
						subs[i].TBegin = 600 - subs[i].TBegin - int32(subs[i].Len)
					}
					slices.Reverse(subs)
				}
				var want []denseChainSummary
				for _, sub := range subs {
					n := int(sub.Len)
					if n >= tc.minScore && n >= tc.minLength {
						want = append(want, denseChainSummary{1, n, n, n,
							int(sub.QBegin), int(sub.QBegin) + n - 1, int(sub.TBegin), int(sub.TBegin) + n - 1, 100})
					}
				}
				chains, _, _, _, _, _, _, _ := chainer.Chain(&subs)
				got := denseChain2Summary(chains)
				if chains != nil {
					RecycleChaining2Result(chains)
				}
				if !reflect.DeepEqual(got, want) {
					t.Fatalf("chains = %+v, want %+v", got, want)
				}
			})
		}
	}
}

// TestChainer2SingleAnchorRetainsExtensionCandidate checks the singleton
// contract independently of the reference kernel and exercises pooled reuse.
func TestChainer2SingleAnchorRetainsExtensionCandidate(t *testing.T) {
	opt := DefaultChaining2Options
	opt.MinScore, opt.MinAlignLen = 35, 50
	chainer := NewChainer2(&opt)
	for _, length := range []int{11, 17, 31, 17} {
		subs := []SubstrPair{{QBegin: 20, TBegin: 40, Len: uint8(length)}}
		chains, matches, aq, at, qb, qe, tb, te := chainer.Chain(&subs)
		if chains == nil {
			t.Fatalf("discarded %d-base extension candidate", length)
		}
		want := []denseChainSummary{{anchors: 1, matched: length, alignedQ: length, alignedT: length,
			qb: 20, qe: 19 + length, tb: 40, te: 39 + length, identity: 100}}
		got := denseChain2Summary(chains)
		RecycleChaining2Result(chains)
		if !reflect.DeepEqual(got, want) ||
			[7]int{matches, aq, at, qb, qe, tb, te} != [7]int{length, length, length, 20, 19 + length, 40, 39 + length} {
			t.Fatalf("length %d: got %+v, want %+v", length, got, want)
		}
	}
}

// Per-call extension mode must not leave a length override in a pooled chainer.
func TestChainer2ExtensionLengthOverride(t *testing.T) {
	opt := DefaultChaining2Options
	opt.MinAlignLen = 50
	chainer := NewChainer2(&opt)
	subs := []SubstrPair{{QBegin: 10, TBegin: 10, Len: 31}, {QBegin: 20, TBegin: 20, Len: 31}}
	for _, extension := range []bool{false, true, false} {
		var chains *[]*Chain2Result
		if extension {
			chains, _, _, _, _, _, _, _ = chainer.chainWithMinScore(&subs, opt.MinScore)
		} else {
			chains, _, _, _, _, _, _, _ = chainer.Chain(&subs)
		}
		if (chains != nil) != extension {
			t.Errorf("extension=%t: retained=%t", extension, chains != nil)
		}
		if chains != nil {
			if len(*chains) != 1 || (*chains)[0].QBegin != 10 || (*chains)[0].QEnd != 50 {
				t.Errorf("unexpected extension candidate: %+v", denseChain2Summary(chains))
			}
			RecycleChaining2Result(chains)
		}
		if opt.MinAlignLen != 50 {
			t.Fatal("mutated shared length threshold")
		}
	}
}

// TestChainer2MinimumScore checks the inclusive score boundary independently
// of the reference kernel, including reversed coordinates and gap penalties.
func TestChainer2MinimumScore(t *testing.T) {
	for _, tc := range []struct {
		name       string
		length     uint8
		gap        int32
		wantChains int
	}{
		{"two-11bp-no-gap", 11, 0, 1},
		{"two-11bp-insertion", 11, 1, 0},
		{"two-11bp-deletion", 11, -1, 0},
		{"two-12bp-gap-two", 12, 2, 1},
	} {
		for _, reverse := range []bool{false, true} {
			t.Run(fmt.Sprintf("%s/reverse=%t", tc.name, reverse), func(t *testing.T) {
				opt := DefaultChaining2Options
				subs := []SubstrPair{{Len: tc.length}, {QBegin: 56, TBegin: 56 + tc.gap, Len: tc.length}}
				if reverse {
					for i := range subs {
						subs[i].QBegin = 200 - subs[i].QBegin - int32(subs[i].Len)
						subs[i].TBegin = 200 - subs[i].TBegin - int32(subs[i].Len)
					}
					slices.Reverse(subs)
				}
				chains, _, _, _, _, _, _, _ := NewChainer2(&opt).Chain(&subs)
				got := denseChain2Summary(chains)
				if chains != nil {
					RecycleChaining2Result(chains)
				}
				if len(got) != tc.wantChains || (len(got) > 0 && got[0].anchors != 2) {
					t.Fatalf("chains = %+v, want %d two-anchor chains", got, tc.wantChains)
				}
			})
		}
	}
}

func denseChain3Summary(c *Chain3Result) *denseChainSummary {
	if c == nil {
		return nil
	}
	return &denseChainSummary{c.NAnchors, c.MatchedBases, c.AlignedBasesQ, c.AlignedBasesT,
		c.QBegin, c.QEnd, c.TBegin, c.TEnd, c.PIdent}
}

func TestChaining3DistanceAndGapMatchesFloat(t *testing.T) {
	coords := []int32{math.MinInt32, math.MinInt32 + 1, -100, -1, 0, 1, 100, math.MaxInt32}
	for _, aq := range coords {
		for _, at := range coords {
			for _, bq := range coords {
				for _, bt := range coords {
					a, b := SubstrPair{QBegin: aq, TBegin: at}, SubstrPair{QBegin: bq, TBegin: bt}
					d, g := chaining3DistanceAndGap(aq, at, bq, bt)
					if float64(d) != distance2(a, b) || float64(g) != gap2(a, b) {
						t.Fatalf("a=%v b=%v: distance/gap=%d/%d", a, b, d, g)
					}
				}
			}
		}
	}
}

// TestChaining3FirstAnchorSelection checks optima independently of the reference kernel.
func TestChaining3FirstAnchorSelection(t *testing.T) {
	for _, tc := range []struct {
		name string
		subs []SubstrPair
		want *denseChainSummary
	}{
		{
			name: "singleton",
			subs: []SubstrPair{{Len: 2}},
			want: &denseChainSummary{matched: 2, alignedQ: 2, alignedT: 2, qe: 1, te: 1, identity: 100},
		},
		{
			name: "first sole maximum",
			subs: []SubstrPair{{Len: 4}, {QBegin: 5, TBegin: 5, Len: 2}},
			want: &denseChainSummary{matched: 4, alignedQ: 4, alignedT: 4, qe: 3, te: 3, identity: 100},
		},
		{
			name: "first tied maximum",
			subs: []SubstrPair{{Len: 4}, {QBegin: 4, TBegin: 4, Len: 2}},
			want: &denseChainSummary{matched: 4, alignedQ: 4, alignedT: 4, qe: 3, te: 3, identity: 100},
		},
		{
			name: "later maximum with tied predecessor score",
			// Equal predecessor and standalone scores must still link to the earlier anchor.
			subs: []SubstrPair{{Len: 2}, {QBegin: 1, TBegin: 1, Len: 4}},
			want: &denseChainSummary{matched: 5, alignedQ: 5, alignedT: 5, qe: 4, te: 4, identity: 100},
		},
		{
			name: "negative score rejected",
			subs: []SubstrPair{{QBegin: 3, TBegin: 3, Len: 2}},
		},
	} {
		t.Run(tc.name, func(t *testing.T) {
			got := NewChainer3(&DefaultChaining3Options).Chain(&tc.subs)
			if got != nil {
				defer poolChain3.Put(got)
			}
			if summary := denseChain3Summary(got); !reflect.DeepEqual(summary, tc.want) {
				t.Fatalf("chain = %+v, want %+v", summary, tc.want)
			}
		})
	}
}

func TestChaining23MatchesFloatReference(t *testing.T) {
	rng := rand.New(rand.NewSource(2302))
	fixtures := [][]SubstrPair{
		{{Len: 31}},
		{{QBegin: 100, TBegin: 200, Len: 5}}, // Negative initial Chainer3 score.
		{{Len: 5}, {QBegin: 1, TBegin: 2, Len: 5}, {QBegin: 2, TBegin: 1, Len: 5}, {QBegin: 4, TBegin: 4, Len: 5}},
		{{QBegin: math.MinInt32, TBegin: 0, Len: 31}, {QBegin: 0, TBegin: math.MinInt32, Len: 31}, {QBegin: math.MaxInt32, TBegin: math.MaxInt32, Len: 31}},
		// Preserve the legacy int32 gap arithmetic and scores above MaxInt32.
		{{QBegin: math.MinInt32, TBegin: 0, Len: 32}, {QBegin: 0, TBegin: 1, Len: 32}, {QBegin: math.MaxInt32, TBegin: 2, Len: 32}},
		{{QBegin: 0, TBegin: math.MaxInt32, Len: 1}, {QBegin: 1, TBegin: math.MaxInt32, Len: 32}, {QBegin: 2, TBegin: math.MaxInt32, Len: 32}},
	}
	// Dense collinear and crossed anchors exercise the prefix score bound and ties.
	for _, crossed := range []bool{false, true} {
		subs := make([]SubstrPair, 256)
		for i := range subs {
			target := i
			if crossed {
				target = i % 64
			}
			subs[i] = SubstrPair{QBegin: int32(i), TBegin: int32(target), Len: 13}
		}
		fixtures = append(fixtures, subs)
	}
	for trial := range 500 {
		n := 1 + rng.Intn(192)
		subs := make([]SubstrPair, n)
		q, target := int32(0), int32(0)
		for i := range subs {
			q += int32(rng.Intn(8))
			target += int32(rng.Intn(8))
			subs[i] = SubstrPair{QBegin: q, TBegin: target, Len: uint8(1 + rng.Intn(32))}
			if trial%3 == 0 { // Repeated or crossed target anchors.
				subs[i].TBegin = int32(rng.Intn(200))
			}
		}
		fixtures = append(fixtures, subs)
	}
	for _, settings := range []struct {
		name                       string
		gap, distance, score, band int
	}{
		{"normal", 50, 100, 50, 50},
		{"score-22", 20, 100, 22, 50},
		{"strict", 0, 0, 1, 0},
		{"permissive", 5000, 5000, -1, 20},
		{"wide", math.MaxInt, math.MaxInt, 1, 50},
		{"reject", 50, 100, math.MaxInt, 50},
	} {
		t.Run(settings.name, func(t *testing.T) {
			opt2 := DefaultChaining2Options
			opt2.MaxGap, opt2.MinScore, opt2.MinAlignLen = settings.gap, settings.score, 1
			opt2.BandCount = settings.band
			opt3 := DefaultChaining3Options
			opt3.MaxGap, opt3.MaxDistance = settings.gap, settings.distance
			opt3.MinScore, opt3.MinAlignLen, opt3.BandCount = settings.score, 1, settings.band
			got2, ref2 := NewChainer2(&opt2), NewChainer2(&opt2)
			got3, ref3 := NewChainer3(&opt3), NewChainer3(&opt3)
			for trial, subs := range fixtures {
				g2, gm, gq, gt, gqb, gqe, gtb, gte := got2.Chain(&subs)
				r2, rm, rq, rt, rqb, rqe, rtb, rte := referenceChaining2(ref2, &subs)
				if !slices.Equal(got2.maxscoresIdxs, ref2.maxscoresIdxs) ||
					[7]int{gm, gq, gt, gqb, gqe, gtb, gte} != [7]int{rm, rq, rt, rqb, rqe, rtb, rte} ||
					!reflect.DeepEqual(denseChain2Summary(g2), denseChain2Summary(r2)) {
					t.Fatalf("Chainer2 trial %d (n=%d): matrices equal=%v, got stats=%v want=%v, got paths=%+v want=%+v", trial, len(subs),
						slices.Equal(got2.maxscoresIdxs, ref2.maxscoresIdxs),
						[7]int{gm, gq, gt, gqb, gqe, gtb, gte}, [7]int{rm, rq, rt, rqb, rqe, rtb, rte},
						denseChain2Summary(g2), denseChain2Summary(r2))
				}
				if g2 != nil {
					RecycleChaining2Result(g2)
				}
				if r2 != nil {
					RecycleChaining2Result(r2)
				}
				g3, r3 := got3.Chain(&subs), referenceChaining3(ref3, &subs)
				if !slices.Equal(got3.maxscoresIdxs, ref3.maxscoresIdxs) ||
					!reflect.DeepEqual(denseChain3Summary(g3), denseChain3Summary(r3)) {
					t.Fatalf("Chainer3 trial %d: packed scores, predecessors or results differ", trial)
				}
				if g3 != nil {
					poolChain3.Put(g3)
				}
				if r3 != nil {
					poolChain3.Put(r3)
				}
			}
		})
	}
}

// TestChainer2TargetIndexMatchesReference exercises long crossed runs, equal
// query starts, and both ordering fallbacks for the target group index.
func TestChainer2TargetIndexMatchesReference(t *testing.T) {
	for _, order := range []string{"sorted", "target-unsorted", "query-unsorted"} {
		subs := make([]SubstrPair, 0, 64*64)
		for group := 0; group < 64; group++ {
			for copy := 0; copy < 64; copy++ {
				target := int32(100000 + group*25 + copy)
				if copy == 0 {
					target = int32(group * 25)
				}
				subs = append(subs, SubstrPair{QBegin: int32(group * 25), TBegin: target, Len: 13})
			}
		}
		if order == "target-unsorted" {
			// Retain Q groups but break their target order.
			subs[65], subs[127] = subs[127], subs[65]
		}
		if order == "query-unsorted" {
			// Break query grouping without changing the original scan contract.
			subs[65], subs[128] = subs[128], subs[65]
		}
		opt := DefaultChaining2Options
		opt.MinScore, opt.MinAlignLen = 1, 1
		got, ref := NewChainer2(&opt), NewChainer2(&opt)
		g, gm, gq, gt, gqb, gqe, gtb, gte := got.Chain(&subs)
		r, rm, rq, rt, rqb, rqe, rtb, rte := referenceChaining2(ref, &subs)
		if (len(got.qGroupStarts) > 0) != (order == "sorted") {
			t.Fatalf("order=%s: unexpected target group index state", order)
		}
		if !slices.Equal(got.maxscoresIdxs, ref.maxscoresIdxs) ||
			[7]int{gm, gq, gt, gqb, gqe, gtb, gte} != [7]int{rm, rq, rt, rqb, rqe, rtb, rte} ||
			!reflect.DeepEqual(denseChain2Summary(g), denseChain2Summary(r)) {
			t.Fatalf("order=%s: chaining differs from reference", order)
		}
		if g != nil {
			RecycleChaining2Result(g)
		}
		if r != nil {
			RecycleChaining2Result(r)
		}
	}
}

// Large mixed-length groups exercise diagonal score bounds beyond the small
// target-group index, including hash collisions and reverse-complement order.
func TestChainer2LargeDiagonalBoundsMatchReference(t *testing.T) {
	for _, reverse := range []bool{false, true} {
		for _, order := range []string{"sorted", "target-unsorted", "query-unsorted"} {
			t.Run(fmt.Sprintf("reverse=%t/%s", reverse, order), func(t *testing.T) {
				subs := make([]SubstrPair, 0, thresholdNSubsLong+8)
				// MaxGap=20 gives 32-base bins; offsets separated by 2^21
				// collide in the folded table and must remain conservative.
				offsets := []int32{0, 0, 19, 20, 21, 100000, 1 << 21, 1<<21 + 20}
				for group := 0; len(subs) <= thresholdNSubsLong; group++ {
					for copy, offset := range offsets {
						q := int32(group * 3)
						// Forward diagonals are negative; mirroring makes them
						// positive, exercising signed binning in both directions.
						a := SubstrPair{QBegin: q + 3000000, TBegin: q + offset + int32(group%7), Len: uint8(19 + copy%3*6)}
						if reverse {
							a.QBegin = 3100000 - a.QBegin - int32(a.Len)
							a.TBegin = 3000000 - a.TBegin - int32(a.Len)
						}
						subs = append(subs, a)
					}
				}
				slices.SortFunc(subs, func(a, b SubstrPair) int {
					if a.QBegin != b.QBegin {
						return int(a.QBegin - b.QBegin)
					}
					if a.Len != b.Len {
						return int(b.Len) - int(a.Len)
					}
					return int(a.TBegin - b.TBegin)
				})
				if order == "target-unsorted" {
					for i := 0; i+1 < len(subs); i++ {
						if subs[i].QBegin == subs[i+1].QBegin {
							subs[i], subs[i+1] = subs[i+1], subs[i]
						}
					}
				} else if order == "query-unsorted" {
					subs[0], subs[len(subs)-1] = subs[len(subs)-1], subs[0]
				}
				opt := DefaultChaining2Options
				opt.MaxGap, opt.MinScore, opt.MinAlignLen = 20, 1, 1
				got, ref := NewChainer2(&opt), NewChainer2(&opt)
				g, gm, gq, gt, gqb, gqe, gtb, gte := got.Chain(&subs)
				r, rm, rq, rt, rqb, rqe, rtb, rte := referenceChaining2(ref, &subs)
				if g != nil {
					defer RecycleChaining2Result(g)
				}
				if r != nil {
					defer RecycleChaining2Result(r)
				}
				if !slices.Equal(got.maxscoresIdxs, ref.maxscoresIdxs) ||
					[7]int{gm, gq, gt, gqb, gqe, gtb, gte} != [7]int{rm, rq, rt, rqb, rqe, rtb, rte} ||
					!reflect.DeepEqual(denseChain2Summary(g), denseChain2Summary(r)) {
					t.Fatal("large-input scores, predecessors, or chains differ from reference")
				}
				// Reuse the chainer after its diagonal-bin prefix bounds. The
				// following small call must restore ordinary global bounds.
				small := []SubstrPair{{Len: 31}, {QBegin: 20, TBegin: 20, Len: 19}, {QBegin: 40, TBegin: 42, Len: 25}}
				gc, _, _, _, _, _, _, _ := got.Chain(&small)
				rc, _, _, _, _, _, _, _ := referenceChaining2(ref, &small)
				if gc != nil {
					defer RecycleChaining2Result(gc)
				}
				if rc != nil {
					defer RecycleChaining2Result(rc)
				}
				if !slices.Equal(got.maxscoresIdxs, ref.maxscoresIdxs) || !reflect.DeepEqual(denseChain2Summary(gc), denseChain2Summary(rc)) {
					t.Fatal("small-call reuse differs from reference")
				}
			})
		}
	}
}

func BenchmarkDenseChaining(b *testing.B) {
	for _, n := range []int{32, 256, 2048} {
		subs := make([]SubstrPair, n)
		for i := range subs {
			subs[i] = SubstrPair{QBegin: int32(i * 3), TBegin: int32(i*3 + i/17), Len: 11}
		}
		for _, reference := range []bool{true, false} {
			label := "integer"
			if reference {
				label = "float"
			}
			b.Run(fmt.Sprintf("chainer2/%d/%s", n, label), func(b *testing.B) {
				chainer := NewChainer2(&DefaultChaining2Options)
				b.ReportAllocs()
				b.ResetTimer()
				for i := 0; i < b.N; i++ {
					var result *[]*Chain2Result
					if reference {
						result, _, _, _, _, _, _, _ = referenceChaining2(chainer, &subs)
					} else {
						result, _, _, _, _, _, _, _ = chainer.Chain(&subs)
					}
					if result != nil {
						RecycleChaining2Result(result)
					}
				}
			})
			b.Run(fmt.Sprintf("chainer3/%d/%s", n, label), func(b *testing.B) {
				chainer := NewChainer3(&DefaultChaining3Options)
				b.ReportAllocs()
				b.ResetTimer()
				for i := 0; i < b.N; i++ {
					var result *Chain3Result
					if reference {
						result = referenceChaining3(chainer, &subs)
					} else {
						result = chainer.Chain(&subs)
					}
					if result != nil {
						poolChain3.Put(result)
					}
				}
			})
		}
	}
}
