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
	"math"
	"math/rand"
	"reflect"
	"slices"
	"testing"
)

// Sparse diagonals exercise whole-row skips; nearby diagonals, negative bins,
// and folded bitmap collisions must preserve all scores and tie predecessors.
func TestChainer2DiagonalSkipsMatchReference(t *testing.T) {
	rng := rand.New(rand.NewSource(20261020))
	opt := DefaultChaining2Options
	got, ref := NewChainer2(&opt), NewChainer2(&opt)
	for trial := 0; trial < 240; trial++ {
		subs := make([]SubstrPair, 511+trial%260)
		for i := range subs {
			q := int32(i / 3 * 7)
			target := q + int32(rng.Intn(8000000))
			if i%5 == 0 {
				target = q + int32(rng.Intn(11)) - 5
			}
			if i%7 == 0 { // collide after folding 64-wide bins into 65,536 bits
				target = q + 64*(1<<16)
			}
			subs[i] = SubstrPair{QBegin: q, TBegin: target, Len: 13}
		}
		slices.SortFunc(subs, func(a, b SubstrPair) int {
			if a.QBegin != b.QBegin {
				return int(a.QBegin - b.QBegin)
			}
			return int(a.TBegin - b.TBegin)
		})
		switch trial % 8 {
		case 1:
			for i := range subs { // shifting both axes changes no gaps
				subs[i].QBegin -= 4000000
				subs[i].TBegin -= 4000000
			}
		case 2:
			subs[len(subs)/2].Len = 31
		case 3:
			slices.Reverse(subs)
		case 4:
			subs[len(subs)-1].TBegin = math.MaxInt32
		}
		opt = DefaultChaining2Options
		opt.MinScore, opt.MinAlignLen = 1, 1
		opt.MaxGap = []int{-1, 0, 1, 50, 63, 64, 100, math.MaxInt}[trial/8%8]
		opt.BandBase = []int{-20, 0, 50, 100, 500}[trial%5]
		opt.BandCount = []int{-1, 0, 1, 7, 50}[trial/5%5]
		g, gm, gq, gt, gqb, gqe, gtb, gte := got.Chain(&subs)
		r, rm, rq, rt, rqb, rqe, rtb, rte := referenceChaining2(ref, &subs)
		if !slices.Equal(got.maxscoresIdxs, ref.maxscoresIdxs) ||
			[7]int{gm, gq, gt, gqb, gqe, gtb, gte} != [7]int{rm, rq, rt, rqb, rqe, rtb, rte} ||
			!reflect.DeepEqual(denseChain2Summary(g, len(subs)), denseChain2Summary(r, len(subs))) {
			t.Fatalf("trial %d options %+v: diagonal filtering changed chaining", trial, opt)
		}
		if g != nil {
			RecycleChaining2Result(g)
		}
		if r != nil {
			RecycleChaining2Result(r)
		}
	}
}

// Compare the complete score/predecessor matrix, not just the winning chain:
// skipped gap failures still affect the legacy band-count stopping position.
func TestChainer2GapGroupSkipsMatchReference(t *testing.T) {
	rng := rand.New(rand.NewSource(20261016))
	for trial := 0; trial < 600; trial++ {
		var subs []SubstrPair
		q := int32(0)
		for group := 0; group < 20; group++ {
			q += int32(1 + rng.Intn(50))
			for copy := 0; copy < 1+rng.Intn(35); copy++ {
				subs = append(subs, SubstrPair{QBegin: q, TBegin: int32(rng.Intn(1600)), Len: 13})
			}
		}
		slices.SortFunc(subs, func(a, b SubstrPair) int {
			if a.QBegin != b.QBegin {
				return int(a.QBegin - b.QBegin)
			}
			return int(a.TBegin - b.TBegin)
		})
		// Mixed lengths, unsorted inputs, and wide coordinates must retain
		// their original behavior without the equal-length group shortcut.
		switch trial % 6 {
		case 1:
			subs[len(subs)/2].Len = 31
		case 2:
			subs[0], subs[len(subs)-1] = subs[len(subs)-1], subs[0]
		case 3:
			subs[len(subs)-1].TBegin = math.MaxInt32
		case 4:
			subs[0].QBegin = math.MinInt32
		}
		opt := DefaultChaining2Options
		opt.MinScore, opt.MinAlignLen = 1, 1
		opt.MaxGap = []int{-1, 0, 1, 13, 50, 100, math.MaxInt}[trial%7]
		opt.BandBase = []int{-20, 0, 13, 100, 500}[trial%5]
		opt.BandCount = []int{-1, 0, 1, 7, 50}[trial/5%5]
		got, ref := NewChainer2(&opt), NewChainer2(&opt)
		g, gm, gq, gt, gqb, gqe, gtb, gte := got.Chain(&subs)
		r, rm, rq, rt, rqb, rqe, rtb, rte := referenceChaining2(ref, &subs)
		if !slices.Equal(got.maxscoresIdxs, ref.maxscoresIdxs) ||
			[7]int{gm, gq, gt, gqb, gqe, gtb, gte} != [7]int{rm, rq, rt, rqb, rqe, rtb, rte} ||
			!reflect.DeepEqual(denseChain2Summary(g, len(subs)), denseChain2Summary(r, len(subs))) {
			t.Fatalf("trial %d options %+v: chaining differs from reference", trial, opt)
		}
		if g != nil {
			RecycleChaining2Result(g)
		}
		if r != nil {
			RecycleChaining2Result(r)
		}
	}
}

// Equal-length containment is exact coordinate duplication. Adjacent,
// overlapping, and different-copy anchors must all survive unchanged.
func TestClearEqualLengthSubstrPairsMatchesGeneral(t *testing.T) {
	rng := rand.New(rand.NewSource(29))
	for _, k := range []int{1, 13, 31} {
		for n := 0; n < 400; n++ {
			subs := make([]SubstrPair, n)
			for i := range subs {
				subs[i] = SubstrPair{QBegin: int32(rng.Intn(25)), TBegin: int32(rng.Intn(35)), Len: uint8(k)}
			}
			want := slices.Clone(subs)
			if n > 1 {
				ClearSubstrPairs(&want, k)
			}
			clearEqualLengthSubstrPairs(&subs)
			if !slices.Equal(subs, want) {
				t.Fatalf("k=%d n=%d: anchors differ", k, n)
			}
		}
	}
}
