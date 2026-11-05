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
	"math/rand"
	"reflect"
	"slices"
	"testing"
)

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
