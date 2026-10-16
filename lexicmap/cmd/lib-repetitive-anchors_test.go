package cmd

import (
	"math"
	"math/rand"
	"reflect"
	"slices"
	"testing"
)

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
