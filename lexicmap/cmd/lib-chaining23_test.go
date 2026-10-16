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

func denseChain2Summary(chains *[]*Chain2Result, nInputs int) []denseChainSummary {
	if chains == nil {
		return nil
	}
	out := make([]denseChainSummary, len(*chains))
	for i, c := range *chains {
		out[i] = denseChainSummary{c.NAnchors, c.MatchedBases, c.AlignedBasesQ, c.AlignedBasesT,
			c.QBegin, c.QEnd, c.TBegin, c.TEnd, c.PIdent}
		if nInputs == 1 {
			// The unchanged single-anchor branch does not write this field.
			// Compare computed outputs, not stale values from different pools.
			out[i].alignedT = 0
		}
	}
	return out
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
					!reflect.DeepEqual(denseChain2Summary(g2, len(subs)), denseChain2Summary(r2, len(subs))) {
					t.Fatalf("Chainer2 trial %d (n=%d): matrices equal=%v, got stats=%v want=%v, got paths=%+v want=%+v", trial, len(subs),
						slices.Equal(got2.maxscoresIdxs, ref2.maxscoresIdxs),
						[7]int{gm, gq, gt, gqb, gqe, gtb, gte}, [7]int{rm, rq, rt, rqb, rqe, rtb, rte},
						denseChain2Summary(g2, len(subs)), denseChain2Summary(r2, len(subs)))
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
			!reflect.DeepEqual(denseChain2Summary(g, len(subs)), denseChain2Summary(r, len(subs))) {
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
