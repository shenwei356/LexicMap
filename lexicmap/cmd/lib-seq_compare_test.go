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

	"github.com/shenwei356/wfa"
	"sync"
	"testing"
)

func TestSeqCompare(t *testing.T) {
	// Identities: 271/288(94%)

	// Query  8    AGGTCCTGCCCCGCGACCTGCACGCCGAATACGTAGCGGCGATCGCCTTAGTCGGTACAG  67
	//             |||||||||||||||||||||||||| ||||||||   |||||||||||||||  |||||
	// Sbjct  15   AGGTCCTGCCCCGCGACCTGCACGCC-AATACGTA-TAGCGATCGCCTTAGTC--TACAG  70

	// Query  68   CCCTGGAAAACATGGCCACCGAAGTTCGTTCCCTGCAACGGACCGAAATCCACGAAGTCG  127
	//             |||||||||||||||||||||||||||||| ||||||||||||||||||||||  |||||
	// Sbjct  71   CCCTGGAAAACATGGCCACCGAAGTTCGTT-CCTGCAACGGACCGAAATCCACTGAGTCG  129

	// Query  128  AAGAACACTTTGCTAAGGGCCAAAAGGGCTCGTCAGCCATGCCGCACAAGCGGAACCCAA  187
	//                || |||| |||||||||||||||||||||||||||||||||||||||||||||||||
	// Sbjct  130  --CAATACTTCGCTAAGGGCCAAAAGGGCTCGTCAGCCATGCCGCACAAGCGGAACCCAA  187

	// Query  188  TTGGCTCCGAAAACATCTGCGGCTGTGCCCGGGTCCTGCGGGGCAACGTGGTGACCGCCT  247
	//             ||||||||||||||||||||||||||||||||||||||||||| ||||||||||||||||
	// Sbjct  188  TTGGCTCCGAAAACATCTGCGGCTGTGCCCGGGTCCTGCGGGG-AACGTGGTGACCGCCT  246

	// Query  248  ACGAAGACGTGACCCTCTGGCACGAACGCGACATCTCCCACTCCAGTG  295
	//             ||||||||||||||||  ||||||||||||||||||||||||||||||
	// Sbjct  247  ACGAAGACGTGACCCTTCGGCACGAACGCGACATCTCCCACTCCAGTG  294

	s1 := []byte("GGTTACGTATTGCTAGGTCCTGCCCCGCGACCTGCACGCCAATACGTATAGCGATCGCCTTAGTCTACAGCCCTGGAAAACATGGCCACCGAAGTTCGTTCCTGCAACGGACCGAAATCCACTGAGTCGCAATACTTCGCTAAGGGCCAAAAGGGCTCGTCAGCCATGCCGCACAAGCGGAACCCAATTGGCTCCGAAAACATCTGCGGCTGTGCCCGGGTCCTGCGGGGAACGTGGTGACCGCCTACGAAGACGTGACCCTTCGGCACGAACGCGACATCTCCCACTCCAGTGAGCAATACGTAACTGAACGAAGAACATCCGCAAAAAAAA")
	s2 := []byte("TCCACCCAGGTCCTGCCCCGCGACCTGCACGCCGAATACGTAGCGGCGATCGCCTTAGTCGGTACAGCCCTGGAAAACATGGCCACCGAAGTTCGTTCCCTGCAACGGACCGAAATCCACGAAGTCGAAGAACACTTTGCTAAGGGCCAAAAGGGCTCGTCAGCCATGCCGCACAAGCGGAACCCAATTGGCTCCGAAAACATCTGCGGCTGTGCCCGGGTCCTGCGGGGCAACGTGGTGACCGCCTACGAAGACGTGACCCTCTGGCACGAACGCGACATCTCCCACTCCAGTGCCGAACGGATGATTCTGCCGGACTCCACGGCGCTGTTG")

	// s1 := []byte("GGTTACGTATTGCTAGGTCCTGCCCCGCGACCTGCACGCCAATACGTATAGCGATCGCCTTAGTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTGCAACGGACCGAAATCCACTGAGTCGCAATACTTCGCTAAGGGCCAAAAGGGCTCGTCAGCCATGCCGCACAAGCGGAACCCAATTGGCTCCGAAAACATCTGCGGCTGTGCCCGGGTCCTGCGGGGcccccccccccccccccccccccccccccccTCGGCACGAACGCGACATCTCCCACTCCAGTGAGCAATACGTAACTGAACGAAGAACATCCGCAAAAAAAA")
	// s2 := []byte("TCCACCCAGGTCCTGCCCCGCGACCTGCACGCCGAATACGTAGCGGCGATCGCCTTAGTCGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCTGCAACGGACCGAAATCCACGAAGTCGAAGAACACTTTGCTAAGGGCCAAAAGGGCTCGTCAGCCATGCCGCACAAGCGGAACCCAATTGGCTCCGAAAACATCTGCGGCTGTGCCCGGGTCCTGCGGGGCaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaTCTGGCACGAACGCGACATCTCCCACTCCAGTGCCGAACGGATGATTCTGCCGGACTCCACGGCGCTGTTG")

	// alignment
	// alg := align.NewAligner(&align.AlignOptions{
	// 	MatchScore:     1,
	// 	MisMatchScore:  -1,
	// 	GapScore:       -1,
	// 	SaveAlignments: true,
	// 	SaveMatrix:     false,
	// })
	// r := alg.Global(s2, s1)

	// t.Logf("matches: %d, gaps: %d, len: %d, identity: %.2f%%\n",
	// 	r.Matches, r.Gaps, r.Len, float64(r.Matches)/float64(r.Len)*100)

	// compare

	cpr := NewSeqComparator(&DefaultSeqComparatorOptions, &sync.Pool{New: func() interface{} {
		return NewChainer2(&DefaultChaining2Options)
	}})

	err := cpr.Index(s1)
	if err != nil {
		t.Logf("%s", err)
		return
	}

	cr, err := cpr.Compare(0, uint32(len(s2)-1), s2, len(s2))
	if err != nil {
		t.Logf("%s", err)
		return
	}
	if cr != nil {
		t.Logf("aligned bases: %d\n", cr.AlignedBases)

		RecycleSeqComparatorResult(cr)
	}
}

//

// reorderedChainSequences puts a short early-query match after a long late-query
// match on the subject, so coordinate order differs from chaining-score order.
func reorderedChainSequences(seed int64) ([]byte, []byte) {
	rng := rand.New(rand.NewSource(seed))
	q := make([]byte, 1020)
	spacer := make([]byte, 150)
	for _, s := range [][]byte{q, spacer} {
		for i := range s {
			s[i] = "ACGT"[rng.Intn(4)]
		}
	}
	s := append([]byte(nil), q[250:]...)
	s = append(s, spacer...)
	s = append(s, q[:100]...)
	return q, s
}

// TestSeqCompareBestChainSurvivesSorting checks both ordering contracts and
// verifies that transferring a chain survives recycling the remaining result.
func TestSeqCompareBestChainSurvivesSorting(t *testing.T) {
	idx := compareTestIndex(t)
	cpr := idx.poolSeqComparator.Get().(*SeqComparator)
	defer idx.poolSeqComparator.Put(cpr)
	defer cpr.RecycleIndex()
	q, s := reorderedChainSequences(719)
	if err := cpr.Index(q); err != nil {
		t.Fatal(err)
	}
	r, err := cpr.Compare(0, uint32(len(q)), s, len(q))
	if err != nil || r == nil {
		t.Fatalf("comparison: result=%v err=%v", r, err)
	}
	if len(*r.Chains) < 2 || r.BestChain == (*r.Chains)[0] || r.BestChain.AlignedBasesQ < 750 {
		t.Fatalf("expected a later, longer best chain: %+v", r)
	}
	for i := 1; i < len(*r.Chains); i++ {
		if (*r.Chains)[i-1].QBegin > (*r.Chains)[i].QBegin {
			t.Fatal("Chains lost coordinate order")
		}
	}
	best := r.BestChain
	best.CIGAR = []byte("retained") // recycling a chain clears its output buffers
	qBegin, qEnd := best.QBegin, best.QEnd
	if r.takeBestChain() != best || r.BestChain != nil || r.takeBestChain() != nil {
		t.Fatal("best-chain ownership was not transferred exactly once")
	}
	for _, chain := range *r.Chains {
		if chain == best {
			t.Fatal("transferred chain remains owned by result")
		}
	}
	RecycleSeqComparatorResult(r)
	if r.BestChain != nil || best.QBegin != qBegin || best.QEnd != qEnd || string(best.CIGAR) != "retained" {
		t.Fatal("recycling result changed the transferred best chain")
	}
	recycleChain2(best)

	// A subsequent single-chain comparison must initialize the pointer too.
	r, err = cpr.Compare(0, uint32(len(q)), q, len(q))
	if err != nil || r == nil {
		t.Fatalf("self comparison: result=%v err=%v", r, err)
	}
	defer RecycleSeqComparatorResult(r)
	if len(*r.Chains) != 1 || r.BestChain != (*r.Chains)[0] {
		t.Fatal("single-chain comparison did not expose its best chain")
	}
}

// TestAlignChainUsesBestPseudoChain checks the default genome alignment path:
// a short coordinate-first match must not hide a later match that passes AF.
func TestAlignChainUsesBestPseudoChain(t *testing.T) {
	idx := compareTestIndex(t)
	cpr := idx.poolSeqComparator.Get().(*SeqComparator)
	defer idx.poolSeqComparator.Put(cpr)
	defer cpr.RecycleIndex()
	q, s := reorderedChainSequences(719)
	if err := cpr.Index(q); err != nil {
		t.Fatal(err)
	}
	algn := wfa.New(wfa.DefaultPenalties, &wfa.Options{GlobalAlignment: true})
	defer wfa.RecycleAligner(algn)
	score := scoreAndEvalue(2, -3, 5, 2, len(s), 0.625, 0.41)
	matched, aligned, _, identity, _, ok := alignChain(q, nil, s,
		&Chain2Result{QBegin: 0, QEnd: len(q) - 1, TBegin: 0, TEnd: len(s) - 1},
		&subjectSketch{seqLen: len(s)}, algn, cpr, 0, 50, 70, 30, idx, &score)
	if !ok || matched < 750 || aligned < 750 || identity < 99 {
		t.Fatalf("best chain was not aligned: matched=%d aligned=%d identity=%g ok=%v", matched, aligned, identity, ok)
	}
}

// TestAlignChainRepeatWindow exercises the production pseudo-alignment and
// global-WFA path, including aligner/result reuse on a second invocation.
func TestAlignChainRepeatWindow(t *testing.T) {
	for _, reverse := range []bool{false, true} {
		for _, maxGap := range []int{20, 100} {
			t.Run(fmt.Sprintf("reverse=%t/maxGap=%d", reverse, maxGap), func(t *testing.T) {
				idx := compareTestIndex(t)
				opt := *idx.seqCompareOption
				opt.Chaining2Options.MaxGap = maxGap
				cpr := NewSeqComparator(&opt, &sync.Pool{New: func() interface{} { return NewChainer2(&opt.Chaining2Options) }})
				defer cpr.RecycleIndex()
				q, target := []byte(repeatWindowQuery), []byte(repeatWindowTarget)
				qb, qe := 45, 399
				if reverse {
					RC(q)
					RC(target)
					qb, qe = len(q)-qe, len(q)-qb
				}
				// Poor-quality clipped ends must not enter the region's correction.
				qual := make([]byte, len(q))
				for i := range qual {
					qual[i] = 33
					if i >= qb && i < qe {
						qual[i] = 73
					}
				}
				if err := cpr.Index(q); err != nil {
					t.Fatal(err)
				}
				algn := wfa.New(wfa.DefaultPenalties, &wfa.Options{GlobalAlignment: true})
				algn.AdaptiveReduction(wfa.DefaultAdaptiveOption)
				defer wfa.RecycleAligner(algn)
				score := scoreAndEvalue(2, -3, 5, 2, len(target), 0.625, 0.41)
				for trial := 0; trial < 2; trial++ {
					matched, aligned, gaps, identity, adjusted, ok := alignChain(q, qual, target,
						&Chain2Result{QBegin: 0, QEnd: len(q) - 1, TBegin: 0, TEnd: len(target) - 1},
						&subjectSketch{seqLen: len(target)}, algn, cpr, 0, 100, 70, 30, idx, &score)
					if adjusted != adjustPIdentByQual(identity, qual[qb:qe]) {
						t.Fatal("quality correction used the clipped query flanks")
					}
					if !ok || matched != 344 || aligned != 354 || gaps != 0 {
						t.Fatalf("trial %d: matched=%d aligned=%d gaps=%d ok=%t; want 344/354 without gaps", trial, matched, aligned, gaps, ok)
					}
				}
			})
		}
	}
}

// TestSeqCompareDefaultGapBoundary compares refinement with ordinary WFA on
// windows produced by the real pseudo-alignment path at MaxGap=20. Gaps above
// that anchor-transition limit may split chains; do not require a single HSP.
func TestSeqCompareDefaultGapBoundary(t *testing.T) {
	seenTwenty := make(map[string]bool)
	for _, repeated := range []bool{false, true} {
		base := []byte(indexTestDNA(1200))
		if repeated {
			for i := 400; i < 800; i++ {
				base[i] = "ACGT"[i%4]
			}
		}
		for _, size := range []int{1, 10, 11, 19, 20, 21, 25, 30, 50, 100} {
			insertion := []byte(indexTestDNA(size))
			query := append([]byte(nil), base[:600]...)
			query = append(query, insertion...)
			query = append(query, base[600:]...)
			for _, swap := range []bool{false, true} {
				for _, reverse := range []bool{false, true} {
					name := fmt.Sprintf("repeat=%t/gap=%d/swap=%t/reverse=%t", repeated, size, swap, reverse)
					t.Run(name, func(t *testing.T) {
						q, target := append([]byte(nil), query...), append([]byte(nil), base...)
						if swap {
							q, target = target, q
						}
						if reverse {
							RC(q)
							RC(target)
						}
						opt := SeqComparatorOptions{K: 31, MinPrefix: 11, Chaining2Options: Chaining2Options{
							MaxGap: 20, MinScore: 35, MinAlignLen: 50, MinIdentity: 70, BandBase: 100, BandCount: 50, HeuristicKmerPidentThreshold: 15}, MinIdentity: 70}
						cpr := NewSeqComparator(&opt, &sync.Pool{New: func() interface{} { return NewChainer2(&opt.Chaining2Options) }})
						if e := cpr.Index(q); e != nil {
							t.Fatal(e)
						}
						defer cpr.RecycleIndex()
						cr, e := cpr.Compare(0, uint32(len(q)-1), target, len(q))
						if e != nil || cr == nil {
							t.Fatalf("missing pseudoalignment: %v", e)
						}
						defer RecycleSeqComparatorResult(cr)
						algn := wfa.New(wfa.DefaultPenalties, &wfa.Options{GlobalAlignment: true})
						algn.AdaptiveReduction(wfa.DefaultAdaptiveOption)
						defer wfa.RecycleAligner(algn)
						score := scoreAndEvalue(2, -3, 5, 2, len(target), 0.625, 0.41)
						for _, chain := range *cr.Chains {
							a, b, _, _, _, _, e := extendMatch(q, target, chain.QBegin, chain.QEnd+1, chain.TBegin, chain.TEnd+1, 50, chain.TBegin, len(target)-1-chain.TEnd, false)
							if e != nil {
								t.Fatal(e)
							}
							before, e := algn.Align(a, b)
							if e != nil {
								t.Fatal(e)
							}
							after, e := alignWithRepeatRefinement(algn, a, b, 31, len(q), 0, score)
							if e != nil {
								wfa.RecycleAlignmentResult(before)
								t.Fatal(e)
							}
							same := reflect.DeepEqual(before, after)
							if !repeated && size == 20 && before.Gaps == 20 {
								seenTwenty[fmt.Sprintf("swap=%t/reverse=%t", swap, reverse)] = true
							}
							wfa.RecycleAlignmentResult(before)
							wfa.RecycleAlignmentResult(after)
							if !same {
								t.Fatal("refinement changed an alignment at the default chaining gap limit")
							}
						}
					})
				}
			}
		}
	}
	if len(seenTwenty) != 4 {
		t.Fatalf("20-base internal gaps were not exercised in all orientations: %v", seenTwenty)
	}
}
