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
