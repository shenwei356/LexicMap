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
	"bytes"
	"fmt"
	"math/rand"
	"reflect"
	"slices"
	"strings"

	rtree "github.com/shenwei356/LexicMap/lexicmap/cmd/tree"

	"github.com/shenwei356/wfa"
	"sync"
	"testing"
)

// A DUST score below the cutoff must retain both orientations in the seed tree.
func TestSeqCompareDustOrientation(t *testing.T) {
	query := []byte("AAAATCAAAATAAACAAAACCATTAAACAAA")
	opt := DefaultSeqComparatorOptions
	opt.K = 31
	cpr := NewSeqComparator(&opt, nil)
	var previous []uint64
	for _, reverse := range []bool{false, true} {
		entries, err := cpr.collectIndexEntries(query)
		if err != nil {
			t.Fatal(err)
		}
		if len(entries) != 2 {
			t.Errorf("reverse=%t: got %d entries, want both strands", reverse, len(entries))
		}
		keys := make([]uint64, len(entries))
		for i, entry := range entries {
			keys[i] = entry.Key
		}
		slices.Sort(keys)
		if reverse && !slices.Equal(keys, previous) {
			t.Errorf("reverse complement changed retained seed keys")
		}
		previous = keys
		RC(query)
	}
}

// Terminal repeat trimming must retain the same anchors after reversing both
// sequences, whether the repeated position is on the query or the target.
func TestTrimSubStrPairsReverseComplement(t *testing.T) {
	for _, swap := range []bool{false, true} {
		anchors := []SubstrPair{{QBegin: 0, TBegin: 0, Len: 31}, {QBegin: 0, TBegin: 10, Len: 31}, {QBegin: 200, TBegin: 200, Len: 31}}
		if swap {
			for i := range anchors {
				anchors[i].QBegin, anchors[i].TBegin = anchors[i].TBegin, anchors[i].QBegin
			}
		}
		for _, reverse := range []bool{false, true} {
			subs := slices.Clone(anchors)
			want := slices.Clone(anchors[1:])
			if reverse {
				for _, pairs := range [][]SubstrPair{subs, want} {
					for i := range pairs {
						pairs[i].QBegin = 400 - pairs[i].QBegin - int32(pairs[i].Len)
						pairs[i].TBegin = 400 - pairs[i].TBegin - int32(pairs[i].Len)
					}
					slices.Reverse(pairs)
				}
			}
			capacity := cap(subs)
			TrimSubStrPairs(&subs, 31, 100)
			if !reflect.DeepEqual(subs, want) || cap(subs) != capacity {
				t.Fatalf("swap=%t reverse=%t: got %+v (capacity %d), want %+v (capacity %d)", swap, reverse, subs, cap(subs), want, capacity)
			}
		}
	}
}

func TestSubstrPairOverlap(t *testing.T) {
	for _, tc := range []struct {
		name string
		a, b SubstrPair
		want int32
	}{
		{"identical", SubstrPair{QBegin: 0, TBegin: 0, Len: 31}, SubstrPair{QBegin: 0, TBegin: 0, Len: 31}, 31},
		{"adjacent", SubstrPair{QBegin: 0, TBegin: 0, Len: 31}, SubstrPair{QBegin: 31, TBegin: 31, Len: 31}, 0},
		{"one base", SubstrPair{QBegin: 0, TBegin: 0, Len: 31}, SubstrPair{QBegin: 30, TBegin: 31, Len: 31}, 1},
		{"contained", SubstrPair{QBegin: 0, TBegin: 0, Len: 31}, SubstrPair{QBegin: 5, TBegin: 50, Len: 10}, 10},
		{"crossed target", SubstrPair{QBegin: 0, TBegin: 20, Len: 31}, SubstrPair{QBegin: 50, TBegin: 10, Len: 20}, 10},
		{"empty", SubstrPair{QBegin: 0, TBegin: 0, Len: 31}, SubstrPair{QBegin: 5, TBegin: 5, Len: 0}, 0},
	} {
		t.Run(tc.name, func(t *testing.T) {
			for _, swap := range []bool{false, true} {
				a, b := tc.a, tc.b
				if swap {
					a, b = b, a
				}
				if got := overlap(a, b); got != tc.want {
					t.Errorf("swap=%t: got overlap %d, want %d", swap, got, tc.want)
				}
				// Reflect both half-open intervals into the opposite orientation.
				a.QBegin, a.TBegin = 100-a.QBegin-int32(a.Len), 100-a.TBegin-int32(a.Len)
				b.QBegin, b.TBegin = 100-b.QBegin-int32(b.Len), 100-b.TBegin-int32(b.Len)
				if got := overlap(a, b); got != tc.want {
					t.Errorf("reverse swap=%t: got overlap %d, want %d", swap, got, tc.want)
				}
			}
		})
	}
}

func TestTrimSubStrPairsOverlapThreshold(t *testing.T) {
	for _, tc := range []struct {
		length uint8
		shift  int32
		trim   bool
	}{
		{31, 7, false}, // 24/31 < 80%; the old extra base made this 25/31.
		{30, 6, false}, // Exactly 80% is retained by the strict cutoff.
		{31, 6, true},  // 25/31 > 80% is a repeat candidate.
	} {
		for _, swap := range []bool{false, true} {
			anchors := []SubstrPair{{QBegin: 0, TBegin: 0, Len: tc.length}, {QBegin: tc.shift, TBegin: 20, Len: tc.length}, {QBegin: 200, TBegin: 200, Len: tc.length}}
			if swap {
				for i := range anchors {
					anchors[i].QBegin, anchors[i].TBegin = anchors[i].TBegin, anchors[i].QBegin
				}
			}
			for _, reverse := range []bool{false, true} {
				subs, want := slices.Clone(anchors), slices.Clone(anchors)
				if tc.trim {
					want = want[1:]
				}
				if reverse {
					for _, pairs := range [][]SubstrPair{subs, want} {
						for i := range pairs {
							pairs[i].QBegin = 400 - pairs[i].QBegin - int32(pairs[i].Len)
							pairs[i].TBegin = 400 - pairs[i].TBegin - int32(pairs[i].Len)
						}
						slices.Reverse(pairs)
					}
				}
				TrimSubStrPairs(&subs, 31, 100)
				if !reflect.DeepEqual(subs, want) {
					t.Errorf("length=%d shift=%d swap=%t reverse=%t: got %+v, want %+v", tc.length, tc.shift, swap, reverse, subs, want)
				}
			}
		}
	}
}

// TestSeqCompareQueryWindow checks containment on both anchor orientations and
// both index forms, including matches that only touch a window boundary.
func TestSeqCompareQueryWindow(t *testing.T) {
	for _, compact := range []bool{false, true} {
		for _, reverse := range []bool{false, true} {
			for _, matchLen := range []int{24, 32} {
				t.Run(fmt.Sprintf("compact=%t/reverse=%t/length=%d", compact, reverse, matchLen), func(t *testing.T) {
					opt := DefaultSeqComparatorOptions
					opt.Chaining2Options = DefaultChaining2Options
					opt.MinAlignLen = 20
					query := []byte(indexTestDNA(300))
					indexed := slices.Clone(query)
					if reverse {
						RC(indexed)
					}
					cpr := NewSeqComparator(&opt, &sync.Pool{New: func() interface{} {
						return NewChainer2(&opt.Chaining2Options)
					}})
					if compact {
						entries, err := cpr.collectIndexEntries(indexed)
						if err != nil {
							t.Fatal(err)
						}
						slices.SortFunc(entries, func(a, b rtree.BatchEntry) int {
							if a.Key < b.Key {
								return -1
							}
							if a.Key > b.Key {
								return 1
							}
							return 0
						})
						cpr.compact = rtree.NewCompactTree(opt.K, entries)
					} else if err := cpr.Index(indexed); err != nil {
						t.Fatal(err)
					}
					defer cpr.RecycleIndex()
					for _, start := range []int{100 - matchLen - 1, 100 - matchLen, 99, 100, 200 - matchLen, 201 - matchLen, 200, 201} {
						// The shorter match is a k-mer suffix with mismatched padding.
						padding := 32 - matchLen
						target := slices.Clone(query[start-padding : start+matchLen])
						for i := 0; i < padding; i++ {
							target[i] = "TGCA"[strings.IndexByte("ACGT", target[i])]
						}
						wantStart := start
						if reverse {
							RC(target)
							wantStart = len(query) - start - matchLen
						}
						for _, bounds := range [][2]uint32{{100, 200}, {100, 100}, {200, 100}, {0, 300}} {
							r, err := cpr.Compare(bounds[0], bounds[1], target, len(query))
							if err != nil {
								t.Fatal(err)
							}
							want := wantStart >= int(bounds[0]) && wantStart+matchLen <= int(bounds[1])
							if (r != nil) != want {
								t.Errorf("start=%d bounds=%v: retained=%t, want %t", wantStart, bounds, r != nil, want)
							}
							if r != nil {
								if want && (len(*r.Chains) != 1 || r.BestChain.QBegin != wantStart || r.BestChain.QEnd != wantStart+matchLen-1 || r.BestChain.MatchedBases != matchLen) {
									t.Errorf("start=%d bounds=%v: wrong match: %+v", wantStart, bounds, r.BestChain)
								}
								RecycleSeqComparatorResult(r)
							}
						}
					}
				})
			}
		}
	}
}

// sparsePrefixSequences creates exactly two matching anchors separated by a
// mismatch-rich interval. Homopolymer padding changes only the target length.
func sparsePrefixSequences(anchorLen, gap, targetLen int) ([]byte, []byte) {
	query := []byte(indexTestDNA(200))
	target := bytes.Repeat([]byte{'A'}, targetLen)
	for i, pos := range []int{20, 100} {
		query[pos-1], query[pos+anchorLen] = 'C', 'C'
		copy(target[pos+i*gap:], query[pos:pos+anchorLen])
	}
	return query, target
}

// TestSeqCompareAdaptiveMinScore covers both sides of every prefix heuristic
// boundary. A one-base gap puts two minimal anchors just below 2*MinPrefix.
func TestSeqCompareAdaptiveMinScore(t *testing.T) {
	for _, base := range []uint8{11, 13} {
		for _, boundary := range []struct{ length, extra int }{
			{9999, 0}, {10000, 2}, {49999, 2}, {50000, 4},
			{249999, 4}, {250000, 6}, {999999, 6}, {1000000, 8},
		} {
			for _, reverse := range []bool{false, true} {
				t.Run(fmt.Sprintf("prefix=%d/length=%d/reverse=%t", base, boundary.length, reverse), func(t *testing.T) {
					opt := DefaultSeqComparatorOptions
					opt.K, opt.MinPrefix = 31, base
					opt.Chaining2Options = DefaultChaining2Options
					before := opt
					pool := &sync.Pool{New: func() interface{} { return NewChainer2(&opt.Chaining2Options) }}
					for _, gap := range []int{0, 1} {
						q, target := sparsePrefixSequences(int(base)+boundary.extra, gap, boundary.length)
						if reverse {
							RC(q)
							RC(target)
						}
						cpr := NewSeqComparator(&opt, pool)
						if err := cpr.Index(q); err != nil {
							t.Fatal(err)
						}
						result, err := cpr.Compare(0, uint32(len(q)), target, len(q))
						cpr.RecycleIndex()
						if err != nil {
							t.Fatal(err)
						}
						got := result != nil
						if result != nil {
							if result.BestChain.NAnchors != 2 {
								t.Errorf("gap %d: got %d anchors, want 2", gap, result.BestChain.NAnchors)
							}
							RecycleSeqComparatorResult(result)
						}
						if got != (gap == 0) {
							t.Errorf("gap %d: retained=%t, want %t", gap, got, gap == 0)
						}
					}
					if opt != before {
						t.Fatal("comparison mutated shared options")
					}
				})
			}
		}
	}
}

// TestSeqCompareConcurrentMinScore alternates thresholds on the same read-only
// query index and shared chainer pool, including short calls after long calls.
func TestSeqCompareConcurrentMinScore(t *testing.T) {
	opt := DefaultSeqComparatorOptions
	opt.K = 31
	opt.Chaining2Options = DefaultChaining2Options
	before := opt
	q, short := sparsePrefixSequences(13, 1, 200)
	_, long := sparsePrefixSequences(13, 1, 10000)
	cpr := NewSeqComparator(&opt, &sync.Pool{New: func() interface{} { return NewChainer2(&opt.Chaining2Options) }})
	if err := cpr.Index(q); err != nil {
		t.Fatal(err)
	}
	defer cpr.RecycleIndex()
	var wg sync.WaitGroup
	for range 8 {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for range 4 {
				for _, target := range [][]byte{short, long, short} {
					result, err := cpr.Compare(0, uint32(len(q)), target, len(q))
					got := result != nil
					if result != nil {
						RecycleSeqComparatorResult(result)
					}
					if err != nil || got != (len(target) == len(short)) {
						t.Errorf("target length %d: retained=%t, error=%v", len(target), got, err)
						return
					}
				}
			}
		}()
	}
	wg.Wait()
	if opt != before {
		t.Fatal("comparison mutated shared options")
	}
}

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

	cr, err := cpr.Compare(0, uint32(len(s2)), s2, len(s2))
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
							MaxGap: 20, MinAlignLen: 50, MinIdentity: 70, BandBase: 100, BandCount: 50, HeuristicKmerPidentThreshold: 15}, MinIdentity: 70}
						cpr := NewSeqComparator(&opt, &sync.Pool{New: func() interface{} { return NewChainer2(&opt.Chaining2Options) }})
						if e := cpr.Index(q); e != nil {
							t.Fatal(e)
						}
						defer cpr.RecycleIndex()
						cr, e := cpr.Compare(0, uint32(len(q)), target, len(q))
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

// TestSeqCompareFirstAnchor retains the longest unlinked match on either strand.
func TestSeqCompareFirstAnchor(t *testing.T) {
	for _, reverse := range []bool{false, true} {
		t.Run(fmt.Sprintf("reverse=%t", reverse), func(t *testing.T) {
			opt := DefaultSeqComparatorOptions
			opt.K = 31
			opt.Chaining2Options = DefaultChaining2Options
			opt.Chaining2Options.MinAlignLen = 20
			query := []byte(indexTestDNA(600))
			target := bytes.Repeat([]byte{'A'}, 600)
			// Bound both exact matches and keep their diagonals too far apart to chain.
			query[19], query[51], query[199], query[211] = 'C', 'C', 'C', 'C'
			copy(target[20:], query[20:51])
			copy(target[450:], query[200:211])
			start := 20
			if reverse {
				RC(query)
				RC(target)
				start = 549
			}
			comparator := NewSeqComparator(&opt, &sync.Pool{New: func() interface{} {
				return NewChainer2(&opt.Chaining2Options)
			}})
			if err := comparator.Index(query); err != nil {
				t.Fatal(err)
			}
			defer comparator.RecycleIndex()
			result, err := comparator.Compare(0, 600, target, 600)
			if err != nil || result == nil {
				t.Fatalf("lost the 31-base exact match: %v", err)
			}
			defer RecycleSeqComparatorResult(result)
			want := []denseChainSummary{{1, 31, 31, 31, start, start + 30, start, start + 30, 100}}
			if got := denseChain2Summary(result.Chains); !reflect.DeepEqual(got, want) {
				t.Fatalf("chains = %+v, want %+v", got, want)
			}
		})
	}
}
