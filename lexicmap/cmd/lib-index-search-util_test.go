package cmd

import (
	"bytes"
	"encoding/json"
	"fmt"
	"math/rand"
	"os"
	"reflect"
	"strings"
	"sync"
	"testing"

	"github.com/shenwei356/wfa"
)

// TestExtendRightFirstAnchor retains the sole matching terminal 2-mer.
func TestExtendRightFirstAnchor(t *testing.T) {
	if q, s := _extendRight([]byte("ACA"), []byte("ACC")); q != 2 || s != 2 {
		t.Fatalf("extension = (%d, %d), want (2, 2)", q, s)
	}
}

// TestExtendMatchFirstAnchor checks terminal matches within sequence and contig bounds on both strands.
func TestExtendMatchFirstAnchor(t *testing.T) {
	for _, tc := range []struct {
		name                       string
		seq1, seq2                 string
		start1, end1, start2, end2 int
		leftLimit, rightLimit      int
		want                       string
		ext                        [4]int // Query left/right, subject left/right extension lengths.
	}{
		{
			name: "right",
			seq1: "GGACA", seq2: "GGACC",
			start1: 0, end1: 2, start2: 0, end2: 2,
			rightLimit: 50,
			want:       "GGAC", ext: [4]int{0, 2, 0, 2},
		},
		{
			name: "left",
			seq1: "ACAGG", seq2: "CCAGG",
			start1: 3, end1: 5, start2: 3, end2: 5,
			leftLimit: 50,
			want:      "CAGG", ext: [4]int{2, 0, 2, 0},
		},
		{
			name: "both clipped to sequence",
			seq1: "ACAGGACA", seq2: "CCAGGACC",
			start1: 3, end1: 5, start2: 3, end2: 5,
			leftLimit: 50, rightLimit: 50,
			want: "CAGGAC", ext: [4]int{2, 2, 2, 2},
		},
		{
			name: "right contig boundary",
			seq1: "GGACA", seq2: "GGACC",
			start1: 0, end1: 2, start2: 0, end2: 2,
			rightLimit: 2,
			want:       "GG",
		},
		{
			name: "left contig boundary",
			seq1: "ACAGG", seq2: "CCAGG",
			start1: 3, end1: 5, start2: 3, end2: 5,
			leftLimit: 2,
			want:      "GG",
		},
		{
			name: "two bases remaining",
			seq1: "GGAC", seq2: "GGAC",
			start1: 0, end1: 2, start2: 0, end2: 2,
			leftLimit: 50, rightLimit: 50,
			want: "GG",
		},
		{
			name: "different query and subject coordinates",
			seq1: "TTGGACA", seq2: "GGACC",
			start1: 2, end1: 4, start2: 0, end2: 2,
			rightLimit: 50,
			want:       "GGAC", ext: [4]int{0, 2, 0, 2},
		},
	} {
		for _, rc := range []bool{false, true} {
			t.Run(fmt.Sprintf("%s/rc=%t", tc.name, rc), func(t *testing.T) {
				tBegin, maxExtLen := tc.leftLimit, tc.rightLimit
				if rc { // The subject's contig limits swap in reverse-complement orientation.
					tBegin, maxExtLen = maxExtLen, tBegin
				}
				q, s, qLeft, qRight, sLeft, sRight, err := extendMatch(
					[]byte(tc.seq1), []byte(tc.seq2), tc.start1, tc.end1, tc.start2, tc.end2,
					50, tBegin, maxExtLen, rc,
				)
				if err != nil {
					t.Fatal(err)
				}
				if got := [4]int{qLeft, qRight, sLeft, sRight}; got != tc.ext {
					t.Fatalf("extension = %v, want %v", got, tc.ext)
				}
				if string(q) != tc.want || string(s) != tc.want {
					t.Fatalf("extended sequences = (%q, %q), want (%q, %q)", q, s, tc.want, tc.want)
				}
			})
		}
	}
}

// Repeat-window fixtures reproduce the original genome-comparison failure.
const repeatWindowQuery = "CCTGTTGCTCCTGTTGCTCCTGTTGGACCTGTTGCTCCATCTGCTCCTGTATTTCCTGTTGCTCCTGTTGGACCTGTTGCTCCATCTGCTCCTGTATTTCCTGTTGCTCCTGTTGCTCCTGTTGGACCTGCTACTCCATCTGCTCCTGTATTTCCTGTTGCTCCTGTTGGACCTGTTGCTCCATCTGCTCCTGTTGGACCTGTCGCTCCTGTTGCCCCTGTTGCTCCTGTTATTCCTACTCCTGTTGGGCCTGTTGCTCCTGTTGGACCTGTTGCTCCAGGTGTTCCTGCTGCTCCTGTTGCTCCTGTTGGGCCTACTAATCCATTTGCTCCTGTTGCCCCTGTTGCTCCTGTTGGACCTGTTACTCCTGTTGCTCCGTCTGCTCCTGTTGCCCCTGTTGCTCCTGTTGCTCCTGTT"
const repeatWindowTarget = "CCTGTTGCTCCTGTTGCTCCTGTTGCTCCTGTTGGGCCTGTTGCTCCTGTATTTCCTGTTGCTCCTGTTGCTCCTGTTGGACCTGCTACTCCATCTGCTCCTGTATTTCCTGTTGCTCCTGTTGGACCTGTTGCTCCATCTGCTCCTGTTGGACCTGTCGCTCCTGTTGCCCCTGTTGCTCCTGTTATTCCTACTCCTGTTGGGCCTGTTGCTCCTGTTGGACCTGTTGCTCCAGGTGTTCCTGCTGCTCCTGTTGCTCCTGTTGGGCCTACTAATCCATTTGCTCCTGTTGCCCCTGTTGCTCCTGTTGGACCTGTTACTCCTGTTGCTCCGTCTGCTCCTGTTGCCCCTGTT"

// TestRepeatRefinementCoordinates covers both clipping directions and checks
// that full/trimmed CIGAR rendering still consumes the original inputs.
func TestRepeatRefinementCoordinates(t *testing.T) {
	for _, adaptive := range []bool{false, true} {
		for _, reverse := range []bool{false, true} {
			for _, swap := range []bool{false, true} {
				t.Run(fmt.Sprintf("adaptive=%t/reverse=%t/swap=%t", adaptive, reverse, swap), func(t *testing.T) {
					q, target := []byte(repeatWindowQuery), []byte(repeatWindowTarget)
					if reverse {
						RC(q)
						RC(target)
					}
					if swap {
						q, target = target, q
					}
					algn := wfa.New(wfa.DefaultPenalties, &wfa.Options{GlobalAlignment: true})
					if adaptive {
						algn.AdaptiveReduction(wfa.DefaultAdaptiveOption)
					}
					defer wfa.RecycleAligner(algn)
					score := scoreAndEvalue(2, -3, 5, 2, len(target), 0.625, 0.41)
					baseline, err := algn.Align(q, target)
					if err != nil {
						t.Fatal(err)
					}
					defer wfa.RecycleAlignmentResult(baseline)
					r, err := alignWithRepeatRefinement(algn, q, target, 31, len(q), 0, score)
					if err != nil {
						t.Fatal(err)
					}
					defer wfa.RecycleAlignmentResult(r)
					t.Logf("old=%d/%d gaps=%d new=%d/%d gaps=%d coords=%d..%d/%d..%d", baseline.Matches, baseline.AlignLen, baseline.Gaps, r.Matches, r.AlignLen, r.Gaps, r.QBegin, r.QEnd, r.TBegin, r.TEnd)
					if r.Matches != 344 || r.AlignLen != 354 || r.Gaps != 0 || r.Score != 40 {
						t.Fatalf("unexpected refined result: %+v", r)
					}
					// Reverse both sequences, then swap query/target if requested.
					// Coordinates and operations must describe the same local alignment.
					bounds := [4]int{46, 399, 1, 354}
					cigar := "5M3X17M2X7M2X2M3X313M"
					if reverse {
						bounds[0], bounds[1] = len(repeatWindowQuery)-bounds[1]+1, len(repeatWindowQuery)-bounds[0]+1
						cigar = "313M3X2M2X7M2X17M3X5M"
					}
					if swap {
						bounds = [4]int{bounds[2], bounds[3], bounds[0], bounds[1]}
					}
					if got := [4]int{r.QBegin, r.QEnd, r.TBegin, r.TEnd}; got != bounds || r.CIGAR(true) != cigar {
						t.Fatalf("bounds/CIGAR = %v/%s, want %v/%s", got, r.CIGAR(true), bounds, cigar)
					}
					if reward, _, _ := score(len(q), r); reward != 658 {
						t.Fatalf("alignment reward = %d, want 658", reward)
					}
					for _, trim := range []bool{false, true} {
						Q, A, T := r.AlignmentText(&q, &target, trim)
						qWant, tWant := q, target
						if trim {
							qWant = q[r.QBegin-1 : r.QEnd]
							tWant = target[r.TBegin-1 : r.TEnd]
						}
						if !bytes.Equal(bytes.ReplaceAll(*Q, []byte{'-'}, nil), qWant) || !bytes.Equal(bytes.ReplaceAll(*T, []byte{'-'}, nil), tWant) || len(*Q) != len(*T) || len(*A) != len(*Q) {
							t.Fatalf("trim=%t: CIGAR/text/coordinates disagree", trim)
						}
						wfa.RecycleAlignmentText(Q, A, T)
					}
					if strings.ContainsAny(r.CIGAR(true), "IDH") {
						t.Fatalf("clipping leaked into aligned CIGAR: %s", r.CIGAR(true))
					}
				})
			}
		}
	}
}

// TestRepeatRefinementPreservesSmallIndels compares every field and operation
// with ordinary WFA for known 1-25 bp indels and larger 30/50/100 bp boundary cases in two backgrounds.
func TestRepeatRefinementPreservesSmallIndels(t *testing.T) {
	rng := rand.New(rand.NewSource(281028))
	randomDNA := func(n int) []byte {
		s := make([]byte, n)
		for i := range s {
			s[i] = "ACGT"[rng.Intn(4)]
		}
		return s
	}
	for _, repetitive := range []bool{false, true} {
		base := randomDNA(300)
		if repetitive {
			copy(base[80:220], strings.Repeat("ACGT", 35))
		}
		// Chaining's default 20-base limit is not a WFA gap-length bound.
		// Include every size through 25 and several larger indels.
		gapSizes := make([]int, 25)
		for i := range gapSizes {
			gapSizes[i] = i + 1
		}
		gapSizes = append(gapSizes, 30, 50, 100)
		for _, size := range gapSizes {
			inserted := randomDNA(size)
			for _, position := range []int{40, 150, 260} {
				q := append([]byte(nil), base[:position]...)
				q = append(q, inserted...)
				q = append(q, base[position:]...)
				// Keep a few substitutions so equality is not limited to perfect flanks.
				for _, p := range []int{17, 213} {
					q[p] = "ACGT"[(strings.IndexByte("ACGT", q[p])+1)%4]
				}
				for _, swap := range []bool{false, true} {
					a, b := q, base
					if swap {
						a, b = b, a
					}
					for _, adaptive := range []bool{false, true} {
						t.Run(fmt.Sprintf("repeat=%t/gap=%d/position=%d/swap=%t/adaptive=%t", repetitive, size, position, swap, adaptive), func(t *testing.T) {
							algn := wfa.New(wfa.DefaultPenalties, &wfa.Options{GlobalAlignment: true})
							defer wfa.RecycleAligner(algn)
							if adaptive {
								algn.AdaptiveReduction(wfa.DefaultAdaptiveOption)
							}
							baseline, e := algn.Align(a, b)
							if e != nil {
								t.Fatal(e)
							}
							defer wfa.RecycleAlignmentResult(baseline)
							if baseline.Gaps == 0 {
								t.Fatal("fixture did not retain an internal indel")
							}
							score := scoreAndEvalue(2, -3, 5, 2, len(b), 0.625, 0.41)
							got, e := alignWithRepeatRefinement(algn, a, b, 31, len(a), 0, score)
							if e != nil {
								t.Fatal(e)
							}
							defer wfa.RecycleAlignmentResult(got)
							if !reflect.DeepEqual(got, baseline) {
								t.Fatalf("small indel changed: old=%+v, new=%+v", baseline, got)
							}
						})
					}
				}
			}
		}
	}
}

// TestRepeatRefinementPreservesTerminalIndels retains short indels near either
// end. One-sided clipping can improve reward by trading an indel for a mismatch,
// but reverses to the zero-offset window that refinement deliberately skips.
func TestRepeatRefinementPreservesTerminalIndels(t *testing.T) {
	for _, tc := range []struct {
		query, target string
		gap           uint32
	}{
		{"AAACC", "AAAC", 1},
		{"AGAAAATATATT", "AGAATATATT", 2},
		{"AGAAAAATATATT", "AGAATATATT", 3},
	} {
		for _, reverse := range []bool{false, true} {
			for _, swap := range []bool{false, true} {
				for _, adaptive := range []bool{false, true} {
					t.Run(fmt.Sprintf("gap=%d/reverse=%t/swap=%t/adaptive=%t", tc.gap, reverse, swap, adaptive), func(t *testing.T) {
						q, target := []byte(tc.query+indexTestDNA(100)), []byte(tc.target+indexTestDNA(100))
						if reverse {
							RC(q)
							RC(target)
						}
						if swap {
							q, target = target, q
						}
						algn := wfa.New(wfa.DefaultPenalties, &wfa.Options{GlobalAlignment: true})
						defer wfa.RecycleAligner(algn)
						if adaptive {
							algn.AdaptiveReduction(wfa.DefaultAdaptiveOption)
						}
						baseline, err := algn.Align(q, target)
						if err != nil {
							t.Fatal(err)
						}
						defer wfa.RecycleAlignmentResult(baseline)
						if baseline.Gaps != tc.gap {
							t.Fatalf("fixture gap bases = %d, want %d", baseline.Gaps, tc.gap)
						}
						score := scoreAndEvalue(2, -3, 5, 2, len(target), 0.625, 0.41)
						got, err := alignWithRepeatRefinement(algn, q, target, 31, len(q), 0, score)
						if err != nil {
							t.Fatal(err)
						}
						defer wfa.RecycleAlignmentResult(got)
						if !reflect.DeepEqual(got, baseline) {
							t.Fatalf("terminal indel changed: old=%s, new=%s", baseline.CIGAR(false), got.CIGAR(false))
						}
					})
				}
			}
		}
	}
}

// TestRepeatRefinementCoverageFallback keeps an existing hit if the shorter
// candidate would fail the caller's query-span coverage threshold.
func TestRepeatRefinementCoverageFallback(t *testing.T) {
	q, target := []byte(repeatWindowQuery), []byte(repeatWindowTarget)
	algn := wfa.New(wfa.DefaultPenalties, &wfa.Options{GlobalAlignment: true})
	defer wfa.RecycleAligner(algn)
	algn.AdaptiveReduction(wfa.DefaultAdaptiveOption)
	score := scoreAndEvalue(2, -3, 5, 2, len(target), 0.625, 0.41)
	baseline, err := algn.Align(q, target)
	if err != nil {
		t.Fatal(err)
	}
	defer wfa.RecycleAlignmentResult(baseline)
	got, err := alignWithRepeatRefinement(algn, q, target, 31, len(q), 90, score)
	if err != nil {
		t.Fatal(err)
	}
	defer wfa.RecycleAlignmentResult(got)
	if !reflect.DeepEqual(got, baseline) {
		t.Fatal("refinement discarded a hit satisfying 90% query coverage")
	}
}

// TestRepeatRefinement16S preserves the actual GCF_001027105.1 search window,
// including its small indels and full-length coverage under the old chaining.
func TestRepeatRefinement16S(t *testing.T) {
	data, err := os.ReadFile("testdata/seq-compare-16s.json")
	if err != nil {
		t.Fatal(err)
	}
	var fixture struct{ Query, Target string }
	if err = json.Unmarshal(data, &fixture); err != nil {
		t.Fatal(err)
	}
	q, target := []byte(fixture.Query), []byte(fixture.Target)
	opt := SeqComparatorOptions{K: 31, MinPrefix: 11, Chaining2Options: Chaining2Options{
		MaxGap: 20, MinScore: 35, MinAlignLen: 50, MinIdentity: 70, BandBase: 100, BandCount: 50, HeuristicKmerPidentThreshold: 15}, MinIdentity: 70}
	cpr := NewSeqComparator(&opt, &sync.Pool{New: func() interface{} { return NewChainer2(&opt.Chaining2Options) }})
	if err := cpr.Index(q); err != nil {
		t.Fatal(err)
	}
	defer cpr.RecycleIndex()
	cr, err := cpr.Compare(0, uint32(len(q)-1), target, len(q))
	if err != nil || cr == nil {
		t.Fatalf("comparison: result=%v err=%v", cr, err)
	}
	defer RecycleSeqComparatorResult(cr)
	c := cr.BestChain
	a, b, qLeft, qRight, _, _, err := extendMatch(q, target, c.QBegin, c.QEnd+1, c.TBegin, c.TEnd+1, 50, c.TBegin, len(target)-1-c.TEnd, false)
	if err != nil {
		t.Fatal(err)
	}
	algn := wfa.New(wfa.DefaultPenalties, &wfa.Options{GlobalAlignment: true})
	defer wfa.RecycleAligner(algn)
	algn.AdaptiveReduction(wfa.DefaultAdaptiveOption)
	baseline, err := algn.Align(a, b)
	if err != nil {
		t.Fatal(err)
	}
	defer wfa.RecycleAlignmentResult(baseline)
	score := scoreAndEvalue(2, -3, 5, 2, len(target), 0.625, 0.41)
	got, err := alignWithRepeatRefinement(algn, a, b, 31, len(q), 0, score)
	if err != nil {
		t.Fatal(err)
	}
	defer wfa.RecycleAlignmentResult(got)
	if !reflect.DeepEqual(got, baseline) {
		t.Fatal("16S alignment changed")
	}
	qb, qe := c.QBegin-qLeft+got.QBegin-1, c.QEnd+qRight-(len(a)-got.QEnd)
	if qb != 4 || qe != 1540 || got.AlignLen != 1568 || got.Matches != 1221 {
		t.Fatalf("16S coverage changed: q=%d..%d result=%+v", qb, qe, got)
	}
}
