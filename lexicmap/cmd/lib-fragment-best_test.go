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
	"cmp"
	"fmt"
	"math"
	"math/rand"
	"slices"
	"testing"
)

// referenceFragmentAlignmentOrder is the original reciprocal-best comparator.
func referenceFragmentAlignmentOrder(a, b *Chain2Result) int {
	if d := cmp.Compare(b.Evalue, a.Evalue); d != 0 {
		return d
	}
	if d := cmp.Compare(a.Score, b.Score); d != 0 {
		return d
	}
	return cmp.Compare(a.BitScore, b.BitScore)
}

func TestSelectBestFragmentAlignmentMatchesSort(t *testing.T) {
	rng := rand.New(rand.NewSource(19))
	for _, n := range []int{0, 1, 2, 5, 100, 1000} {
		chains := make([]*Chain2Result, n)
		for i := range chains {
			score := float64(rng.Intn(5))
			if i%17 == 0 {
				score = math.NaN()
			}
			chains[i] = &Chain2Result{Evalue: score, Score: i / 5, BitScore: i % 5}
		}
		for range 20 {
			rng.Shuffle(n, func(i, j int) { chains[i], chains[j] = chains[j], chains[i] })
			want := slices.Clone(chains)
			slices.SortFunc(want, referenceFragmentAlignmentOrder)
			selectBestFragmentAlignment(chains)
			if n > 0 && chains[0] != want[0] {
				t.Fatalf("n=%d: different best fragment", n)
			}
			// Every losing result must still have one owner for pool recycling.
			got := slices.Clone(chains)
			slices.SortFunc(got, referenceFragmentAlignmentOrder)
			if !slices.Equal(got, want) {
				t.Fatalf("n=%d: lost or duplicated a result", n)
			}
		}
	}
}

func TestChain2ResetClearsAlignmentMetadata(t *testing.T) {
	buffer := make([]byte, 10, 64)
	c := &Chain2Result{NAnchors: 9, Score: 11, BitScore: 29, Evalue: 3.5,
		AlignedFraction: 92, AlignedLength: 100, Gaps: 3, PIdentAdjusted: 91,
		QBegin: 7, QEnd: 108, TBegin: 13, TEnd: 114, MaxExtLen: 17,
		CIGAR: buffer, QSeq: buffer, TSeq: buffer, Alignment: buffer}
	c.Reset()
	if c.Score != 0 || c.BitScore != 0 || c.Evalue != 0 || c.AlignedFraction != 0 ||
		c.AlignedLength != 0 || c.Gaps != 0 || c.PIdentAdjusted != 0 ||
		c.QBegin != 0 || c.QEnd != 0 || c.TBegin != 0 || c.TEnd != 0 || c.MaxExtLen != 0 || c.NAnchors != 0 {
		t.Fatal("reset retained another alignment's metadata")
	}
	for _, s := range [][]byte{c.CIGAR, c.QSeq, c.TSeq, c.Alignment} {
		if len(s) != 0 || cap(s) != cap(buffer) || &s[:cap(s)][0] != &buffer[0] {
			t.Fatal("reset lost reusable output storage")
		}
	}
}

func BenchmarkSelectBestFragmentAlignment(b *testing.B) {
	for _, n := range []int{5, 100, 1000} {
		input := make([]*Chain2Result, n)
		for i := range input {
			input[i] = &Chain2Result{Evalue: float64((i * 37) % 113), Score: i / 5, BitScore: i % 5}
		}
		for _, linear := range []bool{false, true} {
			b.Run(fmt.Sprintf("n=%d/linear=%t", n, linear), func(b *testing.B) {
				work := make([]*Chain2Result, n)
				b.ReportAllocs()
				b.ResetTimer()
				for range b.N {
					copy(work, input)
					if linear {
						selectBestFragmentAlignment(work)
					} else {
						slices.SortFunc(work, referenceFragmentAlignmentOrder)
					}
				}
			})
		}
	}
}
