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

package util

import (
	"bytes"
	"fmt"
	"testing"

	"github.com/shenwei356/kmers"
)

func TestIsLowComplexityDust(t *testing.T) {
	mer := []byte("TAAAAATACCTCAAAAAGAATAAAAATCCCG")
	k := len(mer)

	code, err := kmers.Encode(mer)
	if err != nil {
		t.Error(err)
		return
	}

	fmt.Printf("%s, low-complexity: %v\n", mer, IsLowComplexityDust(code, uint8(k)))
}

// TestIsLowComplexityDustMatchesTripletScore covers every supported k-mer
// length, repetitive and random codes, and the pooled implementation.
func TestIsLowComplexityDustMatchesTripletScore(t *testing.T) {
	var code uint64 = 0x9e3779b97f4a7c15
	for k := uint8(3); k <= 32; k++ {
		for sample := 0; sample < 1000; sample++ {
			code ^= code << 13
			code ^= code >> 7
			code ^= code << 17
			value := code
			if sample%10 == 0 {
				value = 0 // homopolymer
			} else if sample%10 == 1 {
				value = 0x5555555555555555 // alternating bases
			}
			var counts [64]int
			// A k-mer contains exactly k-2 complete triplets.
			for i := uint8(0); i < k-2; i++ {
				counts[(value>>(2*i))&63]++
			}
			score := 0
			for _, count := range counts {
				score += count * (count - 1) / 2
			}
			got, want := IsLowComplexityDust(value, k), score > 50
			if got != want {
				t.Fatalf("k=%d sample=%d code=%x: got %t, want %t", k, sample, value, got, want)
			}
			if slow := IsLowComplexityDustSlow(value, k); got != slow {
				t.Fatalf("k=%d sample=%d code=%x: fast=%t, slow=%t", k, sample, value, got, slow)
			}
		}
	}
}

// DUST counts complete triplets only, so reversing and complementing a k-mer
// must preserve its score and filtering decision.
func TestIsLowComplexityDustReverseComplement(t *testing.T) {
	for k := uint8(0); k < 3; k++ {
		if IsLowComplexityDust(0, k) || IsLowComplexityDustSlow(0, k) {
			t.Errorf("length %d has no complete triplets", k)
		}
	}
	for _, sequence := range []string{
		"AAAATCAAAATAAACAAAACCATTAAACAAA", // score 45, but a padded AAA raises it to 54
		"AAAAAAAAAAAA",                    // 10 complete triplets, score 45
		"AAAAAAAAAAAAA",                   // 11 complete triplets, score 55
	} {
		bases := []byte(sequence)
		var counts [64]int
		for i := 0; i+3 <= len(bases); i++ {
			triplet, err := kmers.Encode(bases[i : i+3])
			if err != nil {
				t.Fatal(err)
			}
			counts[triplet]++
		}
		score := 0
		for _, count := range counts {
			score += count * (count - 1) / 2
		}
		for _, reverse := range []bool{false, true} {
			code, err := kmers.Encode(bases)
			if err != nil {
				t.Fatal(err)
			}
			for name, filter := range map[string]func(uint64, uint8) bool{"fast": IsLowComplexityDust, "slow": IsLowComplexityDustSlow} {
				if got := filter(code, uint8(len(bases))); got != (score > 50) {
					t.Errorf("%s reverse=%t length=%d: got %t, complete-triplet score %d", name, reverse, len(bases), got, score)
				}
			}
			// Reverse complement the synthetic sequence for the second pass.
			for left, right := 0, len(bases)-1; left <= right; left, right = left+1, right-1 {
				complement := func(base byte) byte { return "TGCA"[bytes.IndexByte([]byte("ACGT"), base)] }
				bases[left], bases[right] = complement(bases[right]), complement(bases[left])
			}
		}
	}
}

func TestNs(t *testing.T) {
	var k uint8 = 5
	values := []uint64{
		Ns(0b00, k), // A
		Ns(0b01, k), // C
		Ns(0b10, k), // G
		Ns(0b11, k), // T
	}
	for _, v := range values {
		fmt.Printf("%s, %064b\n", kmers.MustDecode(v, int(k)), v)
	}
}
