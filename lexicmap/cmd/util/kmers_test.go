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
// length, repetitive and random codes, and the retained original implementation.
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
			// Preserve the current DUST window, including its final triplet.
			for i := uint8(0); i <= k-2; i++ {
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
