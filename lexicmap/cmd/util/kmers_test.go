// Copyright © 2023-2024 Wei Shen <shenwei356@gmail.com>
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
	mer := []byte("AAAAAAAACCGGGCAATTGCCCGGTGCTGGA")
	k := len(mer)

	code, err := kmers.Encode(mer)
	if err != nil {
		t.Error(err)
		return
	}

	fmt.Printf("%s, low-complexity: %v\n", mer, IsLowComplexityDust(code, uint8(k)))
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
