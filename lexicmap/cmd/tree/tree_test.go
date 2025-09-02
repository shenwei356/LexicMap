// Copyright © 2023-2025 Wei Shen <shenwei356@gmail.com>
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

package tree

import (
	"testing"

	"github.com/shenwei356/kmers"
	"github.com/shenwei356/lexichash"
)

func TestTree(t *testing.T) {
	var k uint8
	var n uint64
	var i uint64
	var v uint32
	var query string
	var code uint64
	var srs *[]*SearchResult

	for j := 0; j < 1000; j++ {
		k = 6
		n = uint64(1 << (k * 2))

		_t := NewTree(uint8(k))

		for i = 0; i < n; i++ {
			v = uint32(i & 3)
			if v == 3 || v == 0 {
				continue
			}
			_t.Insert(i, v)
		}

		query = "ACTGAC"
		code, _ = kmers.Encode([]byte(query))
		// srs, _ := tree.Search(code, uint8(len(query)), 4)
		srs, _ = _t.Search(code, 5)
		t.Logf("query: %s\n", query)
		for _, sr := range *srs {
			t.Logf("  %s, len(prefix): %d, %v\n",
				lexichash.MustDecode(sr.Kmer, k), sr.LenPrefix, sr.Values)
		}
		_t.RecycleSearchResult(srs)

		RecycleTree(_t)
	}

}
