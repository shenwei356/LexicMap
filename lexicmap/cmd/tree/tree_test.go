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

package tree

import (
	"math/bits"
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

func TestWalkGroupsAndPairs(t *testing.T) {
	const k uint8 = 6
	const p uint8 = 4

	_t := NewTree(k)
	defer RecycleTree(_t)

	n := uint64(1) << (k * 2)
	inserted := make([]uint64, 0, n)
	for i := uint64(0); i < n; i++ {
		v := uint32(i & 3)
		if v == 3 || v == 0 {
			continue
		}
		_t.Insert(i, v)
		inserted = append(inserted, i)
	}

	shift := int(k) - 32
	lcp := func(a, b uint64) uint8 {
		return uint8(bits.LeadingZeros64(a^b)>>1 + shift)
	}

	expectedPairs := make(map[[2]uint64]uint8)
	for i := 0; i < len(inserted); i++ {
		for j := i + 1; j < len(inserted); j++ {
			a, b := inserted[i], inserted[j]
			if a > b {
				a, b = b, a
			}
			if l := lcp(a, b); l >= p {
				expectedPairs[[2]uint64{a, b}] = l
			}
		}
	}

	// -------------------------------------------------------------------------------

	gotKeys := make(map[uint64]bool)
	var groupCount int
	_t.WalkGroups(p, func(keys []uint64, vals [][]uint32, lenPrefix uint8) bool {
		groupCount++
		if len(keys) != len(vals) {
			t.Fatalf("group: len(keys)=%d != len(vals)=%d", len(keys), len(vals))
		}
		if lenPrefix < p {
			t.Fatalf("group lenPrefix %d < p %d", lenPrefix, p)
		}
		for i := 0; i < len(keys); i++ {
			for j := i + 1; j < len(keys); j++ {
				if l := lcp(keys[i], keys[j]); l < lenPrefix {
					t.Fatalf("intra-group LCP %d < reported %d for %x vs %x",
						l, lenPrefix, keys[i], keys[j])
				}
			}
			if gotKeys[keys[i]] {
				t.Fatalf("k-mer %x reported in multiple groups", keys[i])
			}
			gotKeys[keys[i]] = true
		}
		return false
	})
	if len(gotKeys) != len(inserted) {
		t.Fatalf("WalkGroups covered %d k-mers, expected %d", len(gotKeys), len(inserted))
	}
	t.Logf("WalkGroups: %d groups for %d k-mers (p=%d)", groupCount, len(inserted), p)

	// -------------------------------------------------------------------------------

	gotPairs := make(map[[2]uint64]uint8)
	_t.WalkPairs(p, func(a, b uint64, va, vb []uint32, lenPrefix uint8) bool {
		if a > b {
			a, b = b, a
		}
		if _, dup := gotPairs[[2]uint64{a, b}]; dup {
			t.Fatalf("pair (%x,%x) reported twice", a, b)
		}
		gotPairs[[2]uint64{a, b}] = lenPrefix
		return false
	})

	if len(gotPairs) != len(expectedPairs) {
		t.Fatalf("WalkPairs got %d pairs, expected %d", len(gotPairs), len(expectedPairs))
	}
	for pair, want := range expectedPairs {
		got, ok := gotPairs[pair]
		if !ok {
			t.Fatalf("missing pair (%x,%x)", pair[0], pair[1])
		}
		if got != want {
			t.Fatalf("pair (%x,%x) lenPrefix: got %d, want %d", pair[0], pair[1], got, want)
		}
	}

	// -------------------------------------------------------------------------------

	var earlyCount int
	aborted := _t.WalkPairs(p, func(a, b uint64, va, vb []uint32, lenPrefix uint8) bool {
		earlyCount++
		return earlyCount >= 3
	})
	if !aborted {
		t.Fatalf("WalkPairs should report aborted=true on early return")
	}
	if earlyCount != 3 {
		t.Fatalf("early-abort count = %d, expected 3", earlyCount)
	}
}
