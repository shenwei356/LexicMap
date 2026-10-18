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
	"math"
	"math/rand"
	"slices"
	"testing"
)

// Compare pointer order, not just IDs: groups follow their first HSP and all
// score ties and distinct HSP objects within a group must keep their order.
func TestSortBySeqIDMatchesReference(t *testing.T) {
	cases := [][][]byte{
		nil, {[]byte("one")},
		{nil, {}, []byte("b"), nil, []byte("a"), []byte("b"), {}},
		{[]byte("b"), []byte("a"), []byte("b"), []byte("c"), []byte("a")},
		{[]byte("b"), []byte("b"), []byte("a"), []byte("a")},
		{[]byte{0, 255}, []byte{0}, []byte{0, 255}, []byte{255}},
	}
	rng := rand.New(rand.NewSource(20261017))
	for trial := 0; trial < 200; trial++ {
		ids := make([][]byte, rng.Intn(1000))
		groups := 1 + rng.Intn(200)
		for i := range ids {
			// Equal IDs deliberately use different backing arrays.
			ids[i] = []byte(fmt.Sprintf("contig_%d", rng.Intn(groups)))
		}
		cases = append(cases, ids)
	}
	for trial, ids := range cases {
		input := make([]*SimilarityDetail, len(ids))
		for i, id := range ids {
			input[i] = &SimilarityDetail{SeqID: slices.Clone(id), SimilarityScore: float64(i % 3)}
		}
		// The original scan also discards nil slots after the first HSP.
		if len(input) > 20 && trial%2 == 0 {
			input[10] = nil
		}
		gotSlice, wantSlice := slices.Clone(input), slices.Clone(input)
		got, want := &SearchResult{SimilarityDetails: &gotSlice}, &SearchResult{SimilarityDetails: &wantSlice}
		got.SortBySeqID()
		referenceSortBySeqID(want)
		if !slices.Equal(*got.SimilarityDetails, *want.SimilarityDetails) {
			t.Fatalf("trial %d: HSP ordering differs", trial)
		}
		for i, sd := range input {
			if sd != nil && (!bytes.Equal(sd.SeqID, ids[i]) || sd.SimilarityScore != float64(i%3)) {
				t.Fatalf("trial %d: HSP %d was modified", trial, i)
			}
		}
		recycleSimilarityDetailSlice(got.SimilarityDetails)
		recycleSimilarityDetailSlice(want.SimilarityDetails)
	}
}

// The reference retains the full sort and containment scan, including its
// choice of surviving strand flags when duplicate coordinates compare equal.
func TestClearSubstrPairsLongestMatchesReference(t *testing.T) {
	rng := rand.New(rand.NewSource(20261018))
	for trial := 0; trial < 2400; trial++ {
		k := []int{1, 13, 31, 255}[trial%4]
		input := make([]SubstrPair, trial%350)
		for i := range input {
			length := k
			if i%3 == 0 {
				length = rng.Intn(k + 1)
			}
			input[i] = SubstrPair{QBegin: int32(rng.Intn(60)), TBegin: int32(rng.Intn(60)), Len: uint8(length), QRC: rng.Intn(2) == 0, TRC: rng.Intn(2) == 0}
			if i > 0 && i%7 == 0 {
				input[i] = input[i-1]
				input[i].QRC = !input[i].QRC
			}
		}
		if len(input) > 3 {
			// Exercise fallback for a violated maximum length, negative starts,
			// overflowing ends, and k values outside the shortcut's range.
			switch trial % 11 {
			case 1:
				input[0].Len = 255
			case 2:
				input[0].QBegin = -1
			case 3:
				input[0].TBegin = -1
			case 4:
				input[0].QBegin = math.MaxInt32
			case 5:
				input[0].TBegin = math.MaxInt32
			case 6:
				k = 0
			case 7:
				k = 256
			case 8:
				k = -1
			}
		}
		checkClearSubstrPairsReference(t, input, k)
	}
	// Containment need not follow the same diagonal. Shifted equal-length
	// anchors, however, must survive on either sequence.
	checkClearSubstrPairsReference(t, []SubstrPair{
		{QBegin: 5, TBegin: 5, Len: 31},
		{QBegin: 10, TBegin: 12, Len: 10},
		{QBegin: 5, TBegin: 6, Len: 31},
		{QBegin: 6, TBegin: 5, Len: 31},
		{QBegin: 5, TBegin: 5, Len: 31, QRC: true},
	}, 31)
}

func checkClearSubstrPairsReference(t *testing.T, input []SubstrPair, k int) {
	t.Helper()
	got, want := slices.Clone(input), slices.Clone(input)
	capacity := cap(got)
	ClearSubstrPairs(&got, k)
	if len(want) > 1 {
		referenceClearSubstrPairs(&want, k)
	}
	if !slices.Equal(got, want) || cap(got) != capacity {
		t.Fatalf("k=%d input=%+v: got %+v, want %+v", k, input, got, want)
	}
}

// Query-ordered inputs take the group sort; reversed inputs take the full
// sort. Both must match the old equal-length cleaner, including backing capacity.
func TestClearEqualLengthSubstrPairsQueryGroups(t *testing.T) {
	rng := rand.New(rand.NewSource(20261019))
	for trial := 0; trial < 600; trial++ {
		var input []SubstrPair
		for q := 0; q < trial%40; q++ {
			for i := 0; i < 1+rng.Intn(100); i++ {
				input = append(input, SubstrPair{QBegin: int32(q * 7), TBegin: int32(rng.Intn(150)), Len: 13, QRC: trial%2 == 0, TRC: trial%3 == 0})
			}
		}
		if trial%3 == 0 {
			slices.Reverse(input)
		}
		got, want := slices.Clone(input), slices.Clone(input)
		capacity := cap(got)
		clearEqualLengthSubstrPairs(&got)
		referenceClearEqualLengthSubstrPairs(&want)
		if !slices.Equal(got, want) || cap(got) != capacity {
			t.Fatalf("trial %d: equal-length cleanup differs", trial)
		}
	}
}

// Frozen pre-optimization implementations are the differential-test oracles.
func referenceClearSubstrPairs(subs *[]SubstrPair, k int) {

	slices.SortFunc(*subs, func(a, b SubstrPair) int {
		if a.QBegin == b.QBegin {
			if a.QBegin+int32(a.Len) == b.QBegin+int32(b.Len) {
				return int(a.TBegin - b.TBegin)
			}
			return int(b.QBegin) + int(b.Len) - (int(a.QBegin) + int(a.Len))
		}
		return int(a.QBegin - b.QBegin)
	})

	var p SubstrPair
	var upbound, vQEnd, vTBegin, vTEnd int32
	var j int
	markers := poolBoolList.Get().(*[]bool)
	*markers = (*markers)[:0]

	n := len(*subs)
	if cap(*markers) >= n {
		*markers = (*markers)[:n]
		clear(*markers)
	} else {
		*markers = (*markers)[:cap(*markers)]
		clear(*markers)
		for len(*markers) < n {
			*markers = append(*markers, false)
		}
	}

	for i, v := range (*subs)[1:] {

		vQEnd = v.QBegin + int32(v.Len)
		upbound = max(vQEnd-int32(k), 0)
		vTBegin = v.TBegin
		vTEnd = v.TBegin + int32(v.Len)

		start, _ := slices.BinarySearchFunc((*subs)[:i+1], upbound, func(s SubstrPair, target int32) int {
			return int(s.QBegin - target)
		})
		for j = start; j <= i; j++ {
			p = (*subs)[j]
			if vQEnd <= p.QBegin+int32(p.Len) &&
				vTBegin >= p.TBegin && vTEnd <= p.TBegin+int32(p.Len) {
				(*markers)[i+1] = true // because: range (*subs)[1:]
				break
			}
		}
	}

	j = 0
	for i, embedded := range *markers {
		if !embedded {
			(*subs)[j] = (*subs)[i]
			j++
		}
	}
	if j > 0 {
		*subs = (*subs)[:j]
	}

	poolBoolList.Put(markers)
}

func referenceSortBySeqID(sr *SearchResult) {
	if len(*sr.SimilarityDetails) <= 1 {
		return
	}

	sds0 := sr.SimilarityDetails

	sds := poolSimilarityDetails.Get().(*[]*SimilarityDetail) // HSPs in a reference

	// first one
	*sds = append(*sds, (*sds0)[0])
	seqid := (*sds0)[0].SeqID
	(*sds0)[0] = nil

	i := 1
	var nextSeqId []byte
	var sds1 []*SimilarityDetail
	var first bool
	var j, nextI int
	var sd *SimilarityDetail
	n := len(*sds0)
	for {
		sds1 = (*sds0)[i:]
		first = true
		for j, sd = range sds1 {
			if sd == nil {
				continue
			}

			if bytes.Equal(sd.SeqID, seqid) { // has the same seqid
				*sds = append(*sds, sd)
				(*sds0)[i+j] = nil
			} else if first { // record another seqid
				first = false
				nextSeqId = sd.SeqID
				nextI = i + j
			}
		}
		if first { // no other seqids
			break
		}

		*sds = append(*sds, (*sds0)[nextI])
		(*sds0)[nextI] = nil
		seqid = nextSeqId
		i = nextI + 1
		if i >= n { // it's alreadly the last one
			break
		}
	}

	recycleSimilarityDetailSlice(sr.SimilarityDetails)
	sr.SimilarityDetails = sds
}

// These benchmarks separate grouping, query-ordered sorting, and the longest
// anchor shortcut. Both sides receive identical fresh copies of their inputs.
func BenchmarkSearchCleanup(b *testing.B) {
	rng := rand.New(rand.NewSource(42))
	for _, kind := range []string{"longest", "mixed", "shorter"} {
		input := make([]SubstrPair, 4096)
		for i := range input {
			length := uint8(31)
			if kind == "shorter" || kind == "mixed" && i%2 == 0 {
				length = 13
			}
			input[i] = SubstrPair{QBegin: int32(rng.Intn(2048)), TBegin: int32(rng.Intn(2048)), Len: length}
		}
		for _, impl := range []struct {
			name string
			fn   func(*[]SubstrPair, int)
		}{{"reference", referenceClearSubstrPairs}, {"optimized", ClearSubstrPairs}} {
			b.Run("containment/"+kind+"/"+impl.name, func(b *testing.B) {
				work := make([]SubstrPair, len(input))
				b.ReportAllocs()
				b.ResetTimer()
				for i := 0; i < b.N; i++ {
					work = work[:len(input)]
					copy(work, input)
					impl.fn(&work, 31)
				}
			})
		}
	}
	input := make([]SubstrPair, 12800)
	for i := range input {
		input[i] = SubstrPair{QBegin: int32(i / 100), TBegin: int32(rng.Intn(2000)), Len: 13}
	}
	for _, impl := range []struct {
		name string
		fn   func(*[]SubstrPair)
	}{{"reference", referenceClearEqualLengthSubstrPairs}, {"optimized", clearEqualLengthSubstrPairs}} {
		b.Run("queryGroups/"+impl.name, func(b *testing.B) {
			work := make([]SubstrPair, len(input))
			b.ReportAllocs()
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				work = work[:len(input)]
				copy(work, input)
				impl.fn(&work)
			}
		})
	}
	for _, groups := range []int{1, 1000} {
		details := make([]*SimilarityDetail, 5000)
		for i := range details {
			details[i] = &SimilarityDetail{SeqID: []byte(fmt.Sprintf("contig_%d", i%groups))}
		}
		for _, impl := range []struct {
			name string
			fn   func(*SearchResult)
		}{{"reference", referenceSortBySeqID}, {"optimized", (*SearchResult).SortBySeqID}} {
			b.Run(fmt.Sprintf("seqIDs/%d/%s", groups, impl.name), func(b *testing.B) {
				b.ReportAllocs()
				for i := 0; i < b.N; i++ {
					work := poolSimilarityDetails.Get().(*[]*SimilarityDetail)
					*work = append((*work)[:0], details...)
					sr := &SearchResult{SimilarityDetails: work}
					impl.fn(sr)
					recycleSimilarityDetailSlice(sr.SimilarityDetails)
				}
			})
		}
	}
}

// Frozen equal-length cleaner from before query-group sorting was introduced.
func referenceClearEqualLengthSubstrPairs(subs *[]SubstrPair) {
	if len(*subs) < 2 {
		return
	}
	slices.SortFunc(*subs, func(a, b SubstrPair) int {
		if a.QBegin != b.QBegin {
			return int(a.QBegin - b.QBegin)
		}
		return int(a.TBegin - b.TBegin)
	})
	*subs = slices.CompactFunc(*subs, func(a, b SubstrPair) bool {
		return a.QBegin == b.QBegin && a.TBegin == b.TBegin
	})
}
