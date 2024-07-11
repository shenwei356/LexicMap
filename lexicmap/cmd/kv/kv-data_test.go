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

package kv

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/shenwei356/lexichash"
)

func TestKVData(t *testing.T) {
	var k uint8 = 5
	nMasks := 3
	// maxMismatch := -1

	// generate data

	data := make([]*map[uint64]*[]uint64, 0, 4)

	var n uint64 = 1 << (k << 1) // all value for k=5
	var i uint64

	for j := 0; j < nMasks; j++ {
		m := make(map[uint64]*[]uint64, n)
		for i = 0; i < n; i++ {
			m[i] = &[]uint64{i + (i << 30)}
		}
		data = append(data, &m)
	}

	// write data

	file := "t.kv"
	_, err := WriteKVData(k, 0, data, file, 0)
	if err != nil {
		t.Errorf("%s", err)
	}

	// -------------------------------------------------------------------
	// reader
	rdr, err := NewReader(file)
	if err != nil {
		t.Errorf("%s", err)
	}
	var ok bool
	var values *[]uint64
	for i := 0; i < nMasks; i++ {
		m1, err := rdr.ReadDataOfAMaskAsMap()
		if err != nil {
			t.Errorf("%s", err)
		}

		m0 := data[i]

		if len(*m0) != len(*m1) {
			t.Errorf("k-mer number mismatch, expected: %d, result: %d", len(*m0), len(*m1))
			return
		}

		var j uint64
		for j = 0; j < n; j++ {
			if values, ok = (*m1)[j]; !ok {
				t.Errorf("k-mer missing: %d: %s", i, lexichash.MustDecode(j, k))
				return
			}
			if len(*values) != 1 {
				t.Errorf("%s: value number mismatch, expected: %d, result: %d, %d",
					lexichash.MustDecode(j, k), 1, len(*values), *values)
				return
			}
			if (*values)[0] != j+(j<<30) {
				t.Errorf("%s: value mismatch, expected: %d, result: %d",
					lexichash.MustDecode(j, k), j+(j<<30), (*values)[0])
				return
			}
		}

		RecycleKmerData(m1)
	}

	// -------------------------------------------------------------------
	// reader
	rdr2, err := NewReader(file)
	if err != nil {
		t.Errorf("%s", err)
	}

	var v uint64
	for i := 0; i < nMasks; i++ {
		m0 := make(map[uint64]uint64) // input
		var _i uint64
		for _i = 0; _i < n; _i++ {
			m0[_i] = _i + (_i << 30)
		}

		m1, err := rdr2.ReadDataOfAMaskAsList()
		if err != nil {
			t.Errorf("%s", err)
		}

		if len(m0) > len(m1) {
			t.Errorf("k-mer number mismatch, expected: >=%d, result: %d", len(m0), len(m1))
			return
		}

		var j uint64
		for j = 0; j < n; j++ {
			if v, ok = m0[j]; !ok {
				t.Errorf("missing k-mer: %d: %s", i, lexichash.MustDecode(j, k))
				return
			}
			if v != j+(j<<30) {
				t.Errorf("%s: value mismatch, expected: %d, result: %d",
					lexichash.MustDecode(j, k), j+(j<<30), (*values)[0])
				return
			}
			delete(m0, j)
		}
		if len(m0) != 0 {
			t.Errorf("wanted %d k-mer", len(m0))
		}
	}

	// -------------------------------------------------------------------
	// searcher

	scr, err := NewSearcher(file)
	if err != nil {
		t.Errorf("%s", err)
	}

	var mPrefix uint8

	// exactly query
	mPrefix = k

	var hit bool
	var nExpectedResults int
	kmers := make([]uint64, nMasks)
	for mPrefix = 1; mPrefix <= k; mPrefix++ {
		for i = 1; i < n-1; i++ { // not from 1 to n, because aaaaa and ttttt is skipped to search
			// t.Logf("q:%s, prefix:%d, maxMismatch:%d", lexichash.MustDecode(i, scr.K), mPrefix, maxMismatch)
			for j := 0; j < nMasks; j++ {
				kmers[j] = i
			}
			// results, err := scr.Search(kmers, mPrefix, maxMismatch)
			results, err := scr.Search(kmers, mPrefix)
			if err != nil {
				t.Errorf("%s", err)
				return
			}

			// for _, r := range *results {
			// 	for _, v := range r.Values {
			// 		t.Logf("  %s, prefix:%d, mismatch:%d %d\n",
			// 			lexichash.MustDecode(r.Kmer, scr.K), r.LenPrefix, r.Mismatch, v)
			// 	}
			// }

			nExpectedResults = nMasks * (1 << ((k - mPrefix) << 1))
			if len(*results) != nExpectedResults {
				t.Logf("query: %s\n", lexichash.MustDecode(i, k))
				t.Errorf("unexpected number of results: %d, expected: %d", len(*results), nExpectedResults)
				for _, r := range *results {
					for _, v := range r.Values {
						// t.Logf("  %s, prefix:%d, mismatch:%d %d\n",
						// 	lexichash.MustDecode(r.Kmer, scr.K), r.LenPrefix, r.Mismatch, v)
						// t.Logf("  %s, prefix:%d %d\n",
						// 	lexichash.MustDecode(r.Kmer, scr.K), r.LenPrefix, v)
						t.Logf("  prefix:%d %d\n", r.LenPrefix, v)
					}
				}
			}

			hit = false
			for _, r := range *results {
				// if r.Kmer == i {
				if r.LenPrefix == k {
					hit = true
					if r.Values[0] != i+(i<<30) {
						// t.Errorf("unexpected value: %d, expected: %d", r.Values[0], i)
						// return
						hit = false
					}
				}
			}
			if !hit {
				t.Errorf("query (%s) not found in the results:", lexichash.MustDecode(i, k))
				for _, r := range *results {
					for _, v := range r.Values {
						// t.Logf("  %s, prefix:%d, mismatch:%d %d\n",
						// 	lexichash.MustDecode(r.Kmer, scr.K), r.LenPrefix, r.Mismatch, v)
						// t.Logf("  %s, prefix:%d %d\n",
						// 	lexichash.MustDecode(r.Kmer, scr.K), r.LenPrefix, v)
						t.Logf("  prefix:%d %d\n", r.LenPrefix, v)
					}
				}
				return
			}
		}
	}

	// -------------------------------------------------------------------
	// in memory searcher

	scr2, err := NewInMemomrySearcher(file)
	if err != nil {
		t.Errorf("%s", err)
	}

	// exactly query
	mPrefix = k

	for mPrefix = k; mPrefix <= k; mPrefix++ {
		for i = 1; i < n-1; i++ { // not from 1 to n, because aaaaa and ttttt is skipped to search
			// t.Logf("q:%s, prefix:%d, maxMismatch:%d", lexichash.MustDecode(i, scr.K), mPrefix, maxMismatch)
			for j := 0; j < nMasks; j++ {
				kmers[j] = i
			}
			// results, err := scr2.Search(kmers, mPrefix, maxMismatch)
			results, err := scr2.Search(kmers, mPrefix)
			if err != nil {
				t.Errorf("%s", err)
				return
			}

			// for _, r := range *results {
			// 	for _, v := range r.Values {
			// 		t.Logf("  %s, prefix:%d, mismatch:%d %d\n",
			// 			lexichash.MustDecode(r.Kmer, scr.K), r.LenPrefix, r.Mismatch, v)
			// 	}
			// }

			nExpectedResults = nMasks * (1 << ((k - mPrefix) << 1))
			if len(*results) != nExpectedResults {
				t.Logf("query: %s\n", lexichash.MustDecode(i, k))
				t.Errorf("unexpected number of results: %d, expected: %d", len(*results), nExpectedResults)
				for _, r := range *results {
					for _, v := range r.Values {
						// t.Logf("  %s, prefix:%d, mismatch:%d %d\n",
						// 	lexichash.MustDecode(r.Kmer, scr2.K), r.LenPrefix, r.Mismatch, v)
						// t.Logf("  %s, prefix:%d %d\n",
						// 	lexichash.MustDecode(r.Kmer, scr.K), r.LenPrefix, v)
						t.Logf("  prefix:%d %d\n", r.LenPrefix, v)
					}
				}
			}

			hit = false
			for _, r := range *results {
				// if r.Kmer == i {
				if r.LenPrefix == k {
					hit = true
					if r.Values[0] != i+(i<<30) {
						t.Errorf("unexpected value: %d, expected: %d", r.Values[0], i)
						return
					}
				}
			}
			if !hit {
				t.Errorf("query (%s) not found in the results:", lexichash.MustDecode(i, k))
				for _, r := range *results {
					for _, v := range r.Values {
						// t.Logf("  %s, prefix:%d, mismatch:%d %d\n",
						// 	lexichash.MustDecode(r.Kmer, scr2.K), r.LenPrefix, r.Mismatch, v)
						// t.Logf("  %s, prefix:%d %d\n",
						// 	lexichash.MustDecode(r.Kmer, scr.K), r.LenPrefix, v)

						t.Logf("  prefix:%d %d\n", r.LenPrefix, v)
					}
				}
				return
			}
		}
	}

	// -------------------------------------------------------------------

	// clean up
	if os.RemoveAll(file) != nil {
		t.Errorf("failed to remove the kv-data file: %s", file)
		return
	}
	fileIdx := filepath.Clean(file) + KVIndexFileExt
	if os.RemoveAll(fileIdx) != nil {
		t.Errorf("failed to remove the kv-data file: %s", fileIdx)
		return
	}
}
