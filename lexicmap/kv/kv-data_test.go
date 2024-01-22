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
	maxMismatch := -1

	// generate data

	data := make([]map[uint64]*[]uint64, 0, 4)

	var n uint64 = 1 << (k << 1) // all value for k=5
	var i uint64

	for j := 0; j < nMasks; j++ {
		m := make(map[uint64]*[]uint64, n)
		for i = 0; i < n; i++ {
			m[i] = &[]uint64{i}
		}
		data = append(data, m)
	}

	// write data

	file := "t.kv"
	_, err := WriteKVData(k, 0, data, file, 0)
	if err != nil {
		t.Errorf("%s", err)
	}

	// query

	scr, err := NewSearcher(file)
	if err != nil {
		t.Errorf("%s", err)
	}

	var mPrefix uint8

	// exactly query
	mPrefix = k

	var hit bool
	var nExpectedResults int
	for mPrefix = 1; mPrefix <= k; mPrefix++ {
		for i = 0; i < n; i++ {
			// t.Logf("q:%s, prefix:%d, maxMismatch:%d", lexichash.MustDecode(i, scr.K), mPrefix, maxMismatch)
			results, err := scr.Search(i, mPrefix, maxMismatch)
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
						t.Logf("  %s, prefix:%d, mismatch:%d %d\n",
							lexichash.MustDecode(r.Kmer, scr.K), r.LenPrefix, r.Mismatch, v)
					}
				}
			}

			hit = false
			for _, r := range *results {
				if r.Kmer == i {
					hit = true
					if r.Values[0] != i {
						t.Errorf("unexpected value: %d, expected: %d", r.Values[0], i)
						return
					}
				}
			}
			if !hit {
				t.Errorf("query (%s) not found in the results:", lexichash.MustDecode(i, k))
				for _, r := range *results {
					for _, v := range r.Values {
						t.Logf("  %s, prefix:%d, mismatch:%d %d\n",
							lexichash.MustDecode(r.Kmer, scr.K), r.LenPrefix, r.Mismatch, v)
					}
				}
				return
			}
		}
	}

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
