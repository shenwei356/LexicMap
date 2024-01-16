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

package tree

import (
	"math/rand"
	"os"
	"testing"

	"github.com/shenwei356/LexicMap/lexicmap/util"
	"github.com/shenwei356/lexichash"
)

func TestSerialization(t *testing.T) {
	k := 21
	nKmers := 100000 // or the number of sequences
	var seed int64 = 1

	tree := New(uint8(k))

	n := int(float64(nKmers) * 1.2)
	codes := make([]uint64, n)

	r := rand.New(rand.NewSource(seed))
	shift := 64 - k*2
	var kmer uint64
	for i := range codes {
		kmer = util.Hash64(r.Uint64()) >> shift // hash a random int and cut into k*2 bits
		codes[i] = kmer
	}
	util.UniqUint64s(&codes) // remove duplicates
	if len(codes) > nKmers {
		codes = codes[:nKmers]
	}

	for i, kmer := range codes {
		tree.Insert(kmer, uint64(i))
	}

	// ----------------------------------------

	file := "test.ktree.gz"

	N, err := tree.WriteToFile(file)
	if err != nil {
		t.Errorf("writing the tree to file: %s", err)
		return
	}
	t.Logf("%d k-mers are saved to file: %s, number of bytes of uncompressed data: %d",
		tree.NumLeafNodes(), file, N)

	// ----------------------------------------

	tree2, err := NewFromFile(file)
	if err != nil {
		t.Errorf("new tree from file: %s", err)
		return
	}

	if tree.k != tree2.k {
		t.Errorf("Ks unmatched: %d vs %d", tree.k, tree2.k)
		return
	}

	if tree.NumLeafNodes() != tree2.NumLeafNodes() {
		t.Errorf("numbers of k-mers unmatched: %d vs %d", tree.NumLeafNodes(), tree2.NumLeafNodes())
		return
	}

	var v2 []uint64
	var ok bool
	var _v uint64
	var i int
	tree.Walk(func(code uint64, v []uint64) bool {
		v2, ok = tree2.Get(code)
		if !ok {
			t.Errorf("k-mer missing: %s", lexichash.MustDecode(code, uint8(k)))
			return true
		}
		if len(v) != len(v2) {
			t.Errorf("numbers of values unmatched: %d vs %d", len(v), len(v2))
			return true
		}

		for i, _v = range v {
			if _v != v2[i] {
				t.Errorf("values unmatched: %d vs %d", _v, v2[i])
				return true
			}
		}

		return false
	})

	if os.RemoveAll(file) != nil {
		t.Errorf("failed to remove the file: %s", file)
		return
	}
}
