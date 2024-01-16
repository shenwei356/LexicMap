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
	"bytes"
	"testing"

	"github.com/shenwei356/LexicMap/lexicmap/util"
	"github.com/shenwei356/kmers"
)

func parseKmer(s string) ([]byte, uint64, uint8) {
	kmer := []byte(s)
	code, _ := kmers.Encode(kmer)
	k := uint8(len(kmer))
	return kmer, code, k
}

func TestKmerOperations(t *testing.T) {
	kmer, code, k := parseKmer("ACTGACCTGC")

	prefix1, p1, k1 := parseKmer("ACTGCA")
	prefix2, p2, k2 := parseKmer("ACTGC")

	// KmerBaseAt
	var c uint8
	for i, b := range kmer {
		c = util.KmerBaseAt(code, k, uint8(i))
		if bit2base[c] != b {
			t.Errorf("KmerBaseAt error: %d, expected %c, returned %c", i, b, c)
		}
	}

	// KmerPrefix
	var p []byte
	for i := 1; i <= len(kmer); i++ {
		p = kmers.MustDecode(util.KmerPrefix(code, k, uint8(i)), i)
		if !bytes.Equal(p, kmer[:i]) {
			t.Errorf("KmerPrefix error: %d, expected %s, returned %s", i, kmer[:i], p)
		}
	}

	// KmerSuffix
	var s []byte
	for i := 0; i < len(kmer); i++ {
		s = kmers.MustDecode(util.KmerSuffix(code, k, uint8(i)), len(kmer)-i)
		if !bytes.Equal(s, kmer[i:]) {
			t.Errorf("KmerSuffix error: %d, expected %s, returned %s", i, kmer[i:], s)
		}
	}

	// KmerLongestPrefix
	n := util.KmerLongestPrefix(p1, p2, k1, k2)
	if n != 5 {
		t.Errorf("KmerLongestPrefix error: expected %d, returned %d", 5, n)
	}

	n = util.KmerLongestPrefix(p2, p1, k2, k1)
	if n != 5 {
		t.Errorf("KmerLongestPrefix error: expected %d, returned %d", 5, n)
	}

	// KmerHasPrefix
	b := util.KmerHasPrefix(code, p1, k, k1)
	has := bytes.HasPrefix(kmer, prefix1)
	if b != has {
		t.Errorf("KmerHasPrefix error: expected %v, returned %v", has, b)
	}

	b = util.KmerHasPrefix(code, p2, k, k2)
	has = bytes.HasPrefix(kmer, prefix2)
	if b != has {
		t.Errorf("KmerHasPrefix error: expected %v, returned %v", has, b)
	}

}
