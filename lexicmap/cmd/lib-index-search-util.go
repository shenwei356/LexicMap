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

package cmd

import (
	"math"
	"slices"
	"sync"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/tree"
	"github.com/shenwei356/lexichash/iterator"
	"github.com/shenwei356/wfa"
)

// extendMatch an alignment region using a chaining algorithm.
func extendMatch(seq1, seq2 []byte, start1, end1, start2, end2 int, extLen int, tBegin, maxExtLen int, rc bool) ([]byte, []byte, int, int, int, int, error) {
	var m uint8 = 2

	// fmt.Println("before:", start1, end1, start2, end2)

	var _s1, _e1, _s2, _e2 int // extend length
	var _extLen int

	// 3', right
	if end1+int(m) < len(seq1) && end2+int(m) < len(seq2) {
		if rc {
			_extLen = min(extLen, tBegin)
		} else {
			_extLen = min(extLen, maxExtLen)
		}

		if _extLen > 2 {
			e1, e2 := min(end1+_extLen, len(seq1)), min(end2+_extLen, len(seq2))
			_seq1, _seq2 := seq1[end1:e1], seq2[end2:e2]
			// fmt.Printf("seq1: %s\nseq2: %s\n", _seq1, _seq2)

			_e1, _e2 = _extendRight(_seq1, _seq2)
			if _e1 > 0 || _e2 > 0 {
				end1 += _e1
				end2 += _e2
			}
		}
	}

	// 5', left
	if start1 > int(m) && start2 > int(m) {
		if rc {
			_extLen = min(extLen, maxExtLen) // tBegin is 0-based
		} else {
			_extLen = min(extLen, tBegin) // tBegin is 0-based
		}

		if _extLen > 2 {
			s1, s2 := max(start1-_extLen, 0), max(start2-_extLen, 0)
			_seq1, _seq2 := reverseBytes(seq1[s1:start1]), reverseBytes(seq2[s2:start2])
			// fmt.Printf("seq1: %s\nseq2: %s\n", _seq1, _seq2)

			_s1, _s2 = _extendRight(*_seq1, *_seq2)
			if _s1 > 0 || _s2 > 0 {
				start1 -= _s1
				start2 -= _s2
			}
			poolRevBytes.Put(_seq1)
			poolRevBytes.Put(_seq2)
		}
	}

	// fmt.Println("after:", start1, end1, start2, end2)
	return seq1[start1:end1], seq2[start2:end2], _s1, _e1, _s2, _e2, nil
}

func _extendRight(s1, s2 []byte) (int, int) {
	_k := 2
	var m uint8 = 2

	// k-mer iterator
	iter, err := iterator.NewKmerIterator(s1, _k)
	if err != nil {
		return 0, 0
	}

	// index
	t := tree.NewTree(uint8(_k))
	var kmer uint64
	var ok bool
	for {
		kmer, ok, _ = iter.NextPositiveKmer()
		if !ok {
			break
		}
		t.Insert(kmer, uint32(iter.Index()))
	}

	// match
	iter, err = iterator.NewKmerIterator(s2, _k)
	if err != nil {
		return 0, 0
	}

	subs := poolSubs.Get().(*[]*SubstrPair)
	*subs = (*subs)[:0]

	var v, p uint32
	var srs *[]*tree.SearchResult
	var sr *tree.SearchResult

	for {
		kmer, ok, _ = iter.NextPositiveKmer()
		if !ok {
			break
		}

		srs, ok = t.Search(kmer, m)
		if !ok {
			continue
		}

		for _, sr = range *srs {
			// fmt.Printf("%s vs %s, len:%d\n", kmers.MustDecode(kmer, _k), kmers.MustDecode(sr.Kmer, _k), sr.LenPrefix)
			for _, v = range sr.Values {
				p = v

				_sub := poolSub.Get().(*SubstrPair)
				_sub.QBegin = int32(p)
				_sub.TBegin = int32(iter.Index())
				_sub.Len = uint8(sr.LenPrefix)
				_sub.QRC = false
				_sub.TRC = false

				*subs = append(*subs, _sub)
			}
		}
		t.RecycleSearchResult(srs)
	}
	tree.RecycleTree(t)

	if len(*subs) == 0 {
		return 0, 0
	}

	if len(*subs) > 1 {
		// no need to clean as k == min_len
		// ClearSubstrPairs(poolSub, subs, _k)

		slices.SortFunc(*subs, func(a, b *SubstrPair) int {
			if a.QBegin == b.QBegin {
				if a.QBegin+int32(a.Len) == b.QBegin+int32(b.Len) {
					return int(a.TBegin - b.TBegin)
				}
				return int(b.QBegin) + int(b.Len) - (int(a.QBegin) + int(a.Len))
			}
			return int(a.QBegin - b.QBegin)
		})
	}

	// for _, s := range *subs {
	// 	fmt.Println(s)
	// }

	// chaining
	chainer := poolChainers3.Get().(*Chainer3)
	chain := chainer.Chain(subs)

	RecycleSubstrPairs(poolSub, subs)
	poolChainers3.Put(chainer)

	if chain != nil {
		// fmt.Printf("q: %d-%d, t: %d-%d\n", chain.QBegin, chain.QEnd, chain.TBegin, chain.TEnd)
		poolChain3.Put(chain)
		return chain.QEnd + 1, chain.TEnd + 1
	}

	return 0, 0
}

// remember to recycle the result
func reverseBytes(s []byte) *[]byte {
	t := poolRevBytes.Get().(*[]byte)
	if len(s) == len(*t) {

	} else if len(s) < len(*t) {
		*t = (*t)[:len(s)]
	} else {
		n := len(s) - len(*t)
		for i := 0; i < n; i++ {
			*t = append(*t, 0)
		}
	}
	copy(*t, s)

	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		(*t)[i], (*t)[j] = (*t)[j], (*t)[i]
	}

	return t
}

var poolRevBytes = &sync.Pool{New: func() interface{} {
	tmp := make([]byte, 128)
	return &tmp
}}

// ------------------------------------------------------------------------------------------

const OpM = uint64('M')
const OpD = uint64('D')
const OpI = uint64('I')
const OpX = uint64('X')
const OpH = uint64('H')

// trimOps trim ops to keep only aligned region
func trimOps(ops []uint64) []uint64 {
	var start, end int
	start, end = -1, -1
	for i, op := range ops {
		if op>>32 == OpM {
			start = i
			break
		}
	}
	for i := len(ops) - 1; i >= 0; i-- {
		if ops[i]>>32 == OpM {
			end = i
			break
		}
	}
	return ops[start : end+1]
}

func scoreAndEvalue(match, mismatch, gapOpen, gapExt int, totalBase int, lambda, k float64) func(qlen int, cigar *wfa.AlignmentResult) (int, int, float64) {
	// var Kn float64 = float64(k) * float64(totalBase)
	lnK := math.Log(k)

	return func(qlen int, cigar *wfa.AlignmentResult) (int, int, float64) {
		ops := trimOps(cigar.Ops)
		var score, n int
		for _, op := range ops {
			n = int(op & 4294967295)

			// switch op.Op {
			switch op >> 32 {
			// match:
			case OpM:
				score += n * match
			// mismatch
			case OpX:
				score += n * mismatch
			// gap
			case OpI:
				score += gapOpen + n*gapExt
			// case 'D', 'H':
			case OpD, OpH:
				score += gapOpen + n*gapExt
			}
		}

		_score := score

		// from blastn_values_2_3 in ncbi-blast-2.15.0+-src/c++/src/algo/blast/core/blast_stat.c
		// Any odd score must be rounded down to the nearest even number before calculating the e-value
		if _score&1 == 1 {
			_score--
		}

		bitScore := (lambda*float64(_score) - lnK) / math.Ln2

		// evalue := Kn * float64(qlen) * math.Pow(math.E, -lambda*float64(_score))

		evalue := float64(totalBase) * math.Pow(2, -bitScore) * float64(qlen)

		return score, int(math.Ceil(bitScore)), evalue
	}
}
