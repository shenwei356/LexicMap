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
	"sync"
)

// Chaining3Options contains all options in chaining.
type Chaining3Options struct {
	MaxGap      int
	MinScore    int // minimum score of a chain
	MinAlignLen int
	MinIdentity float64
	MaxDistance int
	BandCount   int // only check i in range of  i − A < j < i
	BandBase    int // only check i where i.Qstart+i.Len + A < j.Qstart
}

// DefaultChaining3Options is the defalt vaule of Chaining2Option.
var DefaultChaining3Options = Chaining3Options{
	MaxGap:      5,
	MinScore:    1,
	MinAlignLen: 2,

	MaxDistance: 10,
	BandCount:   20,
	BandBase:    10,
}

// Chainer3 is an object for chaining the anchors in two similar sequences.
// Anchors/seeds/substrings in Chainer3 is denser than those in Chainer,
// and the chaining score function is also much simpler, only considering
// the lengths of anchors and gaps between them.
type Chainer3 struct {
	options *Chaining3Options

	maxscoresIdxs []int64 // pack score (32bit) and index (32bit) to save ram.
}

// NewChainer creates a new chainer.
func NewChainer3(options *Chaining3Options) *Chainer3 {
	c := &Chainer3{
		options: options,

		maxscoresIdxs: make([]int64, 0, 128),
	}
	return c
}

var poolChainers3 = &sync.Pool{New: func() interface{} {
	return NewChainer3(&DefaultChaining3Options)
}}

var poolChain3 = &sync.Pool{New: func() interface{} {
	return &Chain3Result{}
}}

// Chain3Result represents a result of a chain
type Chain3Result struct {
	NAnchors int // The number of substrings

	AlignedFraction float64

	MatchedBases  int     // The number of matched bases.
	AlignedBasesQ int     // The number of aligned bases in Query sequence
	AlignedBasesT int     // The number of aligned bases in Subject sequence
	PIdent        float64 // percentage of identity
	AlignedLength int     // Aligned length, might be longer than AlignedBasesQ or AlignedBasesT
	Gaps          int     // The number of gaps

	QBegin, QEnd int // Query begin/end position (0-based)
	TBegin, TEnd int // Target begin/end position (0-based)
}

// Reset resets a Chain3Result
func (r *Chain3Result) Reset() {
	r.NAnchors = 0
}

var sub0 = &SubstrPair{}

// Chain finds the possible chain path.
// Please remember to call RecycleChaining3Result after using the results.
// Returned results:
//  1. Chain3Result.
//  2. The number of matched bases.
//  3. The number of aligned bases.
//  4. QBegin.
//  5. QEnd.
//  6. TBegin.
//  7. TEnd.
func (ce *Chainer3) Chain(subs *[]*SubstrPair) *Chain3Result {
	n := len(*subs)

	var i, j int
	bandBase := int32(ce.options.BandBase) // band size of banded-DP
	bandCount := ce.options.BandCount
	var _bCount int
	var _bBase int32

	maxscoresIdxs := &ce.maxscoresIdxs
	*maxscoresIdxs = (*maxscoresIdxs)[:0]

	// compute scores
	var s, m, M, g, d float64
	var mj, Mi int
	var a, b *SubstrPair
	maxGap := float64(ce.options.MaxGap)
	maxDistance := float64(ce.options.MaxDistance)

	// initialize
	a = (*subs)[0]
	m = float64(a.Len) - distance2(sub0, a) - gap2(sub0, a)
	*maxscoresIdxs = append(*maxscoresIdxs, int64(m)<<32)

	for i = 1; i < n; i++ {
		a = (*subs)[i] // current seed/anchor

		// just initialize the max score, which comes from the current seed
		m, mj = float64(a.Len)-distance2(sub0, a)-gap2(sub0, a), i

		j = i
		_bCount = 0
		// fmt.Printf("i: %d, %s\n", i, a.String())
		for {
			j--
			if j < 0 {
				break
			}

			b = (*subs)[j] // previous seed/anchor

			// filter out messed/crossed anchors
			if b.QBegin == a.QBegin || b.TBegin > a.TBegin {
				continue
			}

			_bCount++

			_bBase = a.QBegin - b.QBegin - int32(b.Len)
			if !(_bBase <= bandBase || _bCount <= bandCount) {
				break
			}

			d = distance2(a, b)
			if d > maxDistance {
				continue
			}

			g = gap2(a, b)
			if g > maxGap {
				continue
			}

			s = float64((*maxscoresIdxs)[j]>>32) + float64(b.Len) - d - g // compute the score

			if s >= m { // update the max score of current seed/anchor
				m = s
				mj = j
			}
		}

		*maxscoresIdxs = append(*maxscoresIdxs, int64(m)<<32|int64(mj))

		if m > M { // the biggest score in the whole score matrix
			M, Mi = m, i
		}
	}

	// fmt.Printf("M: %f, Mi: %d\n", M, Mi)

	// // print the score matrix
	// fmt.Printf("i\tpair-i\tiMax\tj:scores\n")
	// var _mi int64
	// for i = 0; i < n; i++ {
	// 	_mi = (*maxscoresIdxs)[i]
	// 	fmt.Printf("%d\t%s\t%d:%d", i, (*subs)[i], _mi&4294967295, _mi>>32)
	// 	fmt.Printf("\n")
	// }

	// backtrack -----------------------------------------------------------

	minScore := float64(ce.options.MinScore)
	minAlignLen := ce.options.MinAlignLen

	// check the highest score, for early quit,
	if M < minScore {
		return nil
	}

	var nMatchedBases, nAlignedBasesQ, nAlignedBasesT int

	i = Mi
	var qb, qe, tb, te int32 // the bound (0-based)
	var sub *SubstrPair
	var beginOfNextAnchor int
	var pident float64
	firstAnchorOfAChain := true
	path := poolChain3.Get().(*Chain3Result)
	path.Reset()
	for {
		j = int((*maxscoresIdxs)[i] & 4294967295) // previous seed

		if j < 0 { // the first anchor is not in current region
			break
		}

		sub = (*subs)[i]

		path.NAnchors++

		// fmt.Printf(" AAADDD %d (%s). firstAnchorOfAChain: %v\n", i, *sub, firstAnchorOfAChain)

		if firstAnchorOfAChain {
			// fmt.Printf(" record bound beginning with: %s\n", sub)
			firstAnchorOfAChain = false

			qe = int32(sub.QBegin) + int32(sub.Len) - 1   // end
			te = int32(sub.TBegin) + int32(sub.Len) - 1   // end
			qb, tb = int32(sub.QBegin), int32(sub.TBegin) // in case there's only one anchor

			nMatchedBases += int(sub.Len)
		} else {
			qb, tb = int32(sub.QBegin), int32(sub.TBegin) // begin

			if int(sub.QBegin)+int(sub.Len)-1 >= beginOfNextAnchor {
				nMatchedBases += beginOfNextAnchor - int(sub.QBegin)
			} else {
				nMatchedBases += int(sub.Len)
			}
		}
		beginOfNextAnchor = int(sub.QBegin)

		if i == j { // the path starts here
			if firstAnchorOfAChain { // sadly, there's no anchor added.
				break
			}

			nAlignedBasesQ += int(qe) - int(qb) + 1

			if nAlignedBasesQ < minAlignLen {
				firstAnchorOfAChain = true
				break
			}

			nAlignedBasesT += int(te) - int(tb) + 1

			pident = float64(nMatchedBases) / float64(max(nAlignedBasesQ, nAlignedBasesT)) * 100
			// fmt.Println(nMatchedBases, nAlignedBasesQ, pident)

			// the pident here (pseudo alignment) would be much lower than the real one .
			if pident < 15 {
				firstAnchorOfAChain = true
				break
			}
			if pident > 100 {
				pident = 100
			}

			// reverseInts(path.Chain)
			path.AlignedBasesQ = nAlignedBasesQ
			path.AlignedBasesT = nAlignedBasesT
			path.MatchedBases = nMatchedBases
			path.PIdent = pident
			path.QBegin, path.QEnd = int(qb), int(qe)
			path.TBegin, path.TEnd = int(tb), int(te)

			// fmt.Printf("chain a %d (%d, %d) vs (%d, %d), a:%d, m:%d\n",
			// 	len(*paths), qb, qe, tb, te, nAlignedBasesQ, nMatchedBases)

			firstAnchorOfAChain = true
			break
		}

		i = j
	}

	return path
}
