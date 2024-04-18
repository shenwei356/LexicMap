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
	"sync"
)

// Chaining2Options contains all options in chaining.
type Chaining2Options struct {
	MaxGap      int
	MinScore    int // minimum score of a chain
	MinAlignLen int
	MinIdentity float64
	MaxDistance int
	Band        int // only check i in range of  i − A < j < i
}

// DefaultChaining2Options is the defalt vaule of Chaining2Option.
var DefaultChaining2Options = Chaining2Options{
	MaxGap:      50,
	MinScore:    50,
	MinAlignLen: 50,

	MaxDistance: 50,
	Band:        100,
}

// Chainer2 is an object for chaining the anchors in two similar sequences.
// Anchors/seeds/substrings in Chainer2 is denser than those in Chainer,
// and the chaining score function is also much simpler, only considering
// the lengths of anchors and gaps between them.
type Chainer2 struct {
	options *Chaining2Options

	// scores        []int
	// maxscores     []int
	// maxscoresIdxs []int
	maxscoresIdxs []uint64 // pack score (32bit) and index (32bit) ot save ram.

	bounds []int32 // 4 * chains
}

// NewChainer creates a new chainer.
func NewChainer2(options *Chaining2Options) *Chainer2 {
	c := &Chainer2{
		options: options,

		// scores:        make([]int, 0, 10240),
		// maxscores:     make([]int, 0, 51200),
		maxscoresIdxs: make([]uint64, 0, 51200),
		bounds:        make([]int32, 128),
	}
	return c
}

// RecycleChaining2Result reycles the chaining paths.
// Please remember to call this after using the results.
func RecycleChaining2Result(chains *[]*Chain2Result) {
	for _, chain := range *chains {
		poolChain2.Put(chain)
	}
	poolChains2.Put(chains)
}

var poolChains2 = &sync.Pool{New: func() interface{} {
	tmp := make([]*Chain2Result, 0, 32)
	return &tmp
}}

var poolChain2 = &sync.Pool{New: func() interface{} {
	return &Chain2Result{
		// Chain: make([]int, 0, 128),
	}
}}

// Chain2Result represents a result of a chain
type Chain2Result struct {
	NAnchors int // The number of substrings

	MatchedBases int // The number of matched bases.
	AlignedBases int // The number of aligned bases.
	QBegin, QEnd int // Query begin/end position (0-based)
	TBegin, TEnd int // Target begin/end position (0-based)
}

// Reset resets a Chain2Result
func (r *Chain2Result) Reset() {
	r.NAnchors = 0
}

// Chain finds the possible chain paths.
// Please remember to call RecycleChaining2Result after using the results.
// Returned results:
//  1. Chain2Results.
//  2. The number of matched bases.
//  3. The number of aligned bases.
//  4. QBegin.
//  5. QEnd.
//  6. TBegin.
//  7. TEnd.
func (ce *Chainer2) Chain(subs *[]*SubstrPair) (*[]*Chain2Result, int, int, int, int, int, int) {
	n := len(*subs)

	if n == 1 { // for one seed, just check the seed weight
		paths := poolChains2.Get().(*[]*Chain2Result)
		*paths = (*paths)[:0]

		sub := (*subs)[0]
		slen := int(sub.Len)
		if slen >= ce.options.MinScore && slen >= ce.options.MinAlignLen { // the length of anchor (max 32)
			path := poolChain2.Get().(*Chain2Result)
			path.Reset()

			qe := int(sub.QBegin) + slen - 1           // end
			te := int(sub.TBegin) + slen - 1           // end
			qb, tb := int(sub.QBegin), int(sub.TBegin) // in case there's only one anchor
			path.QBegin, path.QEnd = qb, qe
			path.TBegin, path.TEnd = tb, te
			path.MatchedBases = slen
			path.AlignedBases = slen
			path.NAnchors++
			*paths = append(*paths, path)

			return paths, slen, slen, qb, qe, tb, te
		}

		return paths, 0, 0, 0, 0, 0, 0
	}

	var i, _b, j int
	// var k int
	band := ce.options.Band // band size of banded-DP

	// a list for storing score matrix, the size is band * len(seeds pair)
	// scores := &ce.scores
	// *scores = (*scores)[:0]
	// size := n * (band + 1)
	// for k = 0; k < size; k++ {
	// 	*scores = append(*scores, 0)
	// }

	// reused objects

	// the maximum score for each seed, the size is n
	// maxscores := &ce.maxscores
	// *maxscores = (*maxscores)[:0]
	// index of previous seed, the size is n. pointers for backtracking.
	maxscoresIdxs := &ce.maxscoresIdxs
	*maxscoresIdxs = (*maxscoresIdxs)[:0]

	// initialize
	// *maxscores = append(*maxscores, int((*subs)[0].Len))
	// *maxscoresIdxs = append(*maxscoresIdxs, 0)
	*maxscoresIdxs = append(*maxscoresIdxs, uint64((*subs)[0].Len)<<32)

	// compute scores
	var s, m, M, d, g int
	var mj, Mi int
	var a, b *SubstrPair
	maxGap := ce.options.MaxGap
	maxDistance := ce.options.MaxDistance
	// (*scores)[0] = (*subs)[0].Len
	for i = 1; i < n; i++ {
		a = (*subs)[i] // current seed/anchor
		// k = band * i   // index of current seed in the score matrix

		// just initialize the max score, which comes from the current seed
		m, mj = int(a.Len), i
		// (*scores)[k] = m

		for _b = 1; _b <= band; _b++ { // check previous $band seeds
			j = i - _b // index of the previous seed
			if j < 0 {
				break
			}

			b = (*subs)[j] // previous seed/anchor
			// k++            // index of previous seed in the score matrix

			if b.TBegin > a.TBegin { // filter out messed/crossed anchors
				continue
			}

			d = distance2(a, b)
			if d > maxDistance { // limit the distance. necessary?
				continue
			}

			g = gap2(a, b)
			if g > maxGap { // limit the gap. necessary?
				continue
			}

			// s = (*maxscores)[j] + int(b.Len) - g // compute the score
			s = int((*maxscoresIdxs)[j]>>32) + int(b.Len) - g // compute the score
			// (*scores)[k] = s                // necessary?

			if s >= m { // update the max score of current seed/anchor
				m = s
				mj = j
			}
		}

		// *maxscores = append(*maxscores, m)          // save the max score of the whole
		// *maxscoresIdxs = append(*maxscoresIdxs, mj) // save where the max score comes from
		*maxscoresIdxs = append(*maxscoresIdxs, uint64(m)<<32|uint64(mj))

		if m > M { // the biggest score in the whole score matrix
			M, Mi = m, i
		}
	}

	// fmt.Printf("M: %d, Mi: %d\n", M, Mi)

	// print the score matrix
	// fmt.Printf("i\tpair-i\tiMax\tj:scores\n")
	// var _mi uint64
	// for i = 0; i < n; i++ {
	// 	_mi = (*maxscoresIdxs)[i]
	// 	fmt.Printf("%d\t%s\t%d:%d", i, (*subs)[i], _mi&4294967295, _mi>>32)
	// 	fmt.Printf("\n")
	// }

	// backtrack

	minScore := ce.options.MinScore
	minAlignLen := ce.options.MinAlignLen

	// check the highest score, for early quit,
	if M < minScore {
		return nil, 0, 0, 0, 0, 0, 0
	}

	paths := poolChains2.Get().(*[]*Chain2Result)
	*paths = (*paths)[:0]

	var nMatchedBases, nAlignedBases int
	ce.bounds = ce.bounds[:0]

	_, qB, qE, tB, tE := chainARegion(
		subs,
		// maxscores,
		maxscoresIdxs,
		0,
		minScore,
		minAlignLen,
		ce.options.MinIdentity,
		paths,
		&nMatchedBases,
		&nAlignedBases,
		Mi,
		&ce.bounds,
	)

	return paths, nMatchedBases, nAlignedBases, qB, qE, tB, tE
}

func chainARegion(subs *[]*SubstrPair, // a region of the subs
	// maxscores *[]int, // a region of maxscores
	maxscoresIdxs *[]uint64,
	offset int, // offset of this region of subs
	minScore int, // the threshold
	minAlignLen int,
	minPident float64,
	paths *[]*Chain2Result, // paths
	_nMatchedBases *int,
	_nAlignedBases *int,
	Mi0 int, // found Mi
	bounds *[]int32, // intervals of previous chains
) (
	int, // score
	int, // query begin position (0-based)
	int, // query end position (0-based)
	int, // target begin position (0-based)
	int, // target end position (0-based)
) {
	// fmt.Printf("region: [%d, %d]\n", offset, offset+len(*subs)-1)
	var _mi uint64
	var m, M int
	var i, Mi int
	if Mi0 < 0 { // Mi is not given
		// find the next highest score
		// for i, m = range *maxscores {
		for i, _mi = range *maxscoresIdxs {
			m = int(_mi >> 32)
			if m > M {
				M, Mi = m, i
			}
		}
		if M < minScore { // no valid anchors
			return 0, -1, -1, -1, -1
		}
	} else {
		Mi = Mi0
	}
	// fmt.Printf("  Mi: %d, M: %d\n", Mi, M)

	var nMatchedBases int
	var nAlignedBases int

	i = Mi
	var j int
	var qB, qE, tB, tE int // the bound of the chain (0-based)
	qB, tB = math.MaxInt, math.MaxInt
	var qb, qe, tb, te int32 // the bound (0-based)
	var sub *SubstrPair
	var beginOfNextAnchor int
	var pident float64
	// var overlapped bool
	// var nb, bi, bj int // index of bounds
	firstAnchorOfAChain := true
	path := poolChain2.Get().(*Chain2Result)
	path.Reset()
	for {
		// j = (*maxscoresIdxs)[i] - offset // previous seed
		j = int((*maxscoresIdxs)[i]&4294967295) - offset // previous seed

		if j < 0 { // the first anchor is not in current region
			break
		}

		// check if an anchor overlaps with previous chains
		//
		// Query
		// |        te  / (OK)
		// |        |  /
		// |(NO)/   |____qe
		// |   /   /
		// |qb____/    / (NO)
		// |   /  |   /
		// |OK/   |tb
		// o-------------------- Ref
		//
		sub = (*subs)[i]
		// overlapped = false
		// nb = len(*bounds) >> 2 // len(bounds) / 4
		// for bi = 0; bi < nb; bi++ {
		// 	bj = bi << 2
		// 	if !((sub.QBegin > (*bounds)[bj+1] && sub.TBegin > (*bounds)[bj+3]) || // top right
		// 		(sub.QBegin+int32(sub.Len)-1 < (*bounds)[bj] && sub.TBegin+int32(sub.Len)-1 < (*bounds)[bj+2])) { // bottom left
		// 		overlapped = true
		// 		break
		// 	}
		// }

		// if overlapped {
		// 	// fmt.Printf("  %d (%s) is overlapped previous chain, j=%d\n", i, *sub, j)

		// 	// can not continue here, must check if i == j
		// } else {
		// path.Chain = append(path.Chain, i+offset) // record the seed
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
		// }

		if i == j { // the path starts here
			if firstAnchorOfAChain { // sadly, there's no anchor added.
				break
			}

			nAlignedBases += int(te) - int(tb) + 1

			if nAlignedBases < minAlignLen {
				firstAnchorOfAChain = true
				break
			}

			pident = float64(nMatchedBases) / float64(nAlignedBases) * 100
			if pident < minPident {
				firstAnchorOfAChain = true
				break
			}

			// reverseInts(path.Chain)
			path.AlignedBases = nAlignedBases
			path.MatchedBases = nMatchedBases
			path.QBegin, path.QEnd = int(qb), int(qe)
			path.TBegin, path.TEnd = int(tb), int(te)
			*paths = append(*paths, path)

			*_nAlignedBases += nAlignedBases
			*_nMatchedBases += nMatchedBases

			// fmt.Printf("chain %d (%d, %d) vs (%d, %d), a:%d, m:%d\n",
			// 	len(*paths), qb, qe, tb, te, nAlignedBases, nMatchedBases)

			firstAnchorOfAChain = true
			break
		}

		i = j
	}

	if j < 0 { // the first anchor is not in current region
		// fmt.Printf(" found only part of the chain, nAnchors: %d\n", path.NAnchors)
		// if len(path.Chain) == 0 {
		if path.NAnchors == 0 {
			poolChain2.Put(path)
		} else {
			nAlignedBases += int(qe) - int(qb) + 1

			if nAlignedBases >= minAlignLen {
				pident = float64(nMatchedBases) / float64(nAlignedBases) * 100
				if pident >= minPident {
					// reverseInts(path.Chain)
					path.AlignedBases = nAlignedBases
					path.MatchedBases = nMatchedBases
					path.QBegin, path.QEnd = int(qb), int(qe)
					path.TBegin, path.TEnd = int(tb), int(te)
					*paths = append(*paths, path)

					*_nAlignedBases += nAlignedBases
					*_nMatchedBases += nMatchedBases

					// fmt.Printf("chain %d (%d, %d) vs (%d, %d), a:%d, m:%d\n",
					// 	len(*paths), qb, qe, tb, te, nAlignedBases, nMatchedBases)
				}
			}
		}
	}

	*bounds = append(*bounds, qb)
	*bounds = append(*bounds, qe)
	*bounds = append(*bounds, tb)
	*bounds = append(*bounds, te)

	// initialize the boundary
	qB, qE = int(qb), int(qe)
	tB, tE = int(tb), int(te)

	// fmt.Printf("  i: %d\n", i)

	// the unchecked region on the right
	// if Mi != len(*maxscores)-1 { // Mi is not the last element
	if Mi != len(*maxscoresIdxs)-1 { // Mi is not the last element
		tmp := (*subs)[Mi+1:]
		// tmps := (*maxscores)[Mi+1:]
		tmpsi := (*maxscoresIdxs)[Mi+1:]
		_score, _qB, _qE, _tB, _tE := chainARegion(
			&tmp,
			// &tmps,
			&tmpsi,
			offset+Mi+1,
			minScore,
			minAlignLen,
			minPident,
			paths,
			_nMatchedBases,
			_nAlignedBases,
			-1,
			bounds,
		)
		if _score > 0 {
			if _qB < qB {
				qB = _qB
			}
			if _qE > qE {
				qE = _qE
			}
			if _tB < tB {
				tB = _tB
			}
			if _tE > tE {
				tE = _tE
			}
		}
	}

	// the unchecked region on the left
	if i > 0 { // the first anchor is not the first element
		tmp := (*subs)[:i]
		// tmps := (*maxscores)[:i]
		tmpsi := (*maxscoresIdxs)[:i]
		_score, _qB, _qE, _tB, _tE := chainARegion(
			&tmp,
			// &tmps,
			&tmpsi,
			offset,
			minScore,
			minAlignLen,
			minPident,
			paths,
			_nMatchedBases,
			_nAlignedBases,
			-1,
			bounds,
		)
		if _score > 0 {
			if _qB < qB {
				qB = _qB
			}
			if _qE > qE {
				qE = _qE
			}
			if _tB < tB {
				tB = _tB
			}
			if _tE > tE {
				tE = _tE
			}
		}
	}

	return M, qB, qE, tB, tE
}

func distance2(a, b *SubstrPair) int {
	q := a.QBegin - b.QBegin
	t := a.TBegin - b.TBegin
	if q > t {
		return int(q)
	}
	return int(t)
}

func gap2(a, b *SubstrPair) int {
	g := a.QBegin - b.QBegin - (a.TBegin - b.TBegin)
	if g < 0 {
		return -int(g)
	}
	return int(g)
}
