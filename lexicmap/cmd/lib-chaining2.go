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
	"math"
	"math/bits"
	"sync"
)

// Chaining2Options contains all options in chaining.
type Chaining2Options struct {
	MaxGap      int
	MinScore    int // default for Chain; SeqComparator supplies a per-call score
	MinAlignLen int
	MinIdentity float64
	MaxDistance int
	BandCount   int // only check i in range of  i − A < j < i
	BandBase    int // only check i where i.Qstart+i.Len + A < j.Qstart

	HeuristicKmerPidentThreshold float64 // for filtering less similar chaining results
}

// DefaultChaining2Options contains the default options for Chainer2.
var DefaultChaining2Options = Chaining2Options{
	MaxGap:      50,
	MinScore:    2 * int(defaultSeqComparatorMinPrefix),
	MinAlignLen: 50,

	MaxDistance: 100,
	BandCount:   50,
	BandBase:    100,

	HeuristicKmerPidentThreshold: 15,
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
	maxscoresIdxs []uint64 // pack score (32bit) and index (32bit) to save ram.
	maxPotential  []int64  // packed predecessor score maxima in global or diagonal-bin prefixes
	diagonalSeen  []uint64 // reusable approximate set of earlier anchor diagonals
	qGroupStarts  []int32  // first anchor of each equal-query-position group
	qGroupCursors []int32  // previous target upper bounds, indexed by group start

	// bounds []int32 // 4 * chains, can't limit this or it will miss some cases.
}

// targetUpperBound searches [start, end), reusing the preceding search boundary.
// For sorted T = [10, 20, 100, 200], nearby limits 25 and 30 both return 2.
// Checking the two neighbors avoids another binary search in this common case.
// Otherwise bracket the movement exponentially, then binary-search that interval.
// Both directions are needed across Q groups; the cursor is only a search hint.
func targetUpperBound(subs []SubstrPair, start, end int, limit int32, cursor *int32) int {
	low, high := start, end
	at := min(max(int(*cursor), start), end)
	if at > start && subs[at-1].TBegin > limit {
		high = at - 1
		for step := 1; ; step *= 2 {
			probe := max(start, at-1-step)
			if subs[probe].TBegin <= limit {
				low = probe + 1
				break
			}
			high = probe
			if probe == start {
				low = start
				break
			}
		}
	} else if at < end && subs[at].TBegin <= limit {
		low = at + 1
		for step := 1; ; step *= 2 {
			probe := min(end, at+step)
			if probe == end || subs[probe].TBegin > limit {
				high = probe
				break
			}
			low = probe + 1
		}
	} else {
		return at
	}
	for low < high {
		mid := low + (high-low)/2
		if subs[mid].TBegin <= limit {
			low = mid + 1
		} else {
			high = mid
		}
	}
	*cursor = int32(low)
	return low
}

// NewChainer creates a new chainer.
func NewChainer2(options *Chaining2Options) *Chainer2 {
	c := &Chainer2{
		options: options,

		// scores:        make([]int, 0, 10240),
		// maxscores:     make([]int, 0, 51200),
		maxscoresIdxs: make([]uint64, 0, 102400),
		// bounds:        make([]int32, 128),
	}
	return c
}

// RecycleChaining2Result reycles the chaining paths.
// Please remember to call this after using the results.
func RecycleChaining2Result(chains *[]*Chain2Result) {
	for _, chain := range *chains {
		if chain != nil {
			recycleChain2(chain)
		}
	}
	recycleChaining2ResultSlice(chains)
}

// Recycle only the slice when its HSPs have been transferred to another owner.
func recycleChaining2ResultSlice(chains *[]*Chain2Result) {
	if cap(*chains) > thresholdNSubs {
		*chains = nil
		return
	}
	clear((*chains)[:cap(*chains)])
	*chains = (*chains)[:0]
	poolChains2.Put(chains)
}

const maxPooledAlignmentBytes = 16 << 10 // total capacity of the four output buffers

func recycleChain2(chain *Chain2Result) {
	if cap(chain.CIGAR)+cap(chain.QSeq)+cap(chain.TSeq)+cap(chain.Alignment) > maxPooledAlignmentBytes {
		chain.CIGAR = nil
		chain.QSeq = nil
		chain.TSeq = nil
		chain.Alignment = nil
	} else {
		chain.CIGAR = chain.CIGAR[:0]
		chain.QSeq = chain.QSeq[:0]
		chain.TSeq = chain.TSeq[:0]
		chain.Alignment = chain.Alignment[:0]
	}
	poolChain2.Put(chain)
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

	AlignedFraction float64

	MatchedBases  int     // The number of matched bases.
	AlignedBasesQ int     // The number of aligned bases in Query sequence
	AlignedBasesT int     // The number of aligned bases in Subject sequence
	PIdent        float64 // percentage of identity
	AlignedLength int     // Aligned length, might be longer than AlignedBasesQ or AlignedBasesT
	Gaps          int     // The number of gaps

	PIdentAdjusted float64 // quality-adjusted percentage of identity

	QBegin, QEnd int // Query begin/end position (0-based)
	TBegin, TEnd int // Target begin/end position (0-based)

	MaxExtLen int // max extend length

	tPosOffsetBegin int // start position of the sequence in the concatenated genome

	// for output
	CIGAR     []byte // cigar string
	QSeq      []byte // query seq
	TSeq      []byte // target seq
	Alignment []byte // alignment text

	// statistic
	Score    int
	BitScore int
	Evalue   float64
}

// Reset clears all previous alignment metadata while retaining small output
// buffers for reuse. A new pseudoalignment must not expose another hit's score.
func (r *Chain2Result) Reset() {
	*r = Chain2Result{CIGAR: r.CIGAR[:0], QSeq: r.QSeq[:0], TSeq: r.TSeq[:0], Alignment: r.Alignment[:0]}
}

// Chain finds possible chain paths. Empty anchor sets return no chains.
// Please remember to call RecycleChaining2Result after using the results.
// Returned results:
//  1. Chain2Results.
//  2. The number of matched bases.
//  3. The number of aligned bases.
//  4. QBegin.
//  5. QEnd.
//  6. TBegin.
//  7. TEnd.
func (ce *Chainer2) Chain(subs *[]SubstrPair) (*[]*Chain2Result, int, int, int, int, int, int, int) {
	return ce.chainWithLimits(subs, ce.options.MinScore, ce.options.MinAlignLen)
}

// chainWithMinScore retains extension candidates with a per-call score threshold.
// A chain can be shorter than the final HSP; apply its length limit after alignment.
// Neither override changes shared options or survives reuse of this chainer.
func (ce *Chainer2) chainWithMinScore(subs *[]SubstrPair, minimumScore int) (*[]*Chain2Result, int, int, int, int, int, int, int) {
	return ce.chainWithLimits(subs, minimumScore, 0)
}

// chainWithLimits applies thresholds to the unextended anchor chain.
func (ce *Chainer2) chainWithLimits(subs *[]SubstrPair, minimumScore, minAlignLen int) (*[]*Chain2Result, int, int, int, int, int, int, int) {
	n := len(*subs)
	if n == 0 {
		return nil, 0, 0, 0, 0, 0, 0, 0
	}

	if n == 1 {
		// A lone anchor can extend into a longer alignment. Leave its final
		// length and identity checks to the base-level alignment caller.
		sub := (*subs)[0]
		slen := int(sub.Len)
		paths := poolChains2.Get().(*[]*Chain2Result)

		path := poolChain2.Get().(*Chain2Result)
		path.Reset()

		qe := int(sub.QBegin) + slen - 1           // end
		te := int(sub.TBegin) + slen - 1           // end
		qb, tb := int(sub.QBegin), int(sub.TBegin) // in case there's only one anchor
		path.QBegin, path.QEnd = qb, qe
		path.TBegin, path.TEnd = tb, te
		path.MatchedBases = slen
		path.PIdent = 100
		path.AlignedBasesQ = slen
		path.AlignedBasesT = slen
		path.NAnchors++
		*paths = append(*paths, path)

		return paths, slen, slen, slen, qb, qe, tb, te
	}

	var i, j int
	// var k int
	bandBase := int32(ce.options.BandBase) // band size of banded-DP
	bandCount := ce.options.BandCount
	var _bCount int
	var _bBase int32

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

	usePotential := false
	maxPotential := ce.maxPotential[:0]
	// The legacy int32 gap arithmetic can overflow for widely separated
	// coordinates. Bound scores only when the gap stays nonnegative.
	minQ, maxQ := (*subs)[0].QBegin, (*subs)[0].QBegin
	minT, maxT := (*subs)[0].TBegin, (*subs)[0].TBegin
	querySorted, groupTargetSorted, repeatedQ := true, true, false
	equalLengths := true // one band-base test then applies to a whole Q group
	previousQ, previousT := (*subs)[0].QBegin, (*subs)[0].TBegin
	for index, sub := range *subs {
		if sub.Len != (*subs)[0].Len {
			equalLengths = false
		}
		if sub.QBegin < previousQ {
			querySorted = false
		}
		if index > 0 && sub.QBegin == previousQ && sub.TBegin < previousT {
			groupTargetSorted = false
		}
		if index > 0 && sub.QBegin == previousQ {
			repeatedQ = true
		}
		previousQ = sub.QBegin
		previousT = sub.TBegin
		minQ, maxQ = min(minQ, sub.QBegin), max(maxQ, sub.QBegin)
		minT, maxT = min(minT, sub.TBegin), max(maxT, sub.TBegin)
	}
	usePotential = int64(maxQ)-int64(minQ) <= math.MaxInt32/2 &&
		int64(maxT)-int64(minT) <= math.MaxInt32/2
	// The genome seed path orders equal-Q anchors by target position. Store
	// group boundaries only for such inputs; other callers keep the scan.
	var qGroupStarts []int32
	if querySorted && groupTargetSorted && repeatedQ && n <= thresholdNSubsLong {
		if cap(ce.qGroupStarts) < n {
			ce.qGroupStarts = make([]int32, n)
		}
		qGroupStarts = ce.qGroupStarts[:n]
		if cap(ce.qGroupCursors) < n {
			ce.qGroupCursors = make([]int32, n)
		}
		ce.qGroupCursors = ce.qGroupCursors[:n]
		start := 0
		for index, sub := range *subs {
			if index > 0 && sub.QBegin != (*subs)[index-1].QBegin {
				start = index
			}
			if index == start {
				ce.qGroupCursors[start] = int32(start)
			}
			qGroupStarts[index] = int32(start)
		}
	}
	if usePotential {
		maxPotential = append(maxPotential, int64((*subs)[0].Len))
	}

	// compute scores
	// Scores and gaps are integral. Keep intermediates wider than the packed
	// 32-bit score so comparisons retain their original pre-packing values.
	var s, m, g int64
	// The loop starts at anchor 1; include anchor 0 in the global maximum.
	// Independent tied chains retain their earlier endpoint.
	M := int64((*subs)[0].Len)
	// var d float64
	var mj, Mi int
	// Anchors stay in place during chaining; read by address to avoid value copies.
	var a, b *SubstrPair
	maxGap := int64(ce.options.MaxGap)
	// Only skip gap failures when coordinate arithmetic cannot overflow and
	// every member of a Q group has the same band-base stopping condition.
	skipGapGroups := qGroupStarts != nil && usePotential && equalLengths &&
		maxGap >= 0 && maxGap <= math.MaxInt32
	// maxDistance := float64(ce.options.MaxDistance)
	// (*scores)[0] = (*subs)[0].Len
	var aQBegin, aTBegin, bQBegin, bTBegin int32
	var qDiff, tDiff int32
	var sameQStart int

	// An earlier anchor can pass the gap test only if its T-Q diagonal is
	// within MaxGap of the current diagonal. Folding bins into a bitmap can
	// create false positives, never false negatives. An empty range therefore
	// proves that scanning predecessors cannot change this row's score or tie.
	// For increasing Q and T starts, the gap is exactly |(aT-aQ)-(bT-bQ)|:
	//
	// diagonal:    ... [ d-MaxGap ........ d+MaxGap ] ...
	// earlier bins:      empty    empty    empty
	//                            ^ current d
	//
	// No earlier diagonal in this range means no link can pass the gap test.
	// Small inputs keep the scan to avoid bitmap setup and lookup overhead.
	const diagonalBits = 1 << 16 // 8 KiB per active chainer, allocated on demand
	useDiagonals := n >= 512 && querySorted && equalLengths && usePotential && maxGap >= 0 && maxGap <= math.MaxInt32
	var diagonalSeen []uint64
	var diagonalShift uint
	if useDiagonals {
		if ce.diagonalSeen == nil {
			ce.diagonalSeen = make([]uint64, diagonalBits/64)
		}
		diagonalSeen = ce.diagonalSeen
		clear(diagonalSeen)
		diagonalShift = uint(bits.Len64(uint64(maxGap))) // bin width strictly exceeds MaxGap
		bin := uint64((int64((*subs)[0].TBegin)-int64((*subs)[0].QBegin))>>diagonalShift) & (diagonalBits - 1)
		diagonalSeen[bin>>6] |= uint64(1) << (bin & 63)
	}

	// Large repetitive inputs can have a high-scoring chain on an unrelated
	// diagonal. Its global prefix bound prevents early exits for other chains.
	// Instead, keep prefix maxima along linked lists of folded diagonal bins:
	//
	// bin A:  2 -> 8 -> 13    score maxima for A only
	// bin B:  4 -> 9 -> 15    score maxima for B only
	//
	// A predecessor can link only if its T-Q diagonal is within MaxGap of the
	// current diagonal. At most three bins overlap that interval. Collisions
	// only add impossible candidates, so the bound remains conservative.
	// Lists use original indices, preserving scan order, band counts and ties.
	// Extra storage is 4 bytes/anchor plus 256 KiB, local to this call so large
	// arrays are not retained in pooled chainers. Small/unsorted inputs keep
	// the existing bound; int32 indices and gap arithmetic must not overflow.
	const diagonalPotentialBins = 1 << 16
	useDiagonalPotential := n > thresholdNSubsLong && n <= math.MaxInt32 && querySorted && usePotential && maxGap >= 0 && maxGap <= math.MaxInt32
	var potentialHeads []int32 // latest scored anchor in each folded diagonal bin
	var potentialPrev []int32  // previous anchor in the same bin, or -1
	var potentialShift uint
	if useDiagonalPotential {
		potentialHeads = make([]int32, diagonalPotentialBins)
		for i := range potentialHeads {
			potentialHeads[i] = -1
		}
		potentialPrev = make([]int32, n)
		potentialShift = uint(bits.Len64(uint64(maxGap)))
		bin := uint64((int64((*subs)[0].TBegin)-int64((*subs)[0].QBegin))>>potentialShift) & (diagonalPotentialBins - 1)
		potentialHeads[bin] = 0
		potentialPrev[0] = -1
	}
	for i = 1; i < n; i++ {
		a = &(*subs)[i] // current seed/anchor
		// k = band * i   // index of current seed in the score matrix

		// just initialize the max score, which comes from the current seed
		m, mj = int64(a.Len), i
		// (*scores)[k] = m

		// Cache frequently accessed fields to reduce memory access overhead
		aQBegin, aTBegin = a.QBegin, a.TBegin
		if querySorted && aQBegin != (*subs)[i-1].QBegin {
			sameQStart = i
		}

		// this old one can't handle some frequent k-mers properly, which have many hits in other place.
		// for _b = 1; _b <= band; _b++ { // check previous $band seeds
		j = i
		if querySorted {
			// Equal query starts form one group and never contribute to a chain.
			j = sameQStart
		}

		if useDiagonals {
			diagonal := int64(aTBegin) - int64(aQBegin)
			lower, upper := (diagonal-maxGap)>>diagonalShift, (diagonal+maxGap)>>diagonalShift
			possible := false
			for bin := lower; bin <= upper; bin++ { // at most three adjacent bins
				folded := uint64(bin) & (diagonalBits - 1)
				if diagonalSeen[folded>>6]&(uint64(1)<<(folded&63)) != 0 {
					possible = true
					break
				}
			}
			if !possible {
				j = 0 // no predecessor can link; keep this anchor's own score
			}
			bin := uint64(diagonal>>diagonalShift) & (diagonalBits - 1)
			diagonalSeen[bin>>6] |= uint64(1) << (bin & 63)
		}

		var potentialCursors [3]int32
		nPotentialCursors := 0
		// The bound changes only after the scan passes a list head. Reuse it
		// between heads, but compare every updated score with the cached bound.
		potentialLimit := i
		var potentialBound int64
		if useDiagonalPotential {
			diagonal := int64(aTBegin) - int64(aQBegin)
			lower, upper := (diagonal-maxGap)>>potentialShift, (diagonal+maxGap)>>potentialShift
			for bin := lower; bin <= upper; bin++ {
				potentialCursors[nPotentialCursors] = potentialHeads[uint64(bin)&(diagonalPotentialBins-1)]
				nPotentialCursors++
			}
		}
		_bCount = 0
		// fmt.Printf("i: %d, %s\n", i, a.String())
		for {
			j--
			if j < 0 {
				break
			}
			// Each predecessor t can contribute at most score[t]+a.Len,
			// because the gap penalty is nonnegative. Store prefix maxima of
			// predecessor scores and add the current length at lookup time:
			//
			//   anchor         0      1      2      current
			//   Q/T           0/0   13/13  26/26     39/39
			//   score          13     26     39       ?
			//   current len                         31
			//
			// Linking current to anchor 2 gives m=39+31=70. Before checking
			// anchor 1, maxPotential[1]+31=57, so neither 1 nor 0 can match m.
			// Use strict > to preserve the existing preference for earlier ties.
			if useDiagonalPotential {
				if j < potentialLimit {
					potentialBound = 0
					potentialLimit = -1
					for z := 0; z < nPotentialCursors; z++ {
						head := potentialCursors[z]
						for head > int32(j) {
							head = potentialPrev[head]
						}
						potentialCursors[z] = head
						if head >= 0 {
							potentialBound = max(potentialBound, maxPotential[head])
							potentialLimit = max(potentialLimit, int(head))
						}
					}
				}
				// Strict inequality retains earlier predecessors on equal scores.
				if m > potentialBound+int64(a.Len) {
					break
				}
			} else if usePotential && m > maxPotential[j]+int64(a.Len) {
				break
			}

			b = &(*subs)[j] // previous seed/anchor
			// k++            // index of previous seed in the score matrix

			bQBegin, bTBegin = b.QBegin, b.TBegin

			// no need to compare
			// filter out messed/crossed anchors
			if bQBegin == aQBegin {
				continue
			}
			if bTBegin > aTBegin {
				if qGroupStarts != nil {
					start := int(qGroupStarts[j])
					if j-start >= 8 {
						// T starts are sorted within this Q group. Find its
						// last anchor that can precede the current target.
						j = targetUpperBound(*subs, start, j, aTBegin, &ce.qGroupCursors[start]) // next iteration inspects j-1
					}
				}
				continue
			}

			_bCount++

			_bBase = aQBegin - bQBegin - int32(b.Len)
			// fmt.Printf("  b: %d, j: %d, %s\n", _b, j, b.String())
			if !(_bBase <= bandBase || _bCount <= bandCount) {
				break
			}

			// d = distance2(a, b) + float64(b.Len)
			// if d > maxDistance { // can't do this
			// 	continue
			// }

			// g = gap2(a, b)
			// Inline gap2() using integer arithmetic, convert once at the end
			qDiff = aQBegin - bQBegin
			tDiff = aTBegin - bTBegin
			if qDiff < 0 {
				qDiff = -qDiff
			}
			if tDiff < 0 {
				tDiff = -tDiff
			}
			// Compute absolute difference without converting to float twice
			if qDiff > tDiff {
				g = int64(qDiff - tDiff)
			} else {
				g = int64(tDiff - qDiff)
			}

			if g > maxGap { // limit the gap. necessary?
				// Anchors in one Q group share bQ and have sorted bT values.
				// To connect b to a, the steps along Q and T must differ by
				// at most maxGap. For example, aQ=100, aT=1000, bQ=80,
				// maxGap=5: the Q step is 20, so the T step must be 15..25.
				// This means only bT=975..985 can pass the gap check:
				//
				// Q:       b(80) -------- 20 --------> a(100)
				// T:       b(?) ------- 15..25 ------> a(1000)
				//
				// bT, sorted:  900 970 | 975  980  985 | 990 1000
				//              too low | can pass gap  | too high
				//                      ^               ^
				//                    lower           upper
				//              <--------------------- scan direction
				//
				// On the right, binary-search past all values above upper.
				// Once below lower, skip the rest of this group: every
				// remaining bT is even smaller, making the T step too large.
				// We skip only anchors that would fail the same gap check,
				// so possible links and their score/tie order stay unchanged.
				// Skipped anchors must still count toward BandCount; equal
				// lengths also give this whole Q group the same BandBase test.
				if skipGapGroups {
					start := int(qGroupStarts[j])
					// Keep small groups on the scan to avoid search overhead.
					if j-start >= 8 {
						// The current anchor already passed the crossing and
						// band checks and was counted above.
						upper := int64(aTBegin) - int64(qDiff) + maxGap
						next := start
						if int64(bTBegin) > upper {
							low, high := start, j
							for low < high {
								mid := low + (high-low)/2
								if int64((*subs)[mid].TBegin) <= upper {
									low = mid + 1
								} else {
									high = mid
								}
							}
							next = low
						}
						// Below the lower bound, the rest of this sorted
						// group also fails. Count every skipped anchor and
						// stop here if the original scan would exhaust its band.
						_bCount += j - next
						if _bBase > bandBase && _bCount > bandCount {
							break
						}
						j = next // the next iteration inspects next-1
					}
				}
				continue
			}

			// The predecessor score already includes b. Add a's length so a
			// 13 -> 31 chain scores 44 in either orientation, rather than 26/62.
			s = int64((*maxscoresIdxs)[j]>>32) + int64(a.Len) - g

			if s >= m { // update the max score of current seed/anchor
				m = s
				mj = j
			}
		}

		// *maxscores = append(*maxscores, m)          // save the max score of the whole
		// *maxscoresIdxs = append(*maxscoresIdxs, mj) // save where the max score comes from
		*maxscoresIdxs = append(*maxscoresIdxs, uint64(m)<<32|uint64(mj))
		if usePotential {
			potential := int64(uint32(m))
			if useDiagonalPotential {
				bin := uint64((int64(aTBegin)-int64(aQBegin))>>potentialShift) & (diagonalPotentialBins - 1)
				head := potentialHeads[bin]
				potentialPrev[i] = head
				if head >= 0 {
					potential = max(potential, maxPotential[head])
				}
				maxPotential = append(maxPotential, potential)
				potentialHeads[bin] = int32(i)
			} else {
				maxPotential = append(maxPotential, max(potential, maxPotential[i-1]))
			}
		}

		if m > M || m == M && mj == Mi && extendsChain2End((*subs)[Mi], *a) {
			M, Mi = m, i
		}
	}
	ce.maxPotential = maxPotential

	// fmt.Printf("M: %f, Mi: %d\n", M, Mi)

	// // print the score matrix
	// fmt.Printf("i\tpair-i\tiMax\tj:scores\n")
	// var _mi uint64
	// for i = 0; i < n; i++ {
	// 	_mi = (*maxscoresIdxs)[i]
	// 	fmt.Printf("%d\t%s\t%d:%d", i, (*subs)[i], _mi&4294967295, _mi>>32)
	// 	fmt.Printf("\n")
	// }

	// backtrack

	minScore := int64(minimumScore)

	// check the highest score, for early quit,
	if M < minScore {
		return nil, 0, 0, 0, 0, 0, 0, 0
	}

	paths := poolChains2.Get().(*[]*Chain2Result)

	var nMatchedBases, nAlignedBasesQ, nAlignedBasesT int
	// ce.bounds = ce.bounds[:0]

	_, qB, qE, tB, tE := chainARegion(
		subs,
		// maxscores,
		maxscoresIdxs,
		0,
		float64(minScore),
		minAlignLen,
		ce.options.MinIdentity,
		paths,
		&nMatchedBases,
		&nAlignedBasesQ,
		&nAlignedBasesT,
		Mi,
		nil, // &ce.bounds,
		ce.options.HeuristicKmerPidentThreshold,
	)

	if len(*paths) == 0 {
		recycleChaining2ResultSlice(paths)
		return nil, 0, 0, 0, 0, 0, 0, 0
	}

	return paths, nMatchedBases, nAlignedBasesQ, nAlignedBasesT, qB, qE, tB, tE
}

// extendsChain2End keeps a same-score continuation of the selected endpoint.
// A 31-base anchor followed by 25 bases across a 25-base gap still scores 31:
//
//	query:   [31] -------- [25]
//	target:  [31] ----------- [25]
//
// Keeping only the first anchor loses the terminal region in one orientation.
// Callers require a direct predecessor link, preserving ties between separate
// chains. Neither sequence end may shrink, including on unsorted inputs.
func extendsChain2End(previous, current SubstrPair) bool {
	qEnd := int64(current.QBegin) + int64(current.Len)
	tEnd := int64(current.TBegin) + int64(current.Len)
	previousQEnd := int64(previous.QBegin) + int64(previous.Len)
	previousTEnd := int64(previous.TBegin) + int64(previous.Len)
	return qEnd >= previousQEnd && tEnd >= previousTEnd &&
		(qEnd > previousQEnd || tEnd > previousTEnd)
}

// chain2PeakIndex accelerates repeated maximum-score queries on disjoint regions.
// It is built only after scans exceed four full passes. Most calls finish with
// one or a few chains and allocate no index. Scratch is local to this call.
// A tree over 64-anchor blocks costs at most 1/4 byte per anchor; same-score
// successor links cost 4 bytes per anchor. Boundary blocks are scanned directly.
type chain2PeakIndex struct {
	subs    []SubstrPair // immutable anchors for endpoint-extension ties
	scores  []uint64     // packed scores and global predecessor indexes
	offset  int          // global index of subs[0]
	scanned int          // total entries scanned before constructing the tree
	tree    []uint32     // index+1 of the first maximum in each block-tree node; 0 is empty
	next    []uint32     // index+1 of the first eligible same-score extending child
	size    int          // number of leaves in the block tree
}

const chain2PeakBlockSize = 64 // bound direct boundary scans while keeping the tree small

// better returns the higher score, retaining the earlier index on ties.
func (p *chain2PeakIndex) better(a, b uint32) uint32 {
	if a == 0 {
		return b
	}
	if b == 0 {
		return a
	}
	x, y := p.scores[a-1]>>32, p.scores[b-1]>>32
	if x > y || x == y && a < b {
		return a
	}
	return b
}

// build records first maxima and same-score successors in one local allocation.
func (p *chain2PeakIndex) build() {
	n := len(p.scores)
	blocks := (n + chain2PeakBlockSize - 1) / chain2PeakBlockSize
	p.size = 1 << bits.Len(uint(blocks-1))
	storage := make([]uint32, n+2*p.size)
	p.next, p.tree = storage[:n], storage[n:]
	for i, packed := range p.scores {
		leaf := p.size + i/chain2PeakBlockSize
		p.tree[leaf] = p.better(p.tree[leaf], uint32(i+1))
		previous := int(uint32(packed)) - p.offset
		if previous >= 0 && previous < i && p.next[previous] == 0 &&
			packed>>32 == p.scores[previous]>>32 && extendsChain2End(p.subs[previous], p.subs[i]) {
			p.next[previous] = uint32(i + 1)
		}
	}
	for i := p.size - 1; i > 0; i-- {
		p.tree[i] = p.better(p.tree[i*2], p.tree[i*2+1])
	}
}

// best preserves the scan's tie rule: start at the first maximum and follow
// its first same-score extending child, then that child's first child, etc.
// A later sibling cannot replace the already selected child:
//
//	root -> child A -> grandchild
//	     -> child B              (ignored after A)
//
// These links belong to the selected backtracking path, so following them
// across all disjoint regions visits at most n anchors in total.
func (p *chain2PeakIndex) best(begin, end int) (int, float64) {
	if p.tree == nil && p.scanned >= 4*len(p.scores) && end-begin >= chain2PeakBlockSize && uint64(len(p.scores)) <= math.MaxUint32 {
		p.build()
	}
	if p.tree == nil {
		p.scanned += end - begin
		best, score := begin, uint64(0)
		for i := begin; i < end; i++ {
			packed := p.scores[i]
			m := packed >> 32
			if m > score || m == score && int(uint32(packed)) == p.offset+best && extendsChain2End(p.subs[best], p.subs[i]) {
				best, score = i, m
			}
		}
		return best, float64(score)
	}
	var best uint32
	left := min((begin+chain2PeakBlockSize-1)/chain2PeakBlockSize*chain2PeakBlockSize, end)
	for i := begin; i < left; i++ {
		best = p.better(best, uint32(i+1))
	}
	right := max(left, end/chain2PeakBlockSize*chain2PeakBlockSize)
	for i := right; i < end; i++ {
		best = p.better(best, uint32(i+1))
	}
	for lo, hi := p.size+left/chain2PeakBlockSize, p.size+right/chain2PeakBlockSize; lo < hi; lo, hi = lo/2, hi/2 {
		if lo&1 != 0 {
			best = p.better(best, p.tree[lo])
			lo++
		}
		if hi&1 != 0 {
			hi--
			best = p.better(best, p.tree[hi])
		}
	}
	at := int(best) - 1
	for next := p.next[at]; next != 0 && int(next)-1 < end; next = p.next[at] {
		at = int(next) - 1
	}
	return at, float64(p.scores[at] >> 32)
}

// chainARegion visits the selected chain, then right and left regions, retaining
// the original recursive output order and bounds even for filtered-out chains.
// An explicit stack prevents deep recursion for many independent short chains.
func chainARegion(subs *[]SubstrPair, maxscoresIdxs *[]uint64, offset int,
	minScore float64, minAlignLen int, _ float64, paths *[]*Chain2Result,
	_nMatchedBases, _nAlignedBasesQ, _nAlignedBasesT *int, Mi0 int, _ *[]int32,
	heuristicKmerPidentThreshold float64,
) (float64, int, int, int, int) {
	// Regions use local half-open indexes; packed predecessors remain global.
	type region struct{ begin, end, best int }
	var scratch [32]region
	pending := append(scratch[:0], region{0, len(*subs), Mi0})
	peaks := chain2PeakIndex{subs: *subs, scores: *maxscoresIdxs, offset: offset}
	var rootScore float64
	var qB, qE, tB, tE int
	first := true
	for len(pending) > 0 {
		r := pending[len(pending)-1]
		pending = pending[:len(pending)-1]
		mi, score := r.best, float64(0)
		if mi < 0 {
			mi, score = peaks.best(r.begin, r.end)
			if score < minScore {
				if first {
					return 0, -1, -1, -1, -1
				}
				continue
			}
		}
		if first {
			rootScore = score
		}

		var qb, qe, tb, te int32
		var matched, nAnchors, nextBegin int
		complete := false // the predecessor path starts inside this region
		i := mi
		for {
			previous := int(uint32((*maxscoresIdxs)[i])) - offset
			if previous < r.begin {
				break
			}
			sub := (*subs)[i]
			if nAnchors == 0 {
				qe, te = sub.QBegin+int32(sub.Len)-1, sub.TBegin+int32(sub.Len)-1
				matched = int(sub.Len)
			} else if int(sub.QBegin)+int(sub.Len)-1 >= nextBegin {
				matched += nextBegin - int(sub.QBegin)
			} else {
				matched += int(sub.Len)
			}
			qb, tb = sub.QBegin, sub.TBegin
			nextBegin = int(sub.QBegin)
			nAnchors++
			if i == previous {
				complete = true
				break
			}
			i = previous
		}
		if nAnchors > 0 {
			alignedQ, alignedT := int(qe)-int(qb)+1, int(te)-int(tb)+1
			if alignedQ >= minAlignLen {
				identity := float64(matched) / float64(max(alignedQ, alignedT)) * 100
				keep := identity >= heuristicKmerPidentThreshold
				if complete { // Preserve the original comparison for zero-length paths too.
					keep = !(identity < heuristicKmerPidentThreshold)
				}
				if keep {
					path := poolChain2.Get().(*Chain2Result)
					path.Reset()
					path.NAnchors, path.MatchedBases = nAnchors, matched
					path.AlignedBasesQ, path.AlignedBasesT = alignedQ, alignedT
					path.PIdent = min(identity, 100)
					path.QBegin, path.QEnd, path.TBegin, path.TEnd = int(qb), int(qe), int(tb), int(te)
					*paths = append(*paths, path)
					*_nMatchedBases += matched
					*_nAlignedBasesQ += alignedQ
					*_nAlignedBasesT += alignedT
				}
			}
		}
		if first {
			qB, qE, tB, tE = int(qb), int(qe), int(tb), int(te)
			first = false
		} else if score > 0 {
			qB, qE = min(qB, int(qb)), max(qE, int(qe))
			tB, tE = min(tB, int(tb)), max(tE, int(te))
		}
		// Push left first so the right region is processed first, as before.
		if i > r.begin {
			pending = append(pending, region{r.begin, i, -1})
		}
		if mi+1 < r.end {
			pending = append(pending, region{mi + 1, r.end, -1})
		}
	}
	return rootScore, qB, qE, tB, tE
}

func distance2(a, b SubstrPair) float64 {
	return math.Max(math.Abs(float64(a.QBegin-b.QBegin)), math.Abs(float64(a.TBegin-b.TBegin)))
}

func gap2(a, b SubstrPair) float64 {
	return math.Abs(math.Abs(float64(a.QBegin-b.QBegin)) - math.Abs(float64(a.TBegin-b.TBegin)))
}

func gapScore2(gap float64) float64 {
	// if gap == 0 {
	// 	return 0
	// }
	return 0.1 * gap
}
