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

	rtree "github.com/shenwei356/LexicMap/lexicmap/cmd/tree"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
	"github.com/shenwei356/lexichash/iterator"
)

// SeqComparatorOptions contains options for comparing two sequences.
type SeqComparatorOptions struct {
	// indexing
	K         uint8
	MinPrefix uint8

	// chaining
	Chaining2Options

	// seq similarity
	MinAlignedFraction float64 // minimum query aligned fraction in a HSP

	MinIdentity float64
}

// DefaultSeqComparatorOptions contains the default options for SeqComparatorOptions.
var DefaultSeqComparatorOptions = SeqComparatorOptions{
	K:         32,
	MinPrefix: 11, // can not be too small, or there will be a large number of anchors.

	Chaining2Options: Chaining2Options{
		// should be relative small
		MaxGap: 50,
		// better be larger than MinPrefix
		MinScore: 50,
		// can not be < k
		MaxDistance: 50,
		// can not be two small
		BandBase: 100,
	},

	MinAlignedFraction: 0,
}

// SeqComparator is for fast and accurate similarity estimation of two sequences,
// which are in the same strand (important).
type SeqComparator struct {
	// options
	options *SeqComparatorOptions

	// chainer for chaining anchors,
	// shared variable-length substrings searched by prefix matching.
	// chainer *Chainer2
	poolChainers *sync.Pool

	// each comparator has its own pool
	// poolSub *sync.Pool

	// a prefix tree for matching k-mers
	tree *rtree.Tree

	ccc, ggg, ttt uint64
}

// NewSeqComparator creates a new SeqComparator with given options.
// No options checking now.
func NewSeqComparator(options *SeqComparatorOptions, poolChainers *sync.Pool) *SeqComparator {
	cpr := &SeqComparator{
		options: options,
		// poolChainers: &sync.Pool{New: func() interface{} {
		// 	return NewChainer2(&options.Chaining2Options)
		// }},
		poolChainers: poolChainers,

		// poolSub: &sync.Pool{New: func() interface{} {
		// 	return &SubstrPair{}
		// }},

		ccc: util.Ns(0b01, options.K),
		ggg: util.Ns(0b10, options.K),
		ttt: util.Ns(0b11, options.K),
	}

	return cpr
}

// Index initializes the SeqComparator with the query sequence.
func (cpr *SeqComparator) Index(s []byte) error {
	k := cpr.options.K
	k8 := uint8(k)

	// k-mer iterator
	iter, err := iterator.NewKmerIterator(s, int(k))
	if err != nil {
		return err
	}

	// a reusable Radix tree for searching k-mers sharing at least n-base prefixes.
	t := rtree.NewTree(k)

	var kmer, kmerRC uint64
	var ok bool

	ccc := cpr.ccc
	ggg := cpr.ggg
	ttt := cpr.ttt

	for {
		kmer, kmerRC, ok, _ = iter.NextKmer()
		if !ok {
			break
		}

		// fmt.Printf("%d: %s\n", iter.Index(), lexichash.MustDecode(kmer, k))
		if kmer == 0 || kmer == ccc || kmer == ggg || kmer == ttt ||
			util.IsLowComplexityDust(kmer, k8) {
			continue
		}

		t.Insert(kmer, uint32(iter.Index()<<1))
		t.Insert(kmerRC, uint32(iter.Index()<<1|1))
	}

	cpr.tree = t

	return nil
}

// SeqComparatorResult contains the details of a seq comparison result.
type SeqComparatorResult struct {
	AlignedBases int // The number of aligned bases.

	AlignedFraction float64 // query (original query) coverage per HSP

	MatchedBases int
	PIdent       float64

	QueryLen int // length of the original query, used to compute/update AlignedFraction

	QBegin int
	QEnd   int
	TBegin int
	TEnd   int

	TSeq []byte // target seq

	Chains *[]*Chain2Result
}

// Update updates the data with new chains.
// However it does not considerate gaps.
func (r *SeqComparatorResult) Update(chains *[]*Chain2Result, queryLen int) {
	r.Chains = chains
	r.QueryLen = queryLen

	// alignment regions might have overlap

	r.QBegin, r.TBegin = math.MaxInt, math.MaxInt
	r.QEnd, r.TEnd = -1, -1
	r.MatchedBases = 0

	regions := poolRegions.Get().(*[]*[2]int)
	for _, c := range *r.Chains {
		// fmt.Printf("to merge [%d, %d] vs [%d, %d]\n", c.QBegin, c.QEnd, c.TBegin, c.TEnd)
		region := poolRegion.Get().(*[2]int)
		region[0], region[1] = c.QBegin, c.QEnd
		*regions = append(*regions, region)

		if c.QBegin < r.QBegin {
			r.QBegin = c.QBegin
		}
		if c.TBegin < r.TBegin {
			r.TBegin = c.TBegin
		}
		if c.QEnd > r.QEnd {
			r.QEnd = c.QEnd
		}
		if c.TEnd > r.TEnd {
			r.TEnd = c.TEnd
		}
		r.MatchedBases += c.MatchedBases
	}
	r.AlignedBases = coverageLen(regions)
	recycleRegions(regions)

	r.AlignedFraction = float64(r.AlignedBases) / float64(queryLen) * 100
	r.PIdent = float64(r.MatchedBases) / float64(r.AlignedBases) * 100
	if r.PIdent > 100 {
		r.PIdent = 100
	}
}

// Update2 only compute the aligned fraction for all chains
func (r *SeqComparatorResult) Update2(chains *[]*Chain2Result, queryLen int) {
	r.Chains = chains
	r.QueryLen = queryLen

	// alignment regions might have overlap

	r.QBegin, r.TBegin = math.MaxInt, math.MaxInt
	r.QEnd, r.TEnd = -1, -1
	r.MatchedBases = 0

	regions := poolRegions.Get().(*[]*[2]int)
	for _, c := range *r.Chains {
		// fmt.Printf("to merge [%d, %d] vs [%d, %d]\n", c.QBegin, c.QEnd, c.TBegin, c.TEnd)

		c.AlignedFraction = float64(c.AlignedBasesQ) / float64(queryLen) * 100

		region := poolRegion.Get().(*[2]int)
		region[0], region[1] = c.QBegin, c.QEnd
		*regions = append(*regions, region)

		r.MatchedBases += c.MatchedBases
	}

	r.AlignedBases = coverageLen(regions)
	recycleRegions(regions)

	r.AlignedFraction = float64(r.AlignedBases) / float64(queryLen) * 100
	r.PIdent = float64(r.MatchedBases) / float64(r.AlignedBases) * 100
	if r.PIdent > 100 {
		r.PIdent = 100
	}
}

// -----------------------------------------------------------------------------

func recycleRegions(regions *[]*[2]int) {
	for _, r := range *regions {
		poolRegion.Put(r)
	}
	*regions = (*regions)[:0]
	poolRegions.Put(regions)
}

// coverageLen computes the total covered bases for a list of regions which might have overlaps.
func coverageLen(regions *[]*[2]int) (r int) {
	if len(*regions) == 0 {
		return 0
	}
	if len(*regions) == 1 {
		return (*regions)[0][1] - (*regions)[0][0] + 1
	}

	// sort by the start locations
	// sort.Slice(*regions, func(i, j int) bool {
	// 	return (*regions)[i][0] < (*regions)[j][0]
	// })
	slices.SortFunc(*regions, func(a, b *[2]int) int {
		return a[0] - b[0]
	})

	var region *[2]int // ccurent region
	var start, end int // positions of the a merged region

	region = (*regions)[0] // the first region
	start, end = (*region)[0], (*region)[1]

	for i := 1; i < len(*regions); i++ {
		region = (*regions)[i]
		if region[0] > end { // has no overlap with previous merged region
			r += end - start + 1 // add the length

			start, end = (*region)[0], (*region)[1] // create a new merged region
			continue
		}
		if region[1] <= end { // the current region is in the merged region
			continue
		}
		end = region[1] // merge the current region
	}
	r += end - start + 1 // add the length

	return r
}

var poolRegions = &sync.Pool{New: func() interface{} {
	tmp := make([]*[2]int, 0, 128)
	return &tmp
}}

var poolRegion = &sync.Pool{New: func() interface{} {
	return &[2]int{}
}}

var poolSeqComparatorResult = &sync.Pool{New: func() interface{} {
	return &SeqComparatorResult{}
}}

// RecycleSeqComparatorResult recycles a SeqComparatorResult
func RecycleSeqComparatorResult(r *SeqComparatorResult) {
	RecycleChaining2Result(r.Chains)
	r.Chains = nil
	poolSeqComparatorResult.Put(r)
}

// Compare matchs k-mers for the query sequence (begin: end), chains them up,
// and computes the similarity.
// Please remember to call RecycleSeqComparatorResult() to recycle the result.
func (cpr *SeqComparator) Compare(begin, end uint32, s []byte, queryLen int) (*SeqComparatorResult, error) {
	k8 := cpr.options.K
	k := int(k8)

	m := cpr.options.MinPrefix // 11
	if len(s) >= 1000000 {
		m += 8 // 19
	} else if len(s) >= 250000 {
		m += 6 // 17
	} else if len(s) >= 50000 {
		m += 4 // 15
	} else if len(s) >= 10000 {
		m += 2 // 13
	}

	// --------------------------------------------------------------
	// search on the tree

	// ----------- round 1: k = k ------------------
	iter, err := iterator.NewKmerIterator(s, k)
	if err != nil {
		return nil, err
	}

	t := cpr.tree
	var kmer, kmerRC uint64
	var ok bool
	var v, p uint32
	var srs *[]*rtree.SearchResult
	var sr *rtree.SearchResult

	// poolSub2 := cpr.poolSub
	poolSub2 := poolSub

	// substring pairs/seeds/anchors
	subs := poolSubsLong.Get().(*[]*SubstrPair)
	defer RecycleSubstrPairs(poolSub2, poolSubsLong, subs)

	ccc := cpr.ccc
	ggg := cpr.ggg
	ttt := cpr.ttt

	for {
		kmer, kmerRC, ok, _ = iter.NextKmer()
		if !ok {
			break
		}

		// if kmer == 0 || kmer == ccc || kmer == ggg || kmer == ttt ||
		// 	util.IsLowComplexityDust(kmer, k8) {
		if kmer == 0 || kmer == ccc || kmer == ggg || kmer == ttt {
			continue
		}

		// ------------ positive strand -----------

		srs, ok = t.Search(kmer, m)
		if ok {
			// fmt.Printf("%d: %s\n", iter.Index(), lexichash.MustDecode(kmer, k8))
			for _, sr = range *srs {
				for _, v = range sr.Values {
					p = v >> 1
					// fmt.Printf("  p: %d, len: %d\n", p, sr.LenPrefix)
					if v&1 == 1 || p < begin || p+uint32(sr.LenPrefix) > end { // skip flanking regions
						continue
					}

					_sub2 := poolSub2.Get().(*SubstrPair)
					_sub2.QBegin = int32(p)
					_sub2.TBegin = int32(iter.Index())
					// _sub2.Code = rtree.KmerPrefix(sr.Kmer, k8, sr.LenPrefix)
					_sub2.Len = uint8(sr.LenPrefix)
					_sub2.QRC = false
					_sub2.TRC = false

					*subs = append(*subs, _sub2)
				}
			}

			t.RecycleSearchResult(srs)
		}

		// ------------ negative strand -----------

		srs, ok = t.Search(kmerRC, m)
		if ok {
			// fmt.Printf("%d: %s\n", iter.Index(), lexichash.MustDecode(kmerRC, k8))
			for _, sr = range *srs {
				for _, v = range sr.Values {
					p = v>>1 + uint32(k) - uint32(sr.LenPrefix)
					// fmt.Printf("  p: %d, len: %d\n", p, sr.LenPrefix)
					if v&1 == 0 || p+uint32(sr.LenPrefix) < begin || p > end { // skip flanking regions
						continue
					}

					_sub2 := poolSub2.Get().(*SubstrPair)
					_sub2.QBegin = int32(p)
					_sub2.TBegin = int32(iter.Index() + k - int(sr.LenPrefix))
					// _sub2.Code = rtree.KmerPrefix(sr.Kmer, k8, sr.LenPrefix)
					_sub2.Len = uint8(sr.LenPrefix)
					_sub2.QRC = true
					_sub2.TRC = true

					*subs = append(*subs, _sub2)
				}
			}

			t.RecycleSearchResult(srs)
		}
	}

	if len(*subs) < 1 { // no way, no matches in the pseudo alignment
		return nil, nil
	}

	// --------------------------------------------------------------
	// clear matched substrings

	if len(*subs) > 1 {
		ClearSubstrPairs(poolSub2, subs, k)
	}

	// fmt.Println("----------- cleared anchors ----------")
	// for i, sub := range *subs {
	// 	fmt.Printf("%3d: %s\n", i, sub)
	// }
	// fmt.Println("-------------------------------")

	TrimSubStrPairs(poolSub2, subs, k, 100)
	if len(*subs) == 0 {
		return nil, nil
	}

	// fmt.Println("----------- cleaned anchors ----------")
	// for i, sub := range *subs {
	// 	fmt.Printf("%3d: %s\n", i, sub)
	// }
	// fmt.Println("-------------------------------")

	// --------------------------------------------------------------
	// chaining paired substrings

	chainer := cpr.poolChainers.Get().(*Chainer2)
	chains, _, nAlignedBasesQ, _, _, _, _, _ := chainer.Chain(subs)
	defer func() {
		cpr.poolChainers.Put(chainer)
	}()
	if chains == nil {
		return nil, nil
	}

	af := float64(nAlignedBasesQ) / float64(queryLen) * 100
	// if af < cpr.options.MinAlignedFraction {
	// 	RecycleChaining2Result(chains)
	// 	return nil, nil
	// }

	// result object
	r := poolSeqComparatorResult.Get().(*SeqComparatorResult)
	r.AlignedBases = nAlignedBasesQ
	// r.MatchedBases = nMatchedBases
	r.QueryLen = queryLen
	r.AlignedFraction = af

	// very important
	if len(*chains) > 1 {
		// sort.Slice(*chains, func(i, j int) bool {
		// 	return (*chains)[i].QBegin <= (*chains)[j].QBegin
		// })
		slices.SortFunc(*chains, func(a, b *Chain2Result) int {
			return a.QBegin - b.QBegin
		})
	}

	// fmt.Println("chain2:")
	// for c, chain := range *chains {
	// 	fmt.Printf("  chain2: %d, [%d, %d] vs [%d, %d]\n", c, chain.QBegin, chain.QEnd, chain.TBegin, chain.TEnd)
	// }

	// r.QBegin = qB
	// r.QEnd = qE
	// r.TBegin = tB
	// r.TEnd = tE
	r.Chains = chains

	return r, nil
}

// RecycleIndex recycles the Index (tree data).
// Please call this if you do not need the comparator anymore.
func (cpr *SeqComparator) RecycleIndex() {
	if cpr.tree != nil {
		rtree.RecycleTree(cpr.tree)
	}
}

// TrimSubStrPairs trims anchors for query/subjects with tandem repeats in either end.
//
// case 1: embeded anchor in query/target
//
//	61: 156-186 (+) vs 1163-1193 (+), len:31
//	62: 157-187 (-) vs 1164-1194 (-), len:31
//	63: 158-188 (+) vs 1165-1195 (+), len:31
//	64: 168-195 (-) vs 1168-1195 (-), len:28
//	65: 175-202 (-) vs 1168-1195 (-), len:28 <---
//	66: 182-209 (-) vs 1168-1195 (-), len:28 <---
//	67: 189-216 (-) vs 1168-1195 (-), len:28 <---
//	68: 196-223 (+) vs 1168-1195 (+), len:28 <---
//	69: 203-230 (+) vs 1168-1195 (+), len:28 <---
//	70: 210-237 (-) vs 1168-1195 (-), len:28 <---
//	71: 217-244 (-) vs 1168-1195 (-), len:28 <--- gap=7, overlap=28 (28/28)
//
// case 2: big overlap + big gap
//
//	727: 789-819 (-) vs 789-819 (-), len:31
//	728: 790-820 (-) vs 790-820 (-), len:31
//	729: 804-821 (-) vs 821-838 (-), len:18 <--- gap=17, overlap=17 (17/18)
func TrimSubStrPairs(poolSub *sync.Pool, subs *[]*SubstrPair, k int, minDist float64) {
	if len(*subs) < 2 {
		return
	}

	var _p, p *SubstrPair
	var i int
	last := len(*subs) - 1

	// head
	_p = (*subs)[0]
	start := 0

	for i, p = range (*subs)[1:] {
		if distance(p, _p) < minDist &&
			((p.QBegin == _p.QBegin || p.TBegin == _p.TBegin) || // case 1
				(gap2(_p, p) > 11 && float64(overlap(_p, p))/float64(_p.Len) > 0.8)) { // case 2
			start = i
			_p = p
			continue
		}

		break
	}

	// tail
	_p = (*subs)[last]
	end := last

	for i = len(*subs) - 2; i >= 0; i-- {
		p = (*subs)[i]

		if distance(p, _p) < minDist &&
			((p.QBegin == _p.QBegin || p.TBegin == _p.TBegin) || // case 1
				(gap2(p, _p) > 11 && float64(overlap(p, _p))/float64(_p.Len) > 0.8)) { // case 2
			end = i
			_p = p
			continue
		}

		break
	}

	// fmt.Printf("start: %d, end: %d\n", start, end)

	if start == 0 && end == len(*subs) {
		return
	}

	// remove head overhang
	for i = 0; i < start; i++ {
		poolSub.Put((*subs)[i])
		(*subs)[i] = nil
	}

	// remove tail overhang
	for i = end + 1; i <= last; i++ {
		if (*subs)[i] != nil {
			poolSub.Put((*subs)[i])
			(*subs)[i] = nil
		}
	}

	if start >= end { // all discarded
		*subs = (*subs)[:0]
	} else {
		*subs = (*subs)[start : end+1]
	}
}

// a should be in front of b
func overlap(a, b *SubstrPair) int32 {
	var qo, to int32
	if b.QBegin >= a.QBegin && b.QBegin <= a.QBegin+int32(a.Len) {
		qo = a.QBegin + int32(a.Len) - b.QBegin + 1
	}

	if b.TBegin >= a.TBegin && b.TBegin <= a.TBegin+int32(a.Len) {
		to = a.TBegin + int32(a.Len) - b.TBegin + 1
	}
	return max(qo, to)
}
