// Copyright © 2023-2024 Wei Shen <sheTopLeftei356@gmail.com>
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
	"sort"
	"sync"

	rtree "github.com/shenwei356/LexicMap/lexicmap/cmd/tree"
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
		Band: 100,
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

	// a prefix tree for matching k-mers
	tree *rtree.Tree
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
	}

	return cpr
}

// Index initializes the SeqComparator with the query sequence.
func (cpr *SeqComparator) Index(s []byte) error {
	k := cpr.options.K

	// k-mer iterator
	iter, err := iterator.NewKmerIterator(s, int(k))
	if err != nil {
		return err
	}

	// a reusable Radix tree for searching k-mers sharing at least n-base prefixes.
	t := rtree.NewTree(k)

	// only considering the positive strand
	var kmer uint64
	var ok bool
	ttt := (uint64(1) << (k << 1)) - 1

	for {
		kmer, ok, _ = iter.NextPositiveKmer()
		if !ok {
			break
		}

		// fmt.Printf("%d: %s\n", iter.Index(), lexichash.MustDecode(kmer, k))
		if kmer == 0 || kmer == ttt { // skip AAAAAAAAAA and TTTTTTTTT
			continue
		}

		t.Insert(kmer, uint32(iter.Index()))
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
	*regions = (*regions)[:0]
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
	*regions = (*regions)[:0]
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
	sort.Slice(*regions, func(i, j int) bool {
		return (*regions)[i][0] < (*regions)[j][0]
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
	poolSeqComparatorResult.Put(r)
}

// Compare matchs k-mers for the query sequence (begin: end), chains them up,
// and computes the similarity.
// Please remember to call RecycleSeqComparatorResult() to recycle the result.
func (cpr *SeqComparator) Compare(begin, end uint32, s []byte, queryLen int) (*SeqComparatorResult, error) {
	k8 := cpr.options.K
	k := int(k8)
	m := cpr.options.MinPrefix

	// --------------------------------------------------------------
	// search on the tree

	iter, err := iterator.NewKmerIterator(s, k)
	if err != nil {
		return nil, err
	}

	t := cpr.tree
	var kmer uint64
	var ok bool
	var v uint32
	var srs *[]*rtree.SearchResult
	var sr *rtree.SearchResult

	// substring pairs/seeds/anchors
	subs := poolSubs.Get().(*[]*SubstrPair)
	*subs = (*subs)[:0]

	// only considering k-mers on the positive strand.
	// how can we detect inversion?
	//	-----> <====== ----->
	//	||||||         ||||||
	//	-----> ======> ----->
	ttt := (uint64(1) << (k << 1)) - 1
	for {
		kmer, ok, _ = iter.NextPositiveKmer()
		if !ok {
			break
		}

		if kmer == 0 || kmer == ttt { // skip AAAAAAAAAA and TTTTTTTTT
			continue
		}

		srs, ok = t.Search(kmer, m)
		if !ok {
			continue
		}
		for _, sr = range *srs {
			for _, v = range sr.Values {
				if v+uint32(sr.LenPrefix) < begin || v > end { // skip flanking regions
					continue
				}

				_sub2 := poolSub.Get().(*SubstrPair)
				_sub2.QBegin = int32(v)
				_sub2.TBegin = int32(iter.Index())
				// _sub2.Code = rtree.KmerPrefix(sr.Kmer, k8, sr.LenPrefix)
				_sub2.Len = uint8(sr.LenPrefix)
				_sub2.TRC = false

				*subs = append(*subs, _sub2)
			}
		}
		t.RecycleSearchResult(srs)
	}

	if len(*subs) < 1 { // no way
		return nil, err
	}

	// --------------------------------------------------------------
	// clear matched substrings

	ClearSubstrPairs(subs, k)

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
		RecycleSubstrPairs(subs)
		return nil, nil
	}

	af := float64(nAlignedBasesQ) / float64(queryLen) * 100
	if af < cpr.options.MinAlignedFraction {
		RecycleChaining2Result(chains)
		RecycleSubstrPairs(subs)
		return nil, nil
	}

	// result object
	r := poolSeqComparatorResult.Get().(*SeqComparatorResult)
	r.AlignedBases = nAlignedBasesQ
	// r.MatchedBases = nMatchedBases
	r.QueryLen = queryLen
	r.AlignedFraction = af

	// very important
	sort.Slice(*chains, func(i, j int) bool {
		return (*chains)[i].QBegin <= (*chains)[j].QBegin
	})

	// fmt.Println("chain2:")
	// for c, chain := range *chains {
	// 	fmt.Printf("  chain2: %d, [%d, %d] vs [%d, %d]\n", c, chain.QBegin, chain.QEnd, chain.TBegin, chain.TEnd)
	// }

	// r.QBegin = qB
	// r.QEnd = qE
	// r.TBegin = tB
	// r.TEnd = tE
	r.Chains = chains

	RecycleSubstrPairs(subs)

	return r, nil
}

// RecycleIndex recycles the Index (tree data).
// Please call this if you do not need the comparator anymore.
func (cpr *SeqComparator) RecycleIndex() {
	if cpr.tree != nil {
		rtree.RecycleTree(cpr.tree)
	}
}
