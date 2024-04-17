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
	MinIdentity        float64 // minimum percentage base identity of a segment
	MinSegmentLength   int     // minimum length of a segment
	MinAlignedFraction float64 // minimum query aligned fraction in a chain
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

	MinIdentity:        70,
	MinSegmentLength:   50,
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
	MatchedBases int // The number of matched bases.
	AlignedBases int // The number of aligned bases.

	AlignedFraction float64 // query (original query) coverage per HSP
	PIdentity       float64 // identity (fraction of same bases), percentage

	QueryLen int // length of the original query, used to compute/update AlignedFraction

	// QBegin int
	// QEnd   int
	// TBegin int
	// TEnd   int

	Chains *[]*Chain2Result
}

// Update updates the data with new chains.
func (r *SeqComparatorResult) Update(chains *[]*Chain2Result, queryLen int) {
	r.Chains = chains
	r.QueryLen = queryLen

	r.MatchedBases = 0
	r.AlignedBases = 0
	for _, c := range *r.Chains {
		r.MatchedBases += c.MatchedBases
		r.AlignedBases += c.AlignedBases
	}

	r.PIdentity = float64(r.MatchedBases) / float64(r.AlignedBases) * 100

	af := float64(r.MatchedBases) / float64(r.QueryLen) * 100
	if af > 100 {
		af = 100
	}
	r.AlignedFraction = af
}

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
				_sub2.QBegin = int32(iter.Index())
				_sub2.TBegin = int32(v)
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
	chains, nMatchedBases, nAlignedBases, _, _, _, _ := chainer.Chain(subs)
	defer func() {
		cpr.poolChainers.Put(chainer)
	}()
	if chains == nil {
		RecycleSubstrPairs(subs)
		return nil, nil
	}

	// var i int
	// var sub *SubstrPair
	// for c, chain := range *chains {
	// 	for _, i = range chain.Chain {
	// 		sub = (*subs)[i]
	// 		fmt.Printf("chain: %d, %s\n", c, sub)
	// 	}
	// }
	// fmt.Printf("%d, (%d/%d)\n", len(s), nMatchedBases, nAlignedBases)

	pIdent := float64(nMatchedBases) / float64(nAlignedBases) * 100
	if nAlignedBases < cpr.options.MinSegmentLength || pIdent < cpr.options.MinIdentity {
		RecycleChaining2Result(chains)
		RecycleSubstrPairs(subs)
		return nil, nil
	}

	af := float64(nAlignedBases) / float64(queryLen) * 100
	if af < cpr.options.MinAlignedFraction {
		RecycleChaining2Result(chains)
		RecycleSubstrPairs(subs)
		return nil, nil
	}

	// result object
	r := poolSeqComparatorResult.Get().(*SeqComparatorResult)
	r.AlignedBases = nAlignedBases
	r.MatchedBases = nMatchedBases
	r.QueryLen = queryLen
	r.AlignedFraction = af

	r.PIdentity = pIdent
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
