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

package index

import (
	"sync"

	rtree "github.com/shenwei356/LexicMap/lexicmap/index/recyclable-tree"
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
	MinAlignedFraction float64 // percentage
	MinIdentity        float64 // percentage
}

// DefaultSeqComparatorOptions contains the default options for SeqComparatorOptions.
var DefaultSeqComparatorOptions = SeqComparatorOptions{
	K:         32,
	MinPrefix: 11, // can not be too small, or there will be a large number of anchors.

	Chaining2Options: Chaining2Options{
		// should be relative small
		MaxGap: 32,
		// better be larger than MinPrefix
		MinScore: 20,
		// can not be < k
		MaxDistance: 50,
		// can not be two small
		Band: 20,
	},

	MinAlignedFraction: 70,
	MinIdentity:        70,
}

// SeqComparator is for fast and accurate similarity estimation of two sequences,
// which are in the same strand (important).
type SeqComparator struct {
	// options
	options *SeqComparatorOptions
	// chainer for chaining anchors,
	// shared variable-length substrings searched by prefix matching.
	chainer *Chainer2

	// a prefix tree for matching k-mers
	tree *rtree.Tree
	len  int
}

// NewSeqComparator creates a new SeqComparator with given options.
// No options checking now.
func NewSeqComparator(options *SeqComparatorOptions) *SeqComparator {
	cpr := &SeqComparator{
		options: options,
		chainer: NewChainer2(&options.Chaining2Options),
	}
	return cpr
}

// Index initializes the SeqComparator with a sequence.
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

	for {
		kmer, ok, _ = iter.NextPositiveKmer()
		if !ok {
			break
		}

		t.Insert(kmer, uint32(iter.Index()))
	}

	cpr.tree = t
	cpr.len = len(s)

	return nil
}

// SeqComparatorResult contains the details of a seq comparison result.
type SeqComparatorResult struct {
	MatchedBases int // The number of matched bases.
	AlignedBases int // The number of aligned bases.
	NumChains    int // The number of chains

	AlignedFraction float64 // aligned fraction, percentage
	Identity        float64 // identity (fraction of same bases), percentage
}

var poolSeqComparatorResult = &sync.Pool{New: func() interface{} {
	return &SeqComparatorResult{}
}}

// RecycleSeqComparatorResult recycles a SeqComparatorResult
func RecycleSeqComparatorResult(r *SeqComparatorResult) {
	poolSeqComparatorResult.Put(r)
}

// Compare matchs k-mers for the query sequence, chains them up,
// and computes the similarity.
// Please remember to call RecycleSeqComparatorResult() to recycle the result.
func (cpr *SeqComparator) Compare(s []byte) (*SeqComparatorResult, error) {
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
			for _, v = range sr.Values {
				_sub2 := poolSub.Get().(*SubstrPair)
				_sub2.QBegin = iter.Index()
				_sub2.TBegin = int(v)
				_sub2.Code = rtree.KmerPrefix(sr.Kmer, k8, sr.LenPrefix)
				_sub2.Len = int(sr.LenPrefix)
				_sub2.RC = false

				*subs = append(*subs, _sub2)
			}
		}
		t.RecycleSearchResult(srs)
	}

	if len(*subs) < 1 { // no way, only one match?
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

	chains, nMatchedBases, nAlignedBases := cpr.chainer.Chain(subs)
	if len(*chains) == 0 {
		RecycleChainingResult(chains)
		RecycleSubstrPairs(subs)
		return nil, nil
	}

	// var i int
	// var sub *SubstrPair
	// for c, chain := range *chains {
	// 	for _, i = range *chain {
	// 		sub = (*subs)[i]
	// 		fmt.Printf("chain: %d, %s\n", c, sub)
	// 	}
	// }
	// fmt.Printf("%d, (%d/%d)\n", len(s), nMatchedBases, nAlignedBases)

	// result object
	r := poolSeqComparatorResult.Get().(*SeqComparatorResult)
	r.AlignedBases = nAlignedBases
	r.MatchedBases = nMatchedBases
	r.NumChains = len(*chains)
	r.Identity = float64(nMatchedBases) / float64(nAlignedBases) * 100
	r.AlignedFraction = float64(nAlignedBases) / float64(cpr.len) * 100
	if r.AlignedFraction > 100 {
		r.AlignedFraction = 100
	}

	RecycleChainingResult(chains)
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
