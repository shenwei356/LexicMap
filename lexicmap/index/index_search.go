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
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OFTestSerializationTestSerialization ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

package index

import (
	"fmt"
	"sort"
	"sync"

	"github.com/shenwei356/LexicMap/lexicmap/index/twobit"
	"github.com/shenwei356/LexicMap/lexicmap/tree"
	"github.com/shenwei356/LexicMap/lexicmap/util"
)

// SearchOptions defineds options used in searching.
type SearchOptions struct {
	// basic
	MinPrefix       uint8 // minimum prefix length, e.g, 15
	MinSinglePrefix uint8 // minimum prefix length of the single seed, e.g, 20
	TopN            int   // keep the topN scores

	// chaining
	MaxGap float64
}

// DefaultSearchOptions contains default option values.
var DefaultSearchOptions = SearchOptions{
	MinPrefix:       15,
	MinSinglePrefix: 20,
	TopN:            10,

	MaxGap: 5000,
}

// SetChainingOptions replaces the default chaining option with a new one.
func (idx *Index) SetChainingOptions(co *ChainingOptions) {
	idx.chainingOptions = co
	idx.poolChainers = &sync.Pool{New: func() interface{} {
		return NewChainer(co)
	}}
}

// SetSearchingOptions sets the searching options.
// Note that it overwrites the result of SetChainingOption.
func (idx *Index) SetSearchingOptions(so *SearchOptions) {
	idx.searchOptions = so

	// chaining
	co := &ChainingOptions{
		MaxGap:   so.MaxGap,
		MinScore: seedWeight(float64(so.MinSinglePrefix)),
	}
	idx.chainingOptions = co

	idx.poolChainers = &sync.Pool{New: func() interface{} {
		return NewChainer(co)
	}}
}

// SetAlignOptions sets the alignment options
// func (idx *Index) SetAlignOptions(ao *align.AlignOptions) {
// 	idx.alignOptions = ao
// 	idx.poolAligner = &sync.Pool{New: func() interface{} {
// 		return align.NewAligner(ao)
// 	}}
// }

// SetSeqCompareOptions sets the sequence comparing options
func (idx *Index) SetSeqCompareOptions(sco *SeqComparatorOptions) {
	idx.seqCompareOption = sco
	idx.poolSeqComparator = &sync.Pool{New: func() interface{} {
		return NewSeqComparator(sco)
	}}
}

// --------------------------------------------------------------------------

// SubstrPair represents a pair of found substrings/seeds, it's also called an anchor.
type SubstrPair struct {
	QBegin int    // start position of the substring (0-based) in query
	TBegin int    // start position of the substring (0-based) in reference
	Len    int    // length
	Code   uint64 // k-mer, only for debugging

	RC bool // is the substring from the reference seq on the negative strand.
}

func (s SubstrPair) String() string {
	return fmt.Sprintf("%3d-%3d vs %3d-%3d len:%2d, rc:%v",
		s.QBegin+1, s.QBegin+s.Len, s.TBegin+1, s.TBegin+s.Len, s.Len, s.RC)
}

var poolSub = &sync.Pool{New: func() interface{} {
	return &SubstrPair{}
}}

var poolSubs = &sync.Pool{New: func() interface{} {
	tmp := make([]*SubstrPair, 0, 1024)
	return &tmp
}}

// RecycleSubstrPairs recycles a list of SubstrPairs
func RecycleSubstrPairs(subs *[]*SubstrPair) {
	for _, sub := range *subs {
		poolSub.Put(sub)
	}
	poolSubs.Put(subs)
}

// ClearSubstrPairs removes nested/embedded and same anchors. k is the largest k-mer size.
func ClearSubstrPairs(subs *[]*SubstrPair, k int) {
	if len(*subs) < 2 {
		return
	}

	// sort substrings/seeds in ascending order based on the starting position
	// and in descending order based on the ending position.
	sort.Slice(*subs, func(i, j int) bool {
		a := (*subs)[i]
		b := (*subs)[j]
		if a.QBegin == b.QBegin {
			return a.QBegin+a.Len >= b.QBegin+b.Len
		}
		return a.QBegin < b.QBegin
	})

	var p *SubstrPair
	var upbound, vQEnd, vTEnd int
	var j int
	markers := poolBoolList.Get().(*[]bool)
	*markers = (*markers)[:0]
	for range *subs {
		*markers = append(*markers, false)
	}
	for i, v := range (*subs)[1:] {
		vQEnd = v.QBegin + v.Len
		upbound = vQEnd - k
		vTEnd = v.TBegin + v.Len
		j = i
		for j >= 0 { // have to check previous N seeds
			p = (*subs)[j]
			if p.QBegin < upbound { // no need to check
				break
			}

			// same or nested region
			if vQEnd <= p.QBegin+p.Len &&
				v.TBegin >= p.TBegin && vTEnd <= p.TBegin+p.Len {
				poolSub.Put(v)         // do not forget to recycle the object
				(*markers)[i+1] = true // because of: range (*subs)[1:]
				break
			}

			j--
		}
	}

	j = 0
	for i, embedded := range *markers {
		if !embedded {
			(*subs)[j] = (*subs)[i]
			j++
		}
	}
	if j > 0 {
		*subs = (*subs)[:j]
	}

	poolBoolList.Put(markers)
}

var poolBoolList = &sync.Pool{New: func() interface{} {
	m := make([]bool, 0, 1024)
	return &m
}}

// var poolAlignResults = &sync.Pool{New: func() interface{} {
// 	tmp := make([]*align.AlignResult, 0, 8)
// 	return &tmp
// }}

var poolSimilarityDetail = &sync.Pool{New: func() interface{} {
	return &SimilarityDetail{}
}}

var poolSimilarityDetails = &sync.Pool{New: func() interface{} {
	tmp := make([]*SimilarityDetail, 0, 8)
	return &tmp
}}

var poolSearchResult = &sync.Pool{New: func() interface{} {
	return &SearchResult{}
}}

var poolSearchResults = &sync.Pool{New: func() interface{} {
	tmp := make([]*SearchResult, 0, 16)
	return &tmp
}}

// SearchResult stores a search result for the given query sequence.
type SearchResult struct {
	IdIdx int            // index of the matched reference ID
	Subs  *[]*SubstrPair // matched substring pairs (query,target)

	Score float64 //  score for soring

	Chains *[]*[]int

	// more about the alignment detail
	// AlignResults *[]*align.AlignResult
	SimilarityDetails *[]*SimilarityDetail // sequence comparing
}

type SimilarityDetail struct {
	TBegin int
	TEnd   int
	RC     bool

	SimilarityScore float64
	Similarity      *SeqComparatorResult
	Chain           *[]int
}

func (r *SearchResult) Reset() {
	r.IdIdx = -1
	r.Subs = nil
	r.Score = 0
	r.Chains = nil
	// r.AlignResults = nil
	r.SimilarityDetails = nil
}

// RecycleSearchResults recycles a search result object
func (idx *Index) RecycleSearchResult(r *SearchResult) {
	for _, sub := range *r.Subs {
		poolSub.Put(sub)
	}
	poolSubs.Put(r.Subs)

	if r.Chains != nil {
		for _, chain := range *r.Chains {
			poolChain.Put(chain)
		}
		poolChains.Put(r.Chains)
	}

	// yes, it might be nil for some failed in chaining
	// if r.AlignResults != nil {
	// 	for _, ar := range *r.AlignResults {
	// 		align.RecycleAlignResult(ar)
	// 	}
	// 	poolAlignResults.Put(r.AlignResults)
	// }

	if r.SimilarityDetails != nil {
		for _, sd := range *r.SimilarityDetails {
			RecycleSeqComparatorResult(sd.Similarity)
			poolSimilarityDetail.Put(sd)
		}
		poolSimilarityDetails.Put(r.SimilarityDetails)
	}

	poolSearchResult.Put(r)
}

// RecycleSearchResults recycles search results objects
func (idx *Index) RecycleSearchResults(sr *[]*SearchResult) {
	if sr == nil {
		return
	}

	for _, r := range *sr {
		idx.RecycleSearchResult(r)
	}
	poolSearchResults.Put(sr)
}

var poolSearchResultsMap = &sync.Pool{New: func() interface{} {
	m := make(map[int]*SearchResult, 1024)
	return &m
}}

// Search queries the index with a sequence.
// After using the result, do not forget to call RecycleSearchResult().
func (idx *Index) Search(s []byte) (*[]*SearchResult, error) {
	// ----------------------------------------------------------------
	// mask the query sequence

	_kmers, _locses, err := idx.lh.Mask(s, nil)
	if err != nil {
		return nil, err
	}
	defer idx.lh.RecycleMaskResult(_kmers, _locses)

	// ----------------------------------------------------------------
	// matching the captured k-mers in databases

	// a map for collecting matches for each reference: IdIdx -> result
	m := poolSearchResultsMap.Get().(*map[int]*SearchResult)
	clear(*m) // requires go >= v1.21

	var refpos uint64
	var i int
	var kmer uint64

	// query substring
	var _pos int
	var _begin int
	var _rc bool

	var code uint64
	var K, _k int
	var idIdx, pos, begin int
	var rc bool
	trees := idx.Trees
	K = idx.K()
	K8 := uint8(K)
	var locs []int
	var srs *[]*tree.SearchResult
	var sr *tree.SearchResult
	var ok bool
	minPrefix := idx.searchOptions.MinPrefix
	for i, kmer = range *_kmers { // captured k-mers by the maskes
		srs, ok = trees[i].Search(kmer, minPrefix) // search it on the corresponding tree
		if !ok {                                   // no matcheds
			continue
		}

		locs = (*_locses)[i] // locations in the query

		// multiple locations for each QUERY k-mer,
		// but most of cases, there's only one.
		for _, _pos = range locs {
			_rc = _pos&1 > 0 // if on the reverse complement sequence
			_pos >>= 2

			// different k-mers in subjects,
			// most of cases, there are more than one
			for _, sr = range *srs {
				// matched length
				_k = int(sr.LenPrefix)

				// query
				if _rc { // on the negative strand
					_begin = _pos + K - _k
				} else {
					_begin = _pos
				}

				// matched
				code = util.KmerPrefix(sr.Kmer, K8, sr.LenPrefix)

				// multiple locations for each MATCHED k-mer
				// but most of cases, there's only one.
				for _, refpos = range sr.Values {
					idIdx = int(refpos >> 38)
					pos = int(refpos << 26 >> 28)
					rc = refpos&1 > 0

					if rc {
						begin = pos + K - _k
					} else {
						begin = pos
					}

					_sub2 := poolSub.Get().(*SubstrPair)
					_sub2.QBegin = _begin
					_sub2.TBegin = begin
					_sub2.Code = code
					_sub2.Len = _k
					_sub2.RC = rc

					var r *SearchResult
					if r, ok = (*m)[idIdx]; !ok {
						subs := poolSubs.Get().(*[]*SubstrPair)
						*subs = (*subs)[:0]

						r = poolSearchResult.Get().(*SearchResult)
						r.IdIdx = idIdx
						r.Subs = subs
						r.Score = 0
						r.Chains = nil // important
						// r.AlignResults = nil // important
						r.SimilarityDetails = nil // important

						(*m)[idIdx] = r
					}

					*r.Subs = append(*r.Subs, _sub2)
				}
			}
		}

		trees[i].RecycleSearchResult(srs)
	}

	if len(*m) == 0 { // no results
		poolSearchResultsMap.Put(m)
		return nil, nil
	}

	// ----------------------------------------------------------------
	// chaining matches for all subject sequences

	minChainingScore := idx.chainingOptions.MinScore

	rs := poolSearchResults.Get().(*[]*SearchResult)
	*rs = (*rs)[:0]

	minSinglePrefix := int(idx.searchOptions.MinSinglePrefix)
	for _, r := range *m {
		ClearSubstrPairs(r.Subs, K) // remove duplicates and nested anchors

		// there's no need to chain for a single short seed
		if len(*r.Subs) == 1 && (*r.Subs)[0].Len < minSinglePrefix {
			// do not forget to recycle filtered result
			idx.RecycleSearchResult(r)
			continue
		}

		for _, sub := range *r.Subs {
			r.Score += float64(sub.Len * sub.Len)
		}

		*rs = append(*rs, r)
	}

	// sort subjects in descending order based on the score (simple statistics).
	// just use the standard library for a few seed pairs.
	sort.Slice(*rs, func(i, j int) bool {
		return (*rs)[i].Score > (*rs)[j].Score
	})

	poolSearchResultsMap.Put(m)

	// only keep the top N targets
	topN := idx.searchOptions.TopN
	if topN > 0 && len(*rs) > topN {
		var r *SearchResult
		for i := topN; i < len(*rs); i++ {
			r = (*rs)[i]

			// do not forget to recycle the filtered result
			idx.RecycleSearchResult(r)
		}
		*rs = (*rs)[:topN]
	}

	// chaining
	chainer := idx.poolChainers.Get().(*Chainer)
	j := 0
	for _, r := range *rs {
		r.Chains, r.Score = chainer.Chain(r.Subs)
		if r.Score < minChainingScore {
			idx.RecycleSearchResult(r) // do not forget to recycle unused objects
			continue
		} else {
			(*rs)[j] = r
			j++
		}
	}
	*rs = (*rs)[:j]
	idx.poolChainers.Put(chainer)

	if !idx.saveTwoBit {
		return rs, nil
	}

	// ----------------------------------------------------------------
	// sequence similarity

	var sub *SubstrPair
	qlen := len(s)
	var chain *[]int
	var qs, qe, ts, te, tBegin, tEnd int
	var tSeq *[]byte

	rdr := <-idx.twobitReaders

	// aligner := idx.poolAligner.Get().(*align.Aligner)
	cpr := idx.poolSeqComparator.Get().(*SeqComparator)
	err = cpr.Index(s) // index the query sequence
	if err != nil {
		return nil, err
	}

	minAF := idx.seqCompareOption.MinAlignedFraction
	minIdent := idx.seqCompareOption.MinIdentity

	// check all references
	for _, r := range *rs {
		// ars := poolAlignResults.Get().(*[]*align.AlignResult)
		// *ars = (*ars)[:0]
		sds := poolSimilarityDetails.Get().(*[]*SimilarityDetail)
		*sds = (*sds)[:0]

		// check sequences from all chains
		for i, chain = range *r.Chains {
			// ------------------------------------------------------------------------
			// extract subsequence from the refseq for comparing

			// the first seed pair
			sub = (*r.Subs)[(*chain)[0]]
			qs = sub.QBegin
			ts = sub.TBegin

			// the last seed pair
			sub = (*r.Subs)[(*chain)[len(*chain)-1]]
			qe = sub.QBegin + sub.Len
			te = sub.TBegin + sub.Len

			if len(*r.Subs) == 1 { // if there's only one seed, need to check the strand information
				rc = sub.RC
			} else { // check the strand according to coordinates of seeds
				rc = ts > sub.TBegin
			}

			// estimate the location of target sequence on the reference
			if rc { // reverse complement
				tBegin = sub.TBegin - qlen + qe
				if tBegin < 0 {
					tBegin = 0
				}
				tEnd = ts + sub.Len + qs - 1
			} else {
				tBegin = ts - qs
				if tBegin < 0 {
					tBegin = 0
				}
				tEnd = te + qlen - qe - 1
			}

			// fmt.Printf("subject:%s:%d-%d, rc:%v\n", idx.IDs[r.IdIdx], tBegin+1, tEnd+1, rc)

			// extract target sequence for comparison.
			// Right now, we fetch seq from disk for each seq,
			// Later, we'll buffer frequently accessed references for improving speed.
			tSeq, err = rdr.SubSeq(r.IdIdx, tBegin, tEnd)
			if err != nil {
				return rs, err
			}

			if rc { // reverse complement
				RC(*tSeq)
			}

			// ------------------------------------------------------------------------
			// comparing the two sequences

			// fast filter with sketching comparison

			// costly (pseudo-)alignment

			// ar := aligner.Global(s, *tSeq)
			// *ars = append(*ars, ar)

			cr, err := cpr.Compare(*tSeq)
			if err != nil {
				return nil, err
			}
			if cr == nil {
				twobit.RecycleTwoBit(tSeq)
				continue
			}
			if cr.AlignedFraction < minAF || cr.Identity < minIdent {
				RecycleSeqComparatorResult(cr)
				twobit.RecycleTwoBit(tSeq)
				continue
			}

			sd := poolSimilarityDetail.Get().(*SimilarityDetail)
			sd.TBegin = tBegin
			sd.TEnd = tEnd
			sd.RC = rc
			sd.Chain = (*r.Chains)[i]
			sd.Similarity = cr
			sd.SimilarityScore = cr.AlignedFraction * cr.Identity

			*sds = append(*sds, sd)

			twobit.RecycleSeq(tSeq)
		}
		// r.AlignResults = ars
		sort.Slice(*sds, func(i, j int) bool {
			return (*sds)[i].SimilarityScore > (*sds)[j].SimilarityScore
		})
		r.SimilarityDetails = sds
	}

	// idx.poolAligner.Put(aligner)
	idx.twobitReaders <- rdr

	// recycle the tree data for this query
	cpr.RecycleIndex()
	// recycle this comparator
	idx.poolSeqComparator.Put(cpr)

	return rs, nil
}

// RC computes the reverse complement sequence
func RC(s []byte) []byte {
	n := len(s)
	for i := 0; i < n; i++ {
		s[i] = rcTable[s[i]]
	}
	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		s[i], s[j] = s[j], s[i]
	}
	return s
}

var rcTable = [256]byte{
	0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
	16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
	32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
	48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
	64, 84, 86, 71, 72, 69, 70, 67, 68, 73, 74, 77, 76, 75, 78, 79,
	80, 81, 89, 83, 65, 85, 66, 87, 88, 82, 90, 91, 92, 93, 94, 95,
	96, 116, 118, 103, 104, 101, 102, 99, 100, 105, 106, 109, 108, 107, 110, 111,
	112, 113, 121, 115, 97, 117, 98, 119, 120, 114, 122, 123, 124, 125, 126, 127,
	128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
	144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
	160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
	176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
	192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
	208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
	224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
	240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255,
}
