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

import "sync"

// seedSearchResult holds only state needed before alignment. Each collector
// owns fixed pages so growing the candidate list never moves an active record.
type seedSearchResult struct {
	BatchGenomeIndex uint64
	Subs             []SubstrPair
	Score            float32
	chainRegions     []seedChainRegion
	inlineRegion     [1]seedChainRegion
	inlineSub        [1]SubstrPair
}

const seedSearchResultPageSize = 256

type seedSearchResultPage [seedSearchResultPageSize]seedSearchResult

var poolSeedSearchResultPages = sync.Pool{New: func() any {
	return &seedSearchResultPage{}
}}

var poolSeedSearchResults = sync.Pool{New: func() any {
	results := make([]*seedSearchResult, 0, 16)
	return &results
}}

type seedSearchResultArena struct {
	pages []*seedSearchResultPage
	used  int
}

func (a *seedSearchResultArena) add(id uint64) *seedSearchResult {
	i := a.used % seedSearchResultPageSize
	if i == 0 {
		a.pages = append(a.pages, poolSeedSearchResultPages.Get().(*seedSearchResultPage))
	}
	r := &a.pages[len(a.pages)-1][i]
	r.BatchGenomeIndex = id
	r.Score = 0
	r.resetSeeds()
	r.resetRegions()
	a.used++
	return r
}

func (r *seedSearchResult) resetSeeds() {
	if r.Subs == nil || cap(r.Subs) > thresholdNSubs {
		r.Subs = r.inlineSub[:0]
	} else {
		r.Subs = r.Subs[:0]
	}
}

func (r *seedSearchResult) resetRegions() {
	if r.chainRegions == nil || cap(r.chainRegions) > thresholdNSubs {
		r.chainRegions = r.inlineRegion[:0]
	} else {
		r.chainRegions = r.chainRegions[:0]
	}
}

// alignmentResult copies the bounds into an independently owned result. No
// references to the candidate pages escape to alignment output or its pools.
func (r *seedSearchResult) alignmentResult() *SearchResult {
	result := poolSearchResult.Get().(*SearchResult)
	result.Reset()
	result.BatchGenomeIndex = r.BatchGenomeIndex
	result.GenomeBatch = int(r.BatchGenomeIndex >> BITS_GENOME_IDX)
	result.GenomeIndex = int(r.BatchGenomeIndex & MASK_GENOME_IDX)
	result.Score = r.Score
	result.chainRegions = append(result.chainRegions, r.chainRegions...)
	return result
}

func (a *seedSearchResultArena) recycle() {
	for _, page := range a.pages {
		for i := range page {
			r := &page[i]
			r.BatchGenomeIndex = 0
			r.Score = 0
			r.resetSeeds()
			r.resetRegions()
		}
		poolSeedSearchResultPages.Put(page)
	}
	clear(a.pages)
	a.pages = nil
	a.used = 0
}

func recycleSeedSearchResultSlice(results *[]*seedSearchResult) {
	clear((*results)[:cap(*results)])
	if cap(*results) > maxPooledResults {
		*results = nil
		return
	}
	*results = (*results)[:0]
	poolSeedSearchResults.Put(results)
}
