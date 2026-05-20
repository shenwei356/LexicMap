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
	"cmp"
	"fmt"
	"slices"
	"sync"
	"time"
	"unsafe"

	"github.com/dustin/go-humanize"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
)

// GSearchResultDetail is for storing genome search details
type GSearchResultDetail struct {
	BatchGenomeIndex []uint64 // multiple values belong to the genome chunks of the same genome
	Score            uint64   // score for sorting, total matched bases (masks * unique k-mers * length)
	Hits             []uint32 // count how many many k-mers are matched
}

// RecycleGSearchResultDetailsMap recycles a map of GSearchResultDetail
func (idx *Index) RecycleGSearchResultDetail(r *GSearchResultDetail) {
	r.BatchGenomeIndex = r.BatchGenomeIndex[:0]
	r.Score = 0
	clear(r.Hits)
	idx.poolGSearchResult.Put(r)
}

// RecycleGSearchResultDetailsMap recycles a map of GSearchResultDetail
func (idx *Index) RecycleGSearchResultDetailsMap(m *map[uint64]*GSearchResultDetail) {
	for _, r := range *m {
		r.BatchGenomeIndex = r.BatchGenomeIndex[:0]
		r.Score = 0
		clear(r.Hits)
		idx.poolGSearchResult.Put(r)
	}
	clear(*m)
	idx.poolGSearchResultsMap.Put(m)
}

// RecycleGSearchResultDetails recycles a list of GSearchResultDetail
func (idx *Index) RecycleGSearchResultDetails(rs *[]*GSearchResultDetail) {
	for _, r := range *rs {
		r.BatchGenomeIndex = r.BatchGenomeIndex[:0]
		r.Score = 0
		clear(r.Hits)
		idx.poolGSearchResult.Put(r)
	}
	*rs = (*rs)[:0]
	idx.poolGSearchResults.Put(rs)
}

// RecycleGSearchResult recycles the result of GSearch()
func (idx *Index) RecycleGSearchResult(whiteList *map[uint64]interface{}) {
	clear(*whiteList)
	poolUint64Map.Put(whiteList)
}

// GSearch searchs with a genome and return the list of possible genome internal ids.
func (idx *Index) GSearch(query *GQuery) (*map[uint64]interface{}, error) {
	if idx.opt.Debug {
		startTime0 := time.Now()
		log.Debugf("%s (%s bp): start to search genome", query.id, humanize.Comma(int64(query.genomeSize)))
		defer func() {
			log.Debugf("%s (%s bp): finished searching genome in %.3f seconds",
				query.id, humanize.Comma(int64(query.genomeSize)), time.Since(startTime0).Seconds())
		}()
	}

	// ------------------------------------------------------
	// 1. capture k-mers in overlapped windows

	_kmersW := idx.poolKmers.Get().(*[]*[]uint64)
	// _locsesW := idx.poolLocses.Get().(*[]*[]int)
	defer func() {
		var v *[]uint64
		for _, v = range *_kmersW {
			*v = (*v)[:0]
		}
		idx.poolKmers.Put(_kmersW)

		// var vl *[]int
		// for _, vl = range *_locsesW {
		// 	*vl = (*vl)[:0]
		// }
		// idx.poolLocses.Put(_locsesW)
	}()

	windows := idx.opt.Windows
	lenSeq := len(query.bigSeq)
	step := lenSeq / (windows + 1) // step size
	window := step << 1            // window size

	k := idx.k
	k8 := uint8(idx.lh.K)
	ccc := util.Ns(0b01, k8)
	ggg := util.Ns(0b10, k8)
	ttt := (uint64(1) << (k << 1)) - 1

	var start, end, j int
	var kmer uint64
	for i := 0; i < windows; i++ {
		start = i * step
		if i == windows-1 {
			end = lenSeq
		} else {
			end = start + window
		}
		// fmt.Printf("window #%d: %d-%d\n", i+1, start+1, end)

		_kmers, locses, err := idx.lh.MaskKnownDistinctPrefixes(query.bigSeq[start:end], query.skipRegions, true)
		if err != nil {
			panic(err)
		}

		for j, kmer = range *_kmers {
			if kmer == 0 || kmer == ccc || kmer == ggg || kmer == ttt ||
				util.IsLowComplexityDust(kmer, k8) {
				continue
			}

			*(*_kmersW)[j] = append(*(*_kmersW)[j], kmer)
		}

		if i == windows-1 { // sort k-mers and remove duplicates
			for j = range *_kmers {
				util.UniqUint64s((*_kmersW)[j])
			}
		}

		idx.lh.RecycleMaskResult(_kmers, locses)
	}

	// ------------------------------------------------------
	// 2. search k-mers and return the most similar genomes

	m := idx.poolGSearchResultsMap.Get().(*map[uint64]*GSearchResultDetail)

	inMemorySearch := idx.opt.InMemorySearch
	var searchers []*kv.Searcher
	var searchersIM []*kv.InMemorySearcher
	var nSearchers int
	if inMemorySearch {
		searchersIM = idx.InMemorySearchers
		nSearchers = len(searchersIM)
	} else {
		searchers = idx.Searchers
		nSearchers = len(searchers)
	}

	minPrefix := idx.opt.MinPrefix
	ch := make(chan *[]*kv.SearchResult, nSearchers)
	done := make(chan int) // later, we will reuse this
	var wg sync.WaitGroup
	var beginM, endM int // range of mask of a chunk

	// 2.2) collect search result
	go func() {
		var refpos uint64

		var sr *kv.SearchResult
		var refBatchAndIdxUint64 uint64
		var ok bool

		var filter *map[uint64]bool
		filterByTaxId := idx.filterByTaxId
		filterByPositiveTaxId := idx.filterByPositiveTaxId
		filterByNegativeTaxId := idx.filterByNegativeTaxId
		if filterByTaxId {
			filter = idx.poolTaxIDfilter.Get().(*map[uint64]bool)
		}
		var keepGenome, matchOne bool
		taxon := idx.Taxonomy
		genomeIdx2TaxId := idx.genomeIdx2TaxId
		keepGenomesWithoutTaxId := idx.opt.KeepGenomesWithoutTaxId
		var _taxid, taxid uint32
		taxids := idx.opt.TaxIds
		negativeTaxids := idx.opt.NegativeTaxIds

		for srs := range ch {
			// different k-mers in subjects,
			// most of cases, there are more than one
			for _, sr = range *srs {

				// multiple locations for each MATCHED k-mer
				// but most of cases, there's only one.
				for _, refpos = range sr.Values {
					refBatchAndIdxUint64 = refpos >> BITS_NONE_IDX // batch+refIdx

					// filter by taxid
					if filterByTaxId {
						if keepGenome, ok = (*filter)[refBatchAndIdxUint64]; ok {
							if !keepGenome {
								continue
							}
						} else {

							if taxid, ok = genomeIdx2TaxId[refBatchAndIdxUint64]; ok {
								// black list
								if filterByNegativeTaxId {
									matchOne = false
									for _, _taxid = range negativeTaxids {
										if taxon.LCA(taxid, _taxid) == _taxid {
											matchOne = true
											break
										}
									}
									if matchOne {
										(*filter)[refBatchAndIdxUint64] = false
										continue
									} else if !filterByPositiveTaxId {
										(*filter)[refBatchAndIdxUint64] = true
										continue
									}
								}

								// white list
								if filterByPositiveTaxId {
									matchOne = false
									for _, _taxid = range taxids {
										if taxon.LCA(taxid, _taxid) == _taxid {
											matchOne = true
											break
										}
									}
									if matchOne {
										(*filter)[refBatchAndIdxUint64] = true
									} else {
										(*filter)[refBatchAndIdxUint64] = false
										continue
									}
								}
							} else if !keepGenomesWithoutTaxId {
								(*filter)[refBatchAndIdxUint64] = false
								continue
							}
							(*filter)[refBatchAndIdxUint64] = true
						}
					}

					var r *GSearchResultDetail
					if r, ok = (*m)[refBatchAndIdxUint64]; !ok {
						r = idx.poolGSearchResult.Get().(*GSearchResultDetail)

						r.BatchGenomeIndex = append(r.BatchGenomeIndex, refBatchAndIdxUint64)

						(*m)[refBatchAndIdxUint64] = r
					}

					r.Hits[sr.IQuery]++       // add count to the probe
					r.Score += uint64(sr.Len) // matched length
				}
			}

			kv.RecycleSearchResults(srs)
		}

		if filterByTaxId {
			clear(*filter)
			idx.poolTaxIDfilter.Put(filter)
		}
		done <- 1
	}()

	// 2.1) search with multiple searchers
	// tokensS := make(chan int, idx.opt.MaxSeedingConcurrency)
	for iS := 0; iS < nSearchers; iS++ {
		if inMemorySearch {
			beginM = searchersIM[iS].ChunkIndex
			endM = searchersIM[iS].ChunkIndex + searchersIM[iS].ChunkSize
		} else {
			beginM = searchers[iS].ChunkIndex
			endM = searchers[iS].ChunkIndex + searchers[iS].ChunkSize
		}

		wg.Add(1)
		go func(iS, beginM, endM int) {
			var srs *[]*kv.SearchResult
			var err error
			if inMemorySearch {
				// prefix search
				srs, err = searchersIM[iS].Search2((*_kmersW)[beginM:endM], minPrefix, true, false)
				if err != nil {
					checkError(err)
				}
			} else {
				idx.searcherTokens[iS] <- 1 // get the access to the searcher

				// prefix search
				srs, err = searchers[iS].Search2((*_kmersW)[beginM:endM], minPrefix, true, false)
				if err != nil {
					checkError(err)
				}
			}
			if err != nil {
				checkError(err)
			}

			if len(*srs) == 0 { // no matcheds
				kv.RecycleSearchResults(srs)
			} else {
				ch <- srs // send result
			}

			if !inMemorySearch {
				<-idx.searcherTokens[iS] // return the access
			}

			wg.Done()
		}(iS, beginM, endM)
	}
	wg.Wait()
	close(ch)
	<-done

	if len(*m) == 0 { // no results
		idx.RecycleGSearchResultDetailsMap(m)
		return nil, nil
	}

	// collect and store with a list
	rs := idx.poolGSearchResults.Get().(*[]*GSearchResultDetail)
	for _, r := range *m {
		*rs = append(*rs, r)
	}
	clear(*m)
	idx.RecycleGSearchResultDetailsMap(m)

	// 2.3) handle chunked genomes

	// merge search result from genome chunks, if has split genome
	if idx.hasGenomeChunks {
		var r, rp *GSearchResultDetail
		var i, j, _i int
		var a uint64
		var v uint32
		var ok bool
		var li *[]int
		var ptr uintptr
		gcIdx2List := idx.poolGenomeChunksIdx2List.Get().(*map[uint64]*[]int)
		gcPtr2List := idx.poolGenomeChunksPointer2List.Get().(*map[uintptr]*[]int)

		// collect "i"s belonging to the same genome
		for i, r = range *rs {
			a = r.BatchGenomeIndex[0]
			if li, ok = (*gcIdx2List)[a]; !ok { // not a chunk of some genome
				continue
			}

			*li = append(*li, i) // is is just the index in *rs

			ptr = uintptr(unsafe.Pointer(li))
			if _, ok = (*gcPtr2List)[ptr]; !ok {
				(*gcPtr2List)[ptr] = li
			}
		}

		// merge alignments of the same genome
		for _, li = range *gcPtr2List {
			if len(*li) == 1 { // there's only one genome chunk
				// reset the list
				(*li) = (*li)[:0]
				continue
			}
			i = (*li)[0]
			rp = (*rs)[i]

			for _, j = range (*li)[1:] { // merge j -> i
				r = (*rs)[j]

				// merge r into rp
				rp.Score += r.Score
				for _i, v = range r.Hits {
					rp.Hits[_i] += v
				}

				idx.RecycleGSearchResultDetail(r)
				(*rs)[j] = nil
			}

			// reset the list
			(*li) = (*li)[:0]
		}

		// adjust the result list
		j = 0
		for _, r = range *rs {
			if r == nil {
				continue
			}

			(*rs)[j] = r
			j++
		}
		*rs = (*rs)[:j]

		// recycle datastructure
		clear(*gcPtr2List)
		idx.poolGenomeChunksPointer2List.Put(gcPtr2List)
		// clear(*gcIdx2List) # must not do this !!!
		idx.poolGenomeChunksIdx2List.Put(gcIdx2List)
	}

	// 2.4) sort
	topN := idx.opt.TopN
	slices.SortFunc(*rs, func(a, b *GSearchResultDetail) int {
		return cmp.Compare(b.Score, a.Score)
	})
	if topN > 0 && len(*rs) > topN {
		*rs = (*rs)[:topN]
	}

	whiteList := poolUint64Map.Get().(*map[uint64]interface{})
	clear(*whiteList)

	fmt.Printf("query\tsubject\tscore\thitMasks\thitKmers\thitKmerAvgLen\n")
	var refBatchAndIdxUint64 uint64
	for _, r := range *rs {
		for _, refBatchAndIdxUint64 = range r.BatchGenomeIndex {
			(*whiteList)[refBatchAndIdxUint64] = struct{}{}
		}

		hitKmers := 0
		hitMasks := 0
		for _, v := range r.Hits {
			if v > 0 {
				hitKmers += int(v)
				hitMasks++
			}
		}
		fmt.Printf("%s\t%s\t%d\t%d\t%d\t%.1f\n",
			query.id,
			idx.BatchGenomeIndex2GenomeID[r.BatchGenomeIndex[0]], r.Score,
			hitMasks,
			hitKmers, float64(r.Score)/float64(hitKmers),
		)
	}

	return whiteList, nil
}
