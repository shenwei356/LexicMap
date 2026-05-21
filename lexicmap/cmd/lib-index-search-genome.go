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
	"os"
	"runtime"
	"slices"
	"sync"
	"time"
	"unsafe"

	"github.com/dustin/go-humanize"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
	"github.com/vbauerster/mpb/v8"
	"github.com/vbauerster/mpb/v8/decor"
)

// GSearchResultDetail is for storing genome search details
type GSearchResultDetail struct {
	BatchGenomeIndex []uint64 // multiple values belong to the genome chunks of the same genome
	Score            uint64   // score for sorting, total matched bases (masks * unique k-mers * length)
	Hits             []uint8  // count how many many k-mers are matched for each mask
}

// RecycleGSearchResultDetailsMap recycles a map of GSearchResultDetail
func (idx *Index) RecycleGSearchDetailResult(r *GSearchResultDetail) {
	r.BatchGenomeIndex = r.BatchGenomeIndex[:0]
	r.Score = 0
	clear(r.Hits)
	idx.poolGSearchDetailResult.Put(r)
}

// RecycleGSearchDetailResultsMap recycles a map of GSearchResultDetail
func (idx *Index) RecycleGSearchDetailResultsMap(m *map[uint64]*GSearchResultDetail) {
	for _, r := range *m {
		r.BatchGenomeIndex = r.BatchGenomeIndex[:0]
		r.Score = 0
		clear(r.Hits)
		idx.poolGSearchDetailResult.Put(r)
	}
	clear(*m)
	idx.poolGSearchDetailResultsMap.Put(m)
}

// RecycleGSearchDetailResults recycles a list of GSearchResultDetail
func (idx *Index) RecycleGSearchDetailResults(rs *[]*GSearchResultDetail) {
	for _, r := range *rs {
		r.BatchGenomeIndex = r.BatchGenomeIndex[:0]
		r.Score = 0
		clear(r.Hits)
		idx.poolGSearchDetailResult.Put(r)
	}
	*rs = (*rs)[:0]
	idx.poolGSearchDetailResults.Put(rs)
}

// RecycleGSearchResult recycles the result of GSearch()
func (idx *Index) RecycleGSearchResult(whiteList *map[uint64]interface{}) {
	clear(*whiteList)
	poolUint64Map.Put(whiteList)
}

// GSearchScreen searchs with a genome and return the list of possible genome internal ids.
func (idx *Index) GSearchScreen(query *GQuery, windows int) (*map[uint64]interface{}, error) {
	if windows < 1 {
		return nil, fmt.Errorf("window size needs to be > 0")
	}

	if idx.opt.Debug {
		startTime0 := time.Now()
		log.Debugf("%s (%s bp): start to screen genomes", query.id, humanize.Comma(int64(query.genomeSize)))
		defer func() {
			log.Debugf("%s (%s bp): finished screening genomes in %.3f seconds",
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

	m := idx.poolGSearchDetailResultsMap.Get().(*map[uint64]*GSearchResultDetail)

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
						r = idx.poolGSearchDetailResult.Get().(*GSearchResultDetail)

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
		idx.RecycleGSearchDetailResultsMap(m)
		return nil, nil
	}

	// collect and store with a list
	rs := idx.poolGSearchDetailResults.Get().(*[]*GSearchResultDetail)
	for _, r := range *m {
		*rs = append(*rs, r)
	}
	clear(*m)
	idx.RecycleGSearchDetailResultsMap(m)

	// 2.3) handle chunked genomes

	// merge search result from genome chunks, if has split genome
	if idx.hasGenomeChunks {
		var r, rp *GSearchResultDetail
		var i, j, _i int
		var a uint64
		var v uint8
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
				rp.BatchGenomeIndex = append(rp.BatchGenomeIndex, r.BatchGenomeIndex...)
				rp.Score += r.Score
				for _i, v = range r.Hits {
					rp.Hits[_i] += v
				}

				idx.RecycleGSearchDetailResult(r)
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

	// fmt.Printf("query\tsubject\tscore\thitMasks\thitKmers\thitKmerAvgLen\n")
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
		// fmt.Printf("%s\t%s\t%d\t%d\t%d\t%.1f\n",
		// 	query.id,
		// 	idx.BatchGenomeIndex2GenomeID[r.BatchGenomeIndex[0]], r.Score,
		// 	hitMasks,
		// 	hitKmers, float64(r.Score)/float64(hitKmers),
		// )
	}

	// fmt.Println(whiteList)

	return whiteList, nil
}

// GSearchAlign align fragments of a query to candidates genomes.
func (idx *Index) GSearchAlign(
	query *GQuery, fragLen int, genomeIds *map[uint64]interface{},
	minAF float64,
	maxQueryConcurrency int, gcInterval uint64) error {

	if fragLen < 100 {
		return fmt.Errorf("fragment length is too small")
	}

	debug := idx.opt.Debug

	if debug {
		startTime0 := time.Now()
		log.Debugf("%s (%s bp): start to align fragments", query.id, humanize.Comma(int64(query.genomeSize)))
		defer func() {
			log.Debugf("%s (%s bp): finished aligning fragments in %.3f seconds",
				query.id, humanize.Comma(int64(query.genomeSize)), time.Since(startTime0).Seconds())
		}()
	}

	if maxQueryConcurrency == 0 {
		maxQueryConcurrency = runtime.NumCPU()
	}
	gc := gcInterval > 0
	if gcInterval > 0 {
		gcInterval = uint64(roundup32(uint32(gcInterval)))
		if gcInterval == 1 {
			gcInterval = 2
		}
	}
	gcIntervalMinus1 := gcInterval - 1

	// -----------------------------------------------------------
	// process bar
	var pbs *mpb.Progress
	var bar *mpb.Bar
	var chDuration chan time.Duration
	var doneDuration chan int
	if debug {
		// jobs
		var nJobs int
		var contig *[]byte
		for _, contig = range query.seqs {
			nJobs += len(*contig) / fragLen
			if len(*contig)%fragLen >= 100 {
				nJobs++
			}
		}

		pbs = mpb.New(mpb.WithWidth(40), mpb.WithOutput(os.Stderr))
		bar = pbs.AddBar(int64(nJobs),
			mpb.PrependDecorators(
				decor.Name("checked fragments: ", decor.WC{W: len("checked fragments: "), C: decor.DindentRight}),
				decor.Name("", decor.WCSyncSpaceR),
				decor.CountersNoUnit("%d / %d", decor.WCSyncWidth),
			),
			mpb.AppendDecorators(
				decor.Name("ETA: ", decor.WC{W: len("ETA: ")}),
				decor.EwmaETA(decor.ET_STYLE_GO, 1024),
				decor.OnComplete(decor.Name(""), ". done"),
			),
		)

		chDuration = make(chan time.Duration, idx.opt.NumCPUs)
		doneDuration = make(chan int)
		go func() {
			for t := range chDuration {
				bar.EwmaIncrBy(1, t)
			}
			doneDuration <- 1
		}()
	}

	// -----------------------------------------------------------

	query.result = poolGSearchResults.Get().(*[]*GSearchResult)

	var wg sync.WaitGroup
	tokens := make(chan int, maxQueryConcurrency)

	// 1) collect alignment results
	ch := make(chan *Query, maxQueryConcurrency)
	done := make(chan int)
	var total, matched uint64
	go func() {
		var r *SearchResult
		var sd *SimilarityDetail
		var cr *SeqComparatorResult
		var c *Chain2Result

		m := poolGSearchResultMap.Get().(*map[uint64]*GSearchResult)
		var ok bool
		var gr *GSearchResult

		// 1.1) collect alignment results

		for q := range ch { // each fragment
			total++

			if q.result == nil {
				poolQuery.Put(q)

				if gc && total&gcIntervalMinus1 == 0 {
					runtime.GC()
				}

				continue
			}

			matched++
			// ------------------------------------------------

			for _, r = range *q.result { // each subject genome
				if gr, ok = (*m)[r.BatchGenomeIndex]; !ok {
					gr = poolGSearchResult.Get().(*GSearchResult)
					gr.BatchGenomeIndex = r.BatchGenomeIndex
					gr.GenomeSize = r.GenomeSize
					gr.NumSeqs = r.NumSeqs

					(*m)[r.BatchGenomeIndex] = gr
				}

				for _, sd = range *r.SimilarityDetails { // each chain
					cr = sd.Similarity

					for _, c = range *cr.Chains { // each match
						if c == nil {
							continue
						}

						gr.AlignedFragments++
						gr.AlignedLength += c.AlignedLength
						gr.AlignedMatches += c.MatchedBases
					}
				}
			}

			// ------------------------------------------------

			if gc && total&gcIntervalMinus1 == 0 {
				runtime.GC()
			}
		}

		// ---------------------------------------
		// 1.2 from map to slice
		rs := poolGSearchResults.Get().(*[]*GSearchResult)
		for _, gr = range *m {

			*rs = append(*rs, gr)
		}
		clear(*m)
		poolGSearchResultMap.Put(m)

		// ---------------------------------------
		// 1.3) merge results of genome chunks
		if idx.hasGenomeChunks {
			var r, rp *GSearchResult
			var i, j int
			var a uint64
			var ok bool
			var li *[]int
			var ptr uintptr
			gcIdx2List := idx.poolGenomeChunksIdx2List.Get().(*map[uint64]*[]int)
			gcPtr2List := idx.poolGenomeChunksPointer2List.Get().(*map[uintptr]*[]int)

			// collect "i"s belonging to the same genome
			for i, r = range *rs {
				a = r.BatchGenomeIndex
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
					rp.GenomeSize += r.GenomeSize
					rp.NumSeqs += r.NumSeqs
					rp.AlignedFragments += r.AlignedFragments
					rp.AlignedLength += r.AlignedLength
					rp.AlignedMatches += r.AlignedMatches

					poolGSearchResult.Put(r)
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

		// ---------------------------------------
		// 1.4) compute ani and af
		j := 0
		for _, gr := range *rs {
			gr.ANI = float64(gr.AlignedMatches) / float64(gr.AlignedLength)
			gr.AF = float64(gr.AlignedLength) / float64(query.genomeSize)

			if gr.AF < minAF {
				poolGSearchResult.Put(gr)
				continue
			}

			if gr.AF > 1 {
				gr.AF = 1
			}

			gr.Score = gr.ANI * gr.AF

			(*rs)[j] = gr
			j++

		}
		*rs = (*rs)[:j]

		slices.SortFunc(*rs, func(a, b *GSearchResult) int {
			return cmp.Compare(b.Score, a.Score)
		})

		// ---------------------------------------
		query.result = rs
		done <- 1
	}()

	// 2) alignment fragments
	var i, j, end0, s, e int
	var contig *[]byte
	buf8 := make([]byte, 8)
	fcpus := float64(idx.opt.NumCPUs)
	for i, contig = range query.seqs {
		end0 = len(*contig)
		for j = 0; j < end0; j += fragLen {
			s = j
			e = j + fragLen
			if e > end0 {
				e = end0
				if e-s < 100 { // skip fragments < 100 bp
					continue
				}
			}

			q := poolQuery2.Get().(*Query)
			q.Reset()

			// seq id
			be.PutUint32(buf8[:4], uint32(i))
			be.PutUint32(buf8[4:8], uint32(j))
			q.seqID = append(q.seqID, buf8...)

			// seq
			q.seq = append(q.seq, (*contig)[s:e]...)

			tokens <- 1
			wg.Add(1)
			go func(i, j int, q *Query) {
				timeStart := time.Now()
				defer func() {
					<-tokens
					wg.Done()
					if debug {
						chDuration <- time.Duration(float64(time.Since(timeStart)) / fcpus)
					}
				}()

				var err error
				q.result, err = idx.Search(q, genomeIds, false)
				if err != nil {
					checkError(fmt.Errorf("search contig #%d, fragment %d-%d: %s",
						i+1, j*fragLen+1, (j+1)*fragLen, err))
				}

				ch <- q
			}(i, j, q)
		}
	}

	wg.Wait()
	close(ch)
	<-done

	// process bar
	if debug {
		close(chDuration)
		<-doneDuration
		pbs.Wait()
	}

	return nil
}

var poolQuery2 = &sync.Pool{New: func() interface{} {
	return &Query{
		// 4 bytes for contig index
		// 4 bytes for start position
		seqID: make([]byte, 0, 8),
		seq:   make([]byte, 0, 1<<10), // 2k id enough for the common 1020-bp fragments
	}
}}

// -------------------------------------------------------------------------

// GSearchResult represents a subject genome hit.
type GSearchResult struct {
	BatchGenomeIndex uint64
	GenomeSize       int
	NumSeqs          int

	AlignedFragments int
	AlignedLength    int
	AlignedMatches   int

	ANI   float64
	AF    float64
	Score float64 // for sorting
}

var poolGSearchResult = &sync.Pool{New: func() interface{} {
	return &GSearchResult{}
}}

var poolGSearchResults = &sync.Pool{New: func() interface{} {
	tmp := make([]*GSearchResult, 0, 1024)
	return &tmp
}}

var poolGSearchResultMap = &sync.Pool{New: func() interface{} {
	tmp := make(map[uint64]*GSearchResult, 1024)
	return &tmp
}}

func RecycleGSearchResults(rs *[]*GSearchResult) {
	for _, r := range *rs {
		poolGSearchResult.Put(r)
	}
	*rs = (*rs)[:0]
	poolGSearchResults.Put(*rs)
}

func RecycleGSearchResultMap(rs *map[uint64]*GSearchResult) {
	for _, r := range *rs {
		poolGSearchResult.Put(r)
	}
	clear(*rs)
	poolGSearchResultMap.Put(*rs)
}
