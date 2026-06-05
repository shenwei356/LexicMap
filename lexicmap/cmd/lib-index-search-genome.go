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
	"math"
	"os"
	"slices"
	"sync"
	"time"
	"unsafe"

	"github.com/dustin/go-humanize"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
	"github.com/shenwei356/wfa"
	"github.com/twotwotwo/sorts/sortutil"
	"github.com/vbauerster/mpb/v8"
	"github.com/vbauerster/mpb/v8/decor"
)

// GSearchScreenResultDetail is for storing genome search details
type GSearchScreenResultDetail struct {
	BatchGenomeIndex []uint64 // multiple values belong to the genome chunks of the same genome
	Score            uint64   // score for sorting, total matched bases (masks * unique k-mers * length)
	// Hits             []uint8  // count how many k-mers are matched for each mask
}

// RecycleGSearchResultDetailsMap recycles a map of GSearchResultDetail
func (idx *Index) RecycleGSearchScreenDetailResult(r *GSearchScreenResultDetail) {
	r.BatchGenomeIndex = r.BatchGenomeIndex[:0]
	r.Score = 0
	// clear(r.Hits)
	idx.poolGSearchDetailResult.Put(r)
}

// RecycleGSearchScreenDetailResultsMap recycles a map of GSearchResultDetail
func (idx *Index) RecycleGSearchScreenDetailResultsMap(m *map[uint64]*GSearchScreenResultDetail) {
	for _, r := range *m {
		r.BatchGenomeIndex = r.BatchGenomeIndex[:0]
		r.Score = 0
		// clear(r.Hits)
		idx.poolGSearchDetailResult.Put(r)
	}
	clear(*m)
	idx.poolGSearchDetailResultsMap.Put(m)
}

// RecycleGSearchScreenDetailResults recycles a list of GSearchResultDetail
func (idx *Index) RecycleGSearchScreenDetailResults(rs *[]*GSearchScreenResultDetail) {
	for _, r := range *rs {
		r.BatchGenomeIndex = r.BatchGenomeIndex[:0]
		r.Score = 0
		// clear(r.Hits)
		idx.poolGSearchDetailResult.Put(r)
	}
	*rs = (*rs)[:0]
	idx.poolGSearchDetailResults.Put(rs)
}

var poolUint64ToUint64SliceMap = &sync.Pool{New: func() interface{} {
	tmp := make(map[uint64]*[]uint64, 128)
	return &tmp
}}

// RecycleGSearchScreenResult recycles the result of GSearch()
func (idx *Index) RecycleGSearchScreenResult(whiteList *map[uint64]*[]uint64) {
	clear(*whiteList)
	poolUint64ToUint64SliceMap.Put(whiteList)
}

// GSearchScreen searchs with a genome and return the list of possible genome internal ids.
func (idx *Index) GSearchScreen(query *GQuery, windows int) (*map[uint64]*[]uint64, error) {
	if windows < 1 {
		return nil, fmt.Errorf("window size needs to be > 0")
	}

	whiteList := poolUint64ToUint64SliceMap.Get().(*map[uint64]*[]uint64)
	clear(*whiteList)

	if idx.opt.Debug {
		startTime0 := time.Now()
		log.Debugf("%s (%s bp): start to screen genomes", query.id, humanize.Comma(int64(query.genomeSize)))
		defer func() {
			log.Debugf("%s (%s bp): finished screening genomes in %s",
				query.id, humanize.Comma(int64(query.genomeSize)), time.Since(startTime0))
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
	if windows == 1 {
		window = lenSeq
	}

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

		if i == windows-1 && windows > 1 { // sort k-mers and remove duplicates
			for j = range *_kmers {
				util.UniqUint64s((*_kmersW)[j])
			}
		}

		idx.lh.RecycleMaskResult(_kmers, locses)
	}

	// ------------------------------------------------------
	// 2. search k-mers and return the most similar genomes

	m := idx.poolGSearchDetailResultsMap.Get().(*map[uint64]*GSearchScreenResultDetail)

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

					var r *GSearchScreenResultDetail
					if r, ok = (*m)[refBatchAndIdxUint64]; !ok {
						r = idx.poolGSearchDetailResult.Get().(*GSearchScreenResultDetail)

						r.BatchGenomeIndex = append(r.BatchGenomeIndex, refBatchAndIdxUint64)

						(*m)[refBatchAndIdxUint64] = r
					}

					// r.Hits[sr.IQuery]++       // add count to the probe
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
		idx.RecycleGSearchScreenDetailResultsMap(m)
		return nil, nil
	}

	// collect and store with a list
	rs := idx.poolGSearchDetailResults.Get().(*[]*GSearchScreenResultDetail)
	for _, r := range *m {
		*rs = append(*rs, r)
	}
	clear(*m)
	idx.RecycleGSearchScreenDetailResultsMap(m)

	// 2.3) handle chunked genomes

	// merge search result from genome chunks, if has split genome
	if idx.hasGenomeChunks {
		var r, rp *GSearchScreenResultDetail
		var i, j int
		var a uint64
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
				// for _i, v = range r.Hits {
				// 	rp.Hits[_i] += v
				// }

				idx.RecycleGSearchScreenDetailResult(r)
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
	slices.SortFunc(*rs, func(a, b *GSearchScreenResultDetail) int {
		return cmp.Compare(b.Score, a.Score)
	})
	if topN > 0 && len(*rs) > topN {
		*rs = (*rs)[:topN]
	}

	// fmt.Printf("query\tsubject\tscore\thitMasks\thitKmers\thitKmerAvgLen\n")
	var refBatchAndIdxUint64 uint64
	for _, r := range *rs {
		tmp := make([]uint64, len(r.BatchGenomeIndex))
		copy(tmp, r.BatchGenomeIndex)
		sortutil.Uint64s(tmp)
		for _, refBatchAndIdxUint64 = range r.BatchGenomeIndex {
			(*whiteList)[refBatchAndIdxUint64] = &tmp
		}

		// hitKmers := 0
		// hitMasks := 0
		// for _, v := range r.Hits {
		// 	if v > 0 {
		// 		hitKmers += int(v)
		// 		hitMasks++
		// 	}
		// }
		// fmt.Printf("%s\t%s\t%d\t%d\t%d\t%.1f\n",
		// 	query.id,
		// 	idx.BatchGenomeIndex2GenomeID[r.BatchGenomeIndex[0]], r.Score,
		// 	hitMasks,
		// 	hitKmers, float64(r.Score)/float64(hitKmers),
		// )
	}

	// fmt.Println(whiteList)

	idx.RecycleGSearchScreenDetailResults(rs)

	return whiteList, nil
}

// GSearchAlign2 align fragments of a query to candidates genomes.
// Different from GSearchAlign, this method directly extract candidates genomes for alignment.
func (idx *Index) GSearchAlign2(query *GQuery, fragLen int, minFragLen int, genomeIds *map[uint64]*[]uint64, minAF float64, maxQueryConcurrency int, gcInterval uint64) error {
	debug := idx.opt.Debug

	if debug {
		startTime0 := time.Now()
		log.Debugf("%s (%s bp): start to align query genome fragments", query.id, humanize.Comma(int64(query.genomeSize)))
		defer func() {
			log.Debugf("%s (%s bp): finished aligning query genome fragments in %.3f seconds",
				query.id, humanize.Comma(int64(query.genomeSize)), time.Since(startTime0).Seconds())
		}()
	}

	// --------------------------------------------------------------------------------
	// Step 1. cut query genome into fragments and pre-compute their k-mer entries

	qfrags, qfragLens := seqs2fragments(&query.seqs, fragLen, minFragLen)
	if len(*qfrags) == 0 {
		return fmt.Errorf("no fragments for alignment, are the genome too fragmented with all sequences shorter than the fragment size?")
	}

	// Pre-compute query-side k-mer entries once: the same qfrags is compared
	// against every candidate genome, so collecting its canonical k-mers per
	// goroutine was pure waste. entriesA is read concurrently below; do not
	// mutate it after this point.
	indexer := idx.poolFragmentComparator.Get().(*FragmentComparator)
	entriesA, err := indexer.IndexA(qfrags)
	if err != nil {
		idx.poolFragmentComparator.Put(indexer)
		return err
	}
	idx.poolFragmentComparator.Put(indexer)

	// --------------------------------------------------------------------------------
	// Step 3. collect alignment results

	ch := make(chan *GSearchResult, maxQueryConcurrency)
	done := make(chan int)

	go func() {
		rs := poolGSearchResults.Get().(*[]*GSearchResult)
		*rs = (*rs)[:0]

		for r := range ch {
			*rs = append(*rs, r)
		}

		slices.SortFunc(*rs, func(a, b *GSearchResult) int {
			if d := cmp.Compare(b.ANI, a.ANI); d != 0 {
				return d
			}
			if d := cmp.Compare(b.AFq, a.AFq); d != 0 {
				return d
			}
			return cmp.Compare(b.AFs, a.AFs)
		})

		query.result = rs

		done <- 1
	}()

	// --------------------------------------------------------------------------------
	// Step2. align fragments to candidate genomes

	var wg sync.WaitGroup
	tokens := make(chan int, maxQueryConcurrency)

	alignOption := &wfa.Options{GlobalAlignment: true}

	// only keep one copy of batchIDAndRefIDs for chunks belonging to the same genome
	toDelete := make([]uint64, 0, len(*genomeIds))
	for id, ids := range *genomeIds {
		// fmt.Println(id, *ids)
		if id != (*ids)[0] {
			toDelete = append(toDelete, id)
		}
	}
	for _, id := range toDelete {
		delete(*genomeIds, id)
	}

	// -----------------------------------------------------------
	// process bar
	var pbs *mpb.Progress
	var bar *mpb.Bar
	var chDuration chan time.Duration
	var doneDuration chan int
	if debug {
		pbs = mpb.New(mpb.WithWidth(40), mpb.WithOutput(os.Stderr))
		bar = pbs.AddBar(int64(len(*genomeIds)),
			mpb.PrependDecorators(
				decor.Name("checked subject genomes: ", decor.WC{W: len("checked subject genomes: "), C: decor.DindentRight}),
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

	fcpus := float64(idx.opt.NumCPUs)

	// -----------------------------------------------------------

	for _, batchIDAndRefIDs := range *genomeIds {
		// -------------------------------------------------------------

		tokens <- 1
		wg.Add(1)

		go func(batchIDAndRefIDs *[]uint64) {
			timeStart := time.Now()
			defer func() {
				<-tokens
				wg.Done()
				if debug {
					chDuration <- time.Duration(float64(time.Since(timeStart)) / fcpus)
				}
			}()

			// -------------------------------------------------------------
			// 1. read the genome sequences

			var g *genome.Genome
			genomes := make([]*genome.Genome, len(*batchIDAndRefIDs))
			maxSubjectGenomeSize := idx.opt.MaxSubjectGenomeSize
			for i, batchIDAndRefID := range *batchIDAndRefIDs {
				genomeBatch := int(batchIDAndRefID >> BITS_GENOME_IDX)
				genomeIdx := int(batchIDAndRefID & MASK_GENOME_IDX)

				rdr := <-idx.poolGenomeRdrs[genomeBatch]

				_g, err := rdr.Seqs(genomeIdx)
				if err != nil {
					checkError(fmt.Errorf("fail to read genome sequence for batch %d, genome index %d: %s",
						genomeBatch, genomeIdx, err))
				}

				if i == 0 { // use the first one for later use
					g = _g
				} else { // append sequences to the first genome object
					g.Seqs = append(g.Seqs, _g.Seqs...)
					_g.Seqs = nil
					g.NumSeqs += _g.NumSeqs
					g.GenomeSize += _g.GenomeSize
				}

				if maxSubjectGenomeSize > 0 && g.GenomeSize > maxSubjectGenomeSize {
					log.Warningf("skipped subject genome %s (>= %s bp) which exceeds the maximum allowed size of %s, consider increasing --max-subject-genome-size",
						idx.BatchGenomeIndex2GenomeID[(*batchIDAndRefIDs)[0]],
						humanize.Comma(int64(g.GenomeSize)),
						humanize.Comma(int64(maxSubjectGenomeSize)))

					idx.poolGenomeRdrs[genomeBatch] <- rdr
					for _, gx := range genomes {
						genome.RecycleGenome(gx)
						return
					}
					break
				}

				genomes[i] = _g // recycle them later

				// fmt.Printf("%s, %d seqs, %d bp\n", _g.ID, _g.NumSeqs, _g.GenomeSize)

				idx.poolGenomeRdrs[genomeBatch] <- rdr
			}

			// fmt.Fprintf(os.Stderr, "%s vs %s\n", query.id, g.ID)

			// -------------------------------------------------------------
			// 2. find similar fragment pairs

			sfrags, sfragLens := seqs2fragments(&g.Seqs, fragLen, minFragLen)

			sfragsRC := poolFragments.Get().(*[][]byte)
			n := len(*sfrags)
			if cap(*sfragsRC) >= n {
				*sfragsRC = (*sfragsRC)[:n]
				clear(*sfragsRC)
			} else {
				*sfragsRC = (*sfragsRC)[:cap(*sfragsRC)]
				clear(*sfragsRC)
				for len(*sfragsRC) < n {
					*sfragsRC = append(*sfragsRC, nil)
				}
			}

			fcpr := idx.poolFragmentComparator.Get().(*FragmentComparator)
			pairs, err := fcpr.CompareWithIndexedA(entriesA, sfrags)
			if err != nil {
				checkError(fmt.Errorf("fail to find similar fragments: %s", err))
			}
			sortutil.Uint64s(*pairs)
			// fmt.Fprintf(os.Stderr, "%s vs %s: %d pairs\n", query.id, g.ID, len(*pairs))

			// -------------------------------------------------------------
			// 3. sequence alignment

			cpr := idx.poolSeqComparator.Get().(*SeqComparator)

			algn := wfa.New(wfa.DefaultPenalties, alignOption)
			algn.AdaptiveReduction(wfa.DefaultAdaptiveOption)

			minQcovHSP := idx.seqCompareOption.MinAlignedFraction
			minPIdent := idx.seqCompareOption.MinIdentity

			ma := poolFragAlignResultMap.Get().(*map[uint32]*[]*Chain2Result)
			mb := poolFragAlignResultMap.Get().(*map[uint32]*[]*Chain2Result)

			var a, b, b2 []byte
			var ia, ib uint64
			var cr, cr2 *SeqComparatorResult
			var ls *[]*Chain2Result
			var ok bool
			var c *Chain2Result
			for _, p := range *pairs {
				ia, ib = p>>32, p&4294967295
				a = (*qfrags)[ia]
				b = (*sfrags)[ib]

				// fmt.Printf("%s\t%s\t%d\t%d\t", query.id, g.ID, ia, ib)
				// fmt.Printf("%s\n%s\n", a, b)

				// -----------------------------------------------
				// a) pseudo alignment
				cpr.RecycleIndex()
				err = cpr.Index(a)
				if err != nil {
					checkError(fmt.Errorf("fail to index query fragment: %s", err))
				}

				// positive strand
				cr, err = cpr.Compare(0, uint32(len(a)), b, len(a))
				if err != nil {
					checkError(fmt.Errorf("fail to compare query fragment and subject fragment: %s", err))
				}

				// negative strand
				b2 = (*sfragsRC)[ib]
				if b2 == nil {
					b2 = make([]byte, len(b))
					copy(b2, b)
					RC(b2)
					(*sfragsRC)[ib] = b2
				}
				cr2, err = cpr.Compare(0, uint32(len(a)), b2, len(a))
				if err != nil {
					checkError(fmt.Errorf("fail to compare query fragment and rc subject fragment: %s", err))
				}

				if cr == nil && cr2 == nil {
					// fmt.Printf("no hit\n")
					continue
				}

				if cr != nil && cr2 != nil { // both strands has hits
					// choose the strand with the longer aligned length, if tie, choose the positive strand
					if (*cr.Chains)[0].QEnd-(*cr.Chains)[0].QBegin < (*cr2.Chains)[0].QEnd-(*cr2.Chains)[0].QBegin {
						RecycleSeqComparatorResult(cr)
						cr = cr2
						b = b2
					} else {
						RecycleSeqComparatorResult(cr2)
					}
				} else if cr == nil { // only has hit in the negative strand
					cr = cr2
					b = b2
				} // only has hit in the positive strand

				// -----------------------------------------------
				// b) base-level alignment

				c = (*cr.Chains)[0]   // choose the first chain with the highest chaining score
				(*cr.Chains)[0] = nil // avoid being recycled before we finish processing c

				// fmt.Printf("q[%d,%d]\tt[%d,%d]\t", c.QBegin+1, c.QEnd+1, c.TBegin+1, c.TEnd+1)

				_qseq, _tseq, _, _, _, _, err := extendMatch(a, b, c.QBegin, c.QEnd+1, c.TBegin, c.TEnd+1, idx.opt.ExtendLength2, c.TBegin, idx.opt.ExtendLength2, false)
				if err != nil {
					checkError(fmt.Errorf("fail to extend aligned region: %s", err))
				}

				cigar, err := algn.Align(_qseq, _tseq)
				if err != nil {
					checkError(fmt.Errorf("fail to align sequences: %s", err))
				}

				c.AlignedBasesQ = cigar.QEnd - cigar.QBegin + 1
				c.AlignedLength = int(cigar.AlignLen)
				c.MatchedBases = int(cigar.Matches)
				c.Gaps = int(cigar.Gaps)
				c.AlignedFraction = float64(c.AlignedBasesQ) / float64(cr.QueryLen) * 100
				if c.AlignedFraction > 100 {
					c.AlignedFraction = 100
				}
				c.PIdent = float64(c.MatchedBases) / float64(cigar.AlignLen) * 100

				c.Evalue = c.AlignedFraction * c.PIdent // just for sorting, not the real e-value

				// fmt.Printf("%s\t%s\t%d\t%d\tpident: %.2f, alen: %d, matches: %d, af: %.2f, score: %.2f\n",
				// 	query.id, g.ID, ia, ib, c.PIdent, c.AlignedLength, c.MatchedBases, c.AlignedFraction, c.Evalue)

				// -----------------------------------------------
				// c) filter and store the result

				if c.PIdent >= minPIdent && c.AlignedFraction >= minQcovHSP {
					if ls, ok = (*ma)[uint32(ia)]; !ok {
						ls = poolChains2.Get().(*[]*Chain2Result)
						*ls = (*ls)[:0]
						(*ma)[uint32(ia)] = ls
					}
					*ls = append(*ls, c)

					c.Score = int(ia)    // for finding the corresponding subject fragment in the reciprocal comparison
					c.BitScore = int(ib) // for finding the corresponding subject fragment in the reciprocal comparison

					if ls, ok = (*mb)[uint32(ib)]; !ok {
						ls = poolChains2.Get().(*[]*Chain2Result)
						*ls = (*ls)[:0]
						(*mb)[uint32(ib)] = ls
					}
					*ls = append(*ls, c)
				} else {
					poolChain2.Put(c)
				}

				// -----------------------------------------------
				// clean up
				wfa.RecycleAlignmentResult(cigar)
				RecycleSeqComparatorResult(cr)
			}

			// -------------------------------------------------------------
			// 4. orthologous fragments between two genomes are identified when they showed reciprocal best hit in alignment results.

			fsort := func(a, b *Chain2Result) int {
				if d := cmp.Compare(b.Evalue, a.Evalue); d != 0 {
					return d
				}
				if d := cmp.Compare(a.Score, b.Score); d != 0 {
					return d
				}
				return cmp.Compare(a.BitScore, b.BitScore)
			}

			for _, ls = range *ma {
				if len(*ls) > 1 {
					slices.SortFunc(*ls, fsort)
				}
			}

			for _, ls = range *mb {
				if len(*ls) > 1 {
					slices.SortFunc(*ls, fsort)
				}
			}

			// ANI result
			gr := poolGSearchResult.Get().(*GSearchResult)
			gr.Reset()
			gr.BatchGenomeIndex = (*batchIDAndRefIDs)[0]
			gr.GenomeSize = g.GenomeSize
			gr.NumSeqs = g.NumSeqs

			var _ia, _ib uint32
			var ls2 *[]*Chain2Result
			for _ia, ls = range *ma {
				_ib = uint32((*ls)[0].BitScore)

				if ls2, ok = (*mb)[_ib]; !ok {
					continue
				}
				if (*ls2)[0].Score != int(_ia) {
					continue
				}

				// reciprocal best hit

				c = (*ls)[0]
				// fmt.Printf("%s\t%s\t%d\t%d\tpident: %.2f, alen: %d, matches: %d, af: %.2f, gaps: %d\n",
				// 	query.id, g.ID, _ia, _ib, c.PIdent, c.AlignedLength, c.MatchedBases, c.AlignedFraction, c.Gaps)

				gr.AlignedFragments++
				gr.AlignedLength += c.AlignedLength - c.Gaps
				gr.AlignedMatches += c.MatchedBases
				gr.Pidents = append(gr.Pidents, c.PIdent)
			}

			// gr.ANI = float64(gr.AlignedMatches) / float64(gr.AlignedLength) // shouldn't do this
			sumPident := 0.0
			for _, p := range gr.Pidents {
				sumPident += p
			}
			gr.ANI = sumPident / float64(len(gr.Pidents)) / 100
			gr.AFq = float64(gr.AlignedLength) / float64(qfragLens)
			gr.AFs = float64(gr.AlignedLength) / float64(sfragLens)
			if gr.AFq > 1 {
				gr.AFq = 1
			}
			if gr.AFs > 1 {
				gr.AFs = 1
			}
			gr.Score = gr.ANI // * gr.AF

			if gr.AFq < minAF {
				poolGSearchResult.Put(gr)
			} else {
				ch <- gr
			}

			// manually recycle ma and mb
			for _, ls = range *ma {
				for _, c = range *ls {
					poolChain2.Put(c)
				}
				*ls = (*ls)[:0]
				poolChains2.Put(ls)
			}
			for _, ls = range *mb {
				*ls = (*ls)[:0]
				poolChains2.Put(ls)
			}
			clear(*ma)
			clear(*mb)
			poolFragAlignResultMap.Put(ma)
			poolFragAlignResultMap.Put(mb)

			// -------------------------------------------------------------
			// 5. clean up

			wfa.RecycleAligner(algn)
			idx.poolSeqComparator.Put(cpr)
			RecycleFragmentCompareResult(pairs)
			idx.poolFragmentComparator.Put(fcpr)
			recycleFragments(sfrags)
			recycleFragments(sfragsRC)
			for _, g := range genomes {
				genome.RecycleGenome(g)
			}

		}(batchIDAndRefIDs)
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

	// --------------------------------------------------------------------------------
	// Step4. align fragments to candidate genomes

	recycleFragments(qfrags)
	RecycleResultOfIndexA(entriesA)

	return nil
}

var poolFragAlignResultMap = &sync.Pool{New: func() interface{} {
	// poolChains2 provides pool of *[]*Chain2Result
	tmp := make(map[uint32]*[]*Chain2Result, 4096)
	return &tmp
}}

func RecycleFragAlignResultMap(m *map[uint32]*[]*Chain2Result) {
	for _, rs := range *m {
		RecycleChaining2Result(rs)
	}
	clear(*m)
	poolFragAlignResultMap.Put(m)
}

// -------------------------------------------------------------------------

// GSearchResult represents a subject genome hit.
type GSearchResult struct {
	BatchGenomeIndex uint64
	GenomeSize       int
	NumSeqs          int

	AlignedFragments int
	AlignedLength    int
	AlignedMatches   int
	Pidents          []float64

	ANI   float64
	AFq   float64
	AFs   float64
	Score float64 // for sorting
}

// Reset zeros all fields. Call this when (re)obtaining a GSearchResult from
// the pool, since the pool may return an instance with stale accumulator
// fields (AlignedFragments / AlignedLength / AlignedMatches) that the
// caller increments rather than overwrites.
func (r *GSearchResult) Reset() {
	r.BatchGenomeIndex = math.MaxUint64
	r.GenomeSize = 0
	r.NumSeqs = 0

	r.AlignedFragments = 0
	r.AlignedLength = 0
	r.AlignedMatches = 0
	r.Pidents = r.Pidents[:0]

	r.ANI = 0
	r.AFq = 0
	r.AFs = 0
	r.Score = 0
}

var poolGSearchResult = &sync.Pool{New: func() interface{} {
	return &GSearchResult{
		Pidents: make([]float64, 0, 10240),
	}
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
		r.Reset()
		poolGSearchResult.Put(r)
	}
	*rs = (*rs)[:0]
	poolGSearchResults.Put(rs)
}

func RecycleGSearchResultMap(rs *map[uint64]*GSearchResult) {
	for _, r := range *rs {
		r.Reset()
		poolGSearchResult.Put(r)
	}
	clear(*rs)
	poolGSearchResultMap.Put(rs)
}
