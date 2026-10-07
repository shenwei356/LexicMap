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
	"sync"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
)

const (
	seedAnchorBatchSize    = 512
	seedGenomeSlotPageBits = 10
	seedGenomeSlotPageSize = 1 << seedGenomeSlotPageBits
	seedGenomeSlotPageMask = seedGenomeSlotPageSize - 1
	seedGenomeSlotPages    = 1 << (BITS_GENOME_IDX - seedGenomeSlotPageBits)
)

// seedAnchor is the compact value transferred from seed searchers to the
// collectors. Keeping the transfer bounded avoids retaining another complete
// copy of all anchors while collectors run in parallel.
type seedAnchor struct {
	batchGenomeIndex uint64
	qBegin           int32
	tBegin           int32
	length           uint8
	qrc              bool
	trc              bool
}

type seedGenomeSlotPage [seedGenomeSlotPageSize]*SearchResult
type seedGenomeSlots [seedGenomeSlotPages]*seedGenomeSlotPage

var poolSeedGenomeSlotPages = &sync.Pool{New: func() any {
	return &seedGenomeSlotPage{}
}}

var poolSeedGenomeSlots = &sync.Pool{New: func() any {
	return &seedGenomeSlots{}
}}

var poolSerialSeedSearchResultsMap = &sync.Pool{New: func() any {
	m := make(map[int]*SearchResult, 1024)
	return &m
}}

type seedAnchorCollector struct {
	idx *Index

	channels  []chan []seedAnchor
	results   [][]*SearchResult
	workersWG sync.WaitGroup
	batchPool sync.Pool
}

func newSeedAnchorCollector(idx *Index, nWorkers int) *seedAnchorCollector {
	nWorkers = min(max(nWorkers, 1), idx.info.GenomeBatches)
	c := &seedAnchorCollector{
		idx:      idx,
		channels: make([]chan []seedAnchor, nWorkers),
		results:  make([][]*SearchResult, nWorkers),
	}
	c.batchPool.New = func() any {
		return make([]seedAnchor, 0, seedAnchorBatchSize)
	}
	for i := range nWorkers {
		c.channels[i] = make(chan []seedAnchor, 4)
		c.workersWG.Add(1)
		go c.collect(i)
	}
	return c
}

func (c *seedAnchorCollector) newBuffers() [][]seedAnchor {
	return make([][]seedAnchor, len(c.channels))
}

func (c *seedAnchorCollector) add(buffers [][]seedAnchor, anchor seedAnchor) {
	worker := int(anchor.batchGenomeIndex>>BITS_GENOME_IDX) % len(c.channels)
	batch := buffers[worker]
	if batch == nil {
		batch = c.batchPool.Get().([]seedAnchor)[:0]
	}
	batch = append(batch, anchor)
	if len(batch) == cap(batch) {
		c.channels[worker] <- batch
		buffers[worker] = nil
		return
	}
	buffers[worker] = batch
}

func (c *seedAnchorCollector) flush(buffers [][]seedAnchor) {
	for worker, batch := range buffers {
		if len(batch) > 0 {
			c.channels[worker] <- batch
			buffers[worker] = nil
		}
	}
}

func (c *seedAnchorCollector) finish() {
	for _, ch := range c.channels {
		close(ch)
	}
	c.workersWG.Wait()
}

func (c *seedAnchorCollector) collect(worker int) {
	defer c.workersWG.Done()

	slotsByBatch := make([]*seedGenomeSlots, c.idx.info.GenomeBatches)
	results := make([]*SearchResult, 0, 1024)

	var filter *map[uint64]bool
	if c.idx.filterByTaxId {
		filter = c.idx.poolTaxIDfilter.Get().(*map[uint64]bool)
	}

	for batch := range c.channels[worker] {
		for _, anchor := range batch {
			if filter != nil && !c.idx.keepGenomeByTaxID(filter, anchor.batchGenomeIndex) {
				continue
			}

			genomeBatch := int(anchor.batchGenomeIndex >> BITS_GENOME_IDX)
			genomeIndex := int(anchor.batchGenomeIndex & MASK_GENOME_IDX)
			slots := slotsByBatch[genomeBatch]
			if slots == nil {
				slots = poolSeedGenomeSlots.Get().(*seedGenomeSlots)
				slotsByBatch[genomeBatch] = slots
			}
			pageIndex := genomeIndex >> seedGenomeSlotPageBits
			page := (*slots)[pageIndex]
			if page == nil {
				page = poolSeedGenomeSlotPages.Get().(*seedGenomeSlotPage)
				(*slots)[pageIndex] = page
			}

			pageSlot := genomeIndex & seedGenomeSlotPageMask
			r := (*page)[pageSlot]
			if r == nil {
				subs := poolSubs.Get().(*[]SubstrPair)
				r = poolSearchResult.Get().(*SearchResult)
				r.BatchGenomeIndex = anchor.batchGenomeIndex
				r.GenomeBatch = genomeBatch
				r.GenomeIndex = genomeIndex
				r.GenomeSize = 0
				r.NumSeqs = 0
				r.Subs = subs
				r.Score = 0
				r.resetChainRegions()
				r.SimilarityDetails = nil
				r.AlignedFraction = 0
				(*page)[pageSlot] = r
				results = append(results, r)
			}

			sub := SubstrPair{}
			sub.QBegin = anchor.qBegin
			sub.TBegin = anchor.tBegin
			sub.Len = anchor.length
			sub.QRC = anchor.qrc
			sub.TRC = anchor.trc
			*r.Subs = append(*r.Subs, sub)
		}
		c.batchPool.Put(batch[:0])
	}

	if filter != nil {
		clear(*filter)
		c.idx.poolTaxIDfilter.Put(filter)
	}
	for i, slots := range slotsByBatch {
		if slots != nil {
			for j, page := range slots {
				if page != nil {
					clear(page[:])
					poolSeedGenomeSlotPages.Put(page)
					(*slots)[j] = nil
				}
			}
			poolSeedGenomeSlots.Put(slots)
			slotsByBatch[i] = nil
		}
	}
	c.results[worker] = results
}

func (c *seedAnchorCollector) resultCount() int {
	n := 0
	for _, results := range c.results {
		n += len(results)
	}
	return n
}

func (c *seedAnchorCollector) anchorCount() uint64 {
	var n uint64
	for _, results := range c.results {
		for _, r := range results {
			n += uint64(len(*r.Subs))
		}
	}
	return n
}

func (c *seedAnchorCollector) clearResults() {
	for i, results := range c.results {
		clear(results)
		c.results[i] = nil
	}
}

// collectSeedAnchorsSerial preserves the lower-overhead collector used when
// only one collector worker is available. The batched sharded path pays for
// itself with parallelism, but is slower when all work must remain serial.
func collectSeedAnchorsSerial(
	idx *Index,
	ch <-chan *[]kv.SearchResult,
	locses *[][]int,
	reverseLocses *[][]int,
	genomeIDs *map[uint64]*[]uint64,
) []*SearchResult {
	m := poolSerialSeedSearchResultsMap.Get().(*map[int]*SearchResult)
	var filter *map[uint64]bool
	if idx.filterByTaxId {
		filter = idx.poolTaxIDfilter.Get().(*map[uint64]bool)
	}
	filterByGenomeID := genomeIDs != nil
	K := idx.k

	for srs := range ch {
		for i := range *srs {
			sr := &(*srs)[i]
			kPrefix := int(sr.Len)
			var queryLocs []int
			if !sr.IsSuffix {
				queryLocs = (*locses)[sr.IQuery]
			} else {
				queryLocs = (*locses)[(*reverseLocses)[sr.IQuery][sr.IQuery2]]
			}
			for _, encodedPosQ := range queryLocs {
				rcQ := encodedPosQ&MASK_STRAND > 0
				posQ := encodedPosQ >> BITS_STRAND
				for _, refpos := range sr.Values {
					batchGenomeIndex := refpos >> BITS_NONE_IDX
					if filterByGenomeID {
						if _, ok := (*genomeIDs)[batchGenomeIndex]; !ok {
							continue
						}
					}
					if filter != nil && !idx.keepGenomeByTaxID(filter, batchGenomeIndex) {
						continue
					}

					posT := int(refpos << BITS_IDX >> BITS_IDX_FLAGS)
					rvT := refpos&MASK_REVERSE > 0
					rcT := refpos>>BITS_REVERSE&MASK_STRAND > 0
					var beginQ, beginT int
					if !rvT {
						if rcQ {
							beginQ = posQ + K - kPrefix
						} else {
							beginQ = posQ
						}
						if rcT {
							beginT = posT + K - kPrefix
						} else {
							beginT = posT
						}
					} else {
						if rcQ {
							beginQ = posQ
						} else {
							beginQ = posQ + K - kPrefix
						}
						if rcT {
							beginT = posT
						} else {
							beginT = posT + K - kPrefix
						}
					}

					sub := SubstrPair{}
					sub.QBegin = int32(beginQ)
					sub.TBegin = int32(beginT)
					sub.Len = uint8(kPrefix)
					sub.QRC = rcQ
					sub.TRC = rcT

					key := int(batchGenomeIndex)
					r := (*m)[key]
					if r == nil {
						subs := poolSubs.Get().(*[]SubstrPair)
						r = poolSearchResult.Get().(*SearchResult)
						r.BatchGenomeIndex = batchGenomeIndex
						r.GenomeBatch = key >> BITS_GENOME_IDX
						r.GenomeIndex = key & MASK_GENOME_IDX
						r.GenomeSize = 0
						r.NumSeqs = 0
						r.Subs = subs
						r.Score = 0
						r.resetChainRegions()
						r.SimilarityDetails = nil
						r.AlignedFraction = 0
						(*m)[key] = r
					}
					*r.Subs = append(*r.Subs, sub)
				}
			}
		}
		kv.RecycleSearchResults(srs)
	}

	if filter != nil {
		clear(*filter)
		idx.poolTaxIDfilter.Put(filter)
	}
	results := make([]*SearchResult, 0, len(*m))
	for _, r := range *m {
		results = append(results, r)
	}
	clear(*m)
	poolSerialSeedSearchResultsMap.Put(m)
	return results
}

func clearSeedSearchResults(results [][]*SearchResult) {
	for i, batch := range results {
		clear(batch)
		results[i] = nil
	}
}
