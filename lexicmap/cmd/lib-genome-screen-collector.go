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
	gsearchScreenBatchSize     = 512
	gsearchScreenPageBits      = 10
	gsearchScreenPageSize      = 1 << gsearchScreenPageBits
	gsearchScreenPageMask      = gsearchScreenPageSize - 1
	gsearchScreenPagesPerBatch = 1 << (BITS_GENOME_IDX - gsearchScreenPageBits)
)

type gsearchScreenHit struct {
	batchGenomeIndex uint64
	mask             int32
	length           uint8
}

// The two-level table keeps lookup O(1) without reserving all 2^17 genome
// slots for every batch. Pages are allocated only when a batch has hits in the
// corresponding range of genome indexes.
type gsearchScreenSlotPage [gsearchScreenPageSize]*GSearchScreenResultDetail
type gsearchScreenSlots [gsearchScreenPagesPerBatch]*gsearchScreenSlotPage

var poolGSearchScreenSlotPages = &sync.Pool{New: func() any {
	return &gsearchScreenSlotPage{}
}}

var poolGSearchScreenSlots = &sync.Pool{New: func() any {
	return &gsearchScreenSlots{}
}}

type gsearchScreenCollector struct {
	idx             *Index
	screenMaskCount int

	channels    []chan []gsearchScreenHit
	results     [][]*GSearchScreenResultDetail
	matchCounts []uint64
	workersWG   sync.WaitGroup
	batchPool   sync.Pool
}

func newGSearchScreenCollector(idx *Index, screenMaskCount, nWorkers int) *gsearchScreenCollector {
	nWorkers = min(max(nWorkers, 1), idx.info.GenomeBatches)
	c := &gsearchScreenCollector{
		idx:             idx,
		screenMaskCount: screenMaskCount,
		channels:        make([]chan []gsearchScreenHit, nWorkers),
		results:         make([][]*GSearchScreenResultDetail, nWorkers),
		matchCounts:     make([]uint64, nWorkers),
	}
	c.batchPool.New = func() any {
		return make([]gsearchScreenHit, 0, gsearchScreenBatchSize)
	}
	for i := range nWorkers {
		c.channels[i] = make(chan []gsearchScreenHit, 4)
		c.workersWG.Add(1)
		go c.collect(i)
	}
	return c
}

func (c *gsearchScreenCollector) newBuffers() [][]gsearchScreenHit {
	return make([][]gsearchScreenHit, len(c.channels))
}

func (c *gsearchScreenCollector) add(buffers [][]gsearchScreenHit, hit gsearchScreenHit) {
	worker := int(hit.batchGenomeIndex>>BITS_GENOME_IDX) % len(c.channels)
	batch := buffers[worker]
	if batch == nil {
		batch = c.batchPool.Get().([]gsearchScreenHit)[:0]
	}
	batch = append(batch, hit)
	if len(batch) == cap(batch) {
		c.channels[worker] <- batch
		buffers[worker] = nil
		return
	}
	buffers[worker] = batch
}

func (c *gsearchScreenCollector) flush(buffers [][]gsearchScreenHit) {
	for worker, batch := range buffers {
		if len(batch) > 0 {
			c.channels[worker] <- batch
			buffers[worker] = nil
		}
	}
}

func (c *gsearchScreenCollector) finish() {
	for _, ch := range c.channels {
		close(ch)
	}
	c.workersWG.Wait()
}

func (c *gsearchScreenCollector) collect(worker int) {
	defer c.workersWG.Done()

	slotsByBatch := make([]*gsearchScreenSlots, c.idx.info.GenomeBatches)
	results := make([]*GSearchScreenResultDetail, 0, 1024)
	var nMatches uint64

	var filter *map[uint64]bool
	if c.idx.filterByTaxId {
		filter = c.idx.poolTaxIDfilter.Get().(*map[uint64]bool)
	}

	for batch := range c.channels[worker] {
		for _, hit := range batch {
			if filter != nil && !c.idx.keepGenomeByTaxID(filter, hit.batchGenomeIndex) {
				continue
			}
			nMatches++

			genomeBatch := int(hit.batchGenomeIndex >> BITS_GENOME_IDX)
			genomeIndex := int(hit.batchGenomeIndex & MASK_GENOME_IDX)
			slots := slotsByBatch[genomeBatch]
			if slots == nil {
				slots = poolGSearchScreenSlots.Get().(*gsearchScreenSlots)
				slotsByBatch[genomeBatch] = slots
			}
			pageIndex := genomeIndex >> gsearchScreenPageBits
			page := (*slots)[pageIndex]
			if page == nil {
				page = poolGSearchScreenSlotPages.Get().(*gsearchScreenSlotPage)
				(*slots)[pageIndex] = page
			}

			pageSlot := genomeIndex & gsearchScreenPageMask
			r := (*page)[pageSlot]
			if r == nil {
				r = c.idx.poolGSearchDetailResult.Get().(*GSearchScreenResultDetail)
				r.BatchGenomeIndex = append(r.BatchGenomeIndex, hit.batchGenomeIndex)
				r.SumPrefix = 0
				if r.LongestMatches == nil {
					r.LongestMatches = make([]uint8, c.screenMaskCount)
				}
				(*page)[pageSlot] = r
				results = append(results, r)
			}

			mask := int(hit.mask)
			previous := r.LongestMatches[mask]
			if previous < hit.length {
				r.LongestMatches[mask] = hit.length
				r.SumPrefix += uint64(hit.length - previous)
			}
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
					poolGSearchScreenSlotPages.Put(page)
					(*slots)[j] = nil
				}
			}
			poolGSearchScreenSlots.Put(slots)
			slotsByBatch[i] = nil
		}
	}
	c.results[worker] = results
	c.matchCounts[worker] = nMatches
}

func (c *gsearchScreenCollector) matchCount() uint64 {
	var n uint64
	for _, count := range c.matchCounts {
		n += count
	}
	return n
}

func collectGSearchScreenResultsSerial(
	idx *Index,
	ch <-chan *[]kv.SearchResult,
	screenMaskSlots []int32,
	screenMaskCount int,
) ([]*GSearchScreenResultDetail, uint64) {
	m := idx.poolGSearchDetailResultsMap.Get().(*map[uint64]*GSearchScreenResultDetail)
	var nMatches uint64
	var filter *map[uint64]bool
	if idx.filterByTaxId {
		filter = idx.poolTaxIDfilter.Get().(*map[uint64]bool)
	}

	for srs := range ch {
		for i := range *srs {
			sr := &(*srs)[i]
			iMask := sr.IQuery
			if screenMaskSlots != nil {
				iMask = int(screenMaskSlots[iMask])
				if iMask < 0 {
					continue
				}
			}
			for _, refpos := range sr.Values {
				batchGenomeIndex := refpos >> BITS_NONE_IDX
				if filter != nil && !idx.keepGenomeByTaxID(filter, batchGenomeIndex) {
					continue
				}
				nMatches++

				r := (*m)[batchGenomeIndex]
				if r == nil {
					r = idx.poolGSearchDetailResult.Get().(*GSearchScreenResultDetail)
					r.BatchGenomeIndex = append(r.BatchGenomeIndex, batchGenomeIndex)
					r.SumPrefix = 0
					if r.LongestMatches == nil {
						r.LongestMatches = make([]uint8, screenMaskCount)
					}
					(*m)[batchGenomeIndex] = r
				}

				previous := r.LongestMatches[iMask]
				if previous < sr.Len {
					r.LongestMatches[iMask] = sr.Len
					r.SumPrefix += uint64(sr.Len - previous)
				}
			}
		}
		kv.RecycleSearchResults(srs)
	}

	if filter != nil {
		clear(*filter)
		idx.poolTaxIDfilter.Put(filter)
	}
	results := make([]*GSearchScreenResultDetail, 0, len(*m))
	for _, r := range *m {
		results = append(results, r)
	}
	clear(*m)
	idx.poolGSearchDetailResultsMap.Put(m)
	return results, nMatches
}

func clearGSearchScreenResults(results [][]*GSearchScreenResultDetail) {
	for i, batch := range results {
		clear(batch)
		results[i] = nil
	}
}
