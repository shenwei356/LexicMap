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
	"container/list"
	"fmt"
	"sync"
	"unsafe"

	rtree "github.com/shenwei356/LexicMap/lexicmap/cmd/tree"
	"github.com/twotwotwo/sorts"
)

// preparedCompareGenome owns immutable data shared by concurrent genome pairs.
// It has no alignment results; fragments borrow genome sequences. Cache-owned
// buffers are dropped on eviction instead of being retained in sync.Pool.
// For pairs A-B and A-C, A's query-fragment indexes can be built once and shared.
// Reuse avoids repeated preparation, but retaining several prepared genomes
// can increase peak memory versus uncached pairs, so all data shares a budget.
type preparedCompareGenome struct {
	genome    *GQuery               // sequences and metadata, never a pair's result
	fragments *[][]byte             // contig-aware fragments
	seeds     [][]uint64            // sampled query seeds in ordinary ANI mode
	concat    []byte                // subject forward sequence, gap, and RC
	sketch    *subjectSketch        // ordinary ANI subject lookup
	entries   []rtree.BatchEntry    // sorted canonical entries in OrthoANI mode
	bytes     int64                 // retained charge; lazy growth guarded by cache.mu
	indexes   *queryFragmentIndexes // lazy immutable fine-alignment indexes, owned by this cache entry
}

// prepareCompareGenome prepares both roles once using fixed command options.
func (idx *Index) prepareCompareGenome(q *GQuery, fragLen, minSeedLen int, orthoANI bool) (*preparedCompareGenome, error) {
	if q == nil {
		return nil, nil // empty FASTA: skip its pairs
	}
	p := &preparedCompareGenome{genome: q}
	p.fragments, _ = seqs2fragments(&q.seqs, fragLen, minSeedLen)
	if p.fragments == nil || len(*p.fragments) == 0 {
		return nil, fmt.Errorf("no nonempty genome fragments for alignment")
	}
	if orthoANI {
		cpr := idx.poolFragmentComparator.Get().(*FragmentComparator)
		// Fresh entries avoid retaining an oversized pooled buffer in each cache entry.
		entries, err := cpr.collectEntries(p.fragments, 0, nil)
		idx.poolFragmentComparator.Put(cpr)
		if err != nil {
			return nil, err
		}
		sorts.ByUint64(rtree.BatchEntries(entries))
		p.entries = entries
	} else {
		p.seeds = make([][]uint64, len(*p.fragments))
		for i, frag := range *p.fragments {
			seeds, err := sampleQueryFragment(frag, nil)
			if err != nil {
				return nil, fmt.Errorf("failed to sample query fragment: %w", err)
			}
			p.seeds[i] = seeds
		}
		concat, sketch, err := idx.buildComparisonSubject(q, fragLen, false)
		if err != nil {
			return nil, err
		}
		p.concat, p.sketch = *concat, sketch
	}
	p.bytes = p.retainedBytes()
	return p, nil
}

// retainedBytes charges backing capacities, including a conservative allowance
// for Go map tables. This bounds cache accounting, not total process RSS: active
// uncached pairs, alignment scratch, readers, pools, and runtime memory are extra.
func (p *preparedCompareGenome) retainedBytes() int64 {
	sliceBytes := int64(unsafe.Sizeof([]byte{}))
	// Include cache entry, LRU node, and map bookkeeping in the fixed allowance.
	n := int64(unsafe.Sizeof(*p)+unsafe.Sizeof(GQuery{})) + 512
	q := p.genome
	n += int64(cap(q.id) + cap(q.bigSeq) + cap(q.quals)*int(unsafe.Sizeof((*[]byte)(nil))) + cap(q.skipRegions)*int(unsafe.Sizeof(int(0))))
	n += int64(cap(q.seqs)) * (int64(unsafe.Sizeof((*[]byte)(nil))) + sliceBytes)
	for _, seq := range q.seqs {
		n += int64(cap(*seq))
	}
	for _, qual := range q.quals {
		n += sliceBytes + int64(cap(*qual))
	}
	n += int64(cap(*p.fragments)+cap(p.seeds)) * sliceBytes
	for _, seeds := range p.seeds {
		n += int64(cap(seeds)) * 8
	}
	n += int64(cap(p.entries)) * int64(unsafe.Sizeof(rtree.BatchEntry{}))
	n += int64(cap(p.concat))
	if s := p.sketch; s != nil {
		// Fresh maps start from these capacity hints. Charge above the slot size
		// to cover spare slots, directory/group overhead, and table growth.
		hint := sampledKmerMapCapacity(q.genomeSize, gsa3SamplingScale)
		n += int64(max(hint, len(*s.sampledKmerMap))+max(hint/16, len(*s.repeatedKmerMap))) * 64
		n += int64(cap(*s.repeatedKmerPositions)) * int64(unsafe.Sizeof(repeatedKmerPosition{}))
		n += int64(cap(s.contigBounds)) * int64(unsafe.Sizeof([2]int{}))
		n += int64(unsafe.Sizeof(*s))
	}
	return n
}

// compareGenomeCacheEntry remains pinned through result output. ready publishes
// one load to all concurrent callers, including errors and oversized-genome skips.
type compareGenomeCacheEntry struct {
	key      string
	ready    chan struct{}
	value    *preparedCompareGenome
	err      error
	refs     int           // active pair owners; guarded by cache.mu
	retained bool          // charged against the budget
	idle     *list.Element // idle LRU node; active entries cannot be evicted
}

// compareGenomeCache shares preparation without waiting for pinned entries to
// become evictable. A miss that cannot fit runs uncached, avoiding deadlocks when
// a pair needs two genomes larger than the available budget.
// This is genome compare's --max-genome-cache-memory budget, independent of
// IndexSearchingOptions.MaxGenomeCacheMemory used by genome search.
type compareGenomeCache struct {
	mu                   sync.Mutex
	budget, used, peak   int64
	entries              map[string]*compareGenomeCacheEntry
	idle                 list.List // most recently released first
	load                 func(string) (*preparedCompareGenome, error)
	cacheFragmentIndexes bool // enabled only when more than one pair can reuse fine indexes
	hits, loads          uint64
	evictions            uint64
}

// newCompareGenomeCache scopes keys and fixed preparation options to one command.
func newCompareGenomeCache(budget int64, load func(string) (*preparedCompareGenome, error)) *compareGenomeCache {
	return &compareGenomeCache{budget: budget, entries: make(map[string]*compareGenomeCacheEntry), load: load}
}

// acquire pins an entry; callers must release it even when loading fails.
func (c *compareGenomeCache) acquire(key string) *compareGenomeCacheEntry {
	c.mu.Lock()
	if e := c.entries[key]; e != nil {
		e.refs++
		c.hits++
		if e.idle != nil {
			c.idle.Remove(e.idle)
			e.idle = nil
		}
		c.mu.Unlock()
		<-e.ready
		return e
	}
	e := &compareGenomeCacheEntry{key: key, ready: make(chan struct{}), refs: 1}
	c.entries[key] = e
	c.loads++
	c.mu.Unlock()
	value, err := c.load(key)
	c.mu.Lock()
	e.value, e.err = value, err
	if err == nil && value != nil && value.bytes <= c.budget {
		for c.used > c.budget-value.bytes && c.idle.Len() > 0 {
			old := c.idle.Back().Value.(*compareGenomeCacheEntry)
			c.idle.Remove(old.idle)
			old.idle = nil
			delete(c.entries, old.key)
			c.used -= old.value.bytes
			old.value = nil // do not retain evicted buffers in pools
			c.evictions++
		}
		if c.used <= c.budget-value.bytes {
			e.retained = true
			c.used += value.bytes
			// Slot metadata is bounded and charged too. Only retained genomes
			// cache fine indexes; oversized or pinned-budget misses run uncached.
			if c.cacheFragmentIndexes && value.fragments != nil {
				n := len(*value.fragments)
				charge := queryFragmentSlotBytes(n)
				if charge > 0 && charge <= c.budget-c.used {
					value.indexes = &queryFragmentIndexes{fragments: make([]queryFragmentIndex, n), compareCache: c, compareEntry: e}
					value.bytes += charge
					c.used += charge
				}
			}
			c.peak = max(c.peak, c.used)
		}
	}
	close(e.ready)
	c.mu.Unlock()
	return e
}

// reserveFragmentEntries grows a pinned entry's charge without evicting other
// prepared genomes for optional fine-alignment data. Cache misses keep working
// with private scratch when the remaining budget is too small.
func (c *compareGenomeCache) reserveFragmentEntries(e *compareGenomeCacheEntry, n int64) bool {
	c.mu.Lock()
	defer c.mu.Unlock()
	if !e.retained || e.refs == 0 || n > c.budget-c.used {
		return false
	}
	c.used += n
	e.value.bytes += n
	c.peak = max(c.peak, c.used)
	return true
}

// release makes a successful cached entry evictable after its final pair finishes.
func (c *compareGenomeCache) release(e *compareGenomeCacheEntry) {
	if e == nil {
		return
	}
	c.mu.Lock()
	defer c.mu.Unlock()
	e.refs--
	if e.refs == 0 {
		if e.retained {
			e.idle = c.idle.PushFront(e)
		} else {
			delete(c.entries, e.key)
			e.value = nil
		}
	}
}

// close drops the idle cache after every pair and output owner has released it.
func (c *compareGenomeCache) close() {
	c.mu.Lock()
	defer c.mu.Unlock()
	clear(c.entries)
	c.idle.Init()
	c.used = 0
}
