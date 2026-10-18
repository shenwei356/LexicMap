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
	"sync/atomic"
	"unsafe"

	rtree "github.com/shenwei356/LexicMap/lexicmap/cmd/tree"
	"github.com/twotwotwo/sorts"
)

// The default is configurable through genome search --max-genome-cache-memory.
// It covers complete cached trees and slot metadata across all active queries
// on an Index; private scratch, subjects, pools and runtime memory are extra.
// Genome compare uses its separate compareGenomeCache budget, not this default.
const defaultGenomeSearchFragmentIndexMemory int64 = 4 << 30

// queryFragmentIndex publishes one immutable tree or indexing error. A nil tree
// after once completes selects private indexing. Published trees are never pooled.
type queryFragmentIndex struct {
	once sync.Once
	tree *rtree.CompactTree
	err  error
}

// queryFragmentIndexes prepares only fragments that reach fine alignment.
// Keys, values, and prefix directories are shared until all subjects finish;
// alignment scratch stays private to each comparator.
// In genome search, each query owns its indexes: they are never shared between
// queries and are released after its subjects finish. Laziness avoids building
// fragments that never reach fine alignment; reuse needs no LRU or eviction.
// One query's indexes do not multiply with its subject count, but concurrent
// queries retain separate sets. The shared byte budget bounds their combined
// retention; it does not preallocate memory or limit total process RSS.
// Ordinary search already builds one private query tree outside its subject
// loop, so it never uses this cache. The former entries-only cache saved k-mer
// collection/sorting, but rebuilt the tree and prefix directory for each subject
// even on a cache hit. Sharing CompactTree removes those repeated builds.
type queryFragmentIndexes struct {
	idx       *Index
	fragments []queryFragmentIndex
	reserved  atomic.Int64 // shared tree and slot bytes released after all subjects finish
	// Genome comparisons charge lazy trees to the pinned prepared genome.
	// Its owner releases the charge on eviction, rather than after one pair.
	compareCache *compareGenomeCache
	compareEntry *compareGenomeCacheEntry
}

// reserve charges exact-size arrays before allocation. A full budget falls
// back to private indexing, so pinned genomes never wait for cache space.
func (c *queryFragmentIndexes) reserve(n int64) bool {
	if c.compareCache != nil {
		return c.compareCache.reserveFragmentEntries(c.compareEntry, n)
	}
	for {
		used := c.idx.fragmentIndexBytes.Load()
		if n > c.idx.opt.MaxGenomeCacheMemory-used {
			return false
		}
		if c.idx.fragmentIndexBytes.CompareAndSwap(used, used+n) {
			c.reserved.Add(n)
			return true
		}
	}
}

// newQueryFragmentIndexes avoids cache preparation when only one subject exists.
func newQueryFragmentIndexes(idx *Index, fragments, subjects int) *queryFragmentIndexes {
	c := &queryFragmentIndexes{idx: idx}
	if subjects > 1 && fragments > 0 {
		// Charge every slot before allocating; no fixed fragment-count cutoff
		// should prevent a larger configured budget from covering long queries.
		charge := queryFragmentSlotBytes(fragments)
		if charge > 0 && c.reserve(charge) {
			c.fragments = make([]queryFragmentIndex, fragments)
		}
	}
	return c
}

// queryFragmentSlotBytes returns -1 for impossible metadata sizes so callers
// reject them without wrapping the charge or attempting a huge allocation.
func queryFragmentSlotBytes(n int) int64 {
	const maxCharge = int64(1<<63 - 1)
	object := int64(unsafe.Sizeof(queryFragmentIndexes{}))
	slot := int64(unsafe.Sizeof(queryFragmentIndex{}))
	if int64(n) > (maxCharge-object)/slot {
		return -1
	}
	return object + int64(n)*slot
}

// index borrows a complete immutable tree when available. A nil receiver or full
// cache uses the original Index path without waiting for memory. The caller must
// RecycleIndex before reusing the comparator or releasing the cache owner.
func (c *queryFragmentIndexes) index(cpr *SeqComparator, fragment int, sequence []byte) error {
	if c == nil || fragment >= len(c.fragments) {
		return cpr.Index(sequence)
	}
	slot := &c.fragments[fragment]
	built := false
	slot.once.Do(func() {
		entries, err := cpr.collectIndexEntries(sequence)
		slot.err = err
		built = true
		if err != nil {
			return
		}
		sorts.ByUint64(rtree.BatchEntries(entries))
		if c.reserve(rtree.CompactTreeBytes(cpr.options.K, entries)) {
			slot.tree = rtree.NewCompactTree(cpr.options.K, entries)
		}
	})
	if slot.err != nil {
		return slot.err
	}
	if slot.tree != nil {
		cpr.compact = slot.tree
		return nil
	}
	if built {
		// Reuse the first caller's sorted scratch on a cache miss.
		cpr.indexSortedEntries(cpr.entries)
		return nil
	}
	return cpr.Index(sequence)
}

// close drops shared trees after every subject has finished and detached its
// comparator. Never recycle them: pools must not outlive the released charge.
func (c *queryFragmentIndexes) close() {
	for i := range c.fragments {
		c.fragments[i].tree = nil
	}
	c.fragments = nil
	c.idx.fragmentIndexBytes.Add(-c.reserved.Swap(0))
}
