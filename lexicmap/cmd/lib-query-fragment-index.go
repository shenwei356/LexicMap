// Copyright © 2026 Wei Shen <shenwei356@gmail.com>

package cmd

import (
	"sync"
	"sync/atomic"
	"unsafe"

	rtree "github.com/shenwei356/LexicMap/lexicmap/cmd/tree"
	"github.com/twotwotwo/sorts"
)

// This limit covers cached sorted entries across all active genome-search
// queries on an Index. Entries are allocated on demand; at the default fragment
// length, a fully cached 10 Mb query needs about 296 MiB. Tree/chaining scratch,
// subjects, and query-local slot metadata are additional, so this is not an RSS cap.
const maxGenomeSearchFragmentIndexMemory int64 = 1 << 30

// Bound per-query slot metadata too, including chromosome-sized queries. These
// slots cover about 16.7 Mb at the default fragment length, allowing 10 Mb queries
// to cache their final fragments. Fragments beyond the limit use uncached indexing.
const maxGenomeSearchFragmentIndexSlots = 16384

// queryFragmentIndex publishes one immutable sorted entry slice or indexing
// error. A nil slice after once completes selects the uncached Index path.
type queryFragmentIndex struct {
	once    sync.Once
	entries []rtree.BatchEntry
	err     error
}

// queryFragmentIndexes prepares only fragments that reach fine alignment.
// Entries are shared until all subjects finish; each comparator builds its own
// mutable tree. This preserves existing radix lookup and alignment behavior.
type queryFragmentIndexes struct {
	idx       *Index
	fragments []queryFragmentIndex
	reserved  atomic.Int64 // shared entry bytes released after all subjects finish
}

// newQueryFragmentIndexes avoids cache preparation when only one subject exists.
func newQueryFragmentIndexes(idx *Index, fragments, subjects int) *queryFragmentIndexes {
	c := &queryFragmentIndexes{idx: idx}
	if subjects > 1 {
		c.fragments = make([]queryFragmentIndex, min(fragments, maxGenomeSearchFragmentIndexSlots))
	}
	return c
}

// index builds private trees from shared sorted entries when available. A nil
// receiver or full cache uses the original Index path without waiting for memory.
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
		n := int64(len(entries)) * int64(unsafe.Sizeof(rtree.BatchEntry{}))
		for {
			used := c.idx.fragmentIndexBytes.Load()
			if n > maxGenomeSearchFragmentIndexMemory-used {
				return
			}
			if c.idx.fragmentIndexBytes.CompareAndSwap(used, used+n) {
				// Exact-size ownership avoids retaining the comparator's larger capacity.
				slot.entries = make([]rtree.BatchEntry, len(entries))
				copy(slot.entries, entries)
				c.reserved.Add(n)
				return
			}
		}
	})
	if slot.err != nil {
		return slot.err
	}
	if slot.entries != nil {
		cpr.indexSortedEntries(slot.entries)
		return nil
	}
	if built {
		// The cache-full builder already filtered and sorted its private entries.
		cpr.indexSortedEntries(cpr.entries)
		return nil
	}
	return cpr.Index(sequence)
}

// close drops shared arrays after every subject has finished. Cached entries
// are not pooled, so later queries do not keep them alive through sync.Pool.
func (c *queryFragmentIndexes) close() {
	for i := range c.fragments {
		c.fragments[i].entries = nil
	}
	c.fragments = nil
	c.idx.fragmentIndexBytes.Add(-c.reserved.Swap(0))
}
