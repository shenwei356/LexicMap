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

package tree

import (
	"math/bits"
	"unsafe"
)

// CompactTree is an immutable radix index represented by sorted leaves. Internal
// nodes are implicit: each subtree is a contiguous key interval, whose common
// prefix is the LCP of its first and last keys. All backing arrays are pointer-free.
// Memory note (2026-10): extrapolating two 1020-bp samples at k=31 to a 10-Mb
// query gives about 2.5-3.0 GiB for retained pointer Trees versus 320-336 MiB
// for CompactTrees. These are sample-based capacity estimates, not RSS bounds;
// sequence complexity, fragmentation, scratch and runtime overhead also matter.
// The saving comes from the compact layout, not immutability alone. Compared
// with caching only sorted entries, compact indexes can retain slightly more
// data, but also share the lookup structure instead of rebuilding it per subject.
// Performance note (2026-10): sharing these indexes saved about 4% end-to-end
// time versus the entries-only cache on several measured genome workloads;
// misses can be slower than Tree lookups and process RSS can increase. This
// does not guarantee a speedup for other inputs. A shared full pointer tree
// was slower and larger in the tested workloads.
// Ordinary search already reuses its tree: benchmark its hit/miss mix before
// substituting CompactTree there. Faster hit-only lookups are not sufficient.
type CompactTree struct {
	k           uint8    // common k-mer length (1..32)
	prefixBases uint8    // zero selects the full sorted interval
	keys        []uint64 // distinct keys in lexicographic order
	offsets     []uint32 // value boundaries, including the final sentinel
	values      []uint32 // duplicate values retain their original sorted-entry order
	starts      []uint32 // key boundaries for the optional prefix directory
}

// compactShape counts exact allocation sizes; callers supply sorted entries.
func compactShape(k uint8, entries []BatchEntry) (unique int, prefix uint8) {
	for i, e := range entries {
		if i == 0 || e.Key != entries[i-1].Key {
			unique++
		}
	}
	if k < 4 || unique < 256 {
		return
	}
	prefix = 4
	for prefix < 8 && prefix < k && unique >= 1<<(2*(prefix+1)) {
		prefix++
	}
	return
}

// CompactTreeBytes computes the retained charge before allocation. Builder
// scratch and allocator/runtime overhead are not part of this cache charge.
func CompactTreeBytes(k uint8, entries []BatchEntry) int64 {
	unique, prefix := compactShape(k, entries)
	directory := 0
	if prefix != 0 {
		directory = 1<<(2*prefix) + 1
	}
	return int64(unsafe.Sizeof(CompactTree{})) + int64(unique)*8 + int64(unique+1)*4 + int64(len(entries))*4 + int64(directory)*4
}

// NewCompactTree copies entries sorted by Key for k in [1,32]. It retains no
// builder scratch and is read-only after construction, including during searches.
func NewCompactTree(k uint8, entries []BatchEntry) *CompactTree {
	unique, prefix := compactShape(k, entries)
	c := &CompactTree{k: k, prefixBases: prefix, keys: make([]uint64, unique), offsets: make([]uint32, unique+1), values: make([]uint32, len(entries))}
	j := 0
	for i, e := range entries {
		if i == 0 || e.Key != entries[i-1].Key {
			c.keys[j] = e.Key
			c.offsets[j] = uint32(i)
			j++
		}
		c.values[i] = e.Val
	}
	c.offsets[j] = uint32(len(entries))
	if prefix != 0 {
		c.starts = make([]uint32, 1<<(2*prefix)+1)
		j = 0
		for bucket := range c.starts {
			for j < len(c.keys) && c.keys[j]>>(2*(k-prefix)) < uint64(bucket) {
				j++
			}
			c.starts[bucket] = uint32(j)
		}
	}
	return c
}

// RetainedBytes includes every owned backing-array capacity and the tree object.
func (c *CompactTree) RetainedBytes() int64 {
	return int64(unsafe.Sizeof(*c)) + int64(cap(c.keys))*8 + int64(cap(c.offsets))*4 + int64(cap(c.values))*4 + int64(cap(c.starts))*4
}

// lower finds the first key at least target within an already bounded subtree.
func (c *CompactTree) lower(lo, hi int, target uint64) int {
	for lo < hi {
		mid := lo + (hi-lo)/2
		if c.keys[mid] < target {
			lo = mid + 1
		} else {
			hi = mid
		}
	}
	return lo
}

// SearchEach preserves Tree's leaf order, duplicate order and exact LCP values.
// Values are borrowed and must not be modified. Most searches use the prefix
// directory followed by a short linear scan or a binary search.
func (c *CompactTree) SearchEach(key uint64, p uint8, fn func(uint64, uint8, []uint32)) bool {
	p = min(max(p, 1), c.k)
	shift := 2 * (c.k - p)
	prefix := key >> shift
	lo, hi := 0, len(c.keys)
	if c.prefixBases != 0 && p >= c.prefixBases {
		bucket := key >> (2 * (c.k - c.prefixBases))
		lo, hi = int(c.starts[bucket]), int(c.starts[bucket+1])
	}
	// Small directory buckets are cheaper to scan than to binary-search.
	if hi-lo <= 8 {
		for lo < hi && c.keys[lo]>>shift < prefix {
			lo++
		}
	} else {
		lo = c.lower(lo, hi, prefix<<shift)
	}
	end := lo
	for end < hi && c.keys[end]>>shift == prefix {
		end++
	}
	hi = end
	if lo == hi && p >= 2 && prefix&15 == 0 {
		// An exact prefix always follows the same path in the original tree.
		// Only misses ending in AA can trigger its historical short-edge match.
		lo, hi = c.legacyRange(key, p)
	}

	for i := lo; i < hi; i++ {
		fn(c.keys[i], uint8(bits.LeadingZeros64(key^c.keys[i])>>1+int(c.k)-32), c.values[c.offsets[i]:c.offsets[i+1]])
	}
	return lo < hi
}

// legacyRange follows the same compressed edges as Tree.searchNode without
// storing nodes. For sorted keys, each edge ends at LCP(first key, last key).
//
//	[ACGA..., ACGT...] -> edge ACG -> two contiguous child intervals
//
// This path preserves historical short-edge matches, including their ordering.
func (c *CompactTree) legacyRange(key uint64, p uint8) (int, int) {
	lo, hi, depth := 0, len(c.keys), uint8(0)
	for lo < hi && depth < c.k {
		shift := 2 * (c.k - depth - 1)
		prefix := key >> shift
		lo = c.lower(lo, hi, prefix<<shift)
		end := lo
		// Upper-bound the selected child's prefix; avoid overflow at k=32.
		for end < hi {
			mid := end + (hi-end)/2
			if c.keys[mid]>>shift <= prefix {
				end = mid + 1
			} else {
				hi = mid
			}
		}
		hi = end
		if lo == hi {
			return 0, 0
		}
		next := uint8(bits.LeadingZeros64(c.keys[lo]^c.keys[hi-1])>>1 + int(c.k) - 32)
		edge, k := next-depth, c.k-depth
		search := key & (^uint64(0) >> (64 - 2*k))
		edgePrefix := (c.keys[lo] >> (2 * (c.k - next))) & (^uint64(0) >> (64 - 2*edge))
		if search>>(2*(k-edge)) == edgePrefix {
			if next >= p {
				return lo, hi
			}
			depth = next
			continue
		}
		atleast := p - depth
		// Keep uint8 subtraction, as in the original compressed-edge comparison.
		if search>>((k-atleast)<<1) == edgePrefix>>((edge-atleast)<<1) {
			return lo, hi
		}
		return 0, 0
	}
	return 0, 0
}

// Walk visits all leaves in the same order as Tree.Walk.
func (c *CompactTree) Walk(fn WalkFn) {
	for i, key := range c.keys {
		if fn(key, c.values[c.offsets[i]:c.offsets[i+1]]) {
			return
		}
	}
}
