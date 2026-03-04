// Copyright © 2026-2026 Wei Shen <shenwei356@gmail.com>
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

package rangeindex

import (
	"slices"
	"sync"
)

// RangeIndex stores entries as:
// high 32 bits: location
// low 32 bits: value
type RangeIndex struct {
	data   []uint64
	sorted bool
}

var rangeIndexPool = sync.Pool{
	New: func() any {
		return &RangeIndex{
			data:   make([]uint64, 0, 4096),
			sorted: true,
		}
	},
}

// NewRangeIndex retrieves an object from the pool and resets it.
func NewRangeIndex() *RangeIndex {
	ri := rangeIndexPool.Get().(*RangeIndex)
	ri.data = ri.data[:0]
	ri.sorted = true
	return ri
}

// Release returns the object back to the pool.
// Large backing arrays are dropped to avoid memory bloat.
func (ri *RangeIndex) Release() {
	const maxCap = 1 << 20 // 1M entries safety limit

	if cap(ri.data) > maxCap {
		ri.data = make([]uint64, 0, 4096)
	}
	rangeIndexPool.Put(ri)
}

// Add inserts a (loc, value) pair.
// loc and value are both uint32.
func (ri *RangeIndex) Add(loc, value uint32) {
	ri.data = append(ri.data, (uint64(loc)<<32)|uint64(value))
	ri.sorted = false
}

// Query returns all entries whose loc is within the given range.
// Attention: the length of result might be 0.
// Note that: each returned slice element contain both location and value.
// You can get the location by v >> 32, and value v & 4294967295
func (ri *RangeIndex) Query(left, right uint32) []uint64 {
	if !ri.sorted {
		// sort.Slice(ri.data, func(i, j int) bool {
		// 	return ri.data[i] < ri.data[j]
		// })
		slices.Sort(ri.data) // simpler than code above
		// sortutil.Uint64s(ri.data) // consume extra space

		ri.sorted = true
	}

	data := ri.data
	var target uint64
	var lo, hi, mid int

	// binary search: first index with loc >= left

	target = uint64(left)
	lo, hi = 0, len(data)
	for lo < hi {
		mid = (lo + hi) >> 1
		if data[mid]>>32 < target {
			lo = mid + 1
		} else {
			hi = mid
		}
	}
	start := lo

	// binary search: first index with loc > target

	target = uint64(right)
	hi = len(data)
	for lo < hi {
		mid := (lo + hi) >> 1
		if data[mid]>>32 <= target {
			lo = mid + 1
		} else {
			hi = mid
		}
	}
	end := lo

	return ri.data[start:end]
}

// Unpack extracts location and its value from one element of query result.
func Unpack(v uint64) (uint32, uint32) {
	return uint32(v >> 32), uint32(v & 4294967295)
}
