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

import "sort"

// Small references are cheaper to scan than to allocate an offset directory.
const minIndexedContigs = 32

// contigLookup belongs to one reference alignment. Its ends use the existing
// concatenated coordinates, including the artificial intervals between contigs.
type contigLookup struct {
	ends   []int
	linear bool // malformed/non-monotone layouts retain the original scan
	used   bool // a reference with only one lookup needs no offset directory
}

// locate preserves the first contig satisfying both original K-tolerant bounds.
// It indexes only fragmented references and builds their directory on reuse.
func (c *contigLookup) locate(sizes []int, interval, begin, end, k int) (int, int, int) {
	if len(sizes) < minIndexedContigs || interval < 0 || k < 1 || c.linear {
		return scanContig(sizes, interval, begin, end, k)
	}
	if !c.used {
		c.used = true
		return scanContig(sizes, interval, begin, end, k)
	}
	if c.ends == nil {
		c.ends = make([]int, len(sizes))
		pos := 0
		for i, size := range sizes {
			last := pos + size - 1
			if size < 1 || last < pos {
				c.ends, c.linear = nil, true
				return scanContig(sizes, interval, begin, end, k)
			}
			c.ends[i] = last
			if i+1 < len(sizes) {
				pos = last + interval + 1
				if pos < last {
					c.ends, c.linear = nil, true
					return scanContig(sizes, interval, begin, end, k)
				}
			}
		}
	}
	// Earlier ends cannot contain end-k. If this first candidate starts after
	// begin+k, all later starts do too, so no later contig can satisfy the scan.
	i := sort.Search(len(c.ends), func(i int) bool { return c.ends[i] >= end-k })
	if i == len(c.ends) {
		return -1, 0, 0
	}
	start := 0
	if i > 0 {
		start = c.ends[i-1] + interval + 1
	}
	if begin+k < start {
		return -1, 0, 0
	}
	return i, start, c.ends[i]
}

// scanContig is the original lookup, including its tolerance and early exit.
func scanContig(sizes []int, interval, begin, end, k int) (int, int, int) {
	start, last := 0, 0
	for i, size := range sizes {
		last += size - 1
		if begin+k >= start && end-k <= last {
			return i, start, last
		} else if end < start {
			break
		}
		last += interval + 1
		start = last
	}
	return -1, 0, 0
}
