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

import "sync"

const (
	minGapSize                  = 5
	maxPooledGapRegionsCapacity = 4096
)

var poolGapRegions = &sync.Pool{New: func() any {
	gaps := make([]uint64, 0, 16)
	return &gaps
}}

// findGapRegions returns half-open intervals of N/n runs at least minGapSize
// bases long, matching the former regular expression [Nn]{5,}. Each interval
// packs its start in the high 32 bits and its end in the low 32 bits.
func findGapRegions(seq []byte) *[]uint64 {
	var gaps *[]uint64
	for i := 0; i < len(seq); {
		if seq[i] != 'N' && seq[i] != 'n' {
			i++
			continue
		}

		start := i
		i++
		for i < len(seq) && (seq[i] == 'N' || seq[i] == 'n') {
			i++
		}
		if i-start >= minGapSize {
			if gaps == nil {
				gaps = poolGapRegions.Get().(*[]uint64)
				*gaps = (*gaps)[:0]
			}
			*gaps = append(*gaps, packGapRegion(start, i))
		}
	}
	return gaps
}

func packGapRegion(start, end int) uint64 {
	return uint64(uint32(start))<<32 | uint64(uint32(end))
}

func unpackGapRegion(region uint64) (int, int) {
	return int(region >> 32), int(uint32(region))
}

func recycleGapRegions(gaps *[]uint64) {
	if gaps == nil {
		return
	}
	if cap(*gaps) > maxPooledGapRegionsCapacity {
		return
	}
	*gaps = (*gaps)[:0]
	poolGapRegions.Put(gaps)
}
