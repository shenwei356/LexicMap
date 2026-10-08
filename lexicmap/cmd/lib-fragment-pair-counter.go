// Copyright © 2026 Wei Shen <shenwei356@gmail.com>

package cmd

import (
	"math/bits"

	rtree "github.com/shenwei356/LexicMap/lexicmap/cmd/tree"
)

// Bound retained dense uint16 counters plus their presence bitmap per comparator.
// The bitmap distinguishes an untouched pair from a uint16 count wrapping to
// zero: each cell needs 17 bits. Eight bytes cover bitmap rounding; output pair
// buffers and the original hash fallback are additional. Allocate only the rows
// needed by a comparison; this ceiling allows 20 Mb genomes to use dense blocks.
const maxFragmentPairCounterBytes = 64 << 20
const maxFragmentPairCounterCells = (maxFragmentPairCounterBytes - 8) * 8 / 17

// Limit rescans of the sorted inputs. Larger or sparse fragment-ID spaces use
// the original hash counter rather than multiplying scanning work without bound.
const maxFragmentPairCounterBlocks = 16

// fragmentPairCounterShape chooses query-row blocks within the scratch budget.
// Fragment IDs come from Val>>1, so their maxima fit the original uint32 encoding.
func fragmentPairCounterShape(a, b []rtree.BatchEntry) (rows, cols, blockRows int, ok bool) {
	var maxA, maxB uint32
	for _, entry := range a {
		maxA = max(maxA, entry.Val>>1)
	}
	for _, entry := range b {
		maxB = max(maxB, entry.Val>>1)
	}
	// Check dimensions in uint64 before converting to int, including 32-bit builds.
	nRows, nCols := uint64(maxA)+1, uint64(maxB)+1
	if nCols > maxFragmentPairCounterCells {
		return 0, 0, 0, false
	}
	nBlockRows := min(nRows, uint64(maxFragmentPairCounterCells)/nCols)
	if (nRows+nBlockRows-1)/nBlockRows > maxFragmentPairCounterBlocks {
		return 0, 0, 0, false
	}
	return int(nRows), int(nCols), int(nBlockRows), true
}

// scanPairsDense counts one bounded block of query rows at a time, rescanning
// the immutable sorted streams for each block. It enumerates every original
// position cross product exactly once and retains the same uint16 arithmetic.
// A false result selects the hash path without allocating counter scratch.
func (cpr *FragmentComparator) scanPairsDense(a, b []rtree.BatchEntry) (*[]uint64, bool) {
	if len(a) == 0 || len(b) == 0 {
		cpr.pairCounts = cpr.pairCounts[:0]
		return cpr.selectFragmentPairs(), true
	}
	rows, cols, blockRows, ok := fragmentPairCounterShape(a, b)
	if !ok {
		return nil, false
	}
	cells := blockRows * cols
	if cap(cpr.pairCounter) < cells {
		cpr.pairCounter = make([]uint16, cells)
	} else {
		cpr.pairCounter = cpr.pairCounter[:cells]
	}
	words := (cells + 63) / 64
	if cap(cpr.pairVisited) < words {
		cpr.pairVisited = make([]uint64, words)
	} else {
		cpr.pairVisited = cpr.pairVisited[:words]
	}
	counter, visited := cpr.pairCounter, cpr.pairVisited
	counts := cpr.pairCounts[:0]
	threshold := cpr.options.MinSharedKmers
	for firstRow := 0; firstRow < rows; firstRow += blockRows {
		lastRow := min(firstRow+blockRows, rows)
		for ia, ib := 0, 0; ia < len(a) && ib < len(b); {
			ka, kb := a[ia].Key, b[ib].Key
			if ka < kb {
				ia++
				continue
			}
			if kb < ka {
				ib++
				continue
			}
			ja, jb := ia+1, ib+1
			for ja < len(a) && a[ja].Key == ka {
				ja++
			}
			for jb < len(b) && b[jb].Key == ka {
				jb++
			}
			for _, ea := range a[ia:ja] {
				row := int(ea.Val >> 1)
				if row < firstRow || row >= lastRow {
					continue
				}
				offset := (row - firstRow) * cols
				for _, eb := range b[ib:jb] {
					cell := offset + int(eb.Val>>1)
					if counter[cell] == 0 {
						visited[cell>>6] |= uint64(1) << (cell & 63)
					}
					counter[cell]++
				}
			}
			ia, ib = ja, jb
		}
		// Visit and clear only cells used by this block. Presence remains set
		// across uint16 wraparound, including zero counts when the threshold is zero.
		for word, mask := range visited {
			if mask == 0 {
				continue
			}
			for mask != 0 {
				cell := word*64 + bits.TrailingZeros64(mask)
				if count := counter[cell]; count >= threshold {
					pair := uint64(firstRow+cell/cols)<<32 | uint64(cell%cols)
					counts = append(counts, fragmentPairCount{pair: pair, count: count})
				}
				counter[cell] = 0
				mask &= mask - 1
			}
			visited[word] = 0
		}
	}
	cpr.pairCounts = counts
	return cpr.selectFragmentPairs(), true
}
