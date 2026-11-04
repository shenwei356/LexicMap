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
	"bytes"
	"math"
	"slices"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
	"github.com/shenwei356/kmers"
	"github.com/shenwei356/lexichash"
)

const (
	// Query filling uses the index command's defaults, regardless of index settings.
	defaultSeedMaxDesert        = 100
	defaultSeedInDesertDist     = 50
	defaultQueryDesertThreshold = 50000 // default query length threshold in bases
)

// queryDesertSeeds groups additional k-mers by mask. locIndexes maps each
// k-mer to a separate entry after the original mask entries in query locations.
type queryDesertSeeds struct {
	kmers      [][]uint64 // distinct additional k-mers per mask
	locIndexes [][]int    // flat location slot for each additional k-mer
}

// fillQueryDesertsEnabled keeps the threshold exclusive; zero disables filling.
func fillQueryDesertsEnabled(length, threshold int) bool {
	return threshold > 0 && length > threshold
}

// queryMaskSkipRegions matches index masking around N gaps. Call it only for
// queries selected for desert filling, so other queries retain their behavior.
func queryMaskSkipRegions(sequence []byte) []int {
	gaps := findGapRegions(sequence, 5)
	if gaps == nil {
		return nil
	}
	defer recycleGapRegions(gaps)
	regions := make([]int, 0, len(*gaps)*2)
	for _, gap := range *gaps {
		start, end := unpackGapRegion(gap)
		regions = append(regions, start, end-1)
	}
	return regions
}

// fillQueryDeserts mirrors index desert selection: mask each gap with 1 kb
// flanks, then prefer upstream winners within half the expected seed distance.
// Scratch belongs to this call; concurrent queries never share a window masker.
func fillQueryDeserts(lh *lexichash.LexicHash, sequence []byte, seeds []uint64, locations [][]int, threshold int) (*queryDesertSeeds, [][]int, error) {
	if !fillQueryDesertsEnabled(len(sequence), threshold) {
		return nil, locations, nil
	}
	w, err := lh.NewWindowMasker()
	if err != nil {
		return nil, locations, err
	}
	positions := make([]int, 0, len(seeds))
	for i, seed := range seeds {
		if seed == 0 {
			continue
		}
		for _, loc := range locations[i] {
			positions = append(positions, int(uint32(loc)>>1))
		}
	}
	extra := &queryDesertSeeds{kmers: make([][]uint64, len(seeds)), locIndexes: make([][]int, len(seeds))}
	// Do not grow the borrowed LexicHash slice: its pool requires one slot per mask.
	locs := slices.Clone(locations)
	// The same code can be assigned to multiple masks by shorter-prefix
	// fallback, so retain its mask as well as its sequence when deduplicating.
	type seedKey struct {
		mask int
		code uint64
	}
	type positionKey struct {
		mask, loc int
		code      uint64
	}
	byKmer := make(map[seedKey]int)
	seen := make(map[positionKey]struct{})
	for i, seed := range seeds {
		if seed != 0 {
			byKmer[seedKey{i, seed}] = i
			for _, loc := range locations[i] {
				seen[positionKey{i, loc, seed}] = struct{}{}
			}
		}
	}
	// Merge seed occurrences while preserving each mask's prefix lookup.
	add := func(mask int, code uint64, loc int) {
		position := positionKey{mask, loc, code}
		if _, exists := seen[position]; exists {
			return
		}
		seen[position] = struct{}{}
		key := seedKey{mask, code}
		locIndex, found := byKmer[key]
		if !found {
			locIndex = len(locs)
			byKmer[key] = locIndex
			locs = append(locs, nil)
			extra.kmers[mask] = append(extra.kmers[mask], code)
			extra.locIndexes[mask] = append(extra.locIndexes[mask], locIndex)
		}
		locs[locIndex] = append(locs[locIndex], loc)
	}
	k8 := uint8(lh.K)
	ccc, ggg := util.Ns(1, k8), util.Ns(2, k8)
	ttt := uint64(1)<<(k8<<1) - 1
	// Index filling prefers upstream positions and then the forward strand.
	// Include its choices in both orientations so reverse-complement queries
	// can recover the same index seeds. Coordinates always refer to the input.
	for _, reverse := range []bool{false, true} {
		scanPositions := positions
		if reverse {
			sequence = bytes.Clone(sequence)
			RC(sequence)
			// Ambiguous bases are encoded as A during masking, so merely
			// mirroring the first pass's positions is not equivalent to masking
			// the reverse complement. Recompute and retain those global seeds.
			rSeeds, rLocs, err := lh.MaskKnownDistinctPrefixes(sequence, queryMaskSkipRegions(sequence), true)
			if err != nil {
				return nil, locations, err
			}
			defer lh.RecycleMaskResult(rSeeds, rLocs)
			scanPositions = nil
			for mask, code := range *rSeeds {
				if code == 0 || code == ccc || code == ggg || code == ttt || util.IsLowComplexityDust(code, k8) {
					continue
				}
				for _, loc := range (*rLocs)[mask] {
					p := int(uint32(loc) >> 1)
					scanPositions = append(scanPositions, p)
					add(mask, code, (len(sequence)-lh.K-p)<<1|(loc&1^1))
				}
			}
		}
		slices.Sort(scanPositions)
		scanPositions = append(scanPositions, len(sequence)-lh.K)
		gaps := findGapRegions(sequence, 5)
		if gaps != nil {
			defer recycleGapRegions(gaps)
		}
		previous := 0
		firstGap := 0 // windows advance monotonically; never rescan earlier N gaps
		for _, position := range scanPositions {
			if position-previous < defaultSeedMaxDesert {
				previous = position
				continue
			}
			start, end := max(0, previous-1000), min(len(sequence), position+1000+lh.K)
			codes, masks, err := w.Mask(sequence[start:end])
			if err != nil {
				return nil, locations, err
			}
			// Exclude k-mers crossing runs of N, just as index construction does.
			if gaps != nil {
				for firstGap < len(*gaps) {
					_, b := unpackGapRegion((*gaps)[firstGap])
					if b > start {
						break
					}
					firstGap++
				}
				for _, gap := range (*gaps)[firstGap:] {
					a, b := unpackGapRegion(gap)
					if a-lh.K+1 >= end-lh.K+1 {
						break
					}
					for p := max(start, a-lh.K+1); p < min(b, end-lh.K+1); p++ {
						masks[(p-start)*2], masks[(p-start)*2+1] = -1, -1
					}
				}
			}
			choose := func(p int) bool {
				for strand := 0; strand < 2; strand++ {
					i := (p-start)*2 + strand
					code, mask := codes[i], masks[i]
					if mask < 0 || code == 0 || code == ccc || code == ggg || code == ttt || util.IsLowComplexityDust(code, k8) {
						continue
					}
					loc := p<<1 | strand
					if reverse {
						loc = (len(sequence)-lh.K-p)<<1 | (strand ^ 1)
					}
					add(mask, code, loc)
					return true
				}
				return false
			}
			const radius = defaultSeedInDesertDist / 2
			for p := previous + defaultSeedInDesertDist; p < position; {
				downstream := p + 1
				limit := p - radius
				found := false
				for ; p > limit; p-- {
					if choose(p) {
						found = true
						break
					}
				}
				if !found {
					if downstream >= position {
						break
					}
					limit = min(downstream+radius, position-1)
					for p = downstream; p < limit; p++ {
						if choose(p) {
							break
						}
					}
				}
				p += defaultSeedInDesertDist
			}
			previous = position
		}
	}
	return extra, locs, nil
}

// addReverseSeeds adds suffix probes for desert seeds, preserving the flat
// location mapping used by both in-memory and spilling anchor collectors.
func (extra *queryDesertSeeds) addReverseSeeds(lh *lexichash.LexicHash, reversed *[][]uint64, mapping *[][]int) {
	if extra == nil {
		return
	}
	for mask, codes := range extra.kmers {
		for i, code := range codes {
			code = kmers.MustReverse(code, lh.K)
			matches := lh.MaskKmer(code)
			best, hash := 0, uint64(math.MaxUint64)
			for _, m := range *matches {
				if h := lh.Masks[m] ^ code; h < hash {
					best, hash = m, h
				}
			}
			lh.RecycleMaskKmerResult(matches)
			(*reversed)[best] = append((*reversed)[best], code)
			(*mapping)[best] = append((*mapping)[best], extra.locIndexes[mask][i])
		}
	}
}

// stream searches extra prefix probes and translates the per-mask result index
// into the flat query-location index before any collector sees the result.
func (extra *queryDesertSeeds) stream(searcher seedDataStreamer, begin, end int, prefix uint8, consume func(kv.SearchResult) error) error {
	if extra == nil {
		return nil
	}
	return searcher.Search2Stream(extra.kmers[begin:end], prefix, true, false, func(r kv.SearchResult) error {
		r.IQuery = extra.locIndexes[r.IQuery][r.IQuery2]
		return consume(r)
	})
}

// queryDesertStreamer extends the prefix pass while retaining synchronous,
// borrowed-value callbacks and the existing suffix pass.
type queryDesertStreamer struct {
	seedDataStreamer
	extra      *queryDesertSeeds
	begin, end int // mask range owned by this seed searcher
}

func (s queryDesertStreamer) SearchStream(codes []uint64, prefix uint8, checkStrand, reverse bool, consume func(kv.SearchResult) error) error {
	if err := s.seedDataStreamer.SearchStream(codes, prefix, checkStrand, reverse, consume); err != nil {
		return err
	}
	return s.extra.stream(s.seedDataStreamer, s.begin, s.end, prefix, consume)
}
