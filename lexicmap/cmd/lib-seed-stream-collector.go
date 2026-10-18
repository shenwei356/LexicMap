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

import "github.com/shenwei356/LexicMap/lexicmap/cmd/kv"

// streamSeedAnchors expands borrowed decoded positions into the collector's
// bounded transfer batches. No callback retains Values beyond its invocation.
// The caller owns the searcher token and joins every producer before finish.
func streamSeedAnchors(searcher seedDataStreamer, collector *seedAnchorCollector,
	kmers []uint64, reverseKmers [][]uint64, locses, reverseLocses [][]int,
	minPrefix uint8, genomeIDs *map[uint64]*[]uint64,
) (batches, values uint64, err error) {
	buffers := collector.newBuffers()
	// Drain partial batches on errors too; their arenas are recycled by Search.
	defer collector.flush(buffers)
	emit := func(a seedAnchor) error {
		if genomeIDs != nil {
			if _, ok := (*genomeIDs)[a.batchGenomeIndex]; !ok {
				return nil
			}
		}
		collector.add(buffers, a)
		return nil
	}
	consume := func(sr kv.SearchResult) error {
		batches++
		values += uint64(len(sr.Values))
		locs := locses[sr.IQuery]
		if sr.IsSuffix {
			locs = locses[reverseLocses[sr.IQuery][sr.IQuery2]]
		}
		return expandSeedDataBatch(collector.idx.k, sr, locs, emit)
	}
	err = searcher.SearchStream(kmers, minPrefix, true, false, consume)
	if err == nil {
		err = searcher.Search2Stream(reverseKmers, minPrefix, true, true, consume)
	}
	return
}
