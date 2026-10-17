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
