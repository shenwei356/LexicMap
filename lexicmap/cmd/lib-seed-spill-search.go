package cmd

import (
	"fmt"
	"sync"
	"time"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
)

// seedDataStreamer is the common streaming API implemented by kv.Searcher and
// kv.InMemorySearcher. Both invoke callbacks synchronously with borrowed Values,
// so anchor expansion must finish before returning. No []kv.SearchResult list
// or complete per-key location array is retained by this search path.
type seedDataStreamer interface {
	// SearchStream searches one query k-mer per mask (the prefix-search pass).
	// Its IQuery includes the searcher's mask-chunk offset.
	SearchStream([]uint64, uint8, bool, bool, func(kv.SearchResult) error) error
	// Search2Stream searches multiple query k-mers per mask (the reverse-k-mer
	// suffix-search pass). IQuery2 identifies the query within that mask's list.
	Search2Stream([][]uint64, uint8, bool, bool, func(kv.SearchResult) error) error
}

// collectAndChainSpilledSeeds implements the opt-in collection/chaining path for
// one query. Read in this order:
//
//  1. Acquire a collection-budget slot and arrange cleanup for every return path.
//  2. Search all mask partitions for prefix and suffix matches. Each synchronous
//     callback expands its borrowed seed data, applies genome/TaxId filters,
//     and adds anchors to one serialized seedSpillStore. The store handles spills.
//  3. Wait for every partition: a genome's anchors are incomplete before this.
//  4. Consume grouped anchors and chain one complete genome at a time. Retain
//     candidate scores and chain bounds for the caller's global Top-N selection.
//
// kmers/locses are indexed by mask. reverseKmers groups reversed query k-mers by
// their new masks, and reverseLocses maps each one back to the original mask's
// locses entry. These query arrays remain borrowed until all searchers finish.
// genomeIDs is nil when no ID filter is requested. arena is caller-owned and
// must stay alive through alignment because returned candidates live in it.
// The returned pointer slice uses the existing pool/recycling protocol.
//
// No Top-N selection occurs here; Index.Search applies it after all genome
// scores are known. Files and the budget slot are released before returning.
func (idx *Index) collectAndChainSpilledSeeds(query *Query, kmers *[]uint64, reverseKmers *[][]uint64, locses, reverseLocses *[][]int, genomeIDs *map[uint64]*[]uint64, arena *seedSearchResultArena, debug bool) (results *[]*seedSearchResult, err error) {
	// Step 1: wait for a budget share before taking any disk-searcher tokens.
	store := newSeedSpillStore(idx.seedMemoryBudget)
	defer func() {
		closeErr := store.close()
		if err == nil && closeErr != nil {
			if results != nil {
				recycleSeedSearchResultSlice(results)
				results = nil
			}
			err = fmt.Errorf("remove seed spill: %w", closeErr)
		}
	}()
	// mu protects the non-concurrent store, TaxId cache, counters, and firstErr.
	// Holding it through expansion/spilling also provides backpressure: other
	// searchers wait with only their fixed decoded batch instead of queuing anchors.
	var mu sync.Mutex
	var firstErr error
	var batches, values uint64
	var filter *map[uint64]bool
	if idx.filterByTaxId {
		filter = idx.poolTaxIDfilter.Get().(*map[uint64]bool)
		defer func() { clear(*filter); idx.poolTaxIDfilter.Put(filter) }()
	}
	// Step 2: consume borrowed Values before returning to the decoder. Prefix
	// IQuery names the original mask; suffix IQuery/IQuery2 first maps back to it.
	// One location value can expand into many anchors if the query seed occurs
	// at multiple positions, so budget accounting happens after this expansion.
	consume := func(sr kv.SearchResult) error {
		mu.Lock()
		defer mu.Unlock()
		if firstErr != nil {
			return firstErr
		}
		batches++
		values += uint64(len(sr.Values))
		queryLocs := (*locses)[sr.IQuery]
		if sr.IsSuffix {
			queryLocs = (*locses)[(*reverseLocses)[sr.IQuery][sr.IQuery2]]
		}
		firstErr = expandSeedDataBatch(idx.k, sr, queryLocs, func(a seedAnchor) error {
			// Keep the original filters ahead of collection/chaining. Rejected
			// genomes contribute no anchors to the memory budget or spill files.
			if genomeIDs != nil {
				if _, ok := (*genomeIDs)[a.batchGenomeIndex]; !ok {
					return nil
				}
			}
			if filter != nil && !idx.keepGenomeByTaxID(filter, a.batchGenomeIndex) {
				return nil
			}
			return store.add(a)
		})
		return firstErr
	}
	started := time.Now()
	var wg sync.WaitGroup
	count := len(idx.Searchers)
	if idx.opt.InMemorySearch {
		count = len(idx.InMemorySearchers)
	}
	for i := 0; i < count; i++ {
		wg.Add(1)
		go func(i int) {
			defer wg.Done()
			var searcher seedDataStreamer
			var begin, end int
			if idx.opt.InMemorySearch {
				s := idx.InMemorySearchers[i]
				searcher = s
				begin = s.ChunkIndex
				end = begin + s.ChunkSize
			} else {
				// A disk searcher has mutable decoding resources. Hold its token
				// for both prefix/suffix passes, including synchronous callbacks.
				idx.searcherTokens[i] <- 1
				defer func() { <-idx.searcherTokens[i] }()
				s := idx.Searchers[i]
				searcher = s
				begin = s.ChunkIndex
				end = begin + s.ChunkSize
			}
			e := searcher.SearchStream((*kmers)[begin:end], idx.opt.MinPrefix, true, false, consume)
			if e == nil {
				e = searcher.Search2Stream((*reverseKmers)[begin:end], idx.opt.MinPrefix, true, true, consume)
			}
			if e != nil {
				mu.Lock()
				if firstErr == nil {
					firstErr = e
				}
				mu.Unlock()
			}
		}(i)
	}
	// Step 3: only after every producer finishes can we know that each genome
	// is complete. A recorded error aborts chaining; deferred close removes runs.
	wg.Wait()
	if firstErr != nil {
		return nil, fmt.Errorf("collect seed anchors: %w", firstErr)
	}
	if debug {
		log.Debugf("%s: bounded seed collection: seed data batches=%d, values=%d, anchors=%d, peak collection buffer=%d bytes, spilled=%t, runs=%d, elapsed=%s", query.seqID, batches, values, store.count, store.peakBufferBytes, store.spilled, len(store.runs), time.Since(started))
	}
	started = time.Now()
	// Step 4: reconstruct/chain complete genomes; retain only scores and bounds.
	results, err = idx.chainSeedSpillStore(store, arena)
	if debug && err == nil {
		log.Debugf("%s: finished chaining collected seed anchors (%d genome hits) in %s", query.seqID, len(*results), time.Since(started))
	}
	return results, err
}

// expandSeedDataBatch emits the Cartesian product of query occurrences (locs)
// and reference locations (sr.Values), using the same coordinate/strand rules
// as the original collector. locs encodes query position plus a strand bit;
// reference values pack genome identity, position, strand, and reverse-key mode.
// sr.Len is the matching prefix/suffix length, not necessarily the full k.
// emit receives independent scalar anchors synchronously; its first error stops
// expansion. Neither borrowed input slice is retained beyond this call.
func expandSeedDataBatch(k int, sr kv.SearchResult, locs []int, emit func(seedAnchor) error) error {
	length := int(sr.Len)
	for _, encodedQ := range locs {
		qrc := encodedQ&MASK_STRAND > 0
		posQ := encodedQ >> BITS_STRAND
		for _, refpos := range sr.Values {
			posT := int(refpos << BITS_IDX >> BITS_IDX_FLAGS)
			reverse := refpos&MASK_REVERSE > 0
			trc := refpos>>BITS_REVERSE&MASK_STRAND > 0
			beginQ, beginT := posQ, posT
			// When the occurrence's strand differs from the stored key's reverse
			// mode, its matching part starts k-length bases into the original
			// k-mer. Keep forward-coordinate begins and the original strand flags.
			if qrc != reverse {
				beginQ += k - length
			}
			if trc != reverse {
				beginT += k - length
			}
			if err := emit(seedAnchor{batchGenomeIndex: refpos >> BITS_NONE_IDX, qBegin: int32(beginQ), tBegin: int32(beginT), length: sr.Len, qrc: qrc, trc: trc}); err != nil {
				return err
			}
		}
	}
	return nil
}

// chainSeedSpillStore consumes genome-grouped anchors after collection finishes.
// It accumulates just the current genome in subs, then runs the existing
// ClearSubstrPairs/Chainer logic at each genome boundary and at end of input.
// subs and chainer scratch may exceed the collection share for a very large
// genome; they are explicitly outside that budget.
//
// Passing candidates retain their packed genome identity, chaining score, and
// copied seedChainRegion bounds in the caller's arena, never subs/scratch aliases.
// All candidates must survive until Index.Search ranks them globally, including
// cutoff ties. On error the result container is recycled; the caller still owns
// the arena and store cleanup. On success the caller owns result-slice recycling.
func (idx *Index) chainSeedSpillStore(store *seedSpillStore, arena *seedSearchResultArena) (results *[]*seedSearchResult, err error) {
	results = poolSeedSearchResults.Get().(*[]*seedSearchResult)
	*results = (*results)[:0]
	defer func() {
		if err != nil {
			recycleSeedSearchResultSlice(results)
			results = nil
		}
	}()
	chainer := idx.poolChainers.Get().(*Chainer)
	defer func() {
		if chainer != nil {
			RecycleChainer(idx.poolChainers, chainer)
		}
	}()
	// id identifies the current complete-genome group. subs and regions reuse
	// storage across small genomes; inline handles a single chain without a
	// separately allocated initial region array.
	var id uint64
	var subs []SubstrPair
	var regions []seedChainRegion
	var inline [1]seedChainRegion
	regions = inline[:0]
	// Finish the current genome before accepting the next one's anchors. Copy
	// bounds into the candidate before resetting scratch for reuse.
	finish := func() {
		if len(subs) == 0 {
			return
		}
		if len(subs) > 1 {
			ClearSubstrPairs(&subs, idx.k)
		}
		if chainer == nil {
			chainer = idx.poolChainers.Get().(*Chainer)
		}
		var score float32
		regions, score = chainer.Chain(&subs, regions[:0])
		if score >= idx.chainingOptions.MinScore {
			r := arena.add(id)
			r.Score = score
			r.chainRegions = append(r.chainRegions, regions...)
			*results = append(*results, r)
		}
		if len(chainer.visited) > chainerInitSize {
			// Follow the existing pool policy: abandon an oversized chainer
			// rather than retaining it in the pool or for subsequent genomes.
			chainer = nil
		}
		// One complete genome and its chainer scratch may exceed the collection
		// budget. Do not retain an oversized anchor array for later genomes.
		if cap(subs) > thresholdNSubsLong {
			subs = nil
		} else {
			subs = subs[:0]
		}
		if cap(regions) > thresholdNSubs {
			regions = inline[:0]
		}
	}
	// consume guarantees contiguous genome groups, including anchors split
	// across different runs. Changing IDs therefore marks a complete genome.
	err = store.consume(func(a seedAnchor) error {
		if len(subs) > 0 && id != a.batchGenomeIndex {
			finish()
		}
		id = a.batchGenomeIndex
		subs = append(subs, SubstrPair{QBegin: a.qBegin, TBegin: a.tBegin, Len: a.length, QRC: a.qrc, TRC: a.trc})
		return nil
	})
	if err != nil {
		return results, err
	}
	// No following ID exists to trigger the last group's finish.
	finish()
	return results, nil
}
