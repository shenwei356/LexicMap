package cmd

import (
	"fmt"
	"sync"
	"time"
	"unsafe"

	"github.com/dustin/go-humanize"
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
//     callback expands its borrowed seed data into a budgeted producer batch.
//     Submission applies genome/TaxId filters under the store lock and copies
//     accepted anchors into seedSpillStore, which handles sorting and spills.
//  3. Wait for every partition: a genome's anchors are incomplete before this.
//  4. Consume grouped anchors into budgeted batches of complete genomes for
//     parallel chaining, with a serial fallback for tiny shares or one worker.
//     Retain scores and chain bounds for the caller's global Top-N selection.
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
	if idx.opt.Verbose || idx.opt.Log2File {
		store.onFirstSpill = func() {
			log.Infof("%s (%s bp): seed spill started (%s buffered anchors, query buffer budget: %s)",
				query.seqID, humanize.Comma(int64(len(query.seq))), humanize.Comma(int64(len(store.anchors))),
				humanize.IBytes(uint64(store.budget.anchorsPerQuery)*uint64(seedAnchorBytes)))
		}
	}
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
	// Only submission mutates the store, the TaxId cache, and firstErr. Each
	// searcher expands its borrowed decoded values outside this lock into a fixed
	// producer batch. Reserve all batches up front so parallel expansion cannot
	// quietly exceed the collection share, even while producers wait on I/O.
	var mu sync.Mutex
	var firstErr error
	// Streaming batches can split one k-mer match; values counts decoded reference
	// seed locations before genome/TaxId filtering. store.count counts accepted
	// anchors after filtering, before deduplication/chaining.
	var batches, values uint64
	var nSearchersFinished int // protected by mu; used only for debug progress
	var filter *map[uint64]bool
	if idx.filterByTaxId {
		filter = idx.poolTaxIDfilter.Get().(*map[uint64]bool)
		defer func() { clear(*filter); idx.poolTaxIDfilter.Put(filter) }()
	}
	count := len(idx.Searchers)
	if idx.opt.InMemorySearch {
		count = len(idx.InMemorySearchers)
	}
	batchSize := min(seedAnchorBatchSize, store.budget.anchorsPerQuery/max(1, count)/4)
	store.transferReserve = count * batchSize
	store.recordBufferPeak(0)
	// Called with mu held. Filtering keeps arrival order within each accepted
	// batch. The TaxId cache remains shared, rather than duplicated per searcher.
	submit := func(anchors []seedAnchor) error {
		if firstErr != nil {
			return firstErr
		}
		if genomeIDs != nil || filter != nil {
			kept := anchors[:0]
			for _, a := range anchors {
				if genomeIDs != nil {
					if _, ok := (*genomeIDs)[a.batchGenomeIndex]; !ok {
						continue
					}
				}
				if filter != nil && !idx.keepGenomeByTaxID(filter, a.batchGenomeIndex) {
					continue
				}
				kept = append(kept, a)
			}
			anchors = kept
		}
		firstErr = store.addBatch(anchors)
		return firstErr
	}
	started := time.Now()
	var wg sync.WaitGroup
	for i := 0; i < count; i++ {
		wg.Add(1)
		go func(i int) {
			defer wg.Done()
			var localBatches, localValues uint64
			defer func() {
				mu.Lock()
				batches += localBatches
				values += localValues
				mu.Unlock()
			}()
			buffer := make([]seedAnchor, 0, batchSize)
			flush := func() error {
				if len(buffer) == 0 {
					return nil
				}
				mu.Lock()
				err := submit(buffer)
				mu.Unlock()
				buffer = buffer[:0]
				return err
			}
			consume := func(sr kv.SearchResult) error {
				localBatches++
				localValues += uint64(len(sr.Values))
				queryLocs := (*locses)[sr.IQuery]
				if sr.IsSuffix {
					queryLocs = (*locses)[(*reverseLocses)[sr.IQuery][sr.IQuery2]]
				}
				if batchSize == 0 {
					// Tiny budgets cannot reserve a transfer array per producer.
					// Expand under the lock directly, preserving the old tiny-share path.
					mu.Lock()
					defer mu.Unlock()
					return expandSeedDataBatch(idx.k, sr, queryLocs, func(a seedAnchor) error {
						return submit([]seedAnchor{a})
					})
				}
				return expandSeedDataBatch(idx.k, sr, queryLocs, func(a seedAnchor) error {
					buffer = append(buffer, a)
					if len(buffer) == cap(buffer) {
						return flush()
					}
					return nil
				})
			}
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
			// Streaming callbacks also expand/submit anchors and can wait on
			// spill I/O. Label their inclusive time separately from the original
			// non-streaming path's reading/decoding-only measurement.
			var searchStarted time.Time
			if debug {
				searchStarted = time.Now()
			}
			e := searcher.SearchStream((*kmers)[begin:end], idx.opt.MinPrefix, true, false, consume)
			if e == nil {
				e = searcher.Search2Stream((*reverseKmers)[begin:end], idx.opt.MinPrefix, true, true, consume)
			}
			if e == nil {
				e = flush()
			}
			if debug && e == nil {
				searchDuration := time.Since(searchStarted)
				mu.Lock()
				nSearchersFinished++
				log.Debugf("%s (%s bp): seed searcher finished %d/%d (masks [%d, %d)): streaming search/collection took %s; seed data batches=%s, matched k-mer locations=%s",
					query.seqID, humanize.Comma(int64(len(query.seq))), nSearchersFinished, count,
					begin, end, searchDuration, humanize.Comma(int64(localBatches)), humanize.Comma(int64(localValues)))
				mu.Unlock()
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
	// Producer callbacks and their transfer arrays are now unreachable.
	store.transferReserve = 0
	if firstErr != nil {
		return nil, fmt.Errorf("collect seed anchors: %w", firstErr)
	}
	if debug {
		log.Debugf("%s: bounded seed collection: seed data batches=%s, matched k-mer locations=%s, anchors after filtering=%s, peak collection buffer=%s bytes, spilled=%t, runs=%s, elapsed=%s",
			query.seqID, humanize.Comma(int64(batches)), humanize.Comma(int64(values)), humanize.Comma(int64(store.count)),
			humanize.Comma(store.peakBufferBytes), store.spilled, humanize.Comma(int64(len(store.runs))), time.Since(started))
	}
	started = time.Now()
	// Step 4: reconstruct/chain complete genomes; retain only scores and bounds.
	results, err = idx.chainSeedSpillStore(store, arena)
	if debug && err == nil {
		log.Debugf("%s: finished chaining collected seed anchors (%s genome hits, peak accounted anchor buffers=%s bytes) in %s",
			query.seqID, humanize.Comma(int64(len(*results))), humanize.Comma(store.peakBufferBytes), time.Since(started))
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

// seedReplayBatchBytes caps the charged flat replay buffer independently of a
// large collection share. Genome descriptors/candidate pages and worker chaining
// scratch remain outside the anchor-buffer budget, as in ordinary search.
const seedReplayBatchBytes = 8 << 20

// seedReplayBatchGenomes limits temporary candidate metadata and amortizes worker
// startup. A batch can finish earlier when its flat anchor buffer fills.
const seedReplayBatchGenomes = 4096

// chainSeedSpillStore groups complete genomes into bounded, pointer-free replay
// batches, then uses the ordinary parallel chaining implementation. Tiny shares
// and single-thread searches retain one-genome replay. Every genome's anchors
// stay in their original order, and copied candidates remain in genome order;
// global Top-N selection is still performed only by Index.Search afterward.
func (idx *Index) chainSeedSpillStore(store *seedSpillStore, arena *seedSearchResultArena) (results *[]*seedSearchResult, err error) {
	if idx.opt.NumCPUs <= 1 || int64(store.budget.anchorsPerQuery)*seedAnchorBytes < 128<<10 {
		return idx.chainSeedSpillStoreSerial(store, arena)
	}
	results = poolSeedSearchResults.Get().(*[]*seedSearchResult)
	*results = (*results)[:0]
	var temporary seedSearchResultArena
	var batch []*seedSearchResult
	var buffer, subs []SubstrPair
	var id uint64
	var begin int
	var haveGenome, oversized bool
	defer func() {
		// Candidate pages must not keep aliases of the flat replay buffer when
		// returned to their pool, including after a read fails mid-genome.
		for _, r := range batch {
			r.Subs = nil
		}
		temporary.recycle()
		store.replayBufferBytes = 0
		if err != nil {
			recycleSeedSearchResultSlice(results)
			results = nil
		}
	}()
	// process chains only the completed prefix. An incomplete final genome,
	// if any, is moved to the start afterward, once every worker has joined.
	// This needs no second replay array and never hands a partial genome to a worker.
	process := func() {
		if len(batch) == 0 {
			return
		}
		kept := idx.chainSeedResults([][]*seedSearchResult{batch}, len(batch))
		// Walk the original genome-ordered descriptors, not worker completion
		// order. Copy bounds before recycling any of the temporary candidate pages.
		for _, r := range batch {
			if r.Score >= idx.chainingOptions.MinScore {
				target := arena.add(r.BatchGenomeIndex)
				target.Score = r.Score
				target.chainRegions = append(target.chainRegions, r.chainRegions...)
				*results = append(*results, target)
			}
			r.Subs = nil
		}
		recycleSeedSearchResultSlice(kept)
		temporary.recycle()
		clear(batch)
		batch = batch[:0]
		if oversized {
			// The oversized genome is processed alone. Its array is never kept
			// for the next genome; return to the charged, fixed-capacity buffer.
			subs = buffer[:0]
			oversized = false
		} else {
			n := copy(subs, subs[begin:])
			subs = subs[:n]
		}
		begin = 0
	}
	finishGenome := func() {
		r := temporary.add(id)
		// Full slicing prevents any worker from growing into the following
		// genome's anchors; workers sort/deduplicate only their own segment.
		r.Subs = subs[begin:len(subs):len(subs)]
		batch = append(batch, r)
		begin = len(subs)
		if len(batch) == seedReplayBatchGenomes || oversized {
			process()
		}
	}
	err = store.consume(func(a seedAnchor) error {
		if buffer == nil {
			// consume has already released radix scratch and, on the spill path,
			// the collection array. Charge replay only to the remaining headroom.
			available := int64(store.budget.anchorsPerQuery-cap(store.anchors)-cap(store.sortScratch)-store.transferReserve) * seedAnchorBytes
			pairBytes := int64(unsafe.Sizeof(SubstrPair{}))
			capacity := int(min(int64(seedReplayBatchBytes), available) / pairBytes)
			// Collection already knows the total anchor count. Small queries
			// cannot use more than this, so avoid allocating 8 MiB for each one.
			capacity = int(min(uint64(capacity), store.count))
			buffer = make([]SubstrPair, 0, capacity)
			subs = buffer
			store.replayBufferBytes = int64(capacity) * pairBytes
			store.recordBufferPeak(0)
		}
		if haveGenome && id != a.batchGenomeIndex {
			finishGenome()
		}
		id, haveGenome = a.batchGenomeIndex, true
		if len(subs) == cap(subs) {
			process()
			if len(subs) == cap(subs) {
				// Only the incomplete genome remains and it exceeds the batch
				// buffer. Preserve the existing one-genome budget exclusion.
				oversized = true
			}
		}
		subs = append(subs, SubstrPair{QBegin: a.qBegin, TBegin: a.tBegin, Len: a.length, QRC: a.qrc, TRC: a.trc})
		return nil
	})
	if err != nil {
		return results, err
	}
	if haveGenome {
		finishGenome()
		process()
	}
	return results, nil
}

// chainSeedSpillStoreSerial accumulates one complete genome and uses the same
// ClearSubstrPairs/Chainer logic. It is the fallback for tiny budgets and one
// chaining thread, and preserves the existing unbudgeted one-genome replay.
// The serial path consumes genome-grouped anchors after collection finishes.
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
func (idx *Index) chainSeedSpillStoreSerial(store *seedSpillStore, arena *seedSearchResultArena) (results *[]*seedSearchResult, err error) {
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
