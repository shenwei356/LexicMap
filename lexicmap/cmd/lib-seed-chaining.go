package cmd

import (
	"cmp"
	"slices"
	"sync"
)

// trimSeedSearchResults retains all candidates tied with the Nth chaining score.
func trimSeedSearchResults(rs *[]*SearchResult, topN int) {
	if topN <= 0 || len(*rs) <= topN {
		return
	}
	slices.SortFunc(*rs, func(a, b *SearchResult) int {
		return cmp.Compare(b.Score, a.Score)
	})
	cutoff := (*rs)[topN-1].Score
	end := topN
	for end < len(*rs) && (*rs)[end].Score == cutoff {
		end++
	}
	// Release discarded candidates for GC without retaining them in object pools.
	clear((*rs)[end:])
	*rs = (*rs)[:end]
}

// chainSeedResults batches targets across fixed workers. Each worker owns its
// chainer and output slice; merging happens only after all workers finish.
func (idx *Index) chainSeedResults(seedResults [][]*SearchResult, nTargets int) *[]*SearchResult {
	rs := poolSearchResults.Get().(*[]*SearchResult)
	*rs = (*rs)[:0]
	if nTargets == 0 {
		return rs
	}

	nWorkers := min(idx.opt.NumCPUs, nTargets)
	// Keep enough jobs for load balancing, including small searches, while
	// amortizing channel operations for searches with millions of targets.
	batchSize := min(64, max(1, nTargets/(nWorkers*8)))
	jobs := make(chan []*SearchResult, nWorkers)
	results := make([]*[]*SearchResult, nWorkers)
	var wg sync.WaitGroup
	wg.Add(nWorkers)
	for worker := range nWorkers {
		go func() {
			defer wg.Done()
			local := poolSearchResults.Get().(*[]*SearchResult)
			*local = (*local)[:0]
			var chainer *Chainer
			for batch := range jobs {
				for _, r := range batch {
					if len(*r.Subs) > 1 {
						ClearSubstrPairs(poolSub, r.Subs, idx.k)
					}
					if chainer == nil {
						chainer = idx.poolChainers.Get().(*Chainer)
					}
					r.Chains, r.Score = chainer.Chain(r.Subs)
					// Preserve the existing policy of dropping oversized scratch
					// arrays rather than retaining them for subsequent targets.
					if len(chainer.visited) > chainerInitSize {
						chainer = nil
					}
					if r.Score < idx.chainingOptions.MinScore {
						idx.RecycleSearchResult(r)
						continue
					}
					r.prepareAlignmentRegions()
					*local = append(*local, r)
				}
			}
			if chainer != nil {
				RecycleChainer(idx.poolChainers, chainer)
			}
			results[worker] = local
		}()
	}
	for _, targets := range seedResults {
		for len(targets) > 0 {
			n := min(batchSize, len(targets))
			jobs <- targets[:n]
			targets = targets[n:]
		}
	}
	close(jobs)
	wg.Wait()

	nKept := 0
	for _, local := range results {
		nKept += len(*local)
	}
	*rs = slices.Grow(*rs, nKept)
	for _, local := range results {
		*rs = append(*rs, (*local)...)
		recycleSearchResultSlice(local)
	}
	return rs
}
