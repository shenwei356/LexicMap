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
	"cmp"
	"slices"
	"sync"
)

// trimSeedSearchResults retains all candidates tied with the Nth chaining score.
func trimSeedSearchResults(rs *[]*seedSearchResult, topN int) {
	if topN <= 0 || len(*rs) <= topN {
		return
	}
	// A small Top-N needs only a cutoff score from the discarded majority.
	// Keep full sorting for small inputs or large limits, where heap setup and
	// an extra scan can cost more than the ordering work they avoid.
	if len(*rs) >= 4096 && topN < len(*rs)/8 {
		if cutoff, ok := seedScoreCutoff(*rs, topN); ok {
			kept := 0
			for _, r := range *rs {
				if r.Score >= cutoff { // include every tie at the global cutoff
					(*rs)[kept] = r
					kept++
				}
			}
			clear((*rs)[kept:])
			*rs = (*rs)[:kept]
			// Preserve the ranked-score contract for callers and diagnostics.
			slices.SortFunc(*rs, func(a, b *seedSearchResult) int {
				return cmp.Compare(b.Score, a.Score)
			})
			return
		}
	}
	slices.SortFunc(*rs, func(a, b *seedSearchResult) int {
		return cmp.Compare(b.Score, a.Score)
	})
	cutoff := (*rs)[topN-1].Score
	end := topN
	for end < len(*rs) && (*rs)[end].Score == cutoff {
		end++
	}
	// Remove discarded entries from the ranked list; the query owns their pages.
	clear((*rs)[end:])
	*rs = (*rs)[:end]
}

// seedScoreCutoff maintains the largest n scores in a min-heap. NaNs select
// the original sorting path, preserving its comparison and equality semantics.
func seedScoreCutoff(results []*seedSearchResult, n int) (float32, bool) {
	heap := make([]float32, n)
	for i, r := range results[:n] {
		if r.Score != r.Score {
			return 0, false
		}
		heap[i] = r.Score
	}
	slices.Sort(heap) // an ascending array is already a valid min-heap
	for _, r := range results[n:] {
		score := r.Score
		if score != score {
			return 0, false
		}
		if score <= heap[0] {
			continue
		}
		// Sift the replacement down without interface boxing or heap objects.
		i := 0
		for child := 1; child < n; child = i*2 + 1 {
			if child+1 < n && heap[child+1] < heap[child] {
				child++
			}
			if score <= heap[child] {
				break
			}
			heap[i] = heap[child]
			i = child
		}
		heap[i] = score
	}
	return heap[0], true
}

// chainSeedResults batches targets across fixed workers. Each worker owns its
// chainer and output slice; merging happens only after all workers finish.
func (idx *Index) chainSeedResults(seedResults [][]*seedSearchResult, nTargets int) *[]*seedSearchResult {
	rs := poolSeedSearchResults.Get().(*[]*seedSearchResult)
	*rs = (*rs)[:0]
	if nTargets == 0 {
		return rs
	}

	nWorkers := min(idx.opt.NumCPUs, nTargets)
	// Keep enough jobs for load balancing, including small searches, while
	// amortizing channel operations for searches with millions of targets.
	batchSize := min(64, max(1, nTargets/(nWorkers*8)))
	jobs := make(chan []*seedSearchResult, nWorkers)
	results := make([]*[]*seedSearchResult, nWorkers)
	var wg sync.WaitGroup
	wg.Add(nWorkers)
	for worker := range nWorkers {
		go func() {
			defer wg.Done()
			local := poolSeedSearchResults.Get().(*[]*seedSearchResult)
			*local = (*local)[:0]
			var chainer *Chainer
			for batch := range jobs {
				for _, r := range batch {
					if len(r.Subs) > 1 {
						ClearSubstrPairs(&r.Subs, idx.k)
					}
					if chainer == nil {
						chainer = idx.poolChainers.Get().(*Chainer)
					}
					r.resetRegions()
					r.chainRegions, r.Score = chainer.Chain(&r.Subs, r.chainRegions)
					r.resetSeeds()
					// Preserve the existing policy of dropping oversized scratch
					// arrays rather than retaining them for subsequent targets.
					if len(chainer.visited) > chainerInitSize {
						chainer = nil
					}
					if r.Score < idx.chainingOptions.MinScore {
						continue
					}
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
		recycleSeedSearchResultSlice(local)
	}
	return rs
}
