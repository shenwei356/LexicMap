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
	"cmp"
	"fmt"
	"os"
	"regexp"
	"slices"
	"sync"
	"time"

	"github.com/dustin/go-humanize"
	itree "github.com/rdleal/intervalst/interval"
	"github.com/vbauerster/mpb/v8"
	"github.com/vbauerster/mpb/v8/decor"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
	"github.com/shenwei356/lexichash/iterator"
	"github.com/shenwei356/wfa"
)


// In-memory sketch of a subject genome. We only need it long enough to
// process all query fragments against this subject, then it is recycled.
// The seed layout mirrors what buildAnIndex writes to disk: primary k-mer
// per mask plus extra k-mers from desert filling, stored against the same
// per-mask slot. Positions are encoded as (pos << 1) | strand, matching
// the convention used by MaskKnownDistinctPrefixes.
// subjectSketch is an in-memory sketch of a subject genome for sampled k-mer matching.
// It contains the concatenated sequence (forward + RC) and a k-mer map for fast lookup.
type subjectSketch struct {
	seqLen         int
	sampledKmerMap *map[uint64][]uint32 // k-mer -> positions
	contigBounds   [][2]int             // [start, end) of each contig in forward strand
}

var poolSubjectSketch = &sync.Pool{New: func() interface{} {
	return &subjectSketch{}
}}

// poolConcat is for reusing large byte slices for concatenated genome sequences
var poolConcat = &sync.Pool{New: func() interface{} {
	tmp := make([]byte, 0, 20<<20) // 10MB initial capacity
	return &tmp
}}

// poolKmerMap is for reusing k-mer maps in sampled mode
var poolKmerMap = &sync.Pool{New: func() interface{} {
	m := make(map[uint64][]uint32, 10240) // pre-allocate for ~10k k-mers
	return &m
}}

// poolKmerSlice is for reusing position slices to reduce allocations
var poolKmerSlice = &sync.Pool{New: func() interface{} {
	s := make([]uint32, 0, 8) // pre-allocate for typical k-mer frequency
	return &s
}}

// poolQSeeds is for reusing query seed slices
var poolQSeeds = &sync.Pool{New: func() interface{} {
	s := make([]*[]uint64, 0, 10240)
	return &s
}}

// kmerLookup stores the result of a k-mer map lookup for batch processing
type kmerLookup struct {
	qk    uint64
	qloc  int
	slocs []uint32
}

// poolKmerLookups is for reusing kmerLookup slices in batch processing
var poolKmerLookups = &sync.Pool{New: func() interface{} {
	s := make([]kmerLookup, 0, 64) // pre-allocate for batch processing
	return &s
}}

// Sampling parameters for the simplified seeding strategy.
var gsa3SampledK = 13     // fixed k-mer length for sampling
var gsa3SamplingScale = 4 // sampling rate: keep if hash(kmer) % scale == 0

// buildSubjectSketchSampledOptimized is the optimized version that accepts forwardLen and rcStart.
// When forwardLen and rcStart are provided (> 0), it only processes the forward part and
// mirrors k-mers to the RC part, avoiding redundant scanning (approximately 50% speedup).
func (idx *Index) buildSubjectSketchSampledOptimized(seq []byte, skipRegions [][2]int, contigBounds [][2]int, forwardLen int, rcStart int) (*subjectSketch, error) {
	k := gsa3SampledK
	k8 := uint8(k)
	scale := uint64(gsa3SamplingScale)
	scaleM1 := scale - 1

	if len(seq) < k {
		return nil, fmt.Errorf("sequence too short for k=%d", k)
	}

	// Get k-mer map from pool and clear it
	kmerMap := poolKmerMap.Get().(*map[uint64][]uint32)

	// Build interval tree for skip regions
	tree := itree.NewSearchTree[uint8, int](cmpFn)
	for _, r := range skipRegions {
		tree.Insert(r[0]-k+1, r[1], 1)
	}

	// Low-complexity k-mer values
	ccc := util.Ns(0b01, k8)
	ggg := util.Ns(0b10, k8)
	ttt := (uint64(1) << (k << 1)) - 1

	if forwardLen > 0 && rcStart > 0 && rcStart < len(seq) {
		// Optimized path: only process forward part, then mirror to RC part
		iter, err := iterator.NewKmerIterator(seq[:forwardLen], k)
		if err != nil {
			return nil, err
		}

		pos := 0
		for {
			kmer, kmerRC, ok, _ := iter.NextKmer()
			if !ok {
				break
			}

			// Skip if in a skip region
			if _, inGap := tree.AnyIntersection(pos, pos); inGap {
				pos++
				continue
			}

			// Use canonical k-mer for sampling decision
			canonical := kmer
			if kmerRC < kmer {
				canonical = kmerRC
			}

			// Sample using hash modulo
			if util.Hash64(canonical)&scaleM1 != 0 {
				pos++
				continue
			}

			// Skip low-complexity k-mers
			if kmer == ccc || kmer == ggg || kmer == ttt || util.IsLowComplexityDust(kmer, k8) {
				pos++
				continue
			}

			// Store k-mer in forward part
			if _, exists := (*kmerMap)[kmer]; !exists {
				// Pre-allocate with reasonable capacity to reduce reallocations
				(*kmerMap)[kmer] = make([]uint32, 0, 4)
			}
			(*kmerMap)[kmer] = append((*kmerMap)[kmer], uint32(pos))

			// Mirror to RC part: the RC position is rcStart + (forwardLen - pos - k)
			rcPos := rcStart + (forwardLen - pos - k)
			if _, exists := (*kmerMap)[kmerRC]; !exists {
				(*kmerMap)[kmerRC] = make([]uint32, 0, 4)
			}
			(*kmerMap)[kmerRC] = append((*kmerMap)[kmerRC], uint32(rcPos))

			pos++
		}
	} else {
		// Fallback: process entire sequence (original behavior)
		iter, err := iterator.NewKmerIterator(seq, k)
		if err != nil {
			return nil, err
		}

		pos := 0
		for {
			kmer, kmerRC, ok, _ := iter.NextKmer()
			if !ok {
				break
			}

			// Skip if in a skip region
			if _, inGap := tree.AnyIntersection(pos, pos); inGap {
				pos++
				continue
			}

			// Use canonical k-mer for sampling decision
			canonical := kmer
			if kmerRC < kmer {
				canonical = kmerRC
			}

			// Sample using hash modulo
			if util.Hash64(canonical)&scaleM1 != 0 {
				pos++
				continue
			}

			// Skip low-complexity k-mers
			if kmer == ccc || kmer == ggg || kmer == ttt || util.IsLowComplexityDust(kmer, k8) {
				pos++
				continue
			}

			// Store both forward and reverse k-mers
			(*kmerMap)[kmer] = append((*kmerMap)[kmer], uint32(pos<<1))
			(*kmerMap)[kmerRC] = append((*kmerMap)[kmerRC], uint32(pos<<1|1))

			pos++
		}
	}

	s := poolSubjectSketch.Get().(*subjectSketch)
	s.seqLen = len(seq)
	s.sampledKmerMap = kmerMap
	s.contigBounds = contigBounds

	return s, nil
}

// recycleSubjectSketch returns all per-sketch buffers back to their pools.
func (idx *Index) recycleSubjectSketch(s *subjectSketch) {
	if s == nil {
		return
	}
	if s.sampledKmerMap != nil {
		clear(*s.sampledKmerMap)
		poolKmerMap.Put(s.sampledKmerMap)
		s.sampledKmerMap = nil
	}
	poolSubjectSketch.Put(s)
}

// qFragSeedsSampled holds sampled k-mers for a query fragment in the simplified seeding mode.
// sampleQueryFragment samples fixed-length k-mers from a query fragment.
// Modified to only extract forward strand k-mers.
func sampleQueryFragment(frag []byte) (*[]uint64, error) {
	k := gsa3SampledK
	k8 := uint8(k)
	scale := uint64(gsa3SamplingScale)
	scaleM1 := scale - 1

	if len(frag) < k {
		empty := []uint64{}
		return &empty, nil
	}

	ccc := util.Ns(0b01, k8)
	ggg := util.Ns(0b10, k8)
	ttt := (uint64(1) << (k << 1)) - 1

	sampledKmers := poolKmerAndLocs.Get().(*[]uint64)
	*sampledKmers = (*sampledKmers)[:0]

	iter, err := iterator.NewKmerIterator(frag, k)
	if err != nil {
		return nil, err
	}

	pos := 0
	for {
		kmer, kmerRC, ok, _ := iter.NextKmer()
		if !ok {
			break
		}

		// Use canonical k-mer for sampling decision
		canonical := kmer
		if kmerRC < kmer {
			canonical = kmerRC
		}

		// Sample using hash modulo (fast bitwise AND since scale is power of 2)
		if util.Hash64(canonical)&scaleM1 != 0 {
			pos++
			continue
		}

		// Skip low-complexity k-mers
		if kmer == ccc || kmer == ggg || kmer == ttt || util.IsLowComplexityDust(kmer, k8) {
			pos++
			continue
		}

		// Store only forward strand k-mers
		*sampledKmers = append(*sampledKmers, kmer, uint64(pos))

		pos++
	}

	return sampledKmers, nil
}

// addAnchor emits one SubstrPair from a (qLoc, sLoc) match. The two locs are
// (pos << 1) | strand_bit; relative strand tells us whether the alignment
// is forward (both strands the same) or reverse (different strands). For
// reverse anchors we re-express the subject position in RC-subject
// coordinates so the chainer can treat both groups uniformly.
// alignQueryFragToSubjectSampled matches a query fragment against a sampled subject sketch.
// The sKmerMap should be pre-built once for all fragments to avoid repeated construction.
// Modified version: query only has forward strand k-mers, subject is forward + RC concatenated.
func alignQueryFragToSubjectSampled(
	qfrag []byte,
	qSeeds *[]uint64,
	sketch *subjectSketch,
	sKmerMap *map[uint64][]uint32,
	concat []byte,
	chainer *Chainer2,
	algn *wfa.Aligner,
	K int,
	extLen int,
	extLen2 int,
	minPIdent float64,
	minQcov float64,
	idx *Index,
	fScoreAndEvalue *func(qlen int, cigar *wfa.AlignmentResult) (int, int, float64),
) (int, int, int, float64, bool) {
	// Since we only use forward strand query k-mers and subject is a single concatenated
	// sequence (forward + RC), we only need one set of anchors for unified chaining.
	allSubs := poolSubsLong.Get().(*[]*SubstrPair)
	*allSubs = (*allSubs)[:0]
	defer RecycleSubstrPairs(poolSub, poolSubsLong, allSubs)

	if sKmerMap == nil || len(*sKmerMap) == 0 {
		return 0, 0, 0, 0, false
	}

	qKmers := *qSeeds

	// Match query k-mers against subject k-mers
	// Limit matches per k-mer to avoid excessive anchors from repetitive sequences
	// Use batching to improve cache locality and reduce map lookup overhead
	const maxMatchesPerKmer = 100
	const batchSize = 8 // Process multiple k-mers in a batch

	// Get lookups slice from pool and reuse it
	lookups := poolKmerLookups.Get().(*[]kmerLookup)
	defer func() {
		*lookups = (*lookups)[:0]
		poolKmerLookups.Put(lookups)
	}()

	// Pre-fetch map entries in batches to improve cache hit rate
	for batchStart := 0; batchStart+1 < len(qKmers); batchStart += batchSize * 2 {
		batchEnd := batchStart + batchSize*2
		if batchEnd > len(qKmers) {
			batchEnd = len(qKmers)
		}

		// First pass: lookup all k-mers in the batch to warm up cache
		*lookups = (*lookups)[:0]

		for i := batchStart; i+1 < batchEnd; i += 2 {
			qk := qKmers[i]
			qloc := int(qKmers[i+1])

			if slocs, found := (*sKmerMap)[qk]; found {
				*lookups = append(*lookups, kmerLookup{qk, qloc, slocs})
			}
		}

		// Second pass: process all matches in the batch
		for _, lookup := range *lookups {
			maxN := len(lookup.slocs)
			if maxN > maxMatchesPerKmer {
				maxN = maxMatchesPerKmer
			}

			for j := 0; j < maxN; j++ {
				sloc := int(lookup.slocs[j])
				qpos := lookup.qloc
				spos := sloc

				// Create anchor directly without strand consideration
				sub := poolSub.Get().(*SubstrPair)
				sub.Len = uint8(K)
				sub.QBegin = int32(qpos)
				sub.TBegin = int32(spos)
				sub.QRC = false
				sub.TRC = false
				*allSubs = append(*allSubs, sub)
			}
		}
	}

	chains, chainsOk := chainsFromSubs(allSubs, chainer, K)
	if !chainsOk {
		return 0, 0, 0, 0, false
	}

	// Try all chains and pick the best one
	var bestMatched, bestAligned, bestGaps int
	var bestPident float64
	var bestScore int = -1
	topChains := idx.chainingOptions.TopChains
	onlyTopChains := topChains > 0

	// Pre-index qfrag once for all chains
	cpr := idx.poolSeqComparator.Get().(*SeqComparator)
	defer idx.poolSeqComparator.Put(cpr)
	if err := cpr.Index(qfrag); err != nil {
		return 0, 0, 0, 0, false
	}
	defer cpr.RecycleIndex()

	i := 0
	for _, chain := range *chains {
		if chain == nil {
			continue
		}
		i++
		if onlyTopChains && i > topChains {
			break
		}
		matched, aligned, gaps, pident, ok := alignChain(
			qfrag, concat, chain, sketch, false, algn, cpr,
			extLen, extLen2, minPIdent, minQcov, idx,
			fScoreAndEvalue,
		)

		if ok {
			score := matched * aligned
			if score > bestScore {
				bestScore = score
				bestMatched = matched
				bestAligned = aligned
				bestGaps = gaps
				bestPident = pident
			}
		}
	}
	RecycleChaining2Result(chains)

	if bestScore <= 0 {
		return 0, 0, 0, 0, false
	}

	return bestMatched, bestAligned, bestGaps, bestPident, true
}

// chainsFromSubs runs the chaining pipeline and returns all chains.
func chainsFromSubs(subs *[]*SubstrPair, chainer *Chainer2, K int) (*[]*Chain2Result, bool) {
	if len(*subs) == 0 {
		return nil, false
	}

	if len(*subs) > 1 {
		ClearSubstrPairs(poolSub, subs, K)
	}
	TrimSubStrPairs(poolSub, subs, K, 100)
	if len(*subs) == 0 {
		return nil, false
	}

	chains, _, _, _, _, _, _, _ := chainer.Chain(subs)
	if chains == nil || len(*chains) == 0 {
		if chains != nil {
			RecycleChaining2Result(chains)
		}
		return nil, false
	}

	return chains, true
}

// alignChain performs SeqComparator pseudo-alignment and WFA on a single chain.
// Returns (matched, aligned, gaps, true) on success.
func alignChain(
	qfrag []byte,
	subjectSeq []byte,
	chain *Chain2Result,
	sketch *subjectSketch,
	useRC bool,
	algn *wfa.Aligner,
	cpr *SeqComparator, // pre-indexed comparator
	extLen int,
	extLen2 int,
	minPIdent float64,
	minQcov float64,
	idx *Index,
	fScoreAndEvalue *func(qlen int, cigar *wfa.AlignmentResult) (int, int, float64),
) (int, int, int, float64, bool) {
	// Guard against degenerate chains.
	if chain.QEnd < chain.QBegin || chain.TEnd < chain.TBegin {
		return 0, 0, 0, 0, false
	}

	qLen := len(qfrag)
	subjectLen := sketch.seqLen

	// Locate the contig containing the chain.
	var contigStart, contigEnd int
	contigStart, contigEnd = 0, subjectLen
	if bounds := sketch.contigBounds; len(bounds) > 0 {
		if useRC {
			for i := len(bounds) - 1; i >= 0; i-- {
				cs := subjectLen - bounds[i][1]
				ce := subjectLen - bounds[i][0]
				if chain.TBegin >= cs && chain.TBegin < ce {
					contigStart, contigEnd = cs, ce
					break
				}
			}
		} else {
			for _, b := range bounds {
				if chain.TBegin >= b[0] && chain.TBegin < b[1] {
					contigStart, contigEnd = b[0], b[1]
					break
				}
			}
		}
	}

	// Expand the chain region by extLen.
	tExpBegin := max(chain.TBegin-extLen, contigStart)
	tExpEnd := min(chain.TEnd+extLen, contigEnd-1)
	tSubseq := subjectSeq[tExpBegin : tExpEnd+1]

	qExpBegin := max(chain.QBegin-extLen, 0)
	qExpEnd := min(chain.QEnd+extLen, qLen-1)

	// Pseudo-alignment with SeqComparator (already indexed).
	cr, err := cpr.Compare(uint32(qExpBegin), uint32(qExpEnd), tSubseq, qLen)
	if err != nil {
		return 0, 0, 0, 0, false
	}
	if cr == nil {
		return 0, 0, 0, 0, false
	}
	defer RecycleSeqComparatorResult(cr)

	// WFA alignment on each sub-chain.
	var totMatched, totAligned, totGaps int
	maxEvalue := idx.opt.MaxEvalue
	maxTrials := 2

	trials := 0
	for _, c := range *cr.Chains {
		if c.QEnd < c.QBegin || c.TEnd < c.TBegin {
			continue
		}

		trials++
		if trials > maxTrials { // can't find a valid alignment for the best 3 chains, give up
			break
		}

		cTBegin := c.TBegin
		cMaxExtLen := len(tSubseq) - 1 - c.TEnd

		_qseq, _tseq, _, _, _, _, extErr := extendMatch(
			qfrag, tSubseq,
			c.QBegin, c.QEnd+1,
			c.TBegin, c.TEnd+1,
			extLen2, cTBegin, cMaxExtLen, false,
		)
		if extErr != nil {
			continue
		}

		cigar, alignErr := algn.Align(_qseq, _tseq)
		if alignErr != nil {
			continue
		}

		// score and e-value
		_, _, evalue := (*fScoreAndEvalue)(len(_qseq), cigar)
		if evalue > maxEvalue {
			wfa.RecycleAlignmentResult(cigar)
			continue
		}

		totMatched += int(cigar.Matches)
		totAligned += int(cigar.AlignLen)
		totGaps += int(cigar.Gaps)
		wfa.RecycleAlignmentResult(cigar)

		break // keep the best ONE match
	}

	if totAligned <= 0 {
		return 0, 0, 0, 0, false
	}

	pident := float64(totMatched) / float64(totAligned) * 100
	alignedBasesQ := totAligned - totGaps
	af := float64(alignedBasesQ) / float64(qLen) * 100
	if af > 100 {
		af = 100
	}
	if pident < minPIdent || af < minQcov {
		return 0, 0, 0, 0, false
	}

	return totMatched, totAligned, totGaps, pident, true
}

// GSearchAlign3Sampled is a simplified version of GSearchAlign3 that uses
// sampled fixed-length k-mers instead of LexicHash masking.
func (idx *Index) GSearchAlign3Sampled(query *GQuery, fragLen int, minFragLen int, genomeIds *map[uint64]*[]uint64, minAF float64, maxQueryConcurrency int, gcInterval uint64) error {
	debug := idx.opt.Debug

	startTime0 := time.Now()

	if debug {
		log.Debugf("%s (%s bp): start to preprocess query genome fragments", query.id, humanize.Comma(int64(query.genomeSize)))
	}

	// 1) Cut the query into fragments.
	qfrags, qfragLens := seqs2fragments(&query.seqs, fragLen, minFragLen)
	if len(*qfrags) == 0 {
		return fmt.Errorf("no fragments for alignment, are the genome too fragmented with all sequences shorter than the minimum fragment length (%d bp)?", minFragLen)
	}

	// 2) Sample k-mers from each query fragment.
	qSeeds := poolQSeeds.Get().(*[]*[]uint64)
	*qSeeds = (*qSeeds)[:0]
	if cap(*qSeeds) < len(*qfrags) {
		*qSeeds = make([]*[]uint64, 0, len(*qfrags))
	}
	defer func() {
		for _, seeds := range *qSeeds {
			poolKmerAndLocs.Put(seeds)
		}
		*qSeeds = (*qSeeds)[:0]
		poolQSeeds.Put(qSeeds)
	}()

	for _, qfrag := range *qfrags {
		seeds, err := sampleQueryFragment(qfrag)
		if err != nil {
			return fmt.Errorf("failed to sample query fragment: %w", err)
		}
		*qSeeds = append(*qSeeds, seeds)
	}

	if debug {
		log.Debugf("%s (%s bp): finished preprocessing query genome fragments in %.3f seconds",
			query.id, humanize.Comma(int64(query.genomeSize)), time.Since(startTime0).Seconds())
		log.Debugf("%s (%s bp): start to align query genome fragments", query.id, humanize.Comma(int64(query.genomeSize)))
	}

	startTime := time.Now()

	// 3) Prepare result channel and collector.

	toDelete := make([]uint64, 0, len(*genomeIds))
	for id, ids := range *genomeIds {
		if id != (*ids)[0] {
			toDelete = append(toDelete, id)
		}
	}
	for _, id := range toDelete {
		delete(*genomeIds, id)
	}

	// -----------------------------------------------------------
	// process bar
	var pbs *mpb.Progress
	var bar *mpb.Bar
	var chDuration chan time.Duration
	var doneDuration chan int
	if debug {
		pbs = mpb.New(mpb.WithWidth(40), mpb.WithOutput(os.Stderr))
		bar = pbs.AddBar(int64(len(*genomeIds)),
			mpb.PrependDecorators(
				decor.Name("checked subject genomes: ", decor.WC{W: len("checked subject genomes: "), C: decor.DindentRight}),
				decor.Name("", decor.WCSyncSpaceR),
				decor.CountersNoUnit("%d / %d", decor.WCSyncWidth),
			),
			mpb.AppendDecorators(
				decor.Name("ETA: ", decor.WC{W: len("ETA: ")}),
				decor.EwmaETA(decor.ET_STYLE_GO, 1024),
				decor.OnComplete(decor.Name(""), ". done"),
			),
		)

		chDuration = make(chan time.Duration, idx.opt.NumCPUs)
		doneDuration = make(chan int)
		go func() {
			for t := range chDuration {
				bar.EwmaIncrBy(1, t)
			}
			doneDuration <- 1
		}()
	}

	fcpus := float64(idx.opt.NumCPUs)

	// -----------------------------------------------------------

	ch := make(chan *GSearchResult, 2*idx.opt.NumCPUs)
	done := make(chan int)
	go func() {
		rs := poolGSearchResults.Get().(*[]*GSearchResult)
		*rs = (*rs)[:0]
		for r := range ch {
			*rs = append(*rs, r)
		}
		slices.SortFunc(*rs, func(a, b *GSearchResult) int {
			if d := cmp.Compare(b.ANI, a.ANI); d != 0 {
				return d
			}
			if d := cmp.Compare(b.AFq, a.AFq); d != 0 {
				return d
			}
			return cmp.Compare(b.AFs, a.AFs)
		})

		query.result = rs
		done <- 1
	}()

	// 4) read genomes and align

	K := gsa3SampledK
	contigInterval := int(float64(fragLen) * 1.5)
	if contigInterval < K {
		contigInterval = K
	}
	nnn := bytes.Repeat([]byte{'N'}, contigInterval)
	reGaps := regexp.MustCompile(fmt.Sprintf(`[Nn]{%d,}`, 5))

	alignOption := &wfa.Options{GlobalAlignment: true}
	minPIdent := idx.seqCompareOption.MinIdentity
	minQcovHSP := idx.seqCompareOption.MinAlignedFraction
	extLen := fragLen / 2
	extLen2 := idx.opt.ExtendLength2

	var wg sync.WaitGroup
	tokens := make(chan int, maxQueryConcurrency)

	for _, batchIDAndRefIDs := range *genomeIds {
		tokens <- 1
		wg.Add(1)

		go func(batchIDAndRefIDs *[]uint64) {
			timeStart := time.Now()

			defer func() {
				wg.Done()
				<-tokens

				if debug {
					chDuration <- time.Duration(float64(time.Since(timeStart)) / fcpus)
				}
			}()

			var g *genome.Genome
			genomes := make([]*genome.Genome, len(*batchIDAndRefIDs))
			maxSubjectGenomeSize := idx.opt.MaxSubjectGenomeSize

			for i, batchIDAndRefID := range *batchIDAndRefIDs {
				genomeBatch := int(batchIDAndRefID >> BITS_GENOME_IDX)
				genomeIdx := int(batchIDAndRefID & MASK_GENOME_IDX)

				rdr := <-idx.poolGenomeRdrs[genomeBatch]

				_g, err := rdr.Seqs(genomeIdx)
				if err != nil {
					checkError(fmt.Errorf("fail to read genome sequence for batch %d, genome index %d: %s", genomeBatch, genomeIdx, err))
				}
				if i == 0 {
					g = _g
				} else {
					g.Seqs = append(g.Seqs, _g.Seqs...)
					_g.Seqs = nil
					g.NumSeqs += _g.NumSeqs
					g.GenomeSize += _g.GenomeSize
				}

				if maxSubjectGenomeSize > 0 && g.GenomeSize > maxSubjectGenomeSize {
					log.Warningf("%s (size: %s bp) exceeds the maximum subject genome size which exceeds the maximum allowed size of %s, consider increasing --max-subject-genome-size",
						idx.BatchGenomeIndex2GenomeID[(*batchIDAndRefIDs)[0]],
						humanize.Comma(int64(g.GenomeSize)),
						humanize.Comma(int64(maxSubjectGenomeSize)))

					idx.poolGenomeRdrs[genomeBatch] <- rdr
					for _, gx := range genomes {
						genome.RecycleGenome(gx)
						return
					}
					break
				}

				genomes[i] = _g
				idx.poolGenomeRdrs[genomeBatch] <- rdr
			}

			// b) Concatenate contigs with separators.
			concat := poolConcat.Get().(*[]byte)
			*concat = (*concat)[:0]

			// Calculate total size: forward + contig intervals + RC interval + RC
			var forwardSize int
			for _, s := range g.Seqs {
				forwardSize += len(*s)
			}
			forwardSize += contigInterval * (len(g.Seqs) - 1)

			// Total size = forward + 2*fragLen interval + RC (same as forward)
			rcInterval := fragLen << 1
			totalSize := forwardSize<<1 + rcInterval

			// Pre-allocate the full capacity to avoid reallocation
			if cap(*concat) < totalSize {
				*concat = make([]byte, 0, totalSize)
			}

			var skipRegions [][2]int
			contigBounds := make([][2]int, 0, len(g.Seqs)) // Only forward contigs needed
			for i, s := range g.Seqs {
				if i > 0 {
					boundary := len(*concat)
					skipRegions = append(skipRegions, [2]int{boundary, boundary + contigInterval - 1})
					*concat = append(*concat, nnn...)
				}
				cs := len(*concat)
				*concat = append(*concat, (*s)...)
				contigBounds = append(contigBounds, [2]int{cs, len(*concat)})
			}

			// skip gap regions (N's) in forward strand
			gaps := reGaps.FindAllSubmatchIndex(*concat, -1)
			for _, gap := range gaps {
				skipRegions = append(skipRegions, [2]int{gap[0], gap[1] - 1})
			}

			// c) Append 2*fragLen interval and reverse complement strand
			forwardLen := len(*concat)
			nnnRC := bytes.Repeat([]byte{'N'}, rcInterval)

			// Add interval between forward and RC strands
			*concat = append(*concat, nnnRC...)

			// Append reverse complement of the forward strand
			rcStart := len(*concat)
			// Directly append forward part to concat itself, then RC the newly appended portion
			*concat = append(*concat, (*concat)[:forwardLen]...)
			RC((*concat)[rcStart:])

			// Sort skip regions
			slices.SortFunc(skipRegions, func(a, b [2]int) int {
				return a[0] - b[0]
			})

			// d) Build the subject sketch using sampled k-mers on the combined sequence
			// Pass forwardLen and rcStart for optimized k-mer extraction
			sketch, err := idx.buildSubjectSketchSampledOptimized(*concat, skipRegions, contigBounds, forwardLen, rcStart)
			if err != nil {
				checkError(fmt.Errorf("fail to build subject sketch: %s", err))
			}

			// e) Set up per-subject scratch.
			chainer := idx.poolChainers2.Get().(*Chainer2)
			algn := wfa.New(wfa.DefaultPenalties, alignOption)
			algn.AdaptiveReduction(wfa.DefaultAdaptiveOption)

			gr := poolGSearchResult.Get().(*GSearchResult)
			gr.Reset()
			gr.BatchGenomeIndex = (*batchIDAndRefIDs)[0]
			gr.GenomeSize = g.GenomeSize
			gr.NumSeqs = g.NumSeqs

			// f) Align each query fragment using the pre-built k-mer map
			fScoreAndEvalue := scoreAndEvalue(2, -3, 5, 2, int(g.GenomeSize), 0.625, 0.41)

			for i, qfrag := range *qfrags {
				matched, alignedLen, gaps, pident, ok := alignQueryFragToSubjectSampled(
					qfrag, (*qSeeds)[i], sketch, sketch.sampledKmerMap, (*concat),
					chainer, algn, K, extLen, extLen2,
					minPIdent, minQcovHSP, idx,
					&fScoreAndEvalue,
				)
				if !ok {
					// fmt.Printf("fail to align fragment %d: %s\n", i+1, qfrag)
					continue
				}
				gr.AlignedFragments++
				gr.AlignedLength += alignedLen - gaps
				gr.AlignedMatches += matched
				gr.Pidents = append(gr.Pidents, pident)
			}

			// g) ANI / AF on the accumulated alignment.
			if gr.AlignedLength > 0 {
				sumPident := 0.0
				for _, p := range gr.Pidents {
					sumPident += p
				}
				gr.ANI = sumPident / float64(len(gr.Pidents)) / 100
			}
			gr.AFq = float64(gr.AlignedLength) / float64(qfragLens)
			gr.AFs = float64(gr.AlignedLength) / float64(gr.GenomeSize)
			if gr.AFq > 1 {
				gr.AFq = 1
			}
			if gr.AFs > 1 {
				gr.AFs = 1
			}
			gr.Score = gr.ANI

			if gr.AFq < minAF {
				poolGSearchResult.Put(gr)
			} else {
				ch <- gr
			}

			// h) Cleanup.
			wfa.RecycleAligner(algn)
			idx.poolChainers2.Put(chainer)
			idx.recycleSubjectSketch(sketch)
			for _, gx := range genomes {
				genome.RecycleGenome(gx)
			}

			*concat = (*concat)[:0]
			poolConcat.Put(concat)

		}(batchIDAndRefIDs)
	}

	wg.Wait()
	close(ch)
	<-done

	if debug {
		close(chDuration)
		<-doneDuration
		pbs.Wait()
		log.Debugf("%s (%s bp): finished aligning query genome fragments in %.3f seconds",
			query.id, humanize.Comma(int64(query.genomeSize)), time.Since(startTime).Seconds())
	}

	return nil
}
