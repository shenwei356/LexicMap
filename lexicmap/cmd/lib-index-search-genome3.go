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
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"regexp"
	"slices"
	"sort"
	"sync"
	"time"

	"github.com/dustin/go-humanize"
	"github.com/vbauerster/mpb/v8"
	"github.com/vbauerster/mpb/v8/decor"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/lexichash/iterator"
	"github.com/shenwei356/wfa"
)

var errGenomeTooLarge = errors.New("genome exceeds the configured size limit")

// subjectSketch is an in-memory sketch of a subject genome for sampled k-mer matching.
// It is recycled after all query fragments have been processed against the subject.
type subjectSketch struct {
	seqLen     int
	forwardLen int
	rcStart    int

	// Most sampled k-mers occur only once. Store that common case directly in
	// the main map, and keep repeated positions in one shared arena. Values in
	// sampledKmerMap are position+1, leaving zero available as the map's missing
	// value.
	sampledKmerMap        *map[uint64]uint32
	repeatedKmerMap       *map[uint64]uint64
	repeatedKmerPositions *[]repeatedKmerPosition
	contigBounds          [][2]int // [start, end) of each contig in forward strand
}

type repeatedKmerPosition struct {
	position uint32
	next     uint32 // index+1; zero marks the end of the list
}

var poolSubjectSketch = &sync.Pool{New: func() interface{} {
	return &subjectSketch{}
}}

// poolConcat is for reusing large byte slices for concatenated genome sequences
var poolConcat = &sync.Pool{New: func() interface{} {
	// The required size is known before use and varies with the subject genome.
	// Start empty so small genomes do not each reserve a fixed 10 MiB buffer.
	tmp := make([]byte, 0)
	return &tmp
}}

// These pools intentionally have no fixed-size New function. When a pool is
// empty, allocate from the known subject size instead of using one capacity for
// both short search candidates and whole bacterial genomes.
var poolKmerMap = &sync.Pool{}
var poolRepeatedKmerMap = &sync.Pool{}
var poolRepeatedKmerPositions = &sync.Pool{}

// poolQSeeds is for reusing query seed slices
var poolQSeeds = &sync.Pool{New: func() interface{} {
	s := make([][]uint64, 0, 256)
	return &s
}}

func recycleQuerySeeds(seeds *[][]uint64) {
	const maxPooledQueryFragments = 10240
	if cap(*seeds) > maxPooledQueryFragments {
		clear((*seeds)[:cap(*seeds)])
		*seeds = nil
		return
	}
	// Slots outside the current length may own buffers from a longer query.
	backing := (*seeds)[:cap(*seeds)]
	for i, words := range backing {
		if cap(words) > thresholdNSubsLong {
			backing[i] = nil
		} else {
			backing[i] = words[:0]
		}
	}
	*seeds = (*seeds)[:0]
	poolQSeeds.Put(seeds)
}

func addSampledKmerPosition(kmers *map[uint64]uint32, repeated *map[uint64]uint64, repeatedPositions *[]repeatedKmerPosition, kmer uint64, pos uint32) {
	_, ok := (*kmers)[kmer]
	if !ok {
		(*kmers)[kmer] = pos + 1
		return
	}

	// The first position is stored inline in kmers. Store subsequent positions
	// in one contiguous arena instead of allocating one slice per repeated k-mer.
	// The map value packs head and tail indexes (+1) to preserve insertion order.
	list := (*repeated)[kmer]
	head, tail := uint32(list>>32), uint32(list)
	*repeatedPositions = append(*repeatedPositions, repeatedKmerPosition{position: pos})
	index := uint32(len(*repeatedPositions))
	if head == 0 {
		head = index
	} else {
		(*repeatedPositions)[tail-1].next = index
	}
	(*repeated)[kmer] = uint64(head)<<32 | uint64(index)
}

func sampledKmerMapCapacity(genomeSize int, scale int) int {
	return (genomeSize + scale - 1) / scale
}

func acquireSampledKmerMaps(capacity int) (*map[uint64]uint32, *map[uint64]uint64, *[]repeatedKmerPosition) {
	var kmers *map[uint64]uint32
	if v := poolKmerMap.Get(); v != nil {
		kmers = v.(*map[uint64]uint32)
	} else {
		m := make(map[uint64]uint32, capacity)
		kmers = &m
	}

	repeatedCapacity := capacity / 16
	var repeated *map[uint64]uint64
	if v := poolRepeatedKmerMap.Get(); v != nil {
		repeated = v.(*map[uint64]uint64)
	} else {
		// Repeated canonical 13-mers are normally a small fraction of all
		// sampled k-mers. This hint avoids early map growth without reserving
		// the full primary-map size twice.
		m := make(map[uint64]uint64, repeatedCapacity)
		repeated = &m
	}

	var repeatedPositions *[]repeatedKmerPosition
	if v := poolRepeatedKmerPositions.Get(); v != nil {
		repeatedPositions = v.(*[]repeatedKmerPosition)
	} else {
		positions := make([]repeatedKmerPosition, 0, repeatedCapacity)
		repeatedPositions = &positions
	}

	return kmers, repeated, repeatedPositions
}

// Sampling parameters for the simplified seeding strategy.
var gsa3SampledK = 13     // fixed k-mer length for sampling
var gsa3SamplingScale = 4 // sampling rate: keep if hash(kmer) % scale == 0

// buildSubjectSketchSampledOptimized scans the forward sequence once and records
// enough strand information to address both the forward and RC concatenated copies.
func (idx *Index) buildSubjectSketchSampledOptimized(seq []byte, skipRegions [][2]int, contigBounds [][2]int, genomeSize int, forwardLen int, rcStart int) (*subjectSketch, error) {
	return idx.buildSubjectSketchSampled(seq, skipRegions, contigBounds, genomeSize, forwardLen, rcStart, true)
}

// buildSubjectSketchSampled can allocate maps owned by the bounded compare cache.
func (idx *Index) buildSubjectSketchSampled(seq []byte, skipRegions [][2]int, contigBounds [][2]int, genomeSize int, forwardLen int, rcStart int, pooled bool) (*subjectSketch, error) {
	k := gsa3SampledK
	k8 := uint8(k)
	scale := uint64(gsa3SamplingScale)
	scaleM1 := scale - 1

	if len(seq) < k || forwardLen < k || rcStart <= forwardLen || rcStart >= len(seq) {
		return nil, fmt.Errorf("invalid forward/RC sequence layout for k=%d", k)
	}

	// Low-complexity k-mer values
	ccc := util.Ns(0b01, k8)
	ggg := util.Ns(0b10, k8)
	ttt := (uint64(1) << (k << 1)) - 1

	// Only scan the forward strand. A sampled occurrence is encoded as
	// (position << 1) | canonical-strand and converted to its forward or RC
	// coordinate during lookup.
	iter, err := iterator.NewKmerIterator(seq[:forwardLen], k)
	if err != nil {
		return nil, err
	}

	mapCapacity := sampledKmerMapCapacity(genomeSize, int(scale))
	var kmerMap *map[uint64]uint32
	var repeatedKmerMap *map[uint64]uint64
	var repeatedKmerPositions *[]repeatedKmerPosition
	if pooled {
		kmerMap, repeatedKmerMap, repeatedKmerPositions = acquireSampledKmerMaps(mapCapacity)
	} else {
		m := make(map[uint64]uint32, mapCapacity)
		r := make(map[uint64]uint64, mapCapacity/16)
		p := make([]repeatedKmerPosition, 0, mapCapacity/16)
		kmerMap, repeatedKmerMap, repeatedKmerPositions = &m, &r, &p
	}

	region := 0
	var canonical uint64
	var canonicalRC uint32
	for pos := 0; ; pos++ {
		kmer, kmerRC, ok, _ := iter.NextKmer()
		if !ok {
			break
		}

		// skipRegions is sorted by start. Expand each interval to the left by
		// k-1 so no retained k-mer crosses a gap.
		for region < len(skipRegions) && pos > skipRegions[region][1] {
			region++
		}
		if region < len(skipRegions) && pos >= skipRegions[region][0]-k+1 {
			continue
		}

		canonical = kmer
		canonicalRC = 0
		if kmerRC < kmer {
			canonical = kmerRC
			canonicalRC = 1
		}

		if util.Hash64(canonical)&scaleM1 != 0 {
			continue
		}

		if kmer == ccc || kmer == ggg || kmer == ttt || util.IsLowComplexityDust(kmer, k8) {
			continue
		}

		addSampledKmerPosition(kmerMap, repeatedKmerMap, repeatedKmerPositions, canonical, uint32(pos)<<1|canonicalRC)
	}

	var s *subjectSketch
	if pooled {
		s = poolSubjectSketch.Get().(*subjectSketch)
	} else {
		s = &subjectSketch{}
	}
	s.seqLen = len(seq)
	s.forwardLen = forwardLen
	s.rcStart = rcStart
	s.sampledKmerMap = kmerMap
	s.repeatedKmerMap = repeatedKmerMap
	s.repeatedKmerPositions = repeatedKmerPositions
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
	if s.repeatedKmerMap != nil {
		clear(*s.repeatedKmerMap)
		poolRepeatedKmerMap.Put(s.repeatedKmerMap)
		s.repeatedKmerMap = nil
	}
	if s.repeatedKmerPositions != nil {
		*s.repeatedKmerPositions = (*s.repeatedKmerPositions)[:0]
		poolRepeatedKmerPositions.Put(s.repeatedKmerPositions)
		s.repeatedKmerPositions = nil
	}
	s.contigBounds = nil
	poolSubjectSketch.Put(s)
}

// sampleQueryFragment samples fixed-length k-mers from a query fragment.
func sampleQueryFragment(frag []byte, sampledKmers []uint64) ([]uint64, error) {
	sampledKmers = sampledKmers[:0]
	k := gsa3SampledK
	k8 := uint8(k)
	scale := uint64(gsa3SamplingScale)
	scaleM1 := scale - 1

	if len(frag) < k {
		return sampledKmers, nil
	}

	ccc := util.Ns(0b01, k8)
	ggg := util.Ns(0b10, k8)
	ttt := (uint64(1) << (k << 1)) - 1

	iter, err := iterator.NewKmerIterator(frag, k)
	if err != nil {
		return sampledKmers, err
	}

	pos := 0
	var canonical uint64
	var canonicalRC uint64
	for {
		kmer, kmerRC, ok, _ := iter.NextKmer()
		if !ok {
			break
		}

		// Use canonical k-mer for sampling and lookup. Store the query
		// canonical strand in the low bit of the position.
		canonical = kmer
		canonicalRC = 0
		if kmerRC < kmer {
			canonical = kmerRC
			canonicalRC = 1
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

		sampledKmers = append(sampledKmers, canonical, uint64(pos)<<1|canonicalRC)

		pos++
	}

	return sampledKmers, nil
}

// alignQueryFragToSubjectSampled matches a query fragment against a sampled subject sketch.
func alignQueryFragToSubjectSampled(
	qfrag []byte,
	qqual []byte,
	qSeeds []uint64,
	sketch *subjectSketch,
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
	queryIndexes *queryFragmentIndexes, // nil skips entry reuse in genome compare
	fragmentIndex int,
) (int, int, int, float64, float64, bool) {
	// Since we only use forward strand query k-mers and subject is a single concatenated
	// sequence (forward + RC), we only need one set of anchors for unified chaining.
	allSubs := poolSubsLong.Get().(*[]SubstrPair)
	*allSubs = (*allSubs)[:0]
	defer RecycleSubstrPairs(poolSubsLong, allSubs)

	if sketch.sampledKmerMap == nil || len(*sketch.sampledKmerMap) == 0 {
		return 0, 0, 0, 0, 0, false
	}

	qKmers := qSeeds
	sKmerMap := sketch.sampledKmerMap
	repeatedKmerMap := sketch.repeatedKmerMap
	repeatedKmerPositions := *sketch.repeatedKmerPositions

	// Match query k-mers against subject k-mers
	// Limit matches per k-mer to avoid excessive anchors from repetitive sequences
	const maxMatchesPerKmer = 100

	for i := 0; i+1 < len(qKmers); i += 2 {
		qk := qKmers[i]
		first, found := (*sKmerMap)[qk]
		if !found {
			continue
		}

		qposAndStrand := qKmers[i+1]
		qpos := int32(qposAndStrand >> 1)
		qCanonicalRC := uint32(qposAndStrand & 1)
		// The first occurrence is stored inline in the primary map.
		sub := SubstrPair{}
		sub.Len = uint8(K)
		sub.QBegin = qpos
		sposAndStrand := first - 1
		spos := int(sposAndStrand >> 1)
		if sposAndStrand&1 != qCanonicalRC {
			spos = sketch.rcStart + sketch.forwardLen - spos - K
		}
		sub.TBegin = int32(spos)
		sub.QRC = false
		sub.TRC = false
		*allSubs = append(*allSubs, sub)

		// Leave one slot for the inline first occurrence. Additional positions
		// are linked in insertion order in the shared position arena.
		positionIndex := uint32((*repeatedKmerMap)[qk] >> 32)
		for matches := 1; positionIndex != 0 && matches < maxMatchesPerKmer; matches++ {
			position := repeatedKmerPositions[positionIndex-1]
			posAndStrand := position.position
			spos := int(posAndStrand >> 1)
			if posAndStrand&1 != qCanonicalRC {
				spos = sketch.rcStart + sketch.forwardLen - spos - K
			}
			sub := SubstrPair{}
			sub.Len = uint8(K)
			sub.QBegin = qpos
			sub.TBegin = int32(spos)
			sub.QRC = false
			sub.TRC = false
			*allSubs = append(*allSubs, sub)
			positionIndex = position.next
		}
	}

	chains, chainsOk := chainsFromSubs(allSubs, chainer, K)
	if !chainsOk {
		return 0, 0, 0, 0, 0, false
	}

	// Try all chains and pick the best one
	var bestMatched, bestAligned, bestGaps int
	var bestPident, bestPidentAdjusted float64
	var bestScore int = -1
	topChains := idx.chainingOptions.TopChains
	onlyTopChains := topChains > 0

	// Pre-index qfrag once for all chains.
	cpr := idx.poolSeqComparator.Get().(*SeqComparator)
	defer idx.poolSeqComparator.Put(cpr)
	if err := queryIndexes.index(cpr, fragmentIndex, qfrag); err != nil {
		return 0, 0, 0, 0, 0, false
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
		matched, aligned, gaps, pident, pidentAdjusted, ok := alignChain(
			qfrag, qqual, concat, chain, sketch, algn, cpr,
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
				bestPidentAdjusted = pidentAdjusted
			}
		}
	}
	RecycleChaining2Result(chains)

	if bestScore <= 0 {
		return 0, 0, 0, 0, 0, false
	}

	return bestMatched, bestAligned, bestGaps, bestPident, bestPidentAdjusted, true
}

func adjustedPIdentForAlignment(pident float64, qqual []byte, extendedQueryStart int, cigar *wfa.AlignmentResult) float64 {
	if len(qqual) == 0 || cigar == nil {
		return pident
	}
	qStart := extendedQueryStart + cigar.QBegin - 1
	qEnd := extendedQueryStart + cigar.QEnd
	if qStart < 0 || qEnd > len(qqual) || qStart >= qEnd {
		return pident
	}
	return adjustPIdentByQual(pident, qqual[qStart:qEnd])
}

// chainsFromSubs runs the chaining pipeline and returns all chains.
func chainsFromSubs(subs *[]SubstrPair, chainer *Chainer2, K int) (*[]*Chain2Result, bool) {
	if len(*subs) == 0 {
		return nil, false
	}

	if len(*subs) > 1 {
		ClearSubstrPairs(subs, K)
	}
	TrimSubStrPairs(subs, K, 100)
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
// Returns raw and quality-adjusted identity together with alignment statistics.
func alignChain(
	qfrag []byte,
	qqual []byte,
	subjectSeq []byte,
	chain *Chain2Result,
	sketch *subjectSketch,
	algn *wfa.Aligner,
	cpr *SeqComparator, // pre-indexed comparator
	extLen int,
	extLen2 int,
	minPIdent float64,
	minQcov float64,
	idx *Index,
	fScoreAndEvalue *func(qlen int, cigar *wfa.AlignmentResult) (int, int, float64),
) (int, int, int, float64, float64, bool) {
	// Guard against degenerate chains.
	if chain.QEnd < chain.QBegin || chain.TEnd < chain.TBegin {
		return 0, 0, 0, 0, 0, false
	}

	qLen := len(qfrag)

	// Locate the contig containing the chain on either concatenated strand.
	contigStart, contigEnd := subjectContigBounds(sketch, chain.TBegin)

	// Expand the chain region by extLen.
	tExpBegin := max(chain.TBegin-extLen, contigStart)
	tExpEnd := min(chain.TEnd+extLen, contigEnd-1)
	tSubseq := subjectSeq[tExpBegin : tExpEnd+1]

	qExpBegin := max(chain.QBegin-extLen, 0)
	qExpEnd := min(chain.QEnd+extLen, qLen-1)

	// Pseudo-alignment with SeqComparator (already indexed).
	cr, err := cpr.Compare(uint32(qExpBegin), uint32(qExpEnd), tSubseq, qLen)
	if err != nil {
		return 0, 0, 0, 0, 0, false
	}
	if cr == nil {
		return 0, 0, 0, 0, 0, false
	}
	defer RecycleSeqComparatorResult(cr)

	// WFA alignment on each sub-chain.
	var totMatched, totAligned, totGaps int
	var pidentAdjusted float64
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

		_qseq, _tseq, qLeftExt, _, _, _, extErr := extendMatch(
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
		pident := float64(cigar.Matches) / float64(cigar.AlignLen) * 100
		pidentAdjusted = adjustedPIdentForAlignment(pident, qqual, c.QBegin-qLeftExt, cigar)
		wfa.RecycleAlignmentResult(cigar)

		break // keep the best ONE match
	}

	if totAligned <= 0 {
		return 0, 0, 0, 0, 0, false
	}

	pident := float64(totMatched) / float64(totAligned) * 100
	alignedBasesQ := totAligned - totGaps
	af := float64(alignedBasesQ) / float64(qLen) * 100
	if af > 100 {
		af = 100
	}
	if pident < minPIdent || af < minQcov {
		return 0, 0, 0, 0, 0, false
	}

	return totMatched, totAligned, totGaps, pident, pidentAdjusted, true
}

func subjectContigBounds(sketch *subjectSketch, position int) (int, int) {
	if len(sketch.contigBounds) == 0 {
		return 0, sketch.seqLen
	}
	forwardPosition := position
	if position >= sketch.rcStart {
		forwardPosition = sketch.forwardLen - 1 - (position - sketch.rcStart)
	}
	i := sort.Search(len(sketch.contigBounds), func(i int) bool {
		return sketch.contigBounds[i][1] > forwardPosition
	})
	if i < len(sketch.contigBounds) && forwardPosition >= sketch.contigBounds[i][0] {
		bounds := sketch.contigBounds[i]
		if position >= sketch.rcStart {
			return sketch.rcStart + sketch.forwardLen - bounds[1], sketch.rcStart + sketch.forwardLen - bounds[0]
		}
		return bounds[0], bounds[1]
	}
	return 0, sketch.seqLen
}

// GSearchAlign3Sampled is a simplified version of GSearchAlign3 that uses
// sampled fixed-length k-mers instead of LexicHash masking.
func (idx *Index) GSearchAlign3Sampled(query *GQuery, fragLen int, minFragLen int, genomeIds *map[uint64]*[]uint64, minAF, minANI float64, maxQueryConcurrency int) error {
	debug := idx.opt.Debug

	startTime0 := time.Now()

	if debug {
		log.Debugf("%s (%s bp): start to preprocess query genome fragments", query.id, humanize.Comma(int64(query.genomeSize)))
	}

	// 1) Cut the query into fragments.
	qfrags, qqualFrags, qfragLens, err := seqs2fragmentsWithQual(&query.seqs, &query.quals, fragLen, minFragLen)
	if err != nil {
		return fmt.Errorf("failed to cut query sequences and qualities into fragments: %w", err)
	}
	defer recycleFragments(qfrags)
	defer recycleFragments(qqualFrags)
	if len(*qfrags) == 0 {
		return fmt.Errorf("no fragments for alignment, are the genome too fragmented with all sequences shorter than the minimum fragment length (%d bp)?", minFragLen)
	}

	// 2) Sample k-mers from each query fragment.
	qSeeds := poolQSeeds.Get().(*[][]uint64)
	*qSeeds = slices.Grow((*qSeeds)[:0], len(*qfrags))[:len(*qfrags)]
	defer recycleQuerySeeds(qSeeds)

	for i, qfrag := range *qfrags {
		seeds, err := sampleQueryFragment(qfrag, (*qSeeds)[i])
		(*qSeeds)[i] = seeds
		if err != nil {
			return fmt.Errorf("failed to sample query fragment: %w", err)
		}
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
	queryIndexes := newQueryFragmentIndexes(idx, len(*qfrags), len(*genomeIds))
	defer queryIndexes.close()

	// -----------------------------------------------------------
	// Progress reporting respects quiet mode, independently of debug logging.
	showProgressBar := debug && idx.opt.Verbose
	var pbs *mpb.Progress
	var bar *mpb.Bar
	var chDuration chan time.Duration
	var doneDuration chan int
	if showProgressBar {
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
		trimGSearchResults(rs, idx.opt.TopN)

		query.result = rs
		done <- 1
	}()

	// 4) read genomes and align

	K := gsa3SampledK
	contigInterval := max(K, int(float64(fragLen)*1.5))
	nnn := bytes.Repeat([]byte{'N'}, contigInterval)

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
				if showProgressBar {
					chDuration <- time.Duration(float64(time.Since(timeStart)) / fcpus)
				}
				<-tokens
				wg.Done()
			}()

			var g *genome.Genome
			genomes := make([]*genome.Genome, len(*batchIDAndRefIDs))
			maxSubjectGenomeSize := idx.opt.MaxSubjectGenomeSize
			for i, batchIDAndRefID := range *batchIDAndRefIDs {
				genomeBatch := int(batchIDAndRefID >> BITS_GENOME_IDX)
				genomeIdx := int(batchIDAndRefID & MASK_GENOME_IDX)

				rdr, err := idx.acquireGenomeReader(genomeBatch)
				if err != nil {
					checkError(err)
				}
				_g, err := rdr.Seqs(genomeIdx)
				if err != nil {
					_ = idx.releaseGenomeReader(genomeBatch, rdr)
					checkError(fmt.Errorf("fail to read genome sequence for batch %d, genome index %d: %s", genomeBatch, genomeIdx, err))
				}
				genomes[i] = _g
				if i == 0 {
					g = _g
				} else {
					g.Seqs = append(g.Seqs, _g.Seqs...)
					_g.Seqs = nil
					g.NumSeqs += _g.NumSeqs
					g.GenomeSize += _g.GenomeSize
				}
				if err := idx.releaseGenomeReader(genomeBatch, rdr); err != nil {
					checkError(fmt.Errorf("failed to close genome reader: %s", err))
				}

				if maxSubjectGenomeSize > 0 && g.GenomeSize > maxSubjectGenomeSize {
					log.Warningf("%s (size: %s bp) exceeds the maximum allowed subject genome size of %s, consider increasing --max-subject-genome-size",
						idx.BatchGenomeIndex2GenomeID[(*batchIDAndRefIDs)[0]],
						humanize.Comma(int64(g.GenomeSize)),
						humanize.Comma(int64(maxSubjectGenomeSize)))
					for _, gx := range genomes {
						if gx != nil {
							genome.RecycleGenome(gx)
						}
					}
					return
				}
			}

			concat := poolConcat.Get().(*[]byte)
			*concat = (*concat)[:0]
			forwardSize := contigInterval * (len(g.Seqs) - 1)
			for _, s := range g.Seqs {
				forwardSize += len(*s)
			}
			rcInterval := fragLen << 1
			totalSize := forwardSize<<1 + rcInterval
			if cap(*concat) < totalSize {
				*concat = make([]byte, 0, totalSize)
			}

			var skipRegions [][2]int
			contigBounds := make([][2]int, 0, len(g.Seqs))
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
			if gaps := findGapRegions(*concat, 5); gaps != nil {
				for _, gap := range *gaps {
					start, end := unpackGapRegion(gap)
					skipRegions = append(skipRegions, [2]int{start, end - 1})
				}
				recycleGapRegions(gaps)
			}

			forwardLen := len(*concat)
			*concat = append(*concat, bytes.Repeat([]byte{'N'}, rcInterval)...)
			rcStart := len(*concat)
			*concat = append(*concat, (*concat)[:forwardLen]...)
			RC((*concat)[rcStart:])
			slices.SortFunc(skipRegions, func(a, b [2]int) int { return a[0] - b[0] })

			sketch, err := idx.buildSubjectSketchSampledOptimized(*concat, skipRegions, contigBounds, g.GenomeSize, forwardLen, rcStart)
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
				var qqual []byte
				if qqualFrags != nil {
					qqual = (*qqualFrags)[i]
				}
				matched, alignedLen, gaps, pident, pidentAdjusted, ok := alignQueryFragToSubjectSampled(
					qfrag, qqual, (*qSeeds)[i], sketch, (*concat),
					chainer, algn, K, extLen, extLen2,
					minPIdent, minQcovHSP, idx,
					&fScoreAndEvalue,
					queryIndexes, i,
				)
				if !ok {
					// fmt.Printf("fail to align fragment %d: %s\n", i+1, qfrag)
					continue
				}
				gr.AlignedFragments++
				gr.AlignedLength += alignedLen - gaps
				gr.AlignedMatches += matched
				gr.PidentsSum += pident
				if qqualFrags != nil {
					gr.PidentsAdjustedSum += pidentAdjusted
				}
			}

			// g) ANI / AF on the accumulated alignment.
			if gr.AlignedFragments > 0 {
				gr.ANI = gr.PidentsSum / float64(gr.AlignedFragments) / 100
				if qqualFrags != nil {
					gr.ANIAdjusted = gr.PidentsAdjustedSum / float64(gr.AlignedFragments) / 100
				}
			}
			gr.AFq = float64(gr.AlignedLength) / float64(qfragLens)
			if gr.AFq > 1 {
				gr.AFq = 1
			}
			gr.AFs = float64(gr.AlignedLength) / float64(gr.GenomeSize)
			if gr.AFs > 1 {
				gr.AFs = 1
			}
			gr.Score = gr.ANI

			if gr.AFq < minAF || gr.ANI < minANI {
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

	if showProgressBar {
		close(chDuration)
		<-doneDuration
		pbs.Wait()
		log.Debugf("%s (%s bp): finished aligning query genome fragments in %.3f seconds",
			query.id, humanize.Comma(int64(query.genomeSize)), time.Since(startTime).Seconds())
	}

	return nil
}

// buildComparisonSubject preserves contig gaps and forward/RC coordinates.
// Cached subjects use fresh buffers so evicted maps are not retained in pools.
func (idx *Index) buildComparisonSubject(subject *GQuery, fragLen int, pooled bool) (*[]byte, *subjectSketch, error) {
	K := gsa3SampledK
	contigInterval := int(float64(fragLen) * 1.5)
	if contigInterval < K {
		contigInterval = K
	}
	nnn := bytes.Repeat([]byte{'N'}, contigInterval)

	var concat *[]byte
	if pooled {
		concat = poolConcat.Get().(*[]byte)
		*concat = (*concat)[:0]
	} else {
		buffer := []byte(nil)
		concat = &buffer
	}

	// Calculate total size: forward + contig intervals + RC interval + RC
	var forwardSize int
	for _, s := range subject.seqs {
		forwardSize += len(*s)
	}
	forwardSize += contigInterval * (len(subject.seqs) - 1)

	// Total size = forward + 2*fragLen interval + RC (same as forward)
	rcInterval := fragLen << 1
	totalSize := forwardSize<<1 + rcInterval

	// Pre-allocate the full capacity to avoid reallocation
	if cap(*concat) < totalSize {
		*concat = make([]byte, 0, totalSize)
	}

	var skipRegions [][2]int
	contigBounds := make([][2]int, 0, len(subject.seqs))
	for i, s := range subject.seqs {
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
	gaps := findGapRegions(*concat, 5)
	if gaps != nil {
		for _, gap := range *gaps {
			start, end := unpackGapRegion(gap)
			skipRegions = append(skipRegions, [2]int{start, end - 1})
		}
		recycleGapRegions(gaps)
	}

	// Append 2*fragLen interval and reverse complement strand
	forwardLen := len(*concat)
	nnnRC := bytes.Repeat([]byte{'N'}, rcInterval)

	// Add interval between forward and RC strands
	*concat = append(*concat, nnnRC...)

	// Append reverse complement of the forward strand
	rcStart := len(*concat)
	*concat = append(*concat, (*concat)[:forwardLen]...)
	RC((*concat)[rcStart:])

	// Sort skip regions
	slices.SortFunc(skipRegions, func(a, b [2]int) int {
		return a[0] - b[0]
	})

	// 4) Build the subject sketch using sampled k-mers
	sketch, err := idx.buildSubjectSketchSampled(*concat, skipRegions, contigBounds, subject.genomeSize, forwardLen, rcStart, pooled)
	if err != nil {
		if pooled {
			poolConcat.Put(concat)
		}
		return nil, nil, fmt.Errorf("fail to build subject sketch: %s", err)
	}
	return concat, sketch, nil
}

// CompareTwoGenomes compares two genomes directly without using an index.
// It's adapted from GSearchAlign3Sampled but compares query vs subject directly.
func (idx *Index) CompareTwoGenomes(query, subject *GQuery, fragLen int, minFragLen int, minAF, minANI float64) error {
	return idx.compareTwoGenomesPrepared(query, subject, nil, nil, fragLen, minFragLen, minAF, minANI)
}

// compareTwoGenomesPrepared borrows cached data while keeping results pair-local.
func (idx *Index) compareTwoGenomesPrepared(query, subject *GQuery, qp, sp *preparedCompareGenome, fragLen, minFragLen int, minAF, minANI float64) error {
	var qfrags *[][]byte
	var qfragLens int
	var qSeeds *[][]uint64
	if qp != nil {
		qfrags, qfragLens, qSeeds = qp.fragments, qp.fragmentBases, &qp.seeds
	} else {
		// 1) Cut the query into fragments.
		qfrags, qfragLens = seqs2fragments(&query.seqs, fragLen, minFragLen)
		defer recycleFragments(qfrags)
		if qfrags == nil || len(*qfrags) == 0 {
			return fmt.Errorf("no fragments for alignment, are the genome too fragmented with all sequences shorter than the minimum fragment length (%d bp)?", minFragLen)
		}

		// 2) Sample k-mers from each query fragment.
		qSeeds = poolQSeeds.Get().(*[][]uint64)
		*qSeeds = slices.Grow((*qSeeds)[:0], len(*qfrags))[:len(*qfrags)]
		defer recycleQuerySeeds(qSeeds)

		for i, qfrag := range *qfrags {
			seeds, err := sampleQueryFragment(qfrag, (*qSeeds)[i])
			(*qSeeds)[i] = seeds
			if err != nil {
				return fmt.Errorf("failed to sample query fragment: %w", err)
			}
		}
	}

	// Subject layout and sampled map can be reused in both comparison directions.
	K := gsa3SampledK
	var concat *[]byte
	var sketch *subjectSketch
	if sp != nil {
		concat, sketch = &sp.concat, sp.sketch
	} else {
		var err error
		concat, sketch, err = idx.buildComparisonSubject(subject, fragLen, true)
		if err != nil {
			return err
		}
		defer idx.recycleSubjectSketch(sketch)
		defer func() { *concat = (*concat)[:0]; poolConcat.Put(concat) }()
	}

	// 5) Set up alignment tools
	alignOption := &wfa.Options{GlobalAlignment: true}
	minPIdent := idx.seqCompareOption.MinIdentity
	minQcovHSP := idx.seqCompareOption.MinAlignedFraction
	extLen := fragLen / 2
	extLen2 := idx.opt.ExtendLength2

	chainer := idx.poolChainers2.Get().(*Chainer2)
	defer idx.poolChainers2.Put(chainer)

	algn := wfa.New(wfa.DefaultPenalties, alignOption)
	algn.AdaptiveReduction(wfa.DefaultAdaptiveOption)
	defer wfa.RecycleAligner(algn)

	gr := poolGSearchResult.Get().(*GSearchResult)
	gr.Reset()
	gr.BatchGenomeIndex = 0 // Not from index
	gr.GenomeSize = subject.genomeSize
	gr.NumSeqs = len(subject.seqs)

	// 6) Align each query fragment to subject
	fScoreAndEvalue := scoreAndEvalue(2, -3, 5, 2, int(subject.genomeSize), 0.625, 0.41)
	for i, qfrag := range *qfrags {
		matched, alignedLen, gaps, pident, _, ok := alignQueryFragToSubjectSampled(
			qfrag, nil, (*qSeeds)[i], sketch, (*concat),
			chainer, algn, K, extLen, extLen2,
			minPIdent, minQcovHSP, idx,
			&fScoreAndEvalue,
			nil, 0,
		)
		if !ok {
			continue
		}
		gr.AlignedFragments++
		gr.AlignedLength += alignedLen - gaps
		gr.AlignedMatches += matched
		gr.PidentsSum += pident
	}

	// 7) Calculate ANI / AF on the accumulated alignment
	if gr.AlignedFragments > 0 {
		gr.ANI = gr.PidentsSum / float64(gr.AlignedFragments) / 100
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

	// 8) Store result
	rs := poolGSearchResults.Get().(*[]*GSearchResult)
	*rs = (*rs)[:0]
	if gr.AFq >= minAF && gr.ANI >= minANI {
		*rs = append(*rs, gr)
	} else {
		poolGSearchResult.Put(gr)
		RecycleGSearchResults(rs)
		return nil
	}

	query.result = rs
	return nil
}

// ReadGenome reads a genome from the index
func (idx *Index) ReadGenome(batchIDAndRefIDs *[]uint64, genomeID string) (*GQuery, error) {

	maxSubjectGenomeSize := idx.opt.MaxSubjectGenomeSize

	q := poolGQuery.Get().(*GQuery)
	q.Reset()

	for _, batchIDAndRefID := range *batchIDAndRefIDs {
		genomeBatch := int(batchIDAndRefID >> BITS_GENOME_IDX)
		genomeIdx := int(batchIDAndRefID & MASK_GENOME_IDX)

		rdr, err := idx.acquireGenomeReader(genomeBatch)
		if err != nil {
			RecycleGQuery(q)
			return nil, err
		}

		g, err := rdr.Seqs(genomeIdx)
		if err != nil {
			RecycleGQuery(q)
			_ = idx.releaseGenomeReader(genomeBatch, rdr)
			return nil, fmt.Errorf("fail to read genome sequence for batch %d, genome index %d: %s", genomeBatch, genomeIdx, err)
		}
		if maxSubjectGenomeSize > 0 && g.GenomeSize > maxSubjectGenomeSize-q.genomeSize {
			log.Warningf("skipped genome %s (size: at least %s bp), which exceeds the maximum allowed size of %s; consider increasing --max-genome-size",
				genomeID,
				humanize.Comma(int64(q.genomeSize)+int64(g.GenomeSize)),
				humanize.Comma(int64(maxSubjectGenomeSize)))

			if err := idx.releaseGenomeReader(genomeBatch, rdr); err != nil {
				RecycleGQuery(q)
				genome.RecycleGenome(g)
				return nil, fmt.Errorf("failed to close genome reader: %w", err)
			}
			genome.RecycleGenome(g)
			RecycleGQuery(q)
			return nil, errGenomeTooLarge
		}

		for _, s1 := range g.Seqs {
			s := poolSeq.Get().(*[]byte)
			*s = (*s)[:0]
			*s = append(*s, *s1...)
			q.seqs = append(q.seqs, s)

			q.genomeSize += len(*s1)
		}

		if err := idx.releaseGenomeReader(genomeBatch, rdr); err != nil {
			RecycleGQuery(q)
			genome.RecycleGenome(g)
			return nil, fmt.Errorf("failed to close genome reader: %w", err)
		}
		genome.RecycleGenome(g)
	}

	return q, nil
}

// CompareTwoGenomesOrthoANI compares two genomes using the OrthoANI algorithm.
// This method cuts both query and subject genomes into fragments and only uses
// orthologous fragment pairs (reciprocal best hits) for ANI/AF calculation.
// Based on GSearchAlign2 from lib-index-search-genome.go.
func (idx *Index) CompareTwoGenomesOrthoANI(query, subject *GQuery, fragLen int, minFragLen int, minAF, minANI float64) error {
	return idx.compareTwoGenomesOrthoANIPrepared(query, subject, nil, nil, fragLen, minFragLen, minAF, minANI)
}

// compareTwoGenomesOrthoANIPrepared borrows sorted entries without changing pair selection.
func (idx *Index) compareTwoGenomesOrthoANIPrepared(query, subject *GQuery, qp, sp *preparedCompareGenome, fragLen, minFragLen int, minAF, minANI float64) error {
	var qfrags, sfrags *[][]byte
	var qfragLens, sfragLens int
	if qp != nil {
		qfrags, qfragLens = qp.fragments, qp.fragmentBases
	} else {
		qfrags, qfragLens = seqs2fragments(&query.seqs, fragLen, minFragLen)
		defer recycleFragments(qfrags)
	}
	if sp != nil {
		sfrags, sfragLens = sp.fragments, sp.fragmentBases
	} else {
		sfrags, sfragLens = seqs2fragments(&subject.seqs, fragLen, minFragLen)
		defer recycleFragments(sfrags)
	}
	if qfrags == nil || len(*qfrags) == 0 || sfrags == nil || len(*sfrags) == 0 {
		return fmt.Errorf("no fragments for alignment (minimum fragment length: %d bp)", minFragLen)
	}

	fcpr := idx.poolFragmentComparator.Get().(*FragmentComparator)
	defer idx.poolFragmentComparator.Put(fcpr)
	var pairs *[]uint64
	var err error
	if qp != nil && sp != nil {
		pairs = fcpr.scanPairsMerged(qp.entries, sp.entries)
	} else {
		pairs, err = fcpr.Compare(qfrags, sfrags)
		if err != nil {
			return fmt.Errorf("fail to find similar fragments: %s", err)
		}
	}
	defer RecycleFragmentCompareResult(pairs)

	// Sort pairs for better cache locality
	slices.Sort(*pairs)

	// 4) Prepare reverse complement fragments for subject
	sfragsRC := poolFragments.Get().(*[][]byte)
	n := len(*sfrags)
	if cap(*sfragsRC) >= n {
		*sfragsRC = (*sfragsRC)[:n]
		clear(*sfragsRC)
	} else {
		*sfragsRC = (*sfragsRC)[:cap(*sfragsRC)]
		clear(*sfragsRC)
		for len(*sfragsRC) < n {
			*sfragsRC = append(*sfragsRC, nil)
		}
	}
	defer recycleFragments(sfragsRC)

	// 5) Set up alignment tools
	alignOption := &wfa.Options{GlobalAlignment: true}
	fScoreAndEvalue := scoreAndEvalue(2, -3, 5, 2, int(subject.genomeSize), 0.625, 0.41)
	maxEvalue := idx.opt.MaxEvalue

	cpr := idx.poolSeqComparator.Get().(*SeqComparator)
	defer idx.poolSeqComparator.Put(cpr)
	defer cpr.RecycleIndex()

	algn := wfa.New(wfa.DefaultPenalties, alignOption)
	algn.AdaptiveReduction(wfa.DefaultAdaptiveOption)
	defer wfa.RecycleAligner(algn)

	minQcovHSP := idx.seqCompareOption.MinAlignedFraction
	minPIdent := idx.seqCompareOption.MinIdentity

	// 6) Maps to store alignment results for each fragment
	ma := poolFragAlignResultMap.Get().(*map[uint32]*[]*Chain2Result)
	mb := poolFragAlignResultMap.Get().(*map[uint32]*[]*Chain2Result)
	defer func() {
		for _, ls := range *ma {
			for _, c := range *ls {
				recycleChain2(c)
			}
			recycleChaining2ResultSlice(ls)
		}
		clear(*ma)
		poolFragAlignResultMap.Put(ma)

		for _, ls := range *mb {
			recycleChaining2ResultSlice(ls)
		}
		clear(*mb)
		poolFragAlignResultMap.Put(mb)
	}()

	// 7) Align fragment pairs
	var a, b, b2 []byte
	var ia, ib uint64
	var cr, cr2 *SeqComparatorResult
	var ls *[]*Chain2Result
	var ok bool
	var c *Chain2Result
	var indexedIA uint64
	var hasIndexedIA bool

	for _, p := range *pairs {
		ia, ib = p>>32, p&4294967295
		a = (*qfrags)[ia]
		b = (*sfrags)[ib]

		// a) pseudo alignment
		if !hasIndexedIA || ia != indexedIA {
			cpr.RecycleIndex()
			err = cpr.Index(a)
			if err != nil {
				return fmt.Errorf("fail to index query fragment: %s", err)
			}
			indexedIA = ia
			hasIndexedIA = true
		}

		// positive strand
		cr, err = cpr.Compare(0, uint32(len(a)), b, len(a))
		if err != nil {
			return fmt.Errorf("fail to compare query fragment and subject fragment: %s", err)
		}

		// negative strand
		b2 = (*sfragsRC)[ib]
		if b2 == nil {
			b2 = make([]byte, len(b))
			copy(b2, b)
			RC(b2)
			(*sfragsRC)[ib] = b2
		}
		cr2, err = cpr.Compare(0, uint32(len(a)), b2, len(a))
		if err != nil {
			return fmt.Errorf("fail to compare query fragment and rc subject fragment: %s", err)
		}

		if cr == nil && cr2 == nil {
			continue
		}

		if cr != nil && cr2 != nil { // both strands have hits
			// choose the strand with the longer aligned length, if tie, choose the positive strand
			if (*cr.Chains)[0].QEnd-(*cr.Chains)[0].QBegin < (*cr2.Chains)[0].QEnd-(*cr2.Chains)[0].QBegin {
				RecycleSeqComparatorResult(cr)
				cr = cr2
				b = b2
			} else {
				RecycleSeqComparatorResult(cr2)
			}
		} else if cr == nil { // only has hit in the negative strand
			cr = cr2
			b = b2
		} // only has hit in the positive strand

		// b) base-level alignment
		c = (*cr.Chains)[0]   // choose the first chain with the highest chaining score
		(*cr.Chains)[0] = nil // avoid being recycled before we finish processing c

		_qseq, _tseq, _, _, _, _, err := extendMatch(a, b, c.QBegin, c.QEnd+1, c.TBegin, c.TEnd+1, idx.opt.ExtendLength2, c.TBegin, idx.opt.ExtendLength2, false)
		if err != nil {
			RecycleSeqComparatorResult(cr)
			return fmt.Errorf("fail to extend aligned region: %s", err)
		}

		cigar, err := algn.Align(_qseq, _tseq)
		if err != nil {
			RecycleSeqComparatorResult(cr)
			return fmt.Errorf("fail to align sequences: %s", err)
		}

		_, _, evalue := fScoreAndEvalue(len(_qseq), cigar)
		if evalue > maxEvalue {
			recycleChain2(c)
			wfa.RecycleAlignmentResult(cigar)
			RecycleSeqComparatorResult(cr)
			continue
		}

		c.AlignedBasesQ = cigar.QEnd - cigar.QBegin + 1
		c.AlignedLength = int(cigar.AlignLen)
		c.MatchedBases = int(cigar.Matches)
		c.Gaps = int(cigar.Gaps)
		c.AlignedFraction = float64(c.AlignedBasesQ) / float64(cr.QueryLen) * 100
		if c.AlignedFraction > 100 {
			c.AlignedFraction = 100
		}
		c.PIdent = float64(c.MatchedBases) / float64(cigar.AlignLen) * 100

		c.Evalue = c.AlignedFraction * c.PIdent // just for sorting, not the real e-value

		// c) filter and store the result
		if c.PIdent >= minPIdent && c.AlignedFraction >= minQcovHSP {
			if ls, ok = (*ma)[uint32(ia)]; !ok {
				ls = poolChains2.Get().(*[]*Chain2Result)
				*ls = (*ls)[:0]
				(*ma)[uint32(ia)] = ls
			}
			*ls = append(*ls, c)

			c.Score = int(ia)    // for finding the corresponding subject fragment in the reciprocal comparison
			c.BitScore = int(ib) // for finding the corresponding subject fragment in the reciprocal comparison

			if ls, ok = (*mb)[uint32(ib)]; !ok {
				ls = poolChains2.Get().(*[]*Chain2Result)
				*ls = (*ls)[:0]
				(*mb)[uint32(ib)] = ls
			}
			*ls = append(*ls, c)
		} else {
			recycleChain2(c)
		}

		wfa.RecycleAlignmentResult(cigar)
		RecycleSeqComparatorResult(cr)
	}

	// 8) Identify orthologous fragments (reciprocal best hits)
	fsort := func(a, b *Chain2Result) int {
		// c.Evalue = c.AlignedFraction * c.PIdent // just for sorting, not the real e-value
		if d := cmp.Compare(b.Evalue, a.Evalue); d != 0 {
			return d
		}
		// c.Score = int(ia)    // for finding the corresponding subject fragment in the reciprocal comparison
		// c.BitScore = int(ib) // for finding the corresponding subject fragment in the reciprocal comparison
		if d := cmp.Compare(a.Score, b.Score); d != 0 {
			return d
		}
		return cmp.Compare(a.BitScore, b.BitScore)
	}

	for _, ls = range *ma {
		if len(*ls) > 1 {
			slices.SortFunc(*ls, fsort)
		}
	}

	for _, ls = range *mb {
		if len(*ls) > 1 {
			slices.SortFunc(*ls, fsort)
		}
	}

	// 9) Calculate ANI result from reciprocal best hits
	gr := poolGSearchResult.Get().(*GSearchResult)
	gr.Reset()
	gr.BatchGenomeIndex = 0 // Not from index
	gr.GenomeSize = subject.genomeSize
	gr.NumSeqs = len(subject.seqs)

	var _ia, _ib uint32
	var ls2 *[]*Chain2Result
	for _ia, ls = range *ma {
		_ib = uint32((*ls)[0].BitScore)

		if ls2, ok = (*mb)[_ib]; !ok {
			continue
		}
		if (*ls2)[0].Score != int(_ia) {
			continue
		}

		// reciprocal best hit
		c = (*ls)[0]

		gr.AlignedFragments++
		gr.AlignedMatches += c.MatchedBases
		gr.PidentsSum += c.PIdent
		gr.AlignedLength += c.AlignedLength - c.Gaps
	}

	// 10) Calculate final ANI / AF
	if gr.AlignedFragments > 0 {
		gr.ANI = gr.PidentsSum / float64(gr.AlignedFragments) / 100
	}
	gr.AFq = float64(gr.AlignedLength) / float64(qfragLens)
	gr.AFs = float64(gr.AlignedLength) / float64(sfragLens)
	if gr.AFq > 1 {
		gr.AFq = 1
	}
	if gr.AFs > 1 {
		gr.AFs = 1
	}
	gr.Score = gr.ANI

	// 11) Store result
	rs := poolGSearchResults.Get().(*[]*GSearchResult)
	*rs = (*rs)[:0]
	if gr.AFq >= minAF && gr.ANI >= minANI {
		*rs = append(*rs, gr)
	} else {
		poolGSearchResult.Put(gr)
		RecycleGSearchResults(rs)
		return nil
	}

	query.result = rs
	return nil
}

// ReadGenome reads a genome from a sequence file
func ReadGenomeFromFile(file string, reRefName *regexp.Regexp, fullPathAsRefName bool) (*GQuery, error) {
	fastxReader, err := fastx.NewDefaultReader(file)
	if err != nil {
		return nil, err
	}
	defer fastxReader.Close()

	q := poolGQuery.Get().(*GQuery)
	q.Reset()

	var record *fastx.Record

	for {
		record, err = fastxReader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			RecycleGQuery(q)
			return nil, fmt.Errorf("read seq in %s: %s", file, err)
		}

		s := poolSeq.Get().(*[]byte)
		*s = (*s)[:0]
		*s = append(*s, record.Seq.Seq...)
		q.seqs = append(q.seqs, s)

		q.genomeSize += len(record.Seq.Seq)
	}

	if q.genomeSize == 0 { // no sequence
		RecycleGQuery(q)
		return nil, nil
	}

	baseFile := filepath.Base(file)
	var genomeID string
	if fullPathAsRefName {
		genomeID = file
	} else if reRefName != nil {
		if reRefName.MatchString(baseFile) {
			genomeID = reRefName.FindAllStringSubmatch(baseFile, 1)[0][1]
		} else {
			genomeID, _, _ = filepathTrimExtension(baseFile, nil)
		}
	} else {
		genomeID, _, _ = filepathTrimExtension(baseFile, nil)
	}
	q.id = append(q.id, []byte(genomeID)...)

	return q, nil
}
