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
	"math"
	"math/bits"
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
	"github.com/shenwei356/kmers"
	"github.com/shenwei356/lexichash"
	"github.com/shenwei356/lexichash/iterator"
	"github.com/shenwei356/wfa"
)

// In-memory sketch of a subject genome. We only need it long enough to
// process all query fragments against this subject, then it is recycled.
// The seed layout mirrors what buildAnIndex writes to disk: primary k-mer
// per mask plus extra k-mers from desert filling, stored against the same
// per-mask slot. Positions are encoded as (pos << 1) | strand, matching
// the convention used by MaskKnownDistinctPrefixes.
type subjectSketch struct {
	seqLen     int
	kmers      *[]uint64    // primary k-mer per mask (0 if none)
	locses     *[][]int     // positions per mask, each entry = pos<<1 | strand
	extraKmers *[]*[]uint64 // extras per mask, intermittent (kmer, pos<<1|strand)

	// suffix-side store: for every captured subject k-mer (primary or extra),
	// we also reverse it and route the reversed k-mer to its MaskKmer-best
	// mask. Probed by reversing query k-mers, mirroring the suffix index
	// written by buildAnIndex. Layout: per-mask intermittent (revKmer, pos<<1|strand).
	revKmers *[]*[]uint64

	// contigBounds gives [start, end) of each contig in the forward concat
	// coordinate frame (boundaries are the contigInterval N-runs). Used by
	// alignQueryFragToSubject to clamp extendMatch's left/right extension to
	// the chain's contig instead of letting it cross N-spacers.
	contigBounds [][2]int
}

var poolSubjectSketch = &sync.Pool{New: func() interface{} {
	return &subjectSketch{}
}}

// poolConcat is for reusing large byte slices for concatenated genome sequences
var poolConcat = &sync.Pool{New: func() interface{} {
	tmp := make([]byte, 0, 10<<20) // 10MB initial capacity
	return &tmp
}}

// Defaults for in-memory desert filling. They match the typical index-build
// defaults (--seed-max-desert / --seed-in-desert-dist) so a query fragment
// sees similar seed density to what the persistent index would have.
var gsa3DesertMaxLen = 60
var gsa3DesertExpectedSeedDist = 30

// buildSubjectSketch masks the concatenated subject with the loaded LexicHash,
// drops low-complexity k-mers, and (optionally) fills sketching deserts so
// that no large region is left without a seed. The returned sketch must be
// freed with recycleSubjectSketch.
func (idx *Index) buildSubjectSketch(seq []byte, skipRegions [][2]int, contigBounds [][2]int, fillDeserts bool) (*subjectSketch, error) {
	lh := idx.lh
	k := lh.K
	k8 := uint8(k)

	_kmers, locses, err := lh.MaskKnownDistinctPrefixes(seq, skipRegions, true)
	if err != nil {
		return nil, err
	}

	// drop low-complexity captures so they don't poison anchor lists.
	ccc := util.Ns(0b01, k8)
	ggg := util.Ns(0b10, k8)
	ttt := (uint64(1) << (k << 1)) - 1
	for i, kmer := range *_kmers {
		if kmer == 0 {
			continue
		}
		if kmer == ccc || kmer == ggg || kmer == ttt ||
			util.IsLowComplexityDust(kmer, k8) {
			(*_kmers)[i] = 0
			(*locses)[i] = (*locses)[i][:0]
		}
	}

	// per-mask extras (desert fill output). Always allocate so downstream
	// code can iterate uniformly even when filling is disabled.
	extras := make([]*[]uint64, len(*_kmers))

	if fillDeserts {
		idx.fillSeedDesertsInMemory(seq, _kmers, locses, &extras, skipRegions)
	}

	// build suffix-side store: reverse every captured k-mer and route it to
	// the MaskKmer-best mask, exactly as buildAnIndex does on disk. This
	// lets the matcher find seed hits where the shared region is at the
	// SUFFIX of the (forward) k-mer rather than the prefix.
	revExtras := make([]*[]uint64, len(*_kmers))
	for i, kmer := range *_kmers {
		if kmer == 0 || len((*locses)[i]) == 0 {
			continue
		}
		rk := kmers.MustReverse(kmer, k)
		newMask := bestMaskForKmer(lh, rk)
		for _, loc := range (*locses)[i] {
			appendExtra(&revExtras, newMask, rk, uint64(loc))
		}
	}
	for i, knl := range extras {
		if knl == nil {
			continue
		}
		_ = i
		for j := 0; j+1 < len(*knl); j += 2 {
			kmer := (*knl)[j]
			if kmer == 0 {
				continue
			}
			rk := kmers.MustReverse(kmer, k)
			newMask := bestMaskForKmer(lh, rk)
			appendExtra(&revExtras, newMask, rk, (*knl)[j+1])
		}
	}

	s := poolSubjectSketch.Get().(*subjectSketch)
	s.seqLen = len(seq)
	s.kmers = _kmers
	s.locses = locses
	s.extraKmers = &extras
	s.revKmers = &revExtras
	s.contigBounds = contigBounds

	// Sort k-mers within each mask for faster prefix matching.
	for _, knl := range extras {
		if knl != nil && len(*knl) > 2 {
			sortKmerPairs(knl)
		}
	}
	for _, knl := range revExtras {
		if knl != nil && len(*knl) > 2 {
			sortKmerPairs(knl)
		}
	}

	return s, nil
}

// bestMaskForKmer returns the index of the mask that captures kmer with the
// smallest XOR (the same selection rule used by buildAnIndex's suffix index).
func bestMaskForKmer(lh *lexichash.LexicHash, kmer uint64) int {
	iMasks := lh.MaskKmer(kmer)
	defer lh.RecycleMaskKmerResult(iMasks)
	var minj int
	var minh uint64 = math.MaxUint64
	for _, j := range *iMasks {
		if h := lh.Masks[j] ^ kmer; h < minh {
			minj, minh = j, h
		}
	}
	return minj
}

// recycleSubjectSketch returns all per-sketch buffers back to their pools.
func (idx *Index) recycleSubjectSketch(s *subjectSketch) {
	if s == nil {
		return
	}
	if s.kmers != nil {
		idx.lh.RecycleMaskResult(s.kmers, s.locses)
		s.kmers = nil
		s.locses = nil
	}
	if s.extraKmers != nil {
		for i, knl := range *s.extraKmers {
			if knl != nil {
				*knl = (*knl)[:0]
				poolKmerAndLocs.Put(knl)
				(*s.extraKmers)[i] = nil
			}
		}
		s.extraKmers = nil
	}
	if s.revKmers != nil {
		for i, knl := range *s.revKmers {
			if knl != nil {
				*knl = (*knl)[:0]
				poolKmerAndLocs.Put(knl)
				(*s.revKmers)[i] = nil
			}
		}
		s.revKmers = nil
	}
	poolSubjectSketch.Put(s)
}

// fillSeedDesertsInMemory walks the (already-computed) per-mask seed
// positions, finds spans longer than gsa3DesertMaxLen, and inserts extra
// (k-mer, pos|strand) pairs into the corresponding masks. This is a
// trimmed-down, in-memory variant of the desert-filling pass in buildAnIndex.
func (idx *Index) fillSeedDesertsInMemory(
	seq []byte,
	_kmers *[]uint64,
	locses *[][]int,
	extras *[]*[]uint64,
	skipRegions [][2]int,
) {
	lh := idx.lh
	k := lh.K
	k8 := uint8(k)
	ccc := util.Ns(0b01, k8)
	ggg := util.Ns(0b10, k8)
	ttt := (uint64(1) << (k << 1)) - 1

	// gather all primary seed positions (pos with strand bit) and sort.
	locs := poolInts.Get().(*[]int)
	*locs = (*locs)[:0]
	defer func() {
		*locs = (*locs)[:0]
		poolInts.Put(locs)
	}()
	for _, ls := range *locses {
		*locs = append(*locs, ls...)
	}
	if len(*locs) == 0 {
		return
	}
	slices.Sort(*locs)

	// interval tree of regions we must not place seeds in (contig joins +
	// gaps), so we don't fabricate anchors that span structural breaks.
	tree := itree.NewSearchTree[uint8, int](cmpFn)
	for _, r := range skipRegions {
		tree.Insert(r[0]-int(k)+1, r[1], 1)
	}

	maxDesert := uint32(gsa3DesertMaxLen)
	seedDist := gsa3DesertExpectedSeedDist
	seedPosR := seedDist / 2

	lenSeq := len(seq)

	// scratch buffers (mirroring buildAnIndex's desert pass)
	loc2maskidx := poolLoc2MaskIdx.Get().(*[]int)
	loc2maskidxRC := poolLoc2MaskIdx.Get().(*[]int)
	kmerList := poolKmerKmerRC.Get().(*[]uint64)
	defer func() {
		poolLoc2MaskIdx.Put(loc2maskidx)
		poolLoc2MaskIdx.Put(loc2maskidxRC)
		poolKmerKmerRC.Put(kmerList)
	}()

	// append a pseudo position at the end so the trailing desert is also
	// processed (same trick as buildAnIndex).
	*locs = append(*locs, (lenSeq-int(k))<<1)

	var pre, pos uint32
	pre = 0
	for _, l := range *locs {
		pos = uint32(l) >> 1
		d := pos - pre

		if d < maxDesert {
			pre = pos
			continue
		}

		// scan window around the desert
		start := int(pre) - 1000
		posOfPre := 1000
		if start < 0 {
			posOfPre += start
			start = 0
		}
		end := int(pos) + 1000 + int(k)
		if end > lenSeq {
			end = lenSeq
		}
		posOfCur := posOfPre + int(d)

		// k-mers across the window
		iter, err := iterator.NewKmerIterator(seq[start:end], int(k))
		if err != nil {
			pre = pos
			continue
		}
		*kmerList = (*kmerList)[:0]
		for {
			kmer, kmerRC, ok, _ := iter.NextKmer()
			if !ok {
				break
			}
			*kmerList = append(*kmerList, kmer)
			*kmerList = append(*kmerList, kmerRC)
		}

		// mask the window so we know which mask captures each position.
		_kmers2, _locses2, _ := lh.MaskKnownDistinctPrefixes(seq[start:end], nil, false)

		// loc -> capturing mask, separately for fwd and rc.
		n := end - start
		if cap(*loc2maskidx) < n {
			*loc2maskidx = make([]int, n)
		} else {
			*loc2maskidx = (*loc2maskidx)[:n]
		}
		if cap(*loc2maskidxRC) < n {
			*loc2maskidxRC = make([]int, n)
		} else {
			*loc2maskidxRC = (*loc2maskidxRC)[:n]
		}
		for i := 0; i < n; i++ {
			(*loc2maskidx)[i] = -1
			(*loc2maskidxRC)[i] = -1
		}
		for im, lsW := range *_locses2 {
			if (*_kmers2)[im] == 0 {
				continue
			}
			for _, lw := range lsW {
				p := lw >> 1
				if p < 0 || p >= n {
					continue
				}
				if lw&1 == 0 {
					(*loc2maskidx)[p] = im
				} else {
					(*loc2maskidxRC)[p] = im
				}
			}
		}
		lh.RecycleMaskResult(_kmers2, _locses2)

		// step from previous seed, planting one extra every seedDist.
		j := posOfPre + seedDist
		for j < posOfCur {
			// upstream scan: [j-seedPosR, j]
			startScan := j + 1
			endScan := j - seedPosR
			ok := false
			var kmer, kmerPos uint64
			var im int
			for ; j > endScan; j-- {
				if j < 0 || j >= n {
					continue
				}
				if _, inGap := tree.AnyIntersection(start+j, start+j); inGap {
					continue
				}
				// +strand
				kmer = (*kmerList)[j<<1]
				if kmer != 0 && kmer != ccc && kmer != ggg && kmer != ttt &&
					!util.IsLowComplexityDust(kmer, k8) {
					if im = (*loc2maskidx)[j]; im >= 0 {
						kmerPos = uint64(start+j) << 1
						ok = true
						break
					}
				}
				// -strand
				kmer = (*kmerList)[(j<<1)+1]
				if kmer != 0 && kmer != ccc && kmer != ggg && kmer != ttt &&
					!util.IsLowComplexityDust(kmer, k8) {
					if im = (*loc2maskidxRC)[j]; im >= 0 {
						kmerPos = uint64(start+j)<<1 | 1
						ok = true
						break
					}
				}
			}
			if ok {
				appendExtra(extras, im, kmer, kmerPos)
				j += seedDist
				continue
			}

			if startScan >= posOfCur {
				break
			}
			// downstream scan: [startScan, startScan+seedPosR]
			endScan = startScan + seedPosR
			if endScan >= posOfCur {
				endScan = posOfCur - 1
			}
			for j = startScan; j < endScan; j++ {
				if j < 0 || j >= n {
					continue
				}
				if _, inGap := tree.AnyIntersection(start+j, start+j); inGap {
					continue
				}
				kmer = (*kmerList)[j<<1]
				if kmer != 0 && kmer != ccc && kmer != ggg && kmer != ttt &&
					!util.IsLowComplexityDust(kmer, k8) {
					if im = (*loc2maskidx)[j]; im >= 0 {
						kmerPos = uint64(start+j) << 1
						ok = true
						break
					}
				}
				kmer = (*kmerList)[(j<<1)+1]
				if kmer != 0 && kmer != ccc && kmer != ggg && kmer != ttt &&
					!util.IsLowComplexityDust(kmer, k8) {
					if im = (*loc2maskidxRC)[j]; im >= 0 {
						kmerPos = uint64(start+j)<<1 | 1
						ok = true
						break
					}
				}
			}
			if ok {
				appendExtra(extras, im, kmer, kmerPos)
				j += seedDist
				continue
			}

			// could not plant a seed here (gap / interval / low-complexity);
			// just advance and try the next slot.
			j += seedDist
		}

		pre = pos
	}

	// remove the pseudo position we appended.
	*locs = (*locs)[:len(*locs)-1]
}

func appendExtra(extras *[]*[]uint64, mask int, kmer, kmerPos uint64) {
	knl := (*extras)[mask]
	if knl == nil {
		knl = poolKmerAndLocs.Get().(*[]uint64)
		*knl = (*knl)[:0]
		(*extras)[mask] = knl
	}
	*knl = append(*knl, kmer, kmerPos)
}

// sortKmerPairs sorts (kmer, pos) pairs by kmer value for faster prefix matching.
// The slice format is [kmer0, pos0, kmer1, pos1, ...].
func sortKmerPairs(knl *[]uint64) {
	n := len(*knl) / 2
	if n <= 1 {
		return
	}
	type kmerPair struct {
		kmer uint64
		pos  uint64
	}
	pairs := make([]kmerPair, n)
	for i := 0; i < n; i++ {
		pairs[i] = kmerPair{(*knl)[i*2], (*knl)[i*2+1]}
	}
	slices.SortFunc(pairs, func(a, b kmerPair) int {
		if a.kmer < b.kmer {
			return -1
		} else if a.kmer > b.kmer {
			return 1
		}
		return 0
	})
	for i := 0; i < n; i++ {
		(*knl)[i*2] = pairs[i].kmer
		(*knl)[i*2+1] = pairs[i].pos
	}
}

// GSearchAlign3 aligns query fragments against the full (concatenated)
// sequence of each candidate subject genome via an in-memory LexicHash
// sketch. Unlike GSearchAlign2 — which fragments BOTH sides and keeps only
// reciprocal-best-hit fragment pairs (OrthoANI style) — GSearchAlign3 only
// fragments the query and looks for the best chain of each fragment against
// the whole subject (skani / ANIb style).
//
// Per candidate subject we:
//   - concatenate its contigs with a 1.5*fragLen run of N's so a single
//     query fragment cannot straddle a contig boundary (the user picked
//     this margin so chain post-filtering can stay simple),
//   - build an in-memory LexicHash sketch (with desert filling) of that
//     concatenation,
//   - for every query fragment: match seeds, chain anchors twice (forward
//     and reverse), extend & WFA-align the best chain, then accumulate
//     matched/aligned-length stats.
//
// ANI = aligned_matches / aligned_length (sum over fragments),
// AF  = aligned_length / total_query_fragment_length.
func (idx *Index) GSearchAlign3(query *GQuery, fragLen int, minFragLen int, genomeIds *map[uint64]*[]uint64, minAF float64, maxQueryConcurrency int, gcInterval uint64) error {
	debug := idx.opt.Debug

	startTime0 := time.Now()

	if debug {
		log.Debugf("%s (%s bp): start to preprocess query genome fragments", query.id, humanize.Comma(int64(query.genomeSize)))
	}

	// 1) Cut the query into fragments.
	qfrags, qfragLens := seqs2fragments(&query.seqs, fragLen, minFragLen)
	if len(*qfrags) == 0 {
		return fmt.Errorf("no fragments for alignment, are the genome too fragmented with all sequences shorter than the fragment size?")
	}
	defer recycleFragments(qfrags)

	// 2) Pre-compute each query fragment's masked k-mers once. They are
	// read-only thereafter and shared across all candidate subjects.

	qSeeds := make([]*qFragSeeds, len(*qfrags))
	{
		var wgPre sync.WaitGroup
		tokensPre := make(chan int, maxQueryConcurrency)
		for i := range *qfrags {
			tokensPre <- 1
			wgPre.Add(1)
			go func(i int) {
				defer func() { <-tokensPre; wgPre.Done() }()
				frag := (*qfrags)[i]
				kmers, locses, err := idx.lh.MaskKnownDistinctPrefixes(frag, nil, true)
				if err != nil {
					checkError(err)
				}

				// strip low-complexity captures.
				k8 := uint8(idx.lh.K)
				ccc := util.Ns(0b01, k8)
				ggg := util.Ns(0b10, k8)
				ttt := (uint64(1) << (idx.lh.K << 1)) - 1
				for j, kmer := range *kmers {
					if kmer == 0 {
						continue
					}
					if kmer == ccc || kmer == ggg || kmer == ttt ||
						util.IsLowComplexityDust(kmer, k8) {
						(*kmers)[j] = 0
						(*locses)[j] = (*locses)[j][:0]
					}
				}
				qSeeds[i] = &qFragSeeds{kmers: kmers, locses: locses}
			}(i)
		}
		wgPre.Wait()
	}
	defer func() {
		for _, s := range qSeeds {
			if s != nil {
				idx.lh.RecycleMaskResult(s.kmers, s.locses)
			}
		}
	}()

	// 3) For each candidate subject, chunked or not, drop secondary chunk
	// keys so we process each genome exactly once (same as GSearchAlign2).

	startTime := time.Now()
	if debug {
		log.Debugf("%s (%s bp): finished preprocessing query genome fragments in %.3f seconds",
			query.id, humanize.Comma(int64(query.genomeSize)), time.Since(startTime0).Seconds())
		log.Debugf("%s (%s bp): start to align query genome fragments", query.id, humanize.Comma(int64(query.genomeSize)))
	}

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

	ch := make(chan *GSearchResult, maxQueryConcurrency)
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

	K := idx.lh.K
	contigInterval := int(float64(fragLen) * 1.5)
	if contigInterval < K {
		contigInterval = K
	}
	nnn := bytes.Repeat([]byte{'N'}, contigInterval)
	reGaps := regexp.MustCompile(fmt.Sprintf(`[Nn]{%d,}`, 5))

	alignOption := &wfa.Options{GlobalAlignment: true}
	var minPrefix uint8 = idx.seqCompareOption.MinPrefix
	minPIdent := idx.seqCompareOption.MinIdentity
	minQcovHSP := idx.seqCompareOption.MinAlignedFraction
	extLen := fragLen / 2 // idx.opt.ExtendLength
	extLen2 := idx.opt.ExtendLength2

	var wg sync.WaitGroup
	tokens := make(chan int, maxQueryConcurrency)

	for _, batchIDAndRefIDs := range *genomeIds {
		tokens <- 1
		wg.Add(1)

		go func(batchIDAndRefIDs *[]uint64) {
			var g *genome.Genome
			timeStart := time.Now()
			defer func() {
				<-tokens
				wg.Done()
				if debug {
					// log.Debugf("%s (%s bp): aligning subject genome %s (%s bp) took %s",
					// 	query.id, humanize.Comma(int64(query.genomeSize)),
					// 	idx.BatchGenomeIndex2GenomeID[(*batchIDAndRefIDs)[0]],
					// 	humanize.Comma(int64(g.GenomeSize)), time.Since(timeStart))
					// log.Debug()
					chDuration <- time.Duration(float64(time.Since(timeStart)) / fcpus)
				}
			}()

			// a) Load all chunks of this genome.
			genomes := make([]*genome.Genome, len(*batchIDAndRefIDs))
			maxSubjectGenomeSize := idx.opt.MaxSubjectGenomeSize
			for i, batchIDAndRefID := range *batchIDAndRefIDs {
				genomeBatch := int(batchIDAndRefID >> BITS_GENOME_IDX)
				genomeIdx := int(batchIDAndRefID & MASK_GENOME_IDX)

				rdr := <-idx.poolGenomeRdrs[genomeBatch]
				_g, err := rdr.Seqs(genomeIdx)
				if err != nil {
					checkError(fmt.Errorf("fail to read genome sequence for batch %d, genome index %d: %s",
						genomeBatch, genomeIdx, err))
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
					log.Warningf("skipped subject genome %s (%s bp) which exceeds the maximum allowed size of %s, consider increasing --max-subject-genome-size",
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

			// b) Concatenate contigs with a 1.5*fragLen separator and
			// record the (per-N-block) skip regions for masking.
			concat := poolConcat.Get().(*[]byte)
			*concat = (*concat)[:0]
			var concatSize int
			for _, s := range g.Seqs {
				concatSize += len(*s)
			}
			concatSize += contigInterval * (len(g.Seqs) - 1)
			if concatSize < 0 {
				concatSize = 0
			}
			if cap(*concat) < concatSize {
				*concat = make([]byte, 0, concatSize)
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

			// skip gap regions (N's)
			gaps := reGaps.FindAllSubmatchIndex(*concat, -1)
			if gaps != nil {
				for _, gap := range gaps {
					skipRegions = append(skipRegions, [2]int{gap[0], gap[1] - 1})
				}
				slices.SortFunc(skipRegions, func(a, b [2]int) int {
					return a[0] - b[0]
				})
			}

			// c) Build the subject sketch.
			sketch, err := idx.buildSubjectSketch(*concat, skipRegions, contigBounds, true)
			if err != nil {
				checkError(fmt.Errorf("fail to build subject sketch: %s", err))
			}

			// d) Reverse complement of the concatenated subject. Built once
			// up-front because most queries will exercise both strands; the
			// extra ~|concat| bytes are cheaper than re-RCing per fragment.
			concatRC := poolConcat.Get().(*[]byte)
			if cap(*concatRC) < len(*concat) {
				*concatRC = make([]byte, len(*concat))
			} else {
				*concatRC = (*concatRC)[:len(*concat)]
			}
			copy(*concatRC, *concat)
			RC(*concatRC)

			// e) Set up per-subject scratch (chainer, WFA aligner).
			chainer := idx.poolChainers2.Get().(*Chainer2)
			algn := wfa.New(wfa.DefaultPenalties, alignOption)
			algn.AdaptiveReduction(wfa.DefaultAdaptiveOption)

			gr := poolGSearchResult.Get().(*GSearchResult)
			gr.Reset()
			gr.BatchGenomeIndex = (*batchIDAndRefIDs)[0]
			gr.GenomeSize = g.GenomeSize
			gr.NumSeqs = g.NumSeqs

			// f) Align each query fragment, keeping only the best chain.
			fScoreAndEvalue := scoreAndEvalue(2, -3, 5, 2, int(g.GenomeSize), 0.625, 0.41)

			for i, qfrag := range *qfrags {
				matched, alignedLen, gaps, pident, ok := alignQueryFragToSubject(
					qfrag, qSeeds[i], sketch, (*concat), (*concatRC),
					chainer, algn, minPrefix, K, extLen, extLen2,
					minPIdent, minQcovHSP, idx,
					&fScoreAndEvalue,
				)
				if !ok {
					continue
				}
				gr.AlignedFragments++
				gr.AlignedLength += alignedLen - gaps
				gr.AlignedMatches += matched
				gr.Pidents = append(gr.Pidents, pident)
			}

			// g) ANI / AF on the accumulated alignment.
			if gr.AlignedLength > 0 {
				// gr.ANI = float64(gr.AlignedMatches) / float64(gr.AlignedLength) // shouldn't do this
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
			// fmt.Printf("aligned fragments: %d, aligned length: %d, total length: %d, matched bases: %d, ANI: %.2f, AF: %.2f\n", gr.AlignedFragments, gr.AlignedLength, qfragLens, gr.AlignedMatches, gr.ANI, gr.AF)

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

			// Recycle concat and concatRC
			*concat = (*concat)[:0]
			poolConcat.Put(concat)
			*concatRC = (*concatRC)[:0]
			poolConcat.Put(concatRC)

		}(batchIDAndRefIDs)
	}

	wg.Wait()
	close(ch)
	<-done

	// process bar
	if debug {
		close(chDuration)
		<-doneDuration
		pbs.Wait()
		log.Debugf("%s (%s bp): finished aligning query genome fragments in %.3f seconds",
			query.id, humanize.Comma(int64(query.genomeSize)), time.Since(startTime).Seconds())
	}

	return nil
}

// qFragSeeds holds the LexicHash output for a single query fragment so we
// can reuse it across every candidate subject without re-masking.
type qFragSeeds struct {
	kmers  *[]uint64
	locses *[][]int
}

// alignQueryFragToSubject finds the single best chain of one query fragment
// against the in-memory subject sketch and runs WFA on it.
//
// Returns (matchedBases, alignedColumns, gaps, true) on a hit that passes
// the identity / qcov filters, otherwise (_, _, _, false).
//
// Anchors are split by relative strand: forward anchors (qrc == src) chain
// in (q-forward, s-forward) coords; reverse anchors (qrc != src) chain in
// (q-forward, s-rc) coords against concatRC. Within each group we pick the
// highest-scoring chain.
func alignQueryFragToSubject(
	qfrag []byte,
	qSeeds *qFragSeeds,
	sketch *subjectSketch,
	concat, concatRC []byte,
	chainer *Chainer2,
	algn *wfa.Aligner,
	minPrefix uint8,
	K int,
	extLen int,
	extLen2 int,
	minPIdent float64,
	minQcov float64,
	idx *Index,
	fScoreAndEvalue *func(qlen int, cigar *wfa.AlignmentResult) (int, int, float64),
) (int, int, int, float64, bool) {
	subjectLen := sketch.seqLen
	K8 := uint8(K)

	fwdSubs := poolSubsLong.Get().(*[]*SubstrPair)
	*fwdSubs = (*fwdSubs)[:0]
	revSubs := poolSubsLong.Get().(*[]*SubstrPair)
	*revSubs = (*revSubs)[:0]
	defer RecycleSubstrPairs(poolSub, poolSubsLong, fwdSubs)
	defer RecycleSubstrPairs(poolSub, poolSubsLong, revSubs)

	qKmers := *qSeeds.kmers
	qLocses := *qSeeds.locses
	sKmers := *sketch.kmers
	sLocses := *sketch.locses
	sExtras := *sketch.extraKmers
	sRevKmers := *sketch.revKmers
	shift := K8 - 32

	// build anchors per mask. Each mask carries one primary subject k-mer
	// (with possibly multiple positions) plus zero or more desert-fill
	// extras. The query side only has its primary.
	for i, qk := range qKmers {
		if qk == 0 {
			continue
		}
		qLocs := qLocses[i]
		if len(qLocs) == 0 {
			continue
		}

		// primary subject k-mer.
		if sk := sKmers[i]; sk != 0 {
			if plen := uint8(bits.LeadingZeros64(qk^sk)>>1) + shift; plen >= minPrefix {
				for _, ql := range qLocs {
					for _, sl := range sLocses[i] {
						addAnchor(fwdSubs, revSubs, ql, sl, plen, K, subjectLen)
					}
				}
			}
		}

		// desert-fill extras for this mask.
		if ek := sExtras[i]; ek != nil {
			hitSome := false
			for j := 0; j+1 < len(*ek); j += 2 {
				xk := (*ek)[j]
				if plen := uint8(bits.LeadingZeros64(qk^xk)>>1) + shift; plen >= minPrefix {
					xloc := int((*ek)[j+1])
					for _, ql := range qLocs {
						addAnchor(fwdSubs, revSubs, ql, xloc, plen, K, subjectLen)
					}

					hitSome = true
				} else if hitSome {
					//fmt.Printf("break at %d/%d\n", j, len(*ek))
					break
				}
			}
		}

		// SUFFIX matching: reverse the query k-mer, find its best mask,
		// and probe sketch.revKmers there. A prefix match between qkRev
		// and a stored revSubjectKmer with length plen means the LAST
		// plen bases of qk and the original sk agree. Reversing both
		// k-mers flips both strand interpretations, so we can re-use
		// addAnchor by XOR-ing the strand bits of both locs.
		qkRev := kmers.MustReverse(qk, K)
		jq := bestMaskForKmer(idx.lh, qkRev)
		if rk := sRevKmers[jq]; rk != nil {
			hitSome := false
			for j := 0; j+1 < len(*rk); j += 2 {
				sRev := (*rk)[j]
				if plen := uint8(bits.LeadingZeros64(qkRev^sRev)>>1) + shift; plen >= minPrefix {
					sloc := int((*rk)[j+1])
					for _, ql := range qLocs {
						addAnchor(fwdSubs, revSubs, ql^1, sloc^1, plen, K, subjectLen)
					}

					hitSome = true
				} else if hitSome {
					// fmt.Printf("break at %d/%d\n", j, len(*rk))
					break
				}
			}
		}
	}

	fwdChains, fwdOk := chainsFromSubs(fwdSubs, chainer, K)
	revChains, revOk := chainsFromSubs(revSubs, chainer, K)

	if !fwdOk && !revOk {
		return 0, 0, 0, 0, false
	}

	// Try all chains (both forward and reverse) and pick the one with
	// the best WFA alignment result (matches * aligned_length). For
	// multi-copy conserved genes, multiple chains may look identical
	// in chain-space but differ at the sequence level.
	var bestMatched, bestAligned, bestGaps int
	var bestPident float64
	var bestScore int = -1
	topChains := idx.chainingOptions.TopChains // only check the best N chains
	onlyTopChains := topChains > 0

	if fwdOk {
		i := 0
		for _, chain := range *fwdChains {
			if chain == nil {
				continue
			}
			i++
			if onlyTopChains && i > topChains {
				break
			}
			matched, aligned, gaps, pident, ok := alignChain(
				qfrag, concat, chain, sketch, false, algn,
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
		RecycleChaining2Result(fwdChains)
	}

	if revOk {
		i := 0
		for _, chain := range *revChains {
			if chain == nil {
				continue
			}
			i++
			if onlyTopChains && i > topChains {
				break
			}
			matched, aligned, gaps, pident, ok := alignChain(
				qfrag, concatRC, chain, sketch, true, algn,
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
		RecycleChaining2Result(revChains)
	}

	if bestScore <= 0 {
		return 0, 0, 0, 0, false
	}

	return bestMatched, bestAligned, bestGaps, bestPident, true
}

// addAnchor emits one SubstrPair from a (qLoc, sLoc) match. The two locs are
// (pos << 1) | strand_bit; relative strand tells us whether the alignment
// is forward (both strands the same) or reverse (different strands). For
// reverse anchors we re-express the subject position in RC-subject
// coordinates so the chainer can treat both groups uniformly.
func addAnchor(fwdSubs, revSubs *[]*SubstrPair, qLoc, sLoc int, plen uint8, K, subjectLen int) {
	qpos := qLoc >> 1
	qrc := qLoc&1 == 1
	spos := sLoc >> 1
	src := sLoc&1 == 1
	plenInt := int(plen)

	sub := poolSub.Get().(*SubstrPair)
	sub.Len = plen
	sub.QRC = qrc
	sub.TRC = src

	if qrc == src {
		// same strand -> forward alignment.
		if !qrc {
			sub.QBegin = int32(qpos)
			sub.TBegin = int32(spos)
		} else {
			// both RC: the shared prefix sits at the tail of the
			// forward k-mer, so shift both starts by K - plen.
			sub.QBegin = int32(qpos + K - plenInt)
			sub.TBegin = int32(spos + K - plenInt)
		}
		*fwdSubs = append(*fwdSubs, sub)
		return
	}

	// mixed strand -> reverse alignment. Re-express the subject position
	// in RC-subject coordinates so the chainer can run forward over both
	// query and RC-subject.
	if !qrc {
		// qrc=0, src=1: forward prefix on q starts at qpos; on RC-subject
		// the same prefix starts at subjectLen - K - spos.
		sub.QBegin = int32(qpos)
		sub.TBegin = int32(subjectLen - K - spos)
	} else {
		// qrc=1, src=0: the shared prefix sits at the tail of the
		// forward q-k-mer (qpos+K-plen) and (in forward subject coords)
		// at the start of the s-k-mer (spos). On RC-subject coords the
		// prefix starts at subjectLen - spos - plen.
		sub.QBegin = int32(qpos + K - plenInt)
		sub.TBegin = int32(subjectLen - spos - plenInt)
	}
	*revSubs = append(*revSubs, sub)
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

	// Pseudo-alignment with SeqComparator.
	cpr := idx.poolSeqComparator.Get().(*SeqComparator)
	defer idx.poolSeqComparator.Put(cpr)

	if err := cpr.Index(qfrag); err != nil {
		return 0, 0, 0, 0, false
	}
	defer cpr.RecycleIndex()

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
