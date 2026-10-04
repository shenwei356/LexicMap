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
	"fmt"
	"io"
	"math"
	"path/filepath"
	"regexp"
	"slices"
	"sync"

	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/lexichash"
	"gonum.org/v1/gonum/stat/distuv"
)

// GQuery represents a genome query
type GQuery struct {
	id          []byte
	bigSeq      []byte
	seqs        []*[]byte
	quals       []*[]byte
	skipRegions []int

	genomeSize int

	result *[]*GSearchResult // fragment alignment results

	screenDetails *[]*GSearchScreenResultDetail
}

var poolGQuery = &sync.Pool{New: func() interface{} {
	return &GQuery{
		id:          make([]byte, 0, 127),
		bigSeq:      make([]byte, 0),
		seqs:        make([]*[]byte, 0, 256),
		quals:       make([]*[]byte, 0, 256),
		skipRegions: make([]int, 0, 512),
		genomeSize:  0,
	}
}}

var poolSeq = &sync.Pool{
	New: func() interface{} {
		tmp := make([]byte, 0, 10<<10) // 10K
		return &tmp
	},
}

var poolQual = &sync.Pool{
	New: func() interface{} {
		tmp := make([]byte, 0, 10<<10) // 10K
		return &tmp
	},
}

func (q *GQuery) Reset() {
	q.id = q.id[:0]
	q.bigSeq = q.bigSeq[:0]
	clear(q.seqs)
	q.seqs = q.seqs[:0]
	clear(q.quals)
	q.quals = q.quals[:0]
	q.skipRegions = q.skipRegions[:0]
	q.genomeSize = 0

	q.result = nil
	q.screenDetails = nil
}

func RecycleGQuery(q *GQuery) {
	q.id = q.id[:0]
	q.bigSeq = q.bigSeq[:0]
	if q.seqs != nil {
		for _, s := range q.seqs {
			*s = (*s)[:0]
			poolSeq.Put(s)
		}
		clear(q.seqs)
		q.seqs = q.seqs[:0]
	}
	if q.quals != nil {
		for _, qual := range q.quals {
			*qual = (*qual)[:0]
			poolQual.Put(qual)
		}
		clear(q.quals)
		q.quals = q.quals[:0]
	}
	if q.skipRegions != nil {
		q.skipRegions = q.skipRegions[:0]
	}
	q.genomeSize = 0

	if q.result != nil {
		RecycleGSearchResults(q.result)
		q.result = nil
	}
	q.screenDetails = nil

	poolGQuery.Put(q)
}

// --------------------------------------------------------------

// --------------------------------------------------------------

// GenomeReader is only for `lexicmap genome search`
type GenomeReader struct {
	k         int            // kmer size
	nnn       []byte         // Ns
	reRefName *regexp.Regexp // for extracting genome id from the file name
}

// NewGenomeReader returns a GenomeReader with given k-mer size
func NewGenomeReader(k int, reRefName *regexp.Regexp) *GenomeReader {
	return &GenomeReader{
		k:         k,
		nnn:       bytes.Repeat([]byte{'N'}, k),
		reRefName: reRefName,
	}
}

// Recycle recyle a GQuery object
func (gr *GenomeReader) Recycle(q *GQuery) {
	RecycleGQuery(q)
}

// Read reads a genome from a file or stdin
func (gr *GenomeReader) Read(file string, convertNtoA bool, softMasking bool) (*GQuery, error) {
	return gr.read(file, convertNtoA, softMasking, false)
}

// ReadWithQual reads a genome and retains its Phred+33 qualities.
func (gr *GenomeReader) ReadWithQual(file string, convertNtoA bool, softMasking bool) (*GQuery, error) {
	return gr.read(file, convertNtoA, softMasking, true)
}

func (gr *GenomeReader) read(file string, convertNtoA bool, softMasking bool, readQual bool) (*GQuery, error) {
	fastxReader, err := fastx.NewDefaultReader(file)
	if err != nil {
		return nil, err
	}
	defer fastxReader.Close()

	q := poolGQuery.Get().(*GQuery)
	q.Reset()

	var record *fastx.Record
	i := 0

	var table [256]byte
	if convertNtoA {
		if softMasking {
			table = baseConvertCaseSensitive
		} else {
			table = baseConvert
		}
	}

	for {
		record, err = fastxReader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			RecycleGQuery(q)
			return nil, fmt.Errorf("read seq %d in %s: %s", i, file, err)
		}
		growBy := len(record.Seq.Seq)
		if i > 0 {
			growBy += len(gr.nnn)
		}
		q.bigSeq = slices.Grow(q.bigSeq, growBy)

		if i > 0 {
			q.skipRegions = append(q.skipRegions, len(q.bigSeq), len(q.bigSeq)+gr.k-1)

			q.bigSeq = append(q.bigSeq, gr.nnn...)
		}

		if convertNtoA {
			convertSeq(record.Seq.Seq, table)
		}
		if readQual {
			if !fastxReader.IsFastq || len(record.Seq.Qual) == 0 {
				RecycleGQuery(q)
				return nil, fmt.Errorf("quality-based ANI adjustment requires FASTQ input: %s", file)
			}
			if len(record.Seq.Qual) != len(record.Seq.Seq) {
				RecycleGQuery(q)
				return nil, fmt.Errorf("sequence and quality lengths differ for sequence %d in %s: %d != %d",
					i+1, file, len(record.Seq.Seq), len(record.Seq.Qual))
			}
			for _, v := range record.Seq.Qual {
				if v < 33 || v > 126 {
					RecycleGQuery(q)
					return nil, fmt.Errorf("invalid Phred+33 quality byte %d for sequence %d in %s", v, i+1, file)
				}
			}
		}

		s := poolSeq.Get().(*[]byte)
		*s = (*s)[:0]
		*s = append(*s, record.Seq.Seq...)
		q.seqs = append(q.seqs, s)
		if readQual {
			qual := poolQual.Get().(*[]byte)
			*qual = (*qual)[:0]
			*qual = append(*qual, record.Seq.Qual...)
			q.quals = append(q.quals, qual)
		}
		q.genomeSize += len(record.Seq.Seq)

		q.bigSeq = append(q.bigSeq, record.Seq.Seq...)
	}

	lenSeq := len(q.bigSeq)
	if lenSeq == 0 {
		RecycleGQuery(q)
		return nil, nil
	}

	gaps := findGapRegions(q.bigSeq, 5)
	if gaps != nil {
		for _, gap := range *gaps {
			start, end := unpackGapRegion(gap)
			q.skipRegions = append(q.skipRegions, start, end-1)
		}
		recycleGapRegions(gaps)

		lexichash.SortSkipRegions(q.skipRegions)
	}

	baseFile := filepath.Base(file)
	var genomeID string
	if gr.reRefName != nil {
		if gr.reRefName.MatchString(baseFile) {
			genomeID = gr.reRefName.FindAllStringSubmatch(baseFile, 1)[0][1]
		} else {
			genomeID, _, _ = filepathTrimExtension(baseFile, nil)
		}
	} else {
		genomeID, _, _ = filepathTrimExtension(baseFile, nil)
	}
	q.id = append(q.id, []byte(genomeID)...)

	return q, nil
}

// --------------------------------------------------------------

func convertSeq(seq []byte, table [256]byte) {
	for i, b := range seq {
		seq[i] = table[b]
	}
}

var baseConvert = [256]byte{
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'C', 'C', 'A', 'A', 'A', 'G', 'A', 'A', 'A', 'G', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'C', 'T', 'T', 'A', 'A', 'A', 'C', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'C', 'C', 'A', 'A', 'A', 'G', 'A', 'A', 'A', 'G', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'C', 'T', 'T', 'A', 'A', 'A', 'C', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
}

// baseConvertCaseSensitive converts all lower-cases to A (soft masking)
var baseConvertCaseSensitive = [256]byte{
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'C', 'C', 'A', 'A', 'A', 'G', 'A', 'A', 'A', 'G', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'C', 'T', 'T', 'A', 'A', 'A', 'C', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', // convert all lower cases to A
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', // convert all lower cases to A
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
}

// --------------------------------------------------------------

// cut genome sequences into non-overlapped fragments

var poolFragments = &sync.Pool{
	New: func() interface{} {
		tmp := make([][]byte, 0, 10240) // for a 10Mb genome
		return &tmp
	},
}

func recycleFragments(frags *[][]byte) {
	if frags != nil {
		// Fragments borrow query sequences; retain only the outer buffer.
		clear(*frags)
		*frags = (*frags)[:0]
		poolFragments.Put(frags)
	}
}

// do not forget to call recycleFragments with the non-nil result
func seqs2fragments(seqs *[]*[]byte, fragLen int, minFragLen int) (*[][]byte, int) {
	if seqs == nil || len(*seqs) == 0 {
		return nil, 0
	}

	frags := poolFragments.Get().(*[][]byte)

	var end, s, e int
	var contig *[]byte

	var n int

	for _, contig = range *seqs {
		end = len(*contig)
		for s = 0; s < end; s += fragLen {
			e = s + fragLen

			if e > end {
				e = end
				if e-s < minFragLen { // skip short fragments
					continue
				}
			}

			*frags = append(*frags, (*contig)[s:e])
			n += e - s
		}
	}

	return frags, n
}

// seqs2fragmentsWithQual cuts sequences and their parallel Phred+33 qualities
// using identical fragment boundaries. qualFrags is nil when quals is empty.
// The caller must recycle both non-nil fragment slices.
func seqs2fragmentsWithQual(seqs, quals *[]*[]byte, fragLen int, minFragLen int) (*[][]byte, *[][]byte, int, error) {
	if quals == nil || len(*quals) == 0 {
		frags, n := seqs2fragments(seqs, fragLen, minFragLen)
		return frags, nil, n, nil
	}
	if seqs == nil || len(*seqs) != len(*quals) {
		nSeqs := 0
		if seqs != nil {
			nSeqs = len(*seqs)
		}
		return nil, nil, 0, fmt.Errorf("sequence and quality record counts differ: %d != %d", nSeqs, len(*quals))
	}

	frags := poolFragments.Get().(*[][]byte)
	qualFrags := poolFragments.Get().(*[][]byte)
	*frags = (*frags)[:0]
	*qualFrags = (*qualFrags)[:0]

	var n int
	for i, contig := range *seqs {
		qual := (*quals)[i]
		if len(*contig) != len(*qual) {
			recycleFragments(frags)
			recycleFragments(qualFrags)
			return nil, nil, 0, fmt.Errorf("sequence and quality lengths differ for record %d: %d != %d", i+1, len(*contig), len(*qual))
		}

		end := len(*contig)
		for s := 0; s < end; s += fragLen {
			e := min(s+fragLen, end)
			if e-s < minFragLen {
				continue
			}
			*frags = append(*frags, (*contig)[s:e])
			*qualFrags = append(*qualFrags, (*qual)[s:e])
			n += e - s
		}
	}

	return frags, qualFrags, n, nil
}

var phredErrorProbabilities = func() [94]float64 {
	var probabilities [94]float64
	for q := range probabilities {
		probabilities[q] = math.Pow(10, -float64(q)/10)
	}
	return probabilities
}()

// adjustPIdentByQual estimates the identity before sequencing errors under an
// independent, symmetric substitution-error model. Given the observed identity
// p and the mean per-base error probability e, the estimate is
//
//	p_adjusted = (p - e/3) / (1 - 4e/3).
//
// Both p and p_adjusted are fractions in this formula, although pident and the
// return value are percentages. Qualities are Phred+33, and e is calculated as
// mean_i(10^(-Q_i/10)) over the aligned query region. This is generally not
// equal to converting the arithmetic mean of the Phred scores because the
// conversion is nonlinear. The estimate is clamped to [0, 100]. It does not
// explicitly model insertion or deletion errors.
func adjustPIdentByQual(pident float64, qual []byte) float64 {
	if len(qual) == 0 {
		return pident
	}

	var errorSum float64
	for _, v := range qual {
		q := int(v) - 33
		if q < 0 || q >= len(phredErrorProbabilities) {
			return pident
		}
		errorSum += phredErrorProbabilities[q]
	}
	meanError := errorSum / float64(len(qual))
	denominator := 1 - 4*meanError/3
	if denominator <= 0 {
		return pident
	}

	adjusted := (pident/100 - meanError/3) / denominator * 100
	return min(100, max(0, adjusted))
}

// --------------------------------------------------------------

// Quantiles for the standard normal — pick by desired true-positive retention.
const (
	ZQuantile95  = 1.645 // 95% sensitivity
	ZQuantile975 = 1.96  // 97.5%
	ZQuantile99  = 2.33  // 99%
)

// MinSharedKmersThreshold returns the recommended MinSharedKmers cutoff for
// a FragmentComparator under the Mash / sourmash model (iid mutations,
// sketched shared-count ~ Poisson(μ)).
//
//	μ = (L - k + 1) * ani^k / scaled
//	T = floor(μ - z * sqrt(μ))
//
// L is the fragment length in bp; ani is the minimum identity to tolerate
// (e.g. 0.80); z is the standard-normal quantile (use ZQuantile95 etc.).
// Result is clamped to [1, math.MaxUint16].
func MinSharedKmersThreshold(L int, k uint8, scaled uint32, ani, z float64) uint16 {
	if scaled == 0 {
		scaled = 1
	}
	nk := L - int(k) + 1
	if nk <= 0 {
		return 1
	}
	mu := float64(nk) * math.Pow(ani, float64(k)) / float64(scaled)
	t := math.Floor(mu - z*math.Sqrt(mu))
	if t < 1 {
		return 1
	}
	if t > float64(math.MaxUint16) {
		return math.MaxUint16
	}
	return uint16(t)
}

// MinSharedKmersThresholdExact returns the largest MinSharedKmers cutoff that
// still retains at least `retention` of true ANI-matching fragment pairs,
// using the exact Binomial distribution. Prefer this over the normal-
// approximation variant when μ is small (short fragments or large scaled).
//
//	X ~ Binomial(n, q),   n = L - k + 1,   q = ani^k / scaled
//	T = Quantile(1 - retention)
//
// retention=0.95 → keep 95% of true positives at the given ANI.
// Result is clamped to [1, math.MaxUint16].
func MinSharedKmersThresholdExact(L int, k uint8, scaled uint32, ani, retention float64) uint16 {
	if scaled == 0 {
		scaled = 1
	}
	nk := L - int(k) + 1
	if nk <= 0 {
		return 1
	}
	q := math.Pow(ani, float64(k)) / float64(scaled)
	if q <= 0 {
		return 1
	}
	if q > 1 {
		q = 1
	}
	b := distuv.Binomial{N: float64(nk), P: q}

	// Binomial has no Quantile in gonum; binary-search the CDF for the
	// smallest T with CDF(T-1) <= 1-retention, i.e. P(X >= T) >= retention.
	alpha := 1 - retention
	lo, hi := 0, nk
	for lo < hi {
		mid := (lo + hi) / 2
		if b.CDF(float64(mid)) > alpha {
			hi = mid
		} else {
			lo = mid + 1
		}
	}
	if lo < 1 {
		return 1
	}
	if lo > math.MaxUint16 {
		return math.MaxUint16
	}
	return uint16(lo)
}
