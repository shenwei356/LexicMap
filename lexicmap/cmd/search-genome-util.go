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
	"gonum.org/v1/gonum/stat/distuv"
)

// GQuery represents a genome query
type GQuery struct {
	id          []byte
	bigSeq      []byte
	seqs        []*[]byte
	skipRegions [][2]int

	genomeSize int

	result *[]*GSearchResult // fragment alignment results

	screenDetails *[]*GSearchScreenResultDetail
}

var poolGQuery = &sync.Pool{New: func() interface{} {
	return &GQuery{
		id:          make([]byte, 0, 127),
		bigSeq:      make([]byte, 0, 10<<20), // 10M
		seqs:        make([]*[]byte, 0, 256),
		skipRegions: make([][2]int, 0, 256),
		genomeSize:  0,
	}
}}

var poolSeq = &sync.Pool{
	New: func() interface{} {
		tmp := make([]byte, 0, 10<<10) // 10K
		return &tmp
	},
}

func (q *GQuery) Reset() {
	q.id = q.id[:0]
	q.bigSeq = q.bigSeq[:0]
	q.seqs = q.seqs[:0]
	q.skipRegions = q.skipRegions[:0]
	q.genomeSize = 0

	q.result = nil
}

func RecycleGQuery(q *GQuery) {
	q.id = q.id[:0]
	q.bigSeq = q.bigSeq[:0]
	if q.seqs != nil {
		for _, s := range q.seqs {
			*s = (*s)[:0]
			poolSeq.Put(s)
		}
		q.seqs = q.seqs[:0]
	}
	if q.skipRegions != nil {
		q.skipRegions = q.skipRegions[:0]
	}
	q.genomeSize = 0

	if q.result != nil {
		RecycleGSearchResults(q.result)
		q.result = nil
	}

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
	fastxReader, err := fastx.NewDefaultReader(file)
	if err != nil {
		return nil, err
	}

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

		if i > 0 {
			q.skipRegions = append(q.skipRegions, [2]int{len(q.bigSeq), len(q.bigSeq) + gr.k - 1})

			q.bigSeq = append(q.bigSeq, gr.nnn...)
		}

		if convertNtoA {
			convertSeq(record.Seq.Seq, table)
		}

		s := poolSeq.Get().(*[]byte)
		*s = (*s)[:0]
		*s = append(*s, record.Seq.Seq...)
		q.seqs = append(q.seqs, s)
		q.genomeSize += len(record.Seq.Seq)

		q.bigSeq = append(q.bigSeq, record.Seq.Seq...)
	}

	lenSeq := len(q.bigSeq)
	if lenSeq == 0 {
		RecycleGQuery(q)
		return nil, nil
	}

	gaps := reGaps.FindAllSubmatchIndex(q.bigSeq, -1)
	if gaps != nil {
		for _, gap := range gaps {
			q.skipRegions = append(q.skipRegions, [2]int{gap[0], gap[1] - 1})
		}

		slices.SortFunc(q.skipRegions, func(a, b [2]int) int {
			return a[0] - b[0]
		})
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
