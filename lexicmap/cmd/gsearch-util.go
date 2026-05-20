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
	"path/filepath"
	"regexp"
	"slices"
	"sync"

	"github.com/shenwei356/bio/seqio/fastx"
)

// GQuery represents a genome query
type GQuery struct {
	id          []byte
	bigSeq      []byte
	seqs        []*[]byte
	skipRegions [][2]int

	genomeSize int

	result interface{} // TODO
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
		tmp := make([]byte, 0, 5<<20) // 5M
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
	q.result = nil
	poolGQuery.Put(q)
}

// --------------------------------------------------------------

type GenomeReader struct {
	k         int            // kmer size
	nnn       []byte         // Ns
	reGaps    *regexp.Regexp // for finding gap regions
	reRefName *regexp.Regexp // for extracting genome id from the file name
}

// NewGenomeReader returns a GenomeReader with given k-mer size
func NewGenomeReader(k int, reRefName *regexp.Regexp) *GenomeReader {
	return &GenomeReader{
		k:         k,
		nnn:       bytes.Repeat([]byte{'N'}, k),
		reGaps:    regexp.MustCompile(fmt.Sprintf(`[Nn]{%d,}`, 5)),
		reRefName: reRefName,
	}
}

// Recycle recyle a GQuery object
func (gr *GenomeReader) Recycle(q *GQuery) {
	RecycleGQuery(q)
}

// Read reads a genome from a file or stdin
func (gr *GenomeReader) Read(file string) (*GQuery, error) {
	fastxReader, err := fastx.NewDefaultReader(file)
	if err != nil {
		return nil, err
	}

	q := poolGQuery.Get().(*GQuery)
	q.Reset()

	var record *fastx.Record
	i := 0
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

		s := poolSeq.Get().(*[]byte)
		*s = (*s)[:0]
		*s = append(*s, record.Seq.Seq...)
		q.seqs = append(q.seqs, s)

		q.bigSeq = append(q.bigSeq, record.Seq.Seq...)
	}

	lenSeq := len(q.bigSeq)
	if lenSeq == 0 {
		RecycleGQuery(q)
		return nil, nil
	}

	gaps := gr.reGaps.FindAllSubmatchIndex(q.bigSeq, -1)
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
