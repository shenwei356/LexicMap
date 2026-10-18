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
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"math"
	"testing"
)

func TestWriteLargeSearchResultRowPreservesPayload(t *testing.T) {
	for _, spilled := range []bool{false, true} {
		t.Run(fmt.Sprintf("spilled=%t", spilled), func(t *testing.T) {
			r := newPayloadTestResult(string(bytes.Repeat([]byte("ACGT"), 20000)), string(bytes.Repeat([]byte("AC-T"), 20000)))
			a := &(*r.SimilarityDetails)[0].Alignments[0]
			a.Alignment = bytes.Repeat([]byte("|| |"), 20000)
			row := []byte("query\t70000\t1\tgenome")
			want := bytes.Join([][]byte{row, a.CIGAR, a.QSeq, a.TSeq, a.Alignment}, []byte{'\t'})
			want = append(want, '\n')
			cigar, query, target, alignment := a.CIGAR, a.QSeq, a.TSeq, a.Alignment
			if spilled {
				budget := newAlignmentPayloadBudget(1)
				budget.tempDir = t.TempDir()
				store := newAlignmentPayloadStore(budget)
				defer store.close()
				previous := []*SearchResult{}
				if err := store.retain(r, &previous); err != nil {
					t.Fatal(err)
				}
				if err := store.finalize(); err != nil {
					t.Fatal(err)
				}
				var err error
				cigar, query, target, alignment, _, err = store.payload(a, nil)
				if err != nil {
					t.Fatal(err)
				}
			}
			var got bytes.Buffer
			w := bufio.NewWriterSize(&got, 4096)
			if err := writeLargeSearchResultRow(w, row, cigar, query, target, alignment); err != nil {
				t.Fatal(err)
			}
			if err := w.Flush(); err != nil {
				t.Fatal(err)
			}
			if !bytes.Equal(got.Bytes(), want) {
				t.Fatal("long TSV row differs from concatenated output")
			}
		})
	}
}

// rowFailureWriter fails at one segment, including short writes without an error.
type rowFailureWriter struct {
	at, calls int
	err       error
}

func (w *rowFailureWriter) Write(p []byte) (int, error) {
	w.calls++
	if w.calls == w.at {
		return len(p) - 1, w.err
	}
	return len(p), nil
}

func TestWriteLargeSearchResultRowReturnsWriteErrors(t *testing.T) {
	sentinel := errors.New("output failed")
	for at := 1; at <= 10; at++ {
		for _, err := range []error{nil, sentinel} {
			w := &rowFailureWriter{at: at, err: err}
			got := writeLargeSearchResultRow(w, []byte("header"), []byte("2M"), []byte("AC"), []byte("AG"), []byte("|."))
			want := err
			if want == nil {
				want = io.ErrShortWrite
			}
			if !errors.Is(got, want) || w.calls != at {
				t.Fatalf("segment %d: got %v after %d writes, want %v", at, got, w.calls, want)
			}
		}
	}
}

func BenchmarkLargeSearchResultRow(b *testing.B) {
	for _, size := range []int{32 << 10, 1 << 20} {
		payload := bytes.Repeat([]byte("A"), size)
		for _, direct := range []bool{false, true} {
			b.Run(fmt.Sprintf("size=%d/direct=%t", size, direct), func(b *testing.B) {
				w := bufio.NewWriterSize(io.Discard, 64<<10)
				header := []byte("query\tgenome\t100.000\t1\t1")
				cigar := []byte("100M")
				var row []byte
				b.SetBytes(int64(size * 3))
				b.ReportAllocs()
				b.ResetTimer()
				for range b.N {
					row = append(row[:0], header...)
					if direct {
						if err := writeLargeSearchResultRow(w, row, cigar, payload, payload, payload); err != nil {
							b.Fatal(err)
						}
					} else {
						row = append(row, '\t', '1', '0', '0', 'M')
						for range 3 {
							row = append(row, '\t')
							row = append(row, payload...)
						}
						row = append(row, '\n')
						if _, err := w.Write(row); err != nil {
							b.Fatal(err)
						}
						row = nil // former policy for rows above 64 KiB
					}
				}
				if err := w.Flush(); err != nil {
					b.Fatal(err)
				}
			})
		}
	}
}

func TestAppendSearchResultRowMatchesPrintf(t *testing.T) {
	for _, pident := range []float64{97.123456, math.Copysign(0, -1), math.NaN(), math.Inf(1), math.Inf(-1)} {
		for flags := range 8 {
			opt := searchOutputOptions{showSseqIdx: flags&1 != 0, showAvgQual: flags&2 != 0, prefixName: flags&4 != 0}
			r := &SearchResult{AlignedFraction: 99.9995}
			sd := &SimilarityDetail{RC: true, SeqLen: 100000, SeqID: []byte("seq:α"), ChunkIdx: 1, NChunks: 3, SeqIdx: 7, NSeqs: 17}
			c := &AlignmentResult{AlignedFraction: 12.345678, AlignedLength: 101, PIdent: pident, Gaps: 2,
				QBegin: 0, QEnd: 99, TBegin: 54321, TEnd: 54422, Evalue: 1.23456789e-150, BitScore: 998}
			seqID, genomeID := string(sd.SeqID), "genome α"
			if opt.showSseqIdx {
				seqID = fmt.Sprintf("c%d/%d:s%d/%d:%s", sd.ChunkIdx+1, sd.NChunks, sd.SeqIdx+1, sd.NSeqs, sd.SeqID)
			}
			length := fmt.Sprintf("%d", c.AlignedLength)
			if opt.showAvgQual {
				length = fmt.Sprintf("%d:%.1f", c.AlignedLength, 23.25)
			}
			if opt.prefixName {
				genomeID = "Escherichia coli:" + genomeID
			}
			want := fmt.Sprintf("%s\t%d\t%d\t%s\t%s\t%.3f\t%d\t%d\t%.3f\t%s\t%.3f\t%d\t%d\t%d\t%d\t%d\t%c\t%d\t%.2e\t%d",
				"query α", 1000, 17, genomeID, seqID, r.AlignedFraction, 2, 3, c.AlignedFraction, length,
				c.PIdent, c.Gaps, c.QBegin+1, c.QEnd+1, c.TBegin+1, c.TEnd+1, '-', sd.SeqLen, c.Evalue, c.BitScore)
			got := appendSearchResultRow(nil, []byte("query α"), 1000, 17, []byte("genome α"), "Escherichia coli", r, sd, c, 2, 3, opt, 23.25)
			if string(got) != want {
				t.Fatalf("flags=%d: got %q, want %q", flags, got, want)
			}
		}
	}
}

func TestAppendSearchResultRowReusesBuffer(t *testing.T) {
	buf := make([]byte, 0, 512)
	r, sd, c := &SearchResult{}, &SimilarityDetail{SeqID: []byte("seq")}, &AlignmentResult{}
	queryID, genomeID := []byte("query"), []byte("genome")
	allocs := testing.AllocsPerRun(100, func() {
		buf = appendSearchResultRow(buf[:0], queryID, 1000, 17, genomeID, "", r, sd, c, 1, 1, searchOutputOptions{}, 0)
	})
	if allocs != 0 {
		t.Fatalf("reusing the output buffer allocated %g objects per row", allocs)
	}
}
