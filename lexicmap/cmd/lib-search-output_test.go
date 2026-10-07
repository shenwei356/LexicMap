package cmd

import (
	"fmt"
	"math"
	"testing"
)

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
