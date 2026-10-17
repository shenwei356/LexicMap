package cmd

import (
	"io"
	"strconv"
)

// maxSearchOutputRowBytes bounds the reusable row buffer. Larger alignment
// payloads are written separately so they are not copied into a second buffer.
const maxSearchOutputRowBytes = 64 << 10

var searchOutputTab = []byte{'\t'}
var searchOutputNewline = []byte{'\n'}

// writeLargeSearchResultRow writes the same TSV fields as the small-row path.
// The caller owns every slice until this synchronous write finishes, including
// payloads borrowed from the spill reader's reusable buffer.
func writeLargeSearchResultRow(w io.Writer, row, cigar, query, target, alignment []byte) error {
	for _, part := range [][]byte{row, searchOutputTab, cigar, searchOutputTab, query,
		searchOutputTab, target, searchOutputTab, alignment, searchOutputNewline} {
		n, err := w.Write(part)
		if err != nil {
			return err
		}
		if n != len(part) {
			return io.ErrShortWrite
		}
	}
	return nil
}

type searchOutputOptions struct {
	showSseqIdx bool
	showAvgQual bool
	prefixName  bool
}

// appendSearchResultRow preserves the TSV's numeric formats without temporary
// strings or interface arguments. dst belongs to the serial output handler.
func appendSearchResultRow(dst []byte, queryID []byte, qlen, targets int, genomeID []byte, genomeName string,
	r *SearchResult, sd *SimilarityDetail, c *AlignmentResult, cluster, hsp int, opt searchOutputOptions, avgQual float64,
) []byte {
	dst = append(dst, queryID...)
	dst = append(dst, '\t')
	dst = strconv.AppendInt(dst, int64(qlen), 10)
	dst = append(dst, '\t')
	dst = strconv.AppendInt(dst, int64(targets), 10)
	dst = append(dst, '\t')
	if opt.prefixName {
		dst = append(dst, genomeName...)
		dst = append(dst, ':')
	}
	dst = append(dst, genomeID...)
	dst = append(dst, '\t')
	if opt.showSseqIdx {
		dst = append(dst, 'c')
		dst = strconv.AppendUint(dst, uint64(sd.ChunkIdx+1), 10)
		dst = append(dst, '/')
		dst = strconv.AppendUint(dst, uint64(sd.NChunks), 10)
		dst = append(dst, ':', 's')
		dst = strconv.AppendUint(dst, uint64(sd.SeqIdx+1), 10)
		dst = append(dst, '/')
		dst = strconv.AppendUint(dst, uint64(sd.NSeqs), 10)
		dst = append(dst, ':')
	}
	dst = append(dst, sd.SeqID...)
	dst = append(dst, '\t')
	dst = strconv.AppendFloat(dst, r.AlignedFraction, 'f', 3, 64)
	dst = append(dst, '\t')
	dst = strconv.AppendInt(dst, int64(cluster), 10)
	dst = append(dst, '\t')
	dst = strconv.AppendInt(dst, int64(hsp), 10)
	dst = append(dst, '\t')
	dst = strconv.AppendFloat(dst, c.AlignedFraction, 'f', 3, 64)
	dst = append(dst, '\t')
	dst = strconv.AppendInt(dst, int64(c.AlignedLength), 10)
	if opt.showAvgQual {
		dst = append(dst, ':')
		dst = strconv.AppendFloat(dst, avgQual, 'f', 1, 64)
	}
	dst = append(dst, '\t')
	dst = strconv.AppendFloat(dst, c.PIdent, 'f', 3, 64)
	dst = append(dst, '\t')
	dst = strconv.AppendInt(dst, int64(c.Gaps), 10)
	dst = append(dst, '\t')
	dst = strconv.AppendInt(dst, int64(c.QBegin+1), 10)
	dst = append(dst, '\t')
	dst = strconv.AppendInt(dst, int64(c.QEnd+1), 10)
	dst = append(dst, '\t')
	dst = strconv.AppendInt(dst, int64(c.TBegin+1), 10)
	dst = append(dst, '\t')
	dst = strconv.AppendInt(dst, int64(c.TEnd+1), 10)
	dst = append(dst, '\t')
	strand := byte('+')
	if sd.RC {
		strand = '-'
	}
	dst = append(dst, strand, '\t')
	dst = strconv.AppendInt(dst, int64(sd.SeqLen), 10)
	dst = append(dst, '\t')
	dst = strconv.AppendFloat(dst, c.Evalue, 'e', 2, 64)
	dst = append(dst, '\t')
	return strconv.AppendInt(dst, int64(c.BitScore), 10)
}
