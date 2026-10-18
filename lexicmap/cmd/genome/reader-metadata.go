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

package genome

import (
	"fmt"
	"io"
)

// SubSeqWithoutMetadata extracts sequence before decoding contig lengths/IDs.
// Reuse g only for the same record. LoadMetadata must succeed before using its
// contig boundaries or names. On an extraction error g is recycled, as in SubSeq3.
func (r *Reader) SubSeqWithoutMetadata(idx, start, end int, g *Genome) (*Genome, error) {
	if idx < 0 || idx >= int(r.nSeqs) {
		return nil, fmt.Errorf("sequence index (%d) out of range: [0, %d]", idx, int(r.nSeqs)-1)
	}
	if g == nil {
		offset, ok, err := r.sequenceOffset(idx)
		if err != nil {
			return nil, err
		}
		if !ok {
			// Use the original parser when the indexed bounds cannot contain the DNA.
			return r.SubSeq3(idx, start, end, nil)
		}
		g = PoolGenome.Get().(*Genome)
		g.Reset()
		g.SeqOffSet = offset
		g.Len = int(r.Index[idx*2+1])
	}
	// Use the same bounds and decoder as eager extraction; g already has its offset.
	return r.SubSeq3(idx, start, end, g)
}

// sequenceOffset locates the length header immediately before the packed DNA.
// Records end with ceil(nbases/4) DNA bytes, with no padding or file footer:
//
//	[variable metadata][nbytes, nbases][packed DNA] | next record / EOF
//
// Trust the indexed lengths to avoid a separate header read. LoadMetadata later
// checks the inferred offset; this shortcut is not a file-integrity check.
func (r *Reader) sequenceOffset(idx int) (int64, bool, error) {
	i := idx << 1
	begin, bases := r.Index[i], r.Index[i+1]
	var end uint64
	if idx+1 < int(r.nSeqs) {
		end = r.Index[i+2]
	} else {
		info, err := r.fhData.Stat()
		if err != nil {
			return 0, false, err
		}
		end = uint64(info.Size())
	}
	packed := (bases + 3) / 4 // bases came from a uint32 on disk
	if end > uint64(1<<63-1) || end < begin || end-begin < packed+8 {
		return 0, false, nil
	}
	offset := int64(end - packed - 8)
	return offset, true, nil
}

// LoadMetadata fills the contig metadata for g's record without changing its
// extracted sequence (which may already be reverse-complemented). Repeated calls
// are no-ops. On error g, including partially read IDs, is recycled; use the return.
func (r *Reader) LoadMetadata(idx int, g *Genome) (*Genome, error) {
	if g == nil {
		return nil, fmt.Errorf("cannot load metadata into a nil genome")
	}
	if idx < 0 || idx >= int(r.nSeqs) {
		RecycleGenome(g)
		return nil, fmt.Errorf("sequence index (%d) out of range: [0, %d]", idx, int(r.nSeqs)-1)
	}
	if g.metadataLoaded {
		return g, nil
	}
	offset, length := g.SeqOffSet, g.Len
	if err := r.readMetadata(idx, g); err != nil {
		RecycleGenome(g)
		return nil, err
	}
	if g.SeqOffSet != offset {
		RecycleGenome(g)
		return nil, ErrBrokenFile
	}
	g.Len = length
	return g, nil
}

// readMetadata fills a fresh Genome or one containing only an extracted sequence.
// The caller owns cleanup on error. Append each borrowed ID before reading it so
// RecycleGenome also returns buffers from partially decoded metadata.
func (r *Reader) readMetadata(idx int, g *Genome) error {
	r.buf = resizeByteSlice(r.buf, 24)
	buf := r.buf
	offset := int64(r.Index[idx*2])
	if _, err := r.fhData.Seek(offset, 0); err != nil {
		return err
	}
	br := r.bufReader
	br.Reset(r.fhData)
	if _, err := io.ReadFull(br, buf[:2]); err != nil {
		return ErrBrokenFile
	}
	idLen := int(be.Uint16(buf[:2]))
	if _, err := br.Discard(idLen); err != nil { // genome ID is not used by SubSeq3
		return ErrBrokenFile
	}
	offset += int64(2 + idLen)
	if _, err := io.ReadFull(br, buf[:12]); err != nil {
		return ErrBrokenFile
	}
	g.GenomeSize = int(be.Uint32(buf[:4]))
	g.Len = int(be.Uint32(buf[4:8]))
	g.NumSeqs = int(be.Uint32(buf[8:12]))
	offset += 12
	for i := 0; i < g.NumSeqs; i++ {
		if _, err := io.ReadFull(br, buf[:6]); err != nil {
			return ErrBrokenFile
		}
		g.SeqSizes = append(g.SeqSizes, int(be.Uint32(buf[:4])))
		idLen = int(be.Uint16(buf[4:6]))
		id := poolID.Get().(*[]byte)
		*id = resizeByteSlice(*id, idLen)
		g.SeqIDs = append(g.SeqIDs, id)
		if _, err := io.ReadFull(br, *id); err != nil {
			return ErrBrokenFile
		}
		offset += int64(6 + idLen)
	}
	g.SeqOffSet = offset
	g.metadataLoaded = true
	return nil
}
