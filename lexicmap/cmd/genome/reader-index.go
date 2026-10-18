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
	"bufio"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
)

// ReaderIndex owns immutable record offsets and lengths for one genome file.
// Its readers have independent file handles and buffers. Keep the genome file
// and its index unchanged until all readers using this directory are closed.
type ReaderIndex struct {
	file  string   // absolute path binds every reader to the same genome data file
	batch uint32   // genome batch number
	nSeqs uint32   // number of genome records, not contigs
	index []uint64 // interleaved offset and base count; never mutated after loading
}

// ReadIndex loads a genome file's record directory once for concurrent readers.
// Pass the genome data filename, not its .idx companion.
func ReadIndex(file string) (index *ReaderIndex, err error) {
	if file == "" {
		return nil, fmt.Errorf("genome data file path is empty")
	}
	if strings.HasSuffix(file, GenomeIndexFileExt) {
		return nil, fmt.Errorf("genome file, not the index file should be given")
	}

	// ------------ genome index file ----------------

	file, err = filepath.Abs(file)
	if err != nil {
		return nil, err
	}
	fileIndex := file + GenomeIndexFileExt
	r := &ReaderIndex{file: file}

	fh, err := os.Open(fileIndex)
	if err != nil {
		return nil, err
	}
	defer func() {
		if closeErr := fh.Close(); err == nil && closeErr != nil {
			index, err = nil, closeErr
		}
	}()
	bfh := bufio.NewReaderSize(fh, 64<<10) // 64KB buffer

	var scratch [24]byte
	buf := scratch[:]

	// check the magic number
	n, err := io.ReadFull(bfh, buf[:8])
	if err != nil {
		return nil, err
	}
	if n < 8 {
		return nil, ErrBrokenFile
	}
	same := true
	for i := 0; i < 8; i++ {
		if MagicIdx[i] != buf[i] {
			same = false
			break
		}
	}
	if !same {
		return nil, ErrInvalidFileFormat
	}

	// read metadata
	n, _ = io.ReadFull(bfh, buf[:8])
	if n < 8 {
		return nil, ErrBrokenFile
	}

	// check compatibility
	if MainVersion != buf[0] {
		return nil, ErrVersionMismatch
	}

	// batch number and the number seqs
	n, _ = io.ReadFull(bfh, buf[:8])
	if n < 8 {
		return nil, ErrBrokenFile
	}

	r.batch = be.Uint32(buf[:4])
	r.nSeqs = be.Uint32(buf[4:8])

	// read all index data, because it's small
	info, err := fh.Stat()
	if err != nil {
		return nil, err
	}
	// Reject truncated counts before allocating the shared directory.
	if int64(r.nSeqs) > (info.Size()-24)/12 || uint64(r.nSeqs) > uint64(int(^uint(0)>>1)/2) {
		return nil, ErrBrokenFile
	}
	r.index = make([]uint64, int(r.nSeqs)*2)
	var i2 int
	for i := 0; i < int(r.nSeqs); i++ {
		// offset in the data file and bases
		n, _ = io.ReadFull(bfh, buf[:12])
		if n < 12 {
			return nil, ErrBrokenFile
		}
		i2 = i << 1
		r.index[i2] = be.Uint64(buf[:8])
		r.index[i2+1] = uint64(be.Uint32(buf[8:12]))
	}
	return r, nil
}

// NewReader opens an independent data handle using the shared immutable directory.
// Closing this reader neither closes other readers nor invalidates the directory.
func (index *ReaderIndex) NewReader() (*Reader, error) {
	fh, err := os.Open(index.file)
	if err != nil {
		return nil, err
	}
	r := poolReader.Get().(*Reader)
	r.batch, r.nSeqs, r.Index = index.batch, index.nSeqs, index.index
	r.buf = resizeByteSlice(r.buf, 24)
	r.fhData = fh
	if r.bufReader == nil {
		r.bufReader = bufio.NewReaderSize(nil, 4096)
	} else {
		r.bufReader.Reset(nil)
	}
	return r, nil
}
