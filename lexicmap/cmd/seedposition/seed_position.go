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

package seedposition

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"sync"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
)

var be = binary.BigEndian

// Magic number for checking file format
var Magic = [8]byte{'.', 's', 'e', 'e', 'd', 'l', 'o', 'c'}

// Magic number for the index file
var MagicIdx = [8]byte{'.', 'l', 'o', 'c', 'i', 'd', 'e', 'x'}

// PositionsIndexFileExt is the file extension of seed position index file.
var PositionsIndexFileExt = ".idx"

// MainVersion is use for checking compatibility
var MainVersion uint8 = 0

// MinorVersion is less important
var MinorVersion uint8 = 1

// BufferSize is size of reading and writing buffer
var BufferSize = 65536 // os.Getpagesize()

// ErrInvalidFileFormat means invalid file format.
var ErrInvalidFileFormat = errors.New("seed position data: invalid binary format")

// ErrBrokenFile means the file is not complete.
var ErrBrokenFile = errors.New("seed position data: broken file")

// ErrVersionMismatch means version mismatch between files and program
var ErrVersionMismatch = errors.New("seed position data: version mismatch")

// Writer saves a list of seed positions to a file.
type Writer struct {
	batch uint32
	file  string
	fh    *os.File
	w     *bufio.Writer

	bBuf   bytes.Buffer
	buf    []byte // 24 bytes buffer
	offset int

	// offsets
	index [][2]int
}

// NewWriter creates a new Writer.
// Batch is the batch id for this data file.
func NewWriter(file string, batch uint32) (*Writer, error) {
	w := &Writer{
		batch: batch,
		file:  file,
		index: make([][2]int, 0, 1024),
	}
	var err error
	w.fh, err = os.Create(file)
	if err != nil {
		return nil, err
	}
	w.w = bufio.NewWriterSize(w.fh, BufferSize)

	w.buf = make([]byte, 24)

	// 8-byte magic number
	err = binary.Write(w.w, be, Magic)
	if err != nil {
		return nil, err
	}
	w.offset += 8

	// 8-byte meta info
	// actually, only 2 bytes used and the left 6 bytes is preserved.
	err = binary.Write(w.w, be, [8]uint8{MainVersion, MinorVersion})
	if err != nil {
		return nil, err
	}
	w.offset += 8
	return w, nil
}

// Write writes a list of SORTED uint32s.
// The data should be sorted, because writing is seriallized,
// while sorting can be asynchronous.
func (w *Writer) Write(locs []uint32) error {
	n := len(locs)
	// collect data for the index file
	w.index = append(w.index, [2]int{w.offset, n})

	// write genome information
	buf := w.buf
	buf0 := w.bBuf
	buf0.Reset()

	// The number of positions/records
	be.PutUint32(buf[:4], uint32(n))
	buf0.Write(buf[:4])

	// ------------------------------------------------
	var pre, v1, v2, v3, v4 uint32
	var ctrlByte byte
	var nBytes int
	var i int

	round := (len(locs) + 3) >> 2
	last := round - 1
	if n > 0 && n&3 == 0 { // 4*x (x>0)
		last = round
	}
	for r := 0; r < last; r++ {
		i = r << 2
		v1 = locs[i] - pre
		v2 = locs[i+1] - locs[i]
		v3 = locs[i+2] - locs[i+1]
		v4 = locs[i+3] - locs[i+2]

		ctrlByte, nBytes = util.PutUint32s(buf[1:], v1, v2, v3, v4)
		buf[0] = ctrlByte
		buf0.Write(buf[:nBytes+1])

		pre = locs[i+3]
	}
	if n&3 > 0 {
		i = (round - 1) << 2
		switch n & 3 {
		case 1:
			v1 = locs[i] - pre
			v2, v3, v4 = 0, 0, 0
		case 2:
			v1 = locs[i] - pre
			v2 = locs[i+1] - locs[i]
			v3, v4 = 0, 0
		case 3:
			v1 = locs[i] - pre
			v2 = locs[i+1] - locs[i]
			v3 = locs[i+2] - locs[i+1]
			v4 = 0
		}
		ctrlByte, nBytes = util.PutUint32s(buf[1:], v1, v2, v3, v4)
		buf[0] = ctrlByte
		buf0.Write(buf[:nBytes+1])
	}

	// ------------------------------------------------
	// write data to file
	_, err := w.w.Write(buf0.Bytes())
	if err != nil {
		return err
	}
	w.offset += buf0.Len()

	return err
}

// Close writes the index file and finishes the writing.
func (w *Writer) Close() error {
	err := w.w.Flush()
	if err != nil {
		return err
	}

	err = w.fh.Close()
	if err != nil {
		return err
	}

	// write the index

	fh, err := os.Create(filepath.Clean(w.file) + PositionsIndexFileExt)
	if err != nil {
		return err
	}
	wtr := bufio.NewWriterSize(fh, BufferSize)

	// magic
	err = binary.Write(wtr, be, MagicIdx)
	if err != nil {
		return err
	}

	// versions
	// actually, only 2 bytes used and the left 6 bytes is preserved.
	err = binary.Write(wtr, be, [8]uint8{MainVersion, MinorVersion})
	if err != nil {
		return err
	}

	buf := w.buf
	buf0 := w.bBuf
	buf0.Reset()

	// batch number
	be.PutUint32(buf[:4], w.batch)
	// the number of records
	be.PutUint32(buf[4:8], uint32(len(w.index)))
	buf0.Write(buf[:8])

	buf = w.buf[:12]
	for _, data := range w.index {
		be.PutUint64(buf[:8], uint64(data[0]))   // offset
		be.PutUint32(buf[8:12], uint32(data[1])) // number of positions
		buf0.Write(buf)
	}

	_, err = wtr.Write(buf0.Bytes())
	if err != nil {
		return err
	}

	err = wtr.Flush()
	if err != nil {
		return err
	}

	return fh.Close()
}

// Reader is for reading the seed position of a genome
type Reader struct {
	batch    uint32
	nRecords uint32

	fh     *os.File
	offset int // offset of the first index record

	buf []byte

	fhData *os.File
}

var poolReader = &sync.Pool{New: func() interface{} {
	return &Reader{
		buf: make([]byte, 32),
	}
}}

// NewReader returns a reader from a seed position file.
// The reader is recycled after calling Close().
func NewReader(file string) (*Reader, error) {
	if strings.HasSuffix(file, PositionsIndexFileExt) {
		return nil, fmt.Errorf("seed position file, not the index file should be given")
	}

	// ------------  index file ----------------

	fileIndex := filepath.Clean(file) + PositionsIndexFileExt
	var err error
	r := poolReader.Get().(*Reader)

	r.fh, err = os.Open(fileIndex)
	if err != nil {
		return nil, err
	}

	buf := r.buf

	// check the magic number
	n, err := io.ReadFull(r.fh, buf[:8])
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
	r.offset += 8

	// read metadata
	n, err = io.ReadFull(r.fh, buf[:8])
	if err != nil {
		return nil, err
	}
	if n < 8 {
		return nil, ErrBrokenFile
	}
	r.offset += 8

	// check compatibility
	if MainVersion != buf[0] {
		return nil, ErrVersionMismatch
	}

	// batch number and the number seqs
	n, err = io.ReadFull(r.fh, buf[:8])
	if err != nil {
		return nil, err
	}
	if n < 8 {
		return nil, ErrBrokenFile
	}
	r.offset += 8

	r.batch = be.Uint32(buf[:4])
	r.nRecords = be.Uint32(buf[4:8])

	// ------------ data file ----------------

	r.fhData, err = os.Open(file)
	if err != nil {
		return nil, err
	}

	return r, nil
}

// Close closes and recycles the reader.
func (r *Reader) Close() error {
	err := r.fh.Close()
	if err != nil {
		poolReader.Put(r)
		return err
	}

	err = r.fhData.Close()
	if err != nil {
		poolReader.Put(r)
		return err
	}

	poolReader.Put(r)
	return nil
}

// SeedPositions returns the seed positions with an index of idx (0-based).
func (r *Reader) SeedPositions(idx int, locs *[]uint32) error {
	if idx < 0 || idx >= int(r.nRecords) {
		return fmt.Errorf("genome index (%d) out of range: [0, %d]", idx, int(r.nRecords)-1)
	}

	buf := r.buf

	// -----------------------------------------------------------
	// read index information
	// 24 + 12 * idx
	r.fh.Seek(int64(r.offset)+int64(idx)<<3+int64(idx)<<2, 0)

	// offset in the data file and bases
	n, err := io.ReadFull(r.fh, buf[:12])
	if err != nil {
		return err
	}
	if n < 12 {
		return ErrBrokenFile
	}
	offset := int64(be.Uint64(buf[:8]))
	nRecords := int(be.Uint32(buf[8:12])) // for check end

	if locs == nil {
		tmp := make([]uint32, 0, nRecords)
		locs = &tmp
	} else {
		*locs = (*locs)[:0]
	}

	// ------------------------------------------------
	var pre, v1, v2, v3, v4 uint32
	var ctrlByte byte
	var nBytes, nReaded, nDecoded int

	r.fhData.Seek(offset+4, 0) // 4 is the 4 bytes for storing length of records

	rounds := (nRecords + 3) >> 2

	for round := 0; round < rounds; round++ {
		// read the control byte
		_, err = io.ReadFull(r.fhData, buf[:1])
		if err != nil {
			return err
		}
		ctrlByte = buf[0]

		// parse the control byte
		nBytes = util.CtrlByte2ByteLengthsUint32(ctrlByte)

		// read encoded bytes
		nReaded, err = io.ReadFull(r.fhData, buf[:nBytes])
		if err != nil {
			return err
		}
		if nReaded < nBytes {
			return ErrBrokenFile
		}

		v1, v2, v3, v4, nDecoded = util.Uint32s(ctrlByte, buf[:nBytes])
		if nDecoded == 0 {
			return ErrBrokenFile
		}

		v1 += pre
		*locs = append(*locs, v1)
		v2 += v1
		*locs = append(*locs, v2)
		v3 += v2
		*locs = append(*locs, v3)
		v4 += v3
		*locs = append(*locs, v4)

		pre = v4
	}
	*locs = (*locs)[:nRecords]

	return nil
}
