// Copyright © 2023-2024 Wei Shen <shenwei356@gmail.com>
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

package twobit

import (
	"bufio"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sync"
)

var be = binary.BigEndian

// Magic number for checking file format
var Magic = [8]byte{'2', 'b', 'i', 't', 's', 'e', 'q', 's'}

// the file extension of the 2bit file
const IndexFileExt = ".idx"

// MainVersion is use for checking compatibility
var MainVersion uint8 = 0

// MinorVersion is less important
var MinorVersion uint8 = 1

// BufferSize is size of reading and writing buffer
var BufferSize = 65536 // os.Getpagesize()

// ErrInvalidFileFormat means invalid file format.
var ErrInvalidFileFormat = errors.New("2bit seqs: invalid binary format")

// ErrEmptySeq means the sequence is empty
var ErrEmptySeq = errors.New("2bit seqs: empty seq")

// ErrInvalidTwoBitData means the length of two bit seq slice does not match the number of bases
var ErrInvalidTwoBitData = errors.New("2bit seqs: invalid two-bit data")

// ErrBrokenFile means the file is not complete.
var ErrBrokenFile = errors.New("2bit seqs: broken file")

// ErrKOverflow means K < 1 or K > 32.

// ErrVersionMismatch means version mismatch between files and program
var ErrVersionMismatch = errors.New("2bit seqs: version mismatch")

// Writer saves a list of DNA sequences into 2bit-encoded format.
// The IDs of sequences are not saved.
type Writer struct {
	file string
	fh   *os.File
	w    *bufio.Writer

	buf    []byte // 24 bytes buffer
	offset int

	// offset, #bytes, #bases,
	index [][3]int
}

// NewWriter creates a new Writer.
func NewWriter(file string) (*Writer, error) {
	w := &Writer{file: file}
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

// WriteSeq writes one sequence
func (w *Writer) WriteSeq(s []byte) error {
	b2 := Seq2TwoBit(s)
	err := w.Write2Bit(*b2, len(s))
	RecycleTwoBit(b2)
	return err
}

// Write writes one converted 2bit sequence.
func (w *Writer) Write2Bit(b2 []byte, bases int) error {
	if len(b2) == 0 {
		return ErrEmptySeq
	}
	// possible bases for b2 of n bytes: [n*4-3, n*4]
	if bases < (len(b2)<<2)-3 || bases > len(b2)<<2 {
		return ErrInvalidTwoBitData
	}

	// the number of bytes and bases
	be.PutUint64(w.buf[:8], uint64(len(b2)))
	be.PutUint64(w.buf[8:16], uint64(bases))
	_, err := w.w.Write(w.buf[:16])
	if err != nil {
		return err
	}

	// write 2bit-packed data
	_, err = w.w.Write(b2)
	if err != nil {
		return err
	}

	// collect data for the index file
	w.index = append(w.index, [3]int{w.offset, len(b2), bases})

	w.offset += 16 + len(b2)
	return nil
}

// Close writes the index file and finish writing.
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

	fh, err := os.Create(filepath.Clean(w.file) + IndexFileExt)
	if err != nil {
		return err
	}
	wtr := bufio.NewWriterSize(fh, BufferSize)
	buf := w.buf[:24]

	// the number of records
	be.PutUint64(buf[:8], uint64(len(w.index)))
	_, err = wtr.Write(buf[:8])
	if err != nil {
		return err
	}

	for _, info := range w.index {
		be.PutUint64(buf[:8], uint64(info[0]))    // offset
		be.PutUint64(buf[8:16], uint64(info[1]))  // bytes
		be.PutUint64(buf[16:24], uint64(info[2])) // bases

		_, err = wtr.Write(buf)
		if err != nil {
			return err
		}
	}
	err = wtr.Flush()
	if err != nil {
		return err
	}

	return fh.Close()
}

// Reader is for fast extracting of subsequence of any sequence
type Reader struct {
	fh     *os.File
	offset int

	buf []byte

	index [][3]int

	// use a buffer for frequencly visited sequences.
	// use a counter for visited sequences, if count > 10, let's buf it.
	// every minute, check the visited count in the buf, delete these less visited.
}

// NewReader returns a reader from a file
func NewReader(file string) (*Reader, error) {
	var err error
	r := &Reader{buf: make([]byte, 24)}

	r.fh, err = os.Open(file)
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
		if Magic[i] != buf[i] {
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

	// ------------ index file ----------------

	// if the index file is not existed
	fileIndex := filepath.Clean(file) + IndexFileExt
	rdr, err := os.Open(fileIndex)
	if err != nil {
		return nil, err
	}

	// the number of records
	n, err = io.ReadFull(rdr, buf[:8])
	if err != nil {
		return nil, err
	}
	if n < 8 {
		return nil, ErrBrokenFile
	}

	r.index = make([][3]int, int(be.Uint64(buf[:8])))
	for i := range r.index {
		n, err = io.ReadFull(rdr, buf[:24])
		if err != nil {
			return nil, err
		}
		if n < 24 {
			return nil, ErrBrokenFile
		}

		r.index[i] = [3]int{
			int(be.Uint64(buf[:8])),
			int(be.Uint64(buf[8:16])),
			int(be.Uint64(buf[16:24])),
		}
	}

	return r, nil
}

// Close the file handler.
func (r *Reader) Close() error {
	err := r.fh.Close()
	if err != nil {
		return err
	}
	return nil
}

// Seq returns the sequence with index of idx (0-based).
func (r *Reader) Seq(idx int) (*[]byte, error) {
	if idx < 0 || idx >= len(r.index) {
		return nil, fmt.Errorf("sequence index (%d) out of range: [0, %d]", idx, len(r.index)-1)
	}
	return r.SubSeq(idx, 0, r.index[idx][2]-1)
}

// SubSeq returns the subsequence of sequence (idx is 0-based),
// from start to end (both are 0-based).
// Please call RecycleSeq() after using the result.
func (r *Reader) SubSeq(idx int, start int, end int) (*[]byte, error) {
	if idx < 0 || idx >= len(r.index) {
		return nil, fmt.Errorf("sequence index (%d) out of range: [0, %d]", idx, len(r.index)-1)
	}
	info := r.index[idx]
	offset := info[0] + 16 // 16 is the bytes of #bytes and #bases
	nBases := info[2]
	if start < 0 {
		start = 0
	}
	if end >= nBases-1 {
		end = nBases - 1
	}
	if end < start {
		end = start
	}

	// start of byte
	offset += start >> 2
	_, err := r.fh.Seek(int64(offset), 0)
	if err != nil {
		return nil, err
	}

	nBytes := end>>2 - start>>2 + 1

	var buf []byte
	if nBytes <= len(r.buf) {
		buf = r.buf[:nBytes]
	} else {
		n := nBytes - len(r.buf)
		for i := 0; i < n; i++ {
			r.buf = append(r.buf, 0)
		}
		buf = r.buf
	}
	n, err := io.ReadFull(r.fh, buf)
	if err != nil {
		return nil, err
	}

	if n < nBytes {
		return nil, ErrBrokenFile
	}

	l := end - start + 1

	// initialize with l+4 blank values, because if there less than 4 bases
	// to extract, code below would panic.
	// s := make([]byte, 4, l+4)
	s := poolSubSeq.Get().(*[]byte)
	*s = (*s)[:4]

	// -- first byte --

	b := buf[0]
	j := start & 3

	switch j {
	case 0:
		(*s)[3] = bit2base[b&3]
		b >>= 2
		(*s)[2] = bit2base[b&3]
		b >>= 2
		(*s)[1] = bit2base[b&3]
		b >>= 2
		(*s)[0] = bit2base[b&3]
	case 1:
		(*s)[2] = bit2base[b&3]
		b >>= 2
		(*s)[1] = bit2base[b&3]
		b >>= 2
		(*s)[0] = bit2base[b&3]
	case 2:
		(*s)[1] = bit2base[b&3]
		b >>= 2
		(*s)[0] = bit2base[b&3]
	case 3:
		(*s)[0] = bit2base[b&3]
	}
	j = 4 - j
	*s = (*s)[:j]
	if j >= l {
		tmp := (*s)[:l]
		return &tmp, nil
	}

	// -- middle byte --
	if nBytes > 2 {
		for _, b = range buf[1 : nBytes-1] {
			*s = append(*s, bit2base[b>>6&3])
			*s = append(*s, bit2base[b>>4&3])
			*s = append(*s, bit2base[b>>2&3])
			*s = append(*s, bit2base[b&3])
		}
	}

	if nBytes > 1 {
		// -- last byte --
		b = buf[nBytes-1]
		j = end & 3
		switch j {
		case 0:
			*s = append(*s, bit2base[b>>6&3])
		case 1:
			*s = append(*s, bit2base[b>>6&3])
			*s = append(*s, bit2base[b>>4&3])
		case 2:
			*s = append(*s, bit2base[b>>6&3])
			*s = append(*s, bit2base[b>>4&3])
			*s = append(*s, bit2base[b>>2&3])
		case 3:
			*s = append(*s, bit2base[b>>6&3])
			*s = append(*s, bit2base[b>>4&3])
			*s = append(*s, bit2base[b>>2&3])
			*s = append(*s, bit2base[b&3])
		}
	}

	tmp := (*s)[:l]
	return &tmp, nil
}

// RecycleSeq recycles the sequence
func RecycleSeq(s *[]byte) {
	poolSubSeq.Put(s)
}

var poolSubSeq = &sync.Pool{New: func() interface{} {
	tmp := make([]byte, 4, 10<<10)
	return &tmp
}}

var base2bit = [256]uint8{
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 1, 1, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0,
	0, 0, 0, 1, 3, 3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
	0, 0, 1, 1, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0,
	0, 0, 0, 1, 3, 3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
}

var bit2base = [4]byte{'A', 'C', 'G', 'T'}

// RecycleSeq recycles the sequence
func RecycleTwoBit(b2 *[]byte) {
	poolTwoBit.Put(b2)
}

var poolTwoBit = &sync.Pool{New: func() interface{} {
	tmp := make([]byte, 0, 1<<20)
	return &tmp
}}

// Seq2TwoBit converts a DNA sequence to 2bit-packed sequence.
func Seq2TwoBit(s []byte) *[]byte {
	if s == nil {
		return nil
	}
	if len(s) == 0 {
		return &[]byte{}
	}

	n := len(s) >> 2
	m := len(s) & 3

	// codes := make([]byte, n+1)
	codes := poolTwoBit.Get().(*[]byte)
	*codes = (*codes)[:0]

	var j int
	for i := 0; i < n; i++ {
		j = i << 2

		*codes = append(*codes, base2bit[s[j]]<<6+base2bit[s[j+1]]<<4+base2bit[s[j+2]]<<2+base2bit[s[j+3]])
	}

	if m == 0 {
		tmp := (*codes)[:n]
		return &tmp
	}

	j = n << 2

	switch m {
	case 3:
		*codes = append(*codes, base2bit[s[j]]<<6+base2bit[s[j+1]]<<4+base2bit[s[j+2]]<<2)
	case 2:
		*codes = append(*codes, base2bit[s[j]]<<6+base2bit[s[j+1]]<<4)
	case 1:
		*codes = append(*codes, base2bit[s[j]]<<6)
	}

	return codes
}

// TwoBit2Seq converts a 2bit-packed sequence to DNA.
func TwoBit2Seq(b2 []byte, bases int) ([]byte, error) {
	// possible bases for b2 of n bytes: [n*4-3, n*4]
	if bases < (len(b2)<<2)-3 || bases > len(b2)<<2 {
		return nil, ErrInvalidTwoBitData
	}

	s := make([]byte, bases)
	n := len(s) >> 2
	m := bases & 3
	var b byte
	var j int
	for i := 0; i < n; i++ {
		b = b2[i]
		j = i << 2

		s[j+3] = bit2base[b&3]
		b >>= 2
		s[j+2] = bit2base[b&3]
		b >>= 2
		s[j+1] = bit2base[b&3]
		b >>= 2
		s[j] = bit2base[b&3]
	}
	if m == 0 {
		return s, nil
	}

	b = b2[n]
	j = n << 2
	switch m {
	case 1:
		s[j] = bit2base[b>>6&3]
	case 2:
		b >>= 4
		s[j+1] = bit2base[b&3]
		b >>= 2
		s[j] = bit2base[b&3]
	case 3:
		b >>= 2
		s[j+2] = bit2base[b&3]
		b >>= 2
		s[j+1] = bit2base[b&3]
		b >>= 2
		s[j] = bit2base[b&3]
	}

	return s, nil
}
