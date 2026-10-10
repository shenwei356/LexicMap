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

package kv

import (
	"encoding/binary"
	"fmt"
	"io"
	"slices"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
)

type maskCount struct {
	offset int64
	count  uint64
}

// mergeCursor decodes pair headers, leaving postings in the reader until the
// current key is selected. Its second key follows the first key's postings.
type mergeCursor struct {
	reader                    *Reader
	key, count                uint64
	secondKey, secondCount    uint64
	offset, expected, decoded uint64
	hasSecond, last           bool
	width                     int
}

func (c *mergeCursor) begin() (bool, error) {
	if _, err := io.ReadFull(c.reader.r, c.reader.buf8); err != nil {
		return false, err
	}
	c.expected = be.Uint64(c.reader.buf8)
	c.offset, c.decoded = 0, 0
	if c.expected == 0 {
		return false, nil
	}
	return true, c.readPair()
}

func (c *mergeCursor) readPair() error {
	r, buf := c.reader.r, c.reader.buf
	if _, err := io.ReadFull(r, buf[:1]); err != nil {
		return err
	}
	ctrl := buf[0]
	c.last, c.hasSecond = ctrl&128 != 0, ctrl&64 == 0
	ctrl &= 63
	n := util.CtrlByte2ByteLengthsUint64(ctrl)
	if _, err := io.ReadFull(r, buf[:n]); err != nil {
		return err
	}
	delta1, delta2, _ := util.Uint64s(ctrl, buf[:n])
	c.key = c.offset + delta1
	c.secondKey = c.key + delta2
	if c.key < c.offset || (c.decoded > 0 && c.key == c.offset) ||
		(c.hasSecond && c.secondKey <= c.key) || (!c.hasSecond && !c.last) {
		return ErrBrokenFile
	}
	c.offset = c.secondKey
	c.decoded++
	if c.hasSecond {
		c.decoded++
	}
	if c.decoded > c.expected || c.last != (c.decoded == c.expected) {
		return fmt.Errorf("number of k-mers mismatch or invalid last-pair flag: expected %d, decoded %d", c.expected, c.decoded)
	}
	if _, err := io.ReadFull(r, buf[:1]); err != nil {
		return err
	}
	ctrl = buf[0]
	n = util.CtrlByte2ByteLengthsUint64(ctrl)
	if _, err := io.ReadFull(r, buf[:n]); err != nil {
		return err
	}
	c.count, c.secondCount, _ = util.Uint64s(ctrl, buf[:n])
	return nil
}

func (c *mergeCursor) advance() (bool, error) {
	if c.hasSecond {
		c.key, c.count = c.secondKey, c.secondCount
		c.hasSecond = false
		return true, nil
	}
	if c.last {
		return false, nil
	}
	return true, c.readPair()
}

func (c *mergeCursor) appendPositions(dst []byte, width int) ([]byte, error) {
	// Check before converting the on-disk uint64 count to an allocation size.
	if c.count > uint64((int(^uint(0)>>1)-len(dst))/width) {
		return dst, ErrBrokenFile
	}
	start := len(dst)
	n := int(c.count) * width
	if c.width == width {
		// Grow only as bytes are consumed, so a corrupt count cannot force a
		// huge allocation before its truncated payload is detected.
		for n > 0 {
			batch := min(n, 32<<10)
			start = len(dst)
			dst = slices.Grow(dst, batch)[:start+batch]
			if _, err := io.ReadFull(c.reader.r, dst[start:]); err != nil {
				return dst, err
			}
			n -= batch
		}
		return dst, nil
	}
	// Preserve the old decode/re-encode behavior for mixed-width inputs.
	buf := c.reader.buf2048
	remaining := c.count
	for remaining > 0 {
		count := min(remaining, uint64(len(buf)/c.width))
		if _, err := io.ReadFull(c.reader.r, buf[:int(count)*c.width]); err != nil {
			return dst, err
		}
		start = len(dst)
		dst = slices.Grow(dst, int(count)*width)[:start+int(count)*width]
		for i := 0; i < int(count); i++ {
			var value uint64
			if c.width == 7 {
				value = Uint64ThreeBytes(buf[i*7:])
			} else {
				value = be.Uint64(buf[i*8:])
			}
			if width == 7 {
				PutUint64ThreeBytes(dst[start:], value)
			} else {
				be.PutUint64(dst[start:], value)
			}
			start += width
		}
		remaining -= count
	}
	return dst, nil
}

// Merger merges sorted masks without materializing a whole mask's postings.
// Reuse it across a chunk's masks. Readers and writer remain caller-owned.
type Merger struct {
	writer  *Writer
	cursors []mergeCursor
	heap    []int // cursor indexes, ordered by key then input file order
	values  [2][]byte
	width   int
}

// NewMerger prepares a reusable merger for compatible KV chunks.
func NewMerger(readers []*Reader, writer *Writer) (*Merger, error) {
	m := &Merger{writer: writer, cursors: make([]mergeCursor, len(readers)),
		heap: make([]int, 0, len(readers)), width: 8}
	if writer.use3BytesForSeedPos {
		m.width = 7
	}
	for i, reader := range readers {
		if reader.K != writer.K || reader.ChunkIndex != writer.ChunkIndex || reader.ChunkSize != writer.ChunkSize {
			return nil, fmt.Errorf("incompatible KV chunk: %s", reader.file)
		}
		m.cursors[i] = mergeCursor{reader: reader, width: 8}
		if reader.Use3BytesForSeedPos {
			m.cursors[i].width = 7
		}
	}
	return m, nil
}

func (m *Merger) less(a, b int) bool {
	ka, kb := m.cursors[a].key, m.cursors[b].key
	return ka < kb || (ka == kb && a < b)
}

func (m *Merger) down(i int) {
	h := m.heap
	for {
		child := i*2 + 1
		if child >= len(h) {
			return
		}
		if child+1 < len(h) && m.less(h[child+1], h[child]) {
			child++
		}
		if !m.less(h[child], h[i]) {
			return
		}
		h[i], h[child] = h[child], h[i]
		i = child
	}
}

// next combines equal keys in input order, including duplicate positions.
func (m *Merger) next(slot int) (uint64, error) {
	key := m.cursors[m.heap[0]].key
	dst := m.values[slot][:0]
	for len(m.heap) > 0 && m.cursors[m.heap[0]].key == key {
		c := &m.cursors[m.heap[0]]
		var err error
		dst, err = c.appendPositions(dst, m.width)
		if err != nil {
			return 0, fmt.Errorf("read postings from %s: %w", c.reader.file, err)
		}
		more, err := c.advance()
		if err != nil {
			return 0, fmt.Errorf("read key from %s: %w", c.reader.file, err)
		}
		if !more {
			m.heap[0] = m.heap[len(m.heap)-1]
			m.heap = m.heap[:len(m.heap)-1]
		}
		m.down(0)
	}
	m.values[slot] = dst
	return key, nil
}

// WriteMask consumes one mask from each reader and writes the merged mask.
func (m *Merger) WriteMask() error {
	m.heap = m.heap[:0]
	for i := range m.cursors {
		more, err := m.cursors[i].begin()
		if err != nil {
			return fmt.Errorf("read mask from %s: %w", m.cursors[i].reader.file, err)
		}
		if more {
			m.heap = append(m.heap, i)
		}
	}
	for i := len(m.heap)/2 - 1; i >= 0; i-- {
		m.down(i)
	}
	wtr := m.writer
	headerOffset := int64(wtr.N)
	if err := binary.Write(wtr.w, be, uint64(0)); err != nil {
		return err
	}
	wtr.N += 8
	if len(m.heap) == 0 {
		return binary.Write(wtr.wi, be, uint64(0))
	}
	p2o := wtr.beginMaskIndex()
	// Return the directory on errors too; writeMaskIndex owns it on success.
	indexWritten := false
	defer func() {
		if !indexWritten && p2o != nil {
			wtr.poolP2O.Put(p2o)
		}
	}()
	var offset, count uint64
	for len(m.heap) > 0 {
		key1, err := m.next(0)
		if err != nil {
			return err
		}
		count++
		hasSecond := len(m.heap) > 0
		var key2 uint64
		m.values[1] = m.values[1][:0]
		if hasSecond {
			key2, err = m.next(1)
			if err != nil {
				return err
			}
			count++
		}
		delta2 := uint64(0)
		if hasSecond {
			delta2 = key2 - key1
		}
		ctrl, n := util.PutUint64s(wtr.buf[1:], key1-offset, delta2)
		if len(m.heap) == 0 {
			ctrl |= 128
		}
		if !hasSecond {
			ctrl |= 64
		}
		wtr.buf[0] = ctrl
		used := n + 1
		ctrl, n = util.PutUint64s(wtr.buf[used+1:], uint64(len(m.values[0])/m.width), uint64(len(m.values[1])/m.width))
		wtr.buf[used] = ctrl
		used += n + 1
		pairOffset := uint64(wtr.N) // offsets refer to the encoded pair header
		if _, err = wtr.w.Write(wtr.buf[:used]); err != nil {
			return err
		}
		wtr.N += used
		for _, values := range m.values {
			n, err = wtr.w.Write(values)
			wtr.N += n
			if err != nil {
				return err
			}
		}
		if err = wtr.indexSeedPair(p2o, key1, key2, hasSecond, pairOffset); err != nil {
			return err
		}
		offset = key2
	}
	wtr.maskCounts = append(wtr.maskCounts, maskCount{offset: headerOffset, count: count})
	indexWritten = true
	return wtr.writeMaskIndex(p2o)
}
