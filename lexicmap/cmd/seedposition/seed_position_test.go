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
	"errors"
	"fmt"
	"math/rand"
	"os"
	"path/filepath"
	"reflect"
	"testing"
)

// Reusing one writer must amortize encoding allocations and preserve positions.
func TestWriterBufferReuse(t *testing.T) {
	file := filepath.Join(t.TempDir(), "positions.bin")
	w, err := NewWriter(file, 0)
	if err != nil {
		t.Fatal(err)
	}
	closed := false
	t.Cleanup(func() {
		if !closed {
			w.fh.Close()
		}
	})
	positions := make([]uint32, 100000)
	for i := range positions {
		positions[i] = uint32(i * 20)
	}
	if allocations := testing.AllocsPerRun(20, func() {
		if err := w.Write(positions); err != nil {
			t.Fatal(err)
		}
	}); allocations != 0 {
		t.Errorf("repeated position encoding allocated %.0f times per record", allocations)
	}
	records := len(w.index)
	closed = true // Close owns the handle even when finalization fails.
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	if w.bBuf.Cap() != 0 {
		t.Error("closed writer retains encoding scratch")
	}
	r, err := NewReader(file)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	var got []uint32
	if err := r.SeedPositions(records-1, &got); err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(got, positions) {
		t.Fatal("reused buffer changed the saved positions")
	}
}

// writerFailure injects a delayed I/O error without relying on a full filesystem.
type writerFailure struct{ err error }

func (w writerFailure) Write([]byte) (int, error) { return 0, w.err }

func TestWriterCloseAfterFlushError(t *testing.T) {
	file := filepath.Join(t.TempDir(), "positions.bin")
	w, err := NewWriter(file, 0)
	if err != nil {
		t.Fatal(err)
	}
	defer w.fh.Close()
	index := file + PositionsIndexFileExt
	if err := os.WriteFile(index, []byte("original index"), 0600); err != nil {
		t.Fatal(err)
	}
	want := errors.New("injected write failure")
	w.w = bufio.NewWriter(writerFailure{want})
	if _, err := w.w.WriteString("buffered record"); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); !errors.Is(err, want) {
		t.Fatalf("lost flush error: %v", err)
	}
	if _, err := w.fh.Stat(); !errors.Is(err, os.ErrClosed) {
		t.Errorf("data file still open after flush failure: %v", err)
	}
	if got, err := os.ReadFile(index); err != nil || string(got) != "original index" {
		t.Fatalf("published index after data failure: %q, %v", got, err)
	}
}

func TestSeedPositions(t *testing.T) {
	tests := [][]uint32{
		{},
		{1},
		{1, 15},
		{1, 15, 300},
		{1, 15, 300, 301},
		{1, 15, 300, 301, 2500},
		{1, 15, 300, 301, 2500, 3100},
		{1, 15, 300, 301, 2500, 3100, 3111},
		{1, 15, 300, 301, 2500, 3100, 3111, 5000},
		{1, 15, 300, 301, 2500, 3100, 3111, 5000, 10000},
	}

	file := filepath.Join(t.TempDir(), "positions.bin")

	// ---------------------------------------

	wtr, err := NewWriter(file, 0)
	if err != nil {
		t.Error(err)
		return
	}
	for i, test := range tests {
		err = wtr.Write(test)
		if err != nil {
			t.Errorf("write #%d data: %s", i+1, err)
			return
		}
	}
	err = wtr.Close()
	if err != nil {
		t.Error(err)
		return
	}

	idxs := make([]int, len(tests))
	for i := range tests {
		idxs[i] = i
	}
	rand.Shuffle(len(tests), func(i, j int) { idxs[i], idxs[j] = idxs[j], idxs[i] })

	// ---------------------------------------

	rdr, err := NewReader(file)
	if err != nil {
		t.Error(err)
		return
	}
	// Close before TempDir cleanup, including when a read assertion fails.
	t.Cleanup(func() {
		if err := rdr.Close(); err != nil {
			t.Error(err)
		}
	})

	locs := make([]uint32, 64)
	var test []uint32
	var v uint32
	var j int

	for _, i := range idxs {
		test = tests[i]
		err = rdr.SeedPositions(i, &locs)
		if err != nil {
			t.Errorf("read #%d data: %s", i, err)
			return
		}

		if len(locs) != len(test) {
			t.Errorf("[#%d] unequal of position numbers, expected: %d, returned %d",
				i, len(test), len(locs))
			return
		}

		for j, v = range locs {
			if v != test[j] {
				t.Errorf("[#%d] unequal of positions, expected: %d, returned %d", i, test[j], v)
				return
			}
		}
	}

}

// Closing and reopening pooled readers must reset the record-directory offset.
func TestSeedPositionReaderReuse(t *testing.T) {
	files := []string{filepath.Join(t.TempDir(), "a.bin"), filepath.Join(t.TempDir(), "b.bin")}
	want := [][]uint32{{4, 8, 100}, {20, 60, 900}}
	for i, file := range files {
		writer, err := NewWriter(file, uint32(i))
		if err != nil {
			t.Fatal(err)
		}
		if err := writer.Write(want[i]); err != nil {
			t.Fatal(err)
		}
		if err := writer.Close(); err != nil {
			t.Fatal(err)
		}
	}
	for n := 0; n < 20; n++ {
		i := n % len(files)
		reader, err := NewReader(files[i])
		if err != nil {
			t.Fatal(err)
		}
		var got []uint32
		err = reader.SeedPositions(0, &got)
		closeErr := reader.Close()
		if err != nil || closeErr != nil {
			t.Fatalf("reopen %d: %v, %v", n, err, closeErr)
		}
		if !reflect.DeepEqual(got, want[i]) {
			t.Fatalf("reopen %d: got %v, want %v", n, got, want[i])
		}
	}
}

// Invalid offsets must not decode bytes at the file's previous cursor position.
func TestSeedPositionsInvalidOffset(t *testing.T) {
	for _, offset := range []uint64{1 << 63, ^uint64(0) - 7, ^uint64(0), 1<<63 - 3} {
		t.Run(fmt.Sprint(offset), func(t *testing.T) {
			file := filepath.Join(t.TempDir(), "positions.bin")
			w, err := NewWriter(file, 0)
			if err != nil {
				t.Fatal(err)
			}
			if err := w.Write([]uint32{1, 15, 300, 301}); err != nil {
				t.Fatal(err)
			}
			if err := w.Close(); err != nil {
				t.Fatal(err)
			}
			index := file + PositionsIndexFileExt
			data, err := os.ReadFile(index)
			if err != nil {
				t.Fatal(err)
			}
			be.PutUint64(data[24:32], offset)
			if err := os.WriteFile(index, data, 0600); err != nil {
				t.Fatal(err)
			}
			r, err := NewReader(file)
			if err != nil {
				t.Fatal(err)
			}
			defer r.Close()
			var positions []uint32
			if err := r.SeedPositions(0, &positions); err == nil {
				t.Fatalf("decoded invalid offset as positions %v", positions)
			}
		})
	}
}
