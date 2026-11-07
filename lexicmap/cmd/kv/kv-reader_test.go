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
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"testing"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
)

// writeReaderCountFixture includes a populated mask followed by an empty mask.
func writeReaderCountFixture(t *testing.T, narrow bool) string {
	t.Helper()
	file := filepath.Join(t.TempDir(), "seeds")
	w, err := NewWriter(5, 0, 2, file, 1, 2, narrow)
	if err != nil {
		t.Fatal(err)
	}
	values := []uint64{1, 2}
	writeErr := w.WriteDataOfAMask(map[uint64]*[]uint64{7: &values})
	if writeErr == nil {
		writeErr = w.WriteDataOfAMask(nil)
	}
	if err := errors.Join(writeErr, w.Close()); err != nil {
		t.Fatal(err)
	}
	return file
}

// Corrupt counts must fail before allocating, for every sequential reader API.
func TestReaderRejectsOversizedCounts(t *testing.T) {
	for _, narrow := range []bool{false, true} {
		for _, offset := range []int{16, 24, 32} {
			for _, count := range []uint64{1 << 20, 1 << 62, 1 << 63, ^uint64(0)} {
				if offset == 16 && count <= uint64(^uint(0)>>1)-2 {
					continue // Representable first-mask offsets are valid library inputs.
				}
				t.Run(fmt.Sprintf("narrow%t/offset%d/count%d", narrow, offset, count), func(t *testing.T) {
					file := writeReaderCountFixture(t, narrow)
					data, err := os.ReadFile(file)
					if err != nil {
						t.Fatal(err)
					}
					be.PutUint64(data[offset:offset+8], count)
					if err := os.WriteFile(file, data, 0600); err != nil {
						t.Fatal(err)
					}
					if _, err := NewInMemomrySearcher(file); !errors.Is(err, ErrBrokenFile) {
						t.Fatalf("in-memory load: got %v, want broken file", err)
					}
					for _, mode := range []string{"map", "append", "list", "indexed-list"} {
						r, err := NewReader(file)
						if err == nil {
							switch mode {
							case "map":
								_, err = r.ReadDataOfAMaskAsMap()
							case "append":
								data := make(map[uint64]*[]uint64)
								err = r.ReadDataOfAMaskAndAppendToMap(&data)
							case "list":
								_, err = r.ReadDataOfAMaskAsList()
							case "indexed-list":
								_, _, _, _, err = r.ReadDataOfAMaskAsListAndCreateIndex()
							}
							if closeErr := r.Close(); closeErr != nil {
								t.Fatal(closeErr)
							}
						}
						if !errors.Is(err, ErrBrokenFile) {
							t.Fatalf("%s: got %v, want broken file", mode, err)
						}
					}
				})
			}
		}
	}
}

// Posting lengths also reach slices.Grow when merging masks into an existing map.
func TestReaderRejectsOversizedPostingCount(t *testing.T) {
	for _, narrow := range []bool{false, true} {
		file := writeReaderCountFixture(t, narrow)
		data, err := os.ReadFile(file)
		if err != nil {
			t.Fatal(err)
		}
		// Header + mask count + key control and its two one-byte deltas.
		const countsOffset = 43
		var encoded [17]byte
		control, n := util.PutUint64s(encoded[1:], 1<<62, 0)
		encoded[0] = control
		bad := append(append(append([]byte(nil), data[:countsOffset]...), encoded[:n+1]...), data[countsOffset+3:]...)
		if err := os.WriteFile(file, bad, 0600); err != nil {
			t.Fatal(err)
		}
		r, err := NewReader(file)
		if err != nil {
			t.Fatal(err)
		}
		m := make(map[uint64]*[]uint64)
		err = r.ReadDataOfAMaskAndAppendToMap(&m)
		if closeErr := r.Close(); closeErr != nil {
			t.Fatal(closeErr)
		}
		if !errors.Is(err, ErrBrokenFile) {
			t.Fatalf("narrow=%t: got %v, want broken file", narrow, err)
		}
	}
}

// The lightweight metadata path used by in-memory loading must bound its table too.
func TestInMemorySearcherRejectsInvalidIndexPrefix(t *testing.T) {
	for _, prefix := range []byte{0, 5, 32, 255} {
		file := writeReaderCountFixture(t, true)
		indexFile := file + KVIndexFileExt
		data, err := os.ReadFile(indexFile)
		if err != nil {
			t.Fatal(err)
		}
		data[12] = prefix
		if err := os.WriteFile(indexFile, data, 0600); err != nil {
			t.Fatal(err)
		}
		if _, err := NewInMemomrySearcher(file); !errors.Is(err, ErrInvalidFileFormat) {
			t.Fatalf("anchor prefix=%d: got %v, want invalid format", prefix, err)
		}
	}
}
