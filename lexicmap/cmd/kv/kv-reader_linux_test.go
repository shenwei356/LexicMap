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
	"fmt"
	"io"
	"os"
	"path/filepath"
	"runtime/debug"
	"testing"
)

// Constructors must close rejected inputs immediately, not through finalizers.
func TestReaderConstructorClosesRejectedFiles(t *testing.T) {
	gcPercent := debug.SetGCPercent(-1)
	defer debug.SetGCPercent(gcPercent)
	countOpen := func(file string) int {
		t.Helper()
		entries, err := os.ReadDir("/proc/self/fd")
		if err != nil {
			t.Fatal(err)
		}
		n := 0
		for _, entry := range entries {
			target, err := os.Readlink(filepath.Join("/proc/self/fd", entry.Name()))
			if err == nil && target == file {
				n++
			}
		}
		return n
	}
	for _, index := range []bool{false, true} {
		name, size := "data", 32
		magic := Magic
		open := func(file string) (io.Closer, error) {
			r, err := NewReader(file)
			if err != nil {
				return nil, err
			}
			return r, nil
		}
		if index {
			name, size, magic = "index", 40, MagicIdx
			open = func(file string) (io.Closer, error) {
				r, err := NewIndexReader(file)
				if err != nil {
					return nil, err
				}
				return r, nil
			}
		}
		header := make([]byte, size)
		copy(header, magic[:])
		header[8], header[10] = MainVersion, 5
		be.PutUint64(header[24:32], 1)
		cases := make(map[string][]byte)
		for n := 0; n < len(header); n++ {
			cases[fmt.Sprintf("truncated-%02d", n)] = header[:n]
		}
		badMagic := append([]byte(nil), header...)
		badMagic[0] ^= 0xff
		cases["magic"] = badMagic
		badVersion := append([]byte(nil), header...)
		badVersion[8] = 255
		cases["version"] = badVersion
		if index {
			badMetadata := append([]byte(nil), header...)
			badMetadata[13] = MaskHasIndex15 // tagged offsets require index v2
			cases["metadata"] = badMetadata
		}
		for kind, data := range cases {
			t.Run(name+"/"+kind, func(t *testing.T) {
				file := filepath.Join(t.TempDir(), "seeds")
				if err := os.WriteFile(file, data, 0600); err != nil {
					t.Fatal(err)
				}
				r, err := open(file)
				if err == nil {
					r.Close()
					t.Fatal("expected invalid header to fail")
				}
				if n := countOpen(file); n != 0 {
					t.Errorf("rejected input remains open: %d descriptors", n)
				}
			})
		}
		t.Run(name+"/valid", func(t *testing.T) {
			file := filepath.Join(t.TempDir(), "seeds")
			if err := os.WriteFile(file, header, 0600); err != nil {
				t.Fatal(err)
			}
			r, err := open(file)
			if err != nil {
				t.Fatal(err)
			}
			defer r.Close()
			if n := countOpen(file); n != 1 {
				t.Fatalf("reader does not own input: %d", n)
			}
			if err := r.Close(); err != nil {
				t.Fatal(err)
			}
			if n := countOpen(file); n != 0 {
				t.Fatalf("closed input remains open: %d", n)
			}
		})
	}
}
