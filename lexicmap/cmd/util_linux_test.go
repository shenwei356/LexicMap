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

package cmd

import (
	"bytes"
	"compress/gzip"
	"os"
	"path/filepath"
	"runtime/debug"
	"strings"
	"testing"
)

// Parsing errors must close the input immediately, without relying on file
// finalizers or leaving compressed-stream read-ahead alive.
func TestReadKVsClosesInputOnError(t *testing.T) {
	gcPercent := debug.SetGCPercent(-1)
	defer debug.SetGCPercent(gcPercent)
	for _, compressed := range []bool{false, true} {
		suffix := "/plain"
		if compressed {
			suffix = "/gzip"
		}
		for _, kind := range []string{"text scan", "integer scan", "integer parse"} {
			t.Run(kind+suffix, func(t *testing.T) {
				data := []byte("key\t" + strings.Repeat("1", 70<<10) + "\n")
				if kind == "integer parse" {
					data = []byte("key\tnot-a-taxid\n")
				}
				if compressed {
					var buf bytes.Buffer
					w := gzip.NewWriter(&buf)
					if _, err := w.Write(data); err != nil {
						t.Fatal(err)
					}
					if err := w.Close(); err != nil {
						t.Fatal(err)
					}
					data = buf.Bytes()
				}
				file := filepath.Join(t.TempDir(), "mapping")
				if err := os.WriteFile(file, data, 0600); err != nil {
					t.Fatal(err)
				}
				var err error
				if kind == "text scan" {
					_, err = readKVs(file, false)
				} else {
					_, err = readKVsUint32(file, false, true)
				}
				if err == nil {
					t.Fatal("expected a mapping read error")
				}
				entries, err := os.ReadDir("/proc/self/fd")
				if err != nil {
					t.Fatal(err)
				}
				for _, entry := range entries {
					target, err := os.Readlink(filepath.Join("/proc/self/fd", entry.Name()))
					if err == nil && target == file {
						t.Errorf("mapping remains open after %s error", kind)
					}
				}
			})
		}
	}
}
