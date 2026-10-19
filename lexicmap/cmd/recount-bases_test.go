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
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// writeRecountFixture supplies tiny independent batch files without building seed indexes.
func writeRecountFixture(t *testing.T, dir string, batches int) *IndexInfo {
	t.Helper()
	file := writeSharedReaderFixture(t, dir)
	for _, suffix := range []string{"", ".idx"} {
		data, err := os.ReadFile(file + suffix)
		if err != nil {
			t.Fatal(err)
		}
		for i := 1; i < batches; i++ {
			path := filepath.Join(dir, DirGenomes, batchDir(i), FileGenomes+suffix)
			if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
				t.Fatal(err)
			}
			if err := os.WriteFile(path, data, 0600); err != nil {
				t.Fatal(err)
			}
		}
	}
	info := &IndexInfo{MainVersion: MainVersion, GenomeBatches: batches, InputBases: 77}
	if err := writeIndexInfo(filepath.Join(dir, FileInfo), info); err != nil {
		t.Fatal(err)
	}
	return info
}

func TestUpdateInputBases(t *testing.T) {
	for _, failure := range []string{"", "data", "index"} {
		broken := failure != ""
		dir := t.TempDir()
		info := writeRecountFixture(t, dir, 4)
		fileInfo := filepath.Join(dir, FileInfo)
		original, err := os.ReadFile(fileInfo)
		if err != nil {
			t.Fatal(err)
		}
		switch failure {
		case "data":
			if err := os.Truncate(filepath.Join(dir, DirGenomes, batchDir(2), FileGenomes), 0); err != nil {
				t.Fatal(err)
			}
		case "index":
			if err := os.Remove(filepath.Join(dir, DirGenomes, batchDir(2), FileGenomes+".idx")); err != nil {
				t.Fatal(err)
			}
		}
		got, err := updateInputBases(info, dir, 3)
		if broken {
			if err == nil || !strings.Contains(err.Error(), batchDir(2)) || got != 0 || info.InputBases != 77 {
				t.Fatalf("failed recount published a total: %d, %+v, %v", got, info, err)
			}
			after, err := os.ReadFile(fileInfo)
			if err != nil || !bytes.Equal(after, original) {
				t.Fatalf("failed recount changed index info: %s, %v", after, err)
			}
		} else {
			if err != nil || got != 32 || info.InputBases != got {
				t.Fatalf("recount: %d, %+v, %v", got, info, err)
			}
			saved, err := readIndexInfo(fileInfo)
			if err != nil || saved.InputBases != got {
				t.Fatalf("saved total: %+v, %v", saved, err)
			}
		}
	}
}
