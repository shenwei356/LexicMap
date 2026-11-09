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
	"path/filepath"
	"testing"
)

// Warn only when the final index is empty, including after merging batches.
func TestIndexWarnsWithoutValidSequences(t *testing.T) {
	for _, tc := range []struct {
		name        string
		inputs      []string
		args        []string
		wantWarning bool
	}{
		{"empty", []string{""}, nil, true},
		{"short", []string{">s\nACGT\n"}, nil, true},
		{"filtered", []string{">s\n" + indexTestDNA(100) + "\n"}, []string{"-l", "200"}, true},
		{"empty_batches", []string{"", ">s\nACGT\n"}, []string{"-b", "1"}, true},
		{"mixed_batches", []string{"", ">s\n" + indexTestDNA(100) + "\n"}, []string{"-b", "1"}, false},
	} {
		t.Run(tc.name, func(t *testing.T) {
			dir := t.TempDir()
			db := filepath.Join(dir, "db")
			args := []string{"index", "-j", "1", "-m", "64", "-c", "2", "--partitions", "4", "-O", db}
			args = append(args, tc.args...)
			for i, data := range tc.inputs {
				path := filepath.Join(dir, string(rune('a'+i))+".fa")
				args = append(args, writeTestInput(t, path, data))
			}
			_, stderr, err := runAuditCommand(t, "", args...)
			if err != nil {
				t.Fatalf("index failed: %v, %s", err, stderr)
			}
			count := bytes.Count(stderr, []byte("no valid sequences were indexed; the index is empty"))
			want := 0
			if tc.wantWarning {
				want = 1
			}
			if count != want {
				t.Fatalf("warning count = %d, want %d: %s", count, want, stderr)
			}
			info, err := readIndexInfo(filepath.Join(db, FileInfo))
			if err != nil {
				t.Fatal(err)
			}
			if (info.Genomes == 0) != tc.wantWarning {
				t.Fatalf("unexpected genome count: %d", info.Genomes)
			}
		})
	}
}
