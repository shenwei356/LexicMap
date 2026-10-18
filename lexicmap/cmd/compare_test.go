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
	"os"
	"path/filepath"
	"reflect"
	"testing"
)

func TestCombinationCount(t *testing.T) {
	if got, err := combinationCount(3, "genome IDs"); err != nil || got != 3 {
		t.Fatalf("count: got %d, err=%v, want 3", got, err)
	}
	if _, err := combinationCount(1, "genome IDs"); err == nil {
		t.Fatal("expected an error for fewer than two genome IDs")
	}
}

func TestReadPairFilesSkipsHeaderInEachFile(t *testing.T) {
	dir := t.TempDir()
	file1 := filepath.Join(dir, "pairs1.tsv")
	file2 := filepath.Join(dir, "pairs2.tsv")
	if err := os.WriteFile(file1, []byte("genome1\tgenome2\nA\tB\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(file2, []byte("genome1\tgenome2\nC\tD\textra\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	got, err := readPairFiles([]string{file1, file2}, true)
	if err != nil {
		t.Fatal(err)
	}
	want := []string{"A", "B", "C", "D"}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("pairs: got %v, want %v", got, want)
	}
}
