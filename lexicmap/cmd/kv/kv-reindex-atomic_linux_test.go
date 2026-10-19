//go:build linux

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
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"syscall"
	"testing"
)

func TestReindexWriteFailurePreservesOldIndexes(t *testing.T) {
	if file := os.Getenv("LEXICMAP_REINDEX_WRITE_FAILURE"); file != "" {
		// Limit only the child: the buffered flush fails after the new file opens.
		if err := syscall.Setrlimit(syscall.RLIMIT_FSIZE, &syscall.Rlimit{Cur: 16, Max: 16}); err != nil {
			t.Fatal(err)
		}
		if err := CreateKVIndex(file, 4096); err == nil {
			t.Fatal("expected output size limit to fail")
		}
		return
	}
	file := reindexFixture(t)
	before, err := os.ReadFile(file + KVIndexFileExt)
	if err != nil {
		t.Fatal(err)
	}
	before15, err := os.ReadFile(file + KVIndex15FileExt)
	if err != nil {
		t.Fatal(err)
	}
	command := exec.Command(os.Args[0], "-test.run=^TestReindexWriteFailurePreservesOldIndexes$")
	command.Env = append(os.Environ(), "LEXICMAP_REINDEX_WRITE_FAILURE="+file)
	if out, err := command.CombinedOutput(); err != nil {
		t.Fatalf("fault injection failed: %v\n%s", err, out)
	}
	after, err := os.ReadFile(file + KVIndexFileExt)
	if err != nil || !bytes.Equal(before, after) {
		t.Fatalf("primary changed: %v", err)
	}
	after15, err := os.ReadFile(file + KVIndex15FileExt)
	if err != nil || !bytes.Equal(before15, after15) {
		t.Fatalf("secondary changed: %v", err)
	}
	files, err := filepath.Glob(filepath.Join(filepath.Dir(file), ".seeds.idx-*"))
	if err != nil || len(files) != 0 {
		t.Fatalf("temporary indexes remain: %v, %v", files, err)
	}
	reader, err := NewSearcher(file, 1)
	if err != nil {
		t.Fatal(err)
	}
	if err = reader.Close(); err != nil {
		t.Fatal(err)
	}
}
