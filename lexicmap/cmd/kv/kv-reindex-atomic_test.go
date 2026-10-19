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
	"path/filepath"
	"testing"
)

// reindexFixture supplies a valid adaptive index whose old files must survive errors.
func reindexFixture(t *testing.T) string {
	t.Helper()
	file := filepath.Join(t.TempDir(), "seeds")
	// Exceed the secondary-block threshold so recovery tests exercise lookups,
	// rather than only validating an otherwise unused sidecar header.
	values := make([]uint64, 1024)
	for i := range values {
		values[i] = uint64(i + 1)
	}
	m := map[uint64]*[]uint64{uint64(54) << 18: &values}
	if _, err := WriteKVData(16, 0, []*map[uint64]*[]uint64{&m}, file, 7, 6, 512, false); err != nil {
		t.Fatal(err)
	}
	if _, err := CreateKVIndex15(file, 4096, MinIndex15Threshold); err != nil {
		t.Fatal(err)
	}
	return file
}

func TestReindexReadFailurePreservesOldIndexes(t *testing.T) {
	file := reindexFixture(t)
	before, err := os.ReadFile(file + KVIndexFileExt)
	if err != nil {
		t.Fatal(err)
	}
	before15, err := os.ReadFile(file + KVIndex15FileExt)
	if err != nil {
		t.Fatal(err)
	}
	data, err := os.ReadFile(file)
	if err != nil {
		t.Fatal(err)
	}
	if err = os.Truncate(file, int64(len(data)-1)); err != nil {
		t.Fatal(err)
	}
	if err = CreateKVIndex(file, 4096); err == nil {
		t.Fatal("expected truncated seed data to fail")
	}
	after, err := os.ReadFile(file + KVIndexFileExt)
	if err != nil || !bytes.Equal(before, after) {
		t.Fatalf("primary index changed: %v", err)
	}
	after15, err := os.ReadFile(file + KVIndex15FileExt)
	if err != nil || !bytes.Equal(before15, after15) {
		t.Fatalf("secondary index changed: %v", err)
	}
	if err = os.WriteFile(file, data, 0600); err != nil {
		t.Fatal(err)
	}
	reader, err := NewSearcher(file, 1)
	if err != nil {
		t.Fatal(err)
	}
	if err = reader.Close(); err != nil {
		t.Fatal(err)
	}
	files, err := filepath.Glob(filepath.Join(filepath.Dir(file), ".seeds.idx-*"))
	if err != nil || len(files) != 0 {
		t.Fatalf("temporary indexes remain: %v, %v", files, err)
	}
	// A successful retry replaces the primary and removes its obsolete secondary.
	if err = CreateKVIndex(file, 4096); err != nil {
		t.Fatal(err)
	}
	if _, err = os.Stat(file + KVIndex15FileExt); !os.IsNotExist(err) {
		t.Fatalf("obsolete secondary remains: %v", err)
	}
	reader, err = NewSearcher(file, 1)
	if err != nil {
		t.Fatal(err)
	}
	if err = reader.Close(); err != nil {
		t.Fatal(err)
	}
}
