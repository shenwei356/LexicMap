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
	"reflect"
	"strconv"
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

// TestPrimaryIndexPartitionsInvalid keeps the installed generation searchable on rejection.
func TestPrimaryIndexPartitionsInvalid(t *testing.T) {
	for _, partitions := range []int{0, 1, 8, 1 << 20} {
		t.Run(strconv.Itoa(partitions), func(t *testing.T) {
			file := reindexFixture(t)
			sidecar, err := os.ReadFile(file + KVIndex15FileExt)
			if err != nil {
				t.Fatal(err)
			}
			if err = os.WriteFile(file+KVIndex15FileExt+index15BackupSuffix, sidecar, 0600); err != nil {
				t.Fatal(err)
			}
			before := make(map[string][]byte)
			for _, ext := range []string{"", KVIndexFileExt, KVIndex15FileExt, KVIndex15FileExt + index15BackupSuffix} {
				before[file+ext], err = os.ReadFile(file + ext)
				if err != nil {
					t.Fatal(err)
				}
			}
			want := publicationSearchResults(t, file)
			if len(want[0]) == 0 {
				t.Fatal("fixture produced no search hits")
			}
			if err = CreateKVIndex(file, partitions); err == nil {
				t.Fatal("expected invalid partitions to fail")
			}
			for path, data := range before {
				after, err := os.ReadFile(path)
				if err != nil || !bytes.Equal(data, after) {
					t.Fatalf("%s changed after rejection: %v", path, err)
				}
			}
			entries, err := os.ReadDir(filepath.Dir(file))
			if err != nil || len(entries) != len(before) {
				t.Fatalf("unexpected files after rejection: %v, %v", entries, err)
			}
			if got := publicationSearchResults(t, file); !reflect.DeepEqual(got, want) {
				t.Fatalf("search results changed after rejection: got %v, want %v", got, want)
			}
		})
	}
}

// TestPrimaryIndexPartitionsValid checks legacy bytes and the inclusive primary-prefix limit.
func TestPrimaryIndexPartitionsValid(t *testing.T) {
	for _, settings := range []struct {
		partitions int
		prefix     uint8
	}{{4, 1}, {4096, 6}, {1 << 18, 9}} {
		t.Run(strconv.Itoa(settings.partitions), func(t *testing.T) {
			file := reindexFixture(t)
			want := publicationSearchResults(t, file)[2] // exact, full-length query
			r, err := NewReader(file)
			if err != nil {
				t.Fatal(err)
			}
			data, err := r.ReadDataOfAMaskAsMap()
			if err != nil {
				t.Fatal(err)
			}
			if err = r.Close(); err != nil {
				t.Fatal(err)
			}
			// The unchanged writer supplies the reference single-level disk format.
			baseline := filepath.Join(t.TempDir(), "seeds")
			if _, err = WriteKVData(16, 0, []*map[uint64]*[]uint64{data}, baseline, 7, settings.prefix, 512, false); err != nil {
				t.Fatal(err)
			}
			sidecar, err := os.ReadFile(file + KVIndex15FileExt)
			if err != nil {
				t.Fatal(err)
			}
			if err = os.WriteFile(file+KVIndex15FileExt+index15BackupSuffix, sidecar, 0600); err != nil {
				t.Fatal(err)
			}
			if err = os.Chmod(file+KVIndexFileExt, 0640); err != nil {
				t.Fatal(err)
			}
			if err = CreateKVIndex(file, settings.partitions); err != nil {
				t.Fatal(err)
			}
			info, err := os.Stat(file + KVIndexFileExt)
			if err != nil || info.Mode().Perm() != 0640 {
				t.Fatalf("primary index permissions changed: %v, %v", info, err)
			}
			entries, err := os.ReadDir(filepath.Dir(file))
			if err != nil || len(entries) != 2 {
				t.Fatalf("obsolete sidecars or temporary files remain: %v, %v", entries, err)
			}
			for _, ext := range []string{"", KVIndexFileExt} {
				gotBytes, err := os.ReadFile(file + ext)
				if err != nil {
					t.Fatal(err)
				}
				wantBytes, err := os.ReadFile(baseline + ext)
				if err != nil || !bytes.Equal(gotBytes, wantBytes) {
					t.Fatalf("%s differs from the legacy writer: %v", ext, err)
				}
			}
			disk, err := NewSearcher(file, 1)
			if err != nil {
				t.Fatal(err)
			}
			defer disk.Close()
			memory, err := NewInMemomrySearcher(file)
			if err != nil {
				t.Fatal(err)
			}
			defer memory.Close()
			for _, search := range []func([]uint64, uint8, bool, bool) (*[]SearchResult, error){disk.Search, memory.Search} {
				results, err := search([]uint64{uint64(54) << 18}, 16, false, false)
				if err != nil {
					t.Fatal(err)
				}
				if got := snapshotIndex15Results(results, false); len(got) == 0 || !reflect.DeepEqual(got, want) {
					t.Fatalf("search results changed after rebuilding: got %v, want %v", got, want)
				}
			}
		})
	}
}

// TestPrimaryIndexPartitionsMaskBeyondK guards the unsigned subtraction in validation.
func TestPrimaryIndexPartitionsMaskBeyondK(t *testing.T) {
	file := reindexFixture(t)
	primary, err := os.ReadFile(file + KVIndexFileExt)
	if err != nil {
		t.Fatal(err)
	}
	primary[11] = 17 // malformed mask prefix exceeds K=16
	if err = os.WriteFile(file+KVIndexFileExt, primary, 0600); err != nil {
		t.Fatal(err)
	}
	if err = CreateKVIndex(file, 4); err == nil {
		t.Fatal("expected a mask prefix beyond K to fail")
	}
	after, err := os.ReadFile(file + KVIndexFileExt)
	if err != nil || !bytes.Equal(primary, after) {
		t.Fatalf("malformed primary index changed: %v", err)
	}
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
