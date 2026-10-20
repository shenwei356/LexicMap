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
	"os/exec"
	"path/filepath"
	"slices"
	"strconv"
	"strings"
	"testing"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
)

// TestReindexSeedsInvalidPartitions preserves both indexes, backups, and CLI metadata.
func TestReindexSeedsInvalidPartitions(t *testing.T) {
	dir := t.TempDir()
	seedsDir := filepath.Join(dir, DirSeeds)
	if err := os.MkdirAll(seedsDir, 0700); err != nil {
		t.Fatal(err)
	}
	file := filepath.Join(seedsDir, chunkFile(0))
	values := make([]uint64, 1024)
	for i := range values {
		values[i] = uint64(i + 1)
	}
	data := map[uint64]*[]uint64{uint64(54) << 18: &values}
	if _, err := kv.WriteKVData(16, 0, []*map[uint64]*[]uint64{&data}, file, 7, 6, 512, false); err != nil {
		t.Fatal(err)
	}
	if _, err := kv.CreateKVIndex15(file, 4096, kv.MinIndex15Threshold); err != nil {
		t.Fatal(err)
	}
	sidecar, err := os.ReadFile(file + kv.KVIndex15FileExt)
	if err != nil {
		t.Fatal(err)
	}
	if err = os.WriteFile(file+kv.KVIndex15FileExt+".bak", sidecar, 0600); err != nil {
		t.Fatal(err)
	}
	infoFile := filepath.Join(dir, FileInfo)
	if err = writeIndexInfo(infoFile, &IndexInfo{MainVersion: MainVersion, K: 16, Masks: 1, Chunks: 1, Partitions: 4096}); err != nil {
		t.Fatal(err)
	}
	before := make(map[string][]byte)
	for _, path := range []string{file, file + kv.KVIndexFileExt, file + kv.KVIndex15FileExt, file + kv.KVIndex15FileExt + ".bak", infoFile} {
		before[path], err = os.ReadFile(path)
		if err != nil {
			t.Fatal(err)
		}
	}
	for _, test := range []struct {
		partitions int
		message    string
	}{
		{0, "value of flag --partitions should be greater than 0"},
		{1, "partitions should be a power of 4 and at least 4"},
		{8, "partitions should be a power of 4"},
		{1 << 20, "mask prefix (7) + anchor prefix (10) should not exceed k (16)"},
	} {
		t.Run(strconv.Itoa(test.partitions), func(t *testing.T) {
			command := exec.Command(os.Args[0], "-test.run=^TestSeedIndexCommand$", "--", "utils", "reindex-seeds", "-j", "1", "-d", dir, "--partitions", strconv.Itoa(test.partitions))
			command.Env = append(os.Environ(), "LEXICMAP_SEED_INDEX_TEST_HELPER=1")
			output, err := command.CombinedOutput()
			if err == nil {
				t.Fatalf("invalid partitions succeeded: %s", output)
			}
			if !strings.Contains(string(output), test.message) {
				t.Fatalf("missing validation error: %v, %s", err, output)
			}
			if strings.Contains(string(output), "finished updating") {
				t.Fatalf("failed rebuild announced success: %s", output)
			}
			for path, want := range before {
				got, err := os.ReadFile(path)
				if err != nil || !bytes.Equal(got, want) {
					t.Fatalf("%s changed after rejection: %v", path, err)
				}
			}
			entries, err := os.ReadDir(seedsDir)
			if err != nil || len(entries) != 4 {
				t.Fatalf("unexpected seed files after rejection: %v, %v", entries, err)
			}
			entries, err = os.ReadDir(dir)
			if err != nil || len(entries) != 2 {
				t.Fatalf("unexpected metadata files after rejection: %v, %v", entries, err)
			}
			searcher, err := kv.NewSearcher(file, 1)
			if err != nil {
				t.Fatal(err)
			}
			defer searcher.Close()
			results, err := searcher.Search([]uint64{uint64(54) << 18}, 16, false, false)
			if err != nil {
				t.Fatal(err)
			}
			defer kv.RecycleSearchResults(results)
			if len(*results) != 1 || !slices.Equal((*results)[0].Values, values) {
				t.Fatal("old seed index no longer returns the original values")
			}
		})
	}
}
