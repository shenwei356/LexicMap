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
	"encoding/binary"
	"io/fs"
	"os"
	"os/exec"
	"path/filepath"
	"reflect"
	"slices"
	"strconv"
	"strings"
	"testing"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
)

// Reject a later chunk's bad header before replacing any earlier chunk index.
func TestReindexSeedsRejectsInconsistentLayout(t *testing.T) {
	for _, command := range []string{"reindex-seeds", "reindex-seeds2"} {
		for _, suffix := range []string{"", kv.KVIndexFileExt} {
			t.Run(command+suffix, func(t *testing.T) {
				dir := t.TempDir()
				input := writeTestInput(t, filepath.Join(dir, "g.fa"), ">s\n"+indexTestDNA(500)+"\n")
				db := filepath.Join(dir, "db")
				runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-c", "2", "--partitions", "4", "-O", db, input)

				file := filepath.Join(db, DirSeeds, chunkFile(1)) + suffix
				original, err := os.ReadFile(file)
				if err != nil {
					t.Fatal(err)
				}
				data := bytes.Clone(original)
				binary.BigEndian.PutUint64(data[16:24], 0) // overlaps the first chunk
				if err := os.WriteFile(file, data, 0600); err != nil {
					t.Fatal(err)
				}
				t.Cleanup(func() {
					if err := os.WriteFile(file, original, 0600); err != nil {
						t.Error(err)
					}
				})
				before := make(map[string][]byte)
				if err := filepath.WalkDir(db, func(path string, entry fs.DirEntry, err error) error {
					if err != nil || entry.IsDir() {
						return err
					}
					before[path], err = os.ReadFile(path)
					return err
				}); err != nil {
					t.Fatal(err)
				}
				_, stderr, err := runAuditCommand(t, "", "utils", command, "--quiet", "-j", "1", "-d", db, "--partitions", "16")
				if err == nil || bytes.Contains(stderr, []byte("panic:")) {
					t.Errorf("expected normal chunk validation error: %v\n%s", err, stderr)
				}
				after := make(map[string][]byte)
				if err := filepath.WalkDir(db, func(path string, entry fs.DirEntry, err error) error {
					if err != nil || entry.IsDir() {
						return err
					}
					after[path], err = os.ReadFile(path)
					return err
				}); err != nil {
					t.Fatal(err)
				}
				if !reflect.DeepEqual(before, after) {
					t.Fatal("bad chunk metadata changed index files or left staging files")
				}
			})
		}
	}
}

// Reject invalid work counts before workers or progress bars can rewrite metadata.
func TestReindexSeedsRejectsInvalidSummary(t *testing.T) {
	dir := t.TempDir()
	input := writeTestInput(t, filepath.Join(dir, "g.fa"), ">s\n"+indexTestDNA(500)+"\n")
	db := filepath.Join(dir, "db")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-k", "16", "-m", "64", "-c", "2", "--partitions", "4", "-O", db, input)
	summary := filepath.Join(db, FileInfo)
	original, err := os.ReadFile(summary)
	if err != nil {
		t.Fatal(err)
	}
	// Include file names as well as contents to detect stray staging/backup files.
	snapshot := func(t *testing.T) map[string][]byte {
		t.Helper()
		files := make(map[string][]byte)
		if err := filepath.WalkDir(db, func(path string, entry fs.DirEntry, err error) error {
			if err != nil || entry.IsDir() {
				return err
			}
			files[path], err = os.ReadFile(path)
			return err
		}); err != nil {
			t.Fatal(err)
		}
		return files
	}
	for _, command := range []string{"reindex-seeds", "reindex-seeds2"} {
		for _, tc := range []struct{ name, old, replacement string }{
			{"zero chunks", "chunks = 2", "chunks = 0"},
			{"negative chunks", "chunks = 2", "chunks = -1"},
			{"missing chunks", "chunks = 2", ""},
			{"chunks exceed masks", "chunks = 2", "chunks = 65"},
			{"chunks exceed format limit", "chunks = 2", "chunks = 129"},
			{"huge chunks", "chunks = 2", "chunks = 2147483647"},
			{"zero masks", "masks = 64", "masks = 0"},
			{"zero k", "max-K = 16", "max-K = 0"},
			{"oversized k", "max-K = 16", "max-K = 33"},
			{"wrong version", "main-version = 3", "main-version = 0"},
		} {
			t.Run(command+"/"+tc.name, func(t *testing.T) {
				data := strings.Replace(string(original), tc.old, tc.replacement, 1)
				if data == string(original) {
					t.Fatalf("missing fixture field %q", tc.old)
				}
				if err := os.WriteFile(summary, []byte(data), 0600); err != nil {
					t.Fatal(err)
				}
				before := snapshot(t)
				out, stderr, err := runAuditCommand(t, "", "utils", command, "--quiet", "-j", "2", "-d", db, "--partitions", "16")
				message := string(out) + string(stderr)
				if err == nil || strings.Contains(message, "panic:") ||
					!(strings.Contains(message, "invalid index summary") || strings.Contains(message, "versions do not match")) {
					t.Errorf("expected a normal summary error: %v, %s", err, message)
				}
				if !reflect.DeepEqual(snapshot(t), before) {
					t.Error("rejected summary changed index files")
				}
			})
		}
	}
}

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
