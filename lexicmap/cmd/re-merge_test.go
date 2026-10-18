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
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
	"github.com/shenwei356/lexichash"
)

// makeRemergeBatches supplies completed batch files with empty masks and large postings.
func makeRemergeBatches(t *testing.T, db string, inputs int) []string {
	t.Helper()
	const chunks, masksPerChunk = 2, 32
	lh, err := lexichash.NewWithSeed(16, chunks*masksPerChunk, 1, 0)
	if err != nil {
		t.Fatal(err)
	}
	paths := make([]string, inputs)
	for i := range paths {
		paths[i] = filepath.Join(db+ExtTmpDir, batchDir(i))
		for _, dir := range []string{filepath.Join(paths[i], DirSeeds), filepath.Join(paths[i], DirGenomes, batchDir(i))} {
			if err := os.MkdirAll(dir, 0755); err != nil {
				t.Fatal(err)
			}
		}
		for _, file := range []string{FileGenomeIndex, FileGenomeChunks, filepath.Join(DirGenomes, batchDir(i), "data")} {
			if err := os.WriteFile(filepath.Join(paths[i], file), []byte(fmt.Sprintf("batch %d\n", i)), 0644); err != nil {
				t.Fatal(err)
			}
		}
		if _, err := lh.WriteToFile(filepath.Join(paths[i], FileMasks)); err != nil {
			t.Fatal(err)
		}
		info := &IndexInfo{K: 16, Masks: chunks * masksPerChunk, Partitions: 4, Chunks: chunks,
			InputGenomes: 1, Genomes: 1, GenomeBatches: 1, InputBases: 100}
		if err := writeIndexInfo(filepath.Join(paths[i], FileInfo), info); err != nil {
			t.Fatal(err)
		}
		for chunk := 0; chunk < chunks; chunk++ {
			writer, err := kv.NewWriter(16, chunk*masksPerChunk, masksPerChunk, filepath.Join(paths[i], DirSeeds, chunkFile(chunk)), 3, 1, true)
			if err != nil {
				t.Fatal(err)
			}
			for mask := 0; mask < masksPerChunk; mask++ {
				data := make(map[uint64]*[]uint64)
				if mask%2 == 0 {
					values := make([]uint64, 400)
					for j := range values {
						values[j] = uint64(i)<<47 | uint64(j)<<1
					}
					data[uint64(mask)<<26] = &values
				}
				if err := writer.WriteDataOfAMask(data); err != nil {
					t.Fatal(err)
				}
			}
			if err := writer.Close(); err != nil {
				t.Fatal(err)
			}
		}
	}
	return paths
}

// TestRemergeSeedIndex2 compares the CLI with index's merge path across rounds and thresholds.
func TestRemergeSeedIndex2(t *testing.T) {
	for _, tc := range []struct {
		threshold string // empty disables adaptive indexes
		budget    int    // descriptor budget for the resumed merge
		inputs    int    // completed batches to resume
	}{{"", 12, 3}, {"4K", 1024, 3}, {"8K", 13, 3}, {"8K", 12, 1}} {
		t.Run(fmt.Sprintf("threshold=%s/budget=%d/inputs=%d", tc.threshold, tc.budget, tc.inputs), func(t *testing.T) {
			root := t.TempDir()
			want, got := filepath.Join(root, "want.lmi"), filepath.Join(root, "got.lmi")
			paths := makeRemergeBatches(t, want, tc.inputs)
			makeRemergeBatches(t, got, tc.inputs)
			opt := &IndexBuildingOptions{MaxOpenFiles: 1024, MergeThreads: 2}
			if tc.threshold != "" {
				var err error
				opt.SeedIndex2Threshold, err = parseIndex15Threshold(tc.threshold)
				if err != nil {
					t.Fatal(err)
				}
			}
			if err := mergeIndexes(nil, 3, 1, opt, 2, want, paths, want+ExtTmpDir, 1); err != nil {
				t.Fatal(err)
			}
			logfile := filepath.Join(root, "remerge.log")
			args := []string{"utils", "remerge", "--quiet", "-j", "2", "-J", "2", "-d", got, "--max-open-files", strconv.Itoa(tc.budget), "--log", logfile}
			if tc.threshold != "" {
				args = append(args, "-2", "--seed-index2-threshold", tc.threshold)
			}
			runSeedIndexCommand(t, args...)
			if _, err := os.Stat(got + ExtTmpDir); !os.IsNotExist(err) {
				t.Fatalf("temporary merge directory remains: %v", err)
			}
			summary := readIndex15Summary(t, logfile)
			if tc.threshold == "" {
				if len(summary) != 0 {
					t.Fatalf("ordinary merge logged adaptive counters: %v", summary)
				}
				for chunk := 0; chunk < 2; chunk++ {
					if _, err := os.Stat(filepath.Join(got, DirSeeds, chunkFile(chunk)) + kv.KVIndex15FileExt); !os.IsNotExist(err) {
						t.Fatalf("ordinary merge created a secondary index: %v", err)
					}
				}
			} else {
				// An independent scan verifies totals and confirms intermediate rounds were excluded.
				reindexLog := filepath.Join(root, "reindex.log")
				runSeedIndexCommand(t, "utils", "reindex-seeds2", "--quiet", "-j", "2", "-d", want, "--partitions", "4", "--threshold", tc.threshold, "--log", reindexLog)
				wantSummary := readIndex15Summary(t, reindexLog)
				if len(summary) < 2 || strings.Join(summary, "\n") != strings.Join(wantSummary, "\n") {
					t.Fatalf("remerge statistics differ from reindexing: got %v, want %v", summary, wantSummary)
				}
			}
			if err := filepath.WalkDir(want, func(path string, entry os.DirEntry, err error) error {
				if err != nil || entry.IsDir() {
					return err
				}
				relative, err := filepath.Rel(want, path)
				if err != nil {
					return err
				}
				read := func(file string) []byte {
					data, err := os.ReadFile(file)
					if err != nil {
						t.Fatal(err)
					}
					// Independently written adaptive headers have different pairing tokens.
					if tc.threshold != "" && strings.HasSuffix(file, kv.KVIndexFileExt) {
						clear(data[14:16])
					} else if strings.HasSuffix(file, kv.KVIndex15FileExt) {
						clear(data[46:48])
					}
					return data
				}
				if !bytes.Equal(read(path), read(filepath.Join(got, relative))) {
					t.Errorf("remerge file differs from index's merge: %s", relative)
				}
				return nil
			}); err != nil {
				t.Fatal(err)
			}
		})
	}
}

// TestRemergeRejectsAdaptiveBudget preserves batches when the third output cannot fit.
func TestRemergeRejectsAdaptiveBudget(t *testing.T) {
	db := filepath.Join(t.TempDir(), "index.lmi")
	paths := makeRemergeBatches(t, db, 3)
	args := []string{"utils", "remerge", "--quiet", "-j", "1", "-d", db, "-2", "--max-open-files", "12"}
	command := exec.Command(os.Args[0], append([]string{"-test.run=^TestSeedIndexCommand$", "--"}, args...)...)
	command.Env = append(os.Environ(), "LEXICMAP_SEED_INDEX_TEST_HELPER=1")
	output, err := command.CombinedOutput()
	if err == nil || !bytes.Contains(output, []byte("should be >= 13")) {
		t.Fatalf("unexpected result for insufficient adaptive budget: %v\n%s", err, output)
	}
	entries, err := os.ReadDir(db + ExtTmpDir)
	if err != nil || len(entries) != len(paths) {
		t.Fatalf("merge modified batches before validating the budget: %v, %v", entries, err)
	}
	for _, path := range paths {
		for _, file := range []string{FileMasks, filepath.Join(DirGenomes, filepath.Base(path), "data"), filepath.Join(DirSeeds, chunkFile(0))} {
			if _, err := os.Stat(filepath.Join(path, file)); err != nil {
				t.Fatal(err)
			}
		}
	}
}
