package cmd

import (
	"bytes"
	"fmt"
	"os"
	"path/filepath"
	"testing"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
)

// TestPlanIndexMerge covers group-size and concurrency boundaries and invalid budgets.
func TestPlanIndexMerge(t *testing.T) {
	for _, tc := range []struct {
		// Input count and budget settings, followed by the expected group size and workers.
		inputs, files, threads, batch, workers int
	}{
		{10, 1024, 8, 10, 8},
		{1014, 1024, 8, 1014, 1},
		{1025, 1024, 8, 1014, 1},
		{11, 1024, 8, 11, 8}, // the last group has its own concurrency budget
		{40, 32, 8, 22, 1},
		{2, 32, 8, 2, 6},
		{40, 12, 8, 2, 1}, // minimum budget still makes progress
		{1, 11, 8, 1, 1},
	} {
		batch, workers, err := planIndexMerge(tc.inputs, tc.files, tc.threads)
		if err != nil || batch != tc.batch || workers != tc.workers {
			t.Fatalf("plan(%d, %d, %d) = (%d, %d, %v), want (%d, %d)", tc.inputs, tc.files, tc.threads, batch, workers, err, tc.batch, tc.workers)
		}
	}
	for _, tc := range [][3]int{{0, 1024, 8}, {10, 1024, 0}, {2, 11, 1}, {2, 1, 1}, {1, 10, 1}} {
		if _, _, err := planIndexMerge(tc[0], tc[1], tc[2]); err == nil {
			t.Fatalf("accepted invalid plan %v", tc)
		}
	}
}

// TestIndexMergeRoundsFitFileBudget checks budget compliance and progress in every round.
func TestIndexMergeRoundsFitFileBudget(t *testing.T) {
	for _, files := range []int{12, 13, 16, 32, 1024} {
		for _, inputs := range []int{2, 3, 40, 1014, 1015, 1025, 100000} {
			for _, threads := range []int{1, 8, 128} {
				for n := inputs; n > 1; {
					batch, _, err := planIndexMerge(n, files, threads)
					if err != nil {
						t.Fatal(err)
					}
					next := (n + batch - 1) / batch
					if next >= n {
						t.Fatalf("merge does not progress: %d -> %d", n, next)
					}
					for start := 0; start < n; start += batch {
						size := min(batch, n-start)
						_, workers, err := planIndexMerge(size, files, threads)
						if err != nil || workers < 1 || workers*(size+2)+indexMergeReservedFiles > files {
							t.Fatalf("group exceeds budget: inputs=%d workers=%d files=%d err=%v", size, workers, files, err)
						}
					}
					n = next
				}
			}
		}
	}
}

// TestBuildIndexRejectsInsufficientMergeBudgetBeforeBuilding prevents wasted batch construction.
func TestBuildIndexRejectsInsufficientMergeBudgetBeforeBuilding(t *testing.T) {
	outdir := filepath.Join(t.TempDir(), "index")
	opt := &IndexBuildingOptions{GenomeBatchSize: 1, MaxOpenFiles: 11, MergeThreads: 1}
	if err := BuildIndex(outdir, []string{"missing1.fna", "missing2.fna"}, opt); err == nil {
		t.Fatal("accepted a budget too small to merge two indexes")
	}
	if _, err := os.Stat(outdir + ExtTmpDir); !os.IsNotExist(err) {
		t.Fatalf("created batches before validating the budget: %v", err)
	}
}

// TestMergeIndexesMultipleRoundsMatchesSingleRound preserves seed and metadata bytes across rounds.
func TestMergeIndexesMultipleRoundsMatchesSingleRound(t *testing.T) {
	for _, narrow := range []bool{false, true} {
		t.Run(fmt.Sprintf("three-byte-positions=%t", narrow), func(t *testing.T) {
			// Small fixtures exercise multiple workers and multiple merge rounds.
			const inputs, chunks, masksPerChunk = 17, 3, 2
			root := t.TempDir()
			outputs := make([]string, 2)
			for run, budget := range []int{1024, 16} {
				tmp := filepath.Join(root, fmt.Sprint(run))
				paths := make([]string, inputs)
				for i := range paths {
					paths[i] = filepath.Join(tmp, batchDir(i))
					genomes := filepath.Join(paths[i], DirGenomes, batchDir(i))
					for _, dir := range []string{genomes, filepath.Join(paths[i], DirSeeds)} {
						if err := os.MkdirAll(dir, 0755); err != nil {
							t.Fatal(err)
						}
					}
					for _, file := range []string{FileGenomeIndex, FileGenomeChunks, FileMasks, filepath.Join(DirGenomes, batchDir(i), "data")} {
						if err := os.WriteFile(filepath.Join(paths[i], file), []byte(fmt.Sprintf("batch %d\n", i)), 0644); err != nil {
							t.Fatal(err)
						}
					}
					info := &IndexInfo{InputGenomes: 1, Genomes: 1, GenomeBatches: 1, InputBases: 100, Chunks: chunks}
					if err := writeIndexInfo(filepath.Join(paths[i], FileInfo), info); err != nil {
						t.Fatal(err)
					}
					for chunk := 0; chunk < chunks; chunk++ {
						w, err := kv.NewWriter(32, chunk*masksPerChunk, masksPerChunk, filepath.Join(paths[i], DirSeeds, chunkFile(chunk)), 1, 2, narrow)
						if err != nil {
							t.Fatal(err)
						}
						for mask := 0; mask < masksPerChunk; mask++ {
							values := []uint64{uint64(i)<<30 | uint64(chunk*10+mask), uint64(i)<<30 | 3}
							data := map[uint64]*[]uint64{0: &values, uint64(i + 1): &[]uint64{uint64(i)}}
							if err := w.WriteDataOfAMask(data); err != nil {
								t.Fatal(err)
							}
						}
						if err := w.Close(); err != nil {
							t.Fatal(err)
						}
					}
				}
				outputs[run] = filepath.Join(root, fmt.Sprintf("output%d", run))
				opt := &IndexBuildingOptions{MaxOpenFiles: budget, MergeThreads: chunks}
				if err := mergeIndexes(nil, 1, 2, opt, chunks, outputs[run], paths, tmp, 1); err != nil {
					t.Fatal(err)
				}
			}
			err := filepath.WalkDir(outputs[0], func(path string, entry os.DirEntry, err error) error {
				if err != nil || entry.IsDir() {
					return err
				}
				relative, err := filepath.Rel(outputs[0], path)
				if err != nil {
					return err
				}
				want, err := os.ReadFile(path)
				if err != nil {
					return err
				}
				got, err := os.ReadFile(filepath.Join(outputs[1], relative))
				if err != nil {
					return err
				}
				if !bytes.Equal(want, got) {
					t.Errorf("file differs across merge rounds: %s", relative)
				}
				return nil
			})
			if err != nil {
				t.Fatal(err)
			}
		})
	}
}

// TestAdaptiveIndexMergeFileBudget checks the third output at minimum and large budgets.
func TestAdaptiveIndexMergeFileBudget(t *testing.T) {
	for _, tc := range [][5]int{{2, 13, 8, 2, 1}, {3, 13, 8, 2, 1}, {10, 32, 8, 10, 1}, {2, 32, 8, 2, 4}, {1, 12, 8, 1, 1}} {
		batch, workers, err := planIndexMergeWithOutputs(tc[0], tc[1], tc[2], 3)
		if err != nil || batch != tc[3] || workers != tc[4] || workers*(batch+3)+indexMergeReservedFiles > tc[1] {
			t.Fatalf("adaptive plan %v: batch=%d workers=%d err=%v", tc, batch, workers, err)
		}
	}
	if _, _, err := planIndexMergeWithOutputs(2, 12, 1, 3); err == nil {
		t.Fatal("accepted a budget too small for two inputs and three outputs")
	}
	outdir := filepath.Join(t.TempDir(), "index")
	opt := &IndexBuildingOptions{GenomeBatchSize: 1, MaxOpenFiles: 12, MergeThreads: 1, SeedIndex2Threshold: kv.MinIndex15Threshold}
	if err := BuildIndex(outdir, []string{"missing1.fa", "missing2.fa"}, opt); err == nil {
		t.Fatal("started adaptive construction before validating the merge budget")
	}
	if _, err := os.Stat(outdir + ExtTmpDir); !os.IsNotExist(err) {
		t.Fatalf("created batches before validating the budget: %v", err)
	}
}
