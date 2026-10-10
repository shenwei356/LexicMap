package cmd

import (
	"bytes"
	"fmt"
	"math/rand"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
)

// TestSeedIndexCommand runs each CLI invocation in its own process to isolate flags and pools.
func TestSeedIndexCommand(t *testing.T) {
	if os.Getenv("LEXICMAP_SEED_INDEX_TEST_HELPER") != "1" {
		return
	}
	for i, arg := range os.Args {
		if arg == "--" {
			RootCmd.SetArgs(os.Args[i+1:])
			if err := RootCmd.Execute(); err != nil {
				os.Exit(1)
			}
			os.Exit(0)
		}
	}
	t.Fatal("missing command arguments")
}

// runSeedIndexCommand exercises the public index, reindex, and search commands.
func runSeedIndexCommand(t *testing.T, args ...string) []byte {
	t.Helper()
	command := exec.Command(os.Args[0], append([]string{"-test.run=^TestSeedIndexCommand$", "--"}, args...)...)
	command.Env = append(os.Environ(), "LEXICMAP_SEED_INDEX_TEST_HELPER=1")
	var output, stderr bytes.Buffer
	command.Stdout, command.Stderr = &output, &stderr
	if err := command.Run(); err != nil {
		t.Fatalf("command %v failed: %v\n%s", args, err, stderr.String())
	}
	return output.Bytes()
}

// TestIndexSeedIndex2 matches integrated builds against a separate reindex-seeds2 run.
func TestIndexSeedIndex2(t *testing.T) {
	root := t.TempDir()
	// Repeated sequence produces blocks large enough to exercise secondary lookups.
	rng := rand.New(rand.NewSource(1))
	motif := make([]byte, 256)
	for i := range motif {
		motif[i] = "ACGT"[rng.Intn(4)]
	}
	sequence := bytes.Repeat(motif, 128)
	files := make([]string, 3)
	for i := range files {
		files[i] = filepath.Join(root, fmt.Sprintf("ref%d.fa", i))
		// Distinct match scores avoid the unspecified ordering of tied genomes.
		referenceMotif := bytes.Clone(motif)
		for mutation := 0; mutation < i; mutation++ {
			pos := 50 + mutation*100
			referenceMotif[pos] = "ACGT"[(strings.IndexByte("ACGT", referenceMotif[pos])+1)%4]
		}
		if err := os.WriteFile(files[i], append([]byte(">seq\n"), bytes.Repeat(referenceMotif, 128)...), 0600); err != nil {
			t.Fatal(err)
		}
	}
	query := filepath.Join(root, "query.fa")
	if err := os.WriteFile(query, append([]byte(">query\n"), sequence[:600]...), 0600); err != nil {
		t.Fatal(err)
	}
	for _, settings := range [][2]int{{1, 1024}, {3, 1024}, {1, 13}} {
		batchSize, maxOpenFiles := settings[0], settings[1]
		t.Run(fmt.Sprintf("batch-size=%d/max-open-files=%d", batchSize, maxOpenFiles), func(t *testing.T) {
			dir := t.TempDir()
			plain, adaptive := filepath.Join(dir, "plain.lmi"), filepath.Join(dir, "adaptive.lmi")
			buildArgs := []string{"index", "--quiet", "-j", "1", "-c", "2", "-m", "64", "-k", "16", "--partitions", "4", "-b", strconv.Itoa(batchSize), "--max-open-files", strconv.Itoa(maxOpenFiles)}
			runSeedIndexCommand(t, append(append([]string{}, buildArgs...), append([]string{"-O", plain}, files...)...)...)
			search := func(db string) []byte {
				return runSeedIndexCommand(t, "search", "--quiet", "-j", "1", "-d", db, "-p", "12", "-P", "12", "-n", "1", "-N", "1", query)
			}
			wantSearch := search(plain)
			if bytes.Count(wantSearch, []byte("\n")) < 2 {
				t.Fatalf("fixture produced no search hits: %s", wantSearch)
			}
			// One case uses the default 8K threshold; the other overrides it to 4K.
			threshold := "8K"
			// Exercise the numeric shorthand in the constrained multi-batch build.
			seedIndex2Flag := "--seed-index2"
			if maxOpenFiles == 13 {
				seedIndex2Flag = "-2"
			}
			adaptiveArgs := append(append([]string{}, buildArgs...), seedIndex2Flag, "-O", adaptive)
			if batchSize == 1 {
				threshold = "4K"
				adaptiveArgs = append(adaptiveArgs, "--seed-index2-threshold", threshold)
			}
			runSeedIndexCommand(t, append(adaptiveArgs, files...)...)
			if _, err := os.Stat(adaptive + ExtTmpDir); batchSize == 1 && !os.IsNotExist(err) {
				t.Fatalf("temporary merge directory remains: %v", err)
			}
			if got := search(adaptive); !bytes.Equal(got, wantSearch) {
				t.Fatalf("adaptive search differs from the default build:\ngot: %s\nwant: %s", got, wantSearch)
			}
			for chunk := 0; chunk < 2; chunk++ {
				if _, err := os.Stat(filepath.Join(plain, DirSeeds, chunkFile(chunk)) + kv.KVIndex15FileExt); !os.IsNotExist(err) {
					t.Fatalf("default build created a secondary index: %v", err)
				}
			}
			runSeedIndexCommand(t, "utils", "reindex-seeds2", "--quiet", "-j", "2", "-d", plain, "--partitions", "4", "--threshold", threshold)
			var secondaryBytes int
			for chunk := 0; chunk < 2; chunk++ {
				for _, ext := range []string{"", kv.KVIndexFileExt, kv.KVIndex15FileExt} {
					read := func(db string) []byte {
						data, err := os.ReadFile(filepath.Join(db, DirSeeds, chunkFile(chunk)) + ext)
						if err != nil {
							t.Fatal(err)
						}
						// Each invocation uses a fresh pairing token in the two headers.
						if ext == kv.KVIndexFileExt {
							clear(data[14:16])
						} else if ext == kv.KVIndex15FileExt {
							clear(data[46:48])
						}
						return data
					}
					got, want := read(adaptive), read(plain)
					if !bytes.Equal(got, want) {
						t.Fatalf("chunk %d%s differs from separate reindexing", chunk, ext)
					}
					if ext == kv.KVIndex15FileExt {
						secondaryBytes += len(got)
					}
				}
			}
			if secondaryBytes <= 2*48 {
				t.Fatal("fixture did not create any secondary lookup blocks")
			}
			// Switching partitions and rebuilding twice covers removal and absent sidecars.
			for i := 0; i < 2; i++ {
				runSeedIndexCommand(t, "utils", "reindex-seeds", "--quiet", "-j", "2", "-d", adaptive, "--partitions", "16")
				for chunk := 0; chunk < 2; chunk++ {
					if _, err := os.Stat(filepath.Join(adaptive, DirSeeds, chunkFile(chunk)) + kv.KVIndex15FileExt); !os.IsNotExist(err) {
						t.Fatalf("obsolete secondary index was not removed: %v", err)
					}
				}
			}
			info, err := readIndexInfo(filepath.Join(adaptive, FileInfo))
			if err != nil || info.Partitions != 16 {
				t.Fatalf("wrong partition count after reindexing: %+v, %v", info, err)
			}
			if got := search(adaptive); !bytes.Equal(got, wantSearch) {
				t.Fatalf("single-level search differs after switching indexes:\ngot: %s\nwant: %s", got, wantSearch)
			}
		})
	}
}

// TestParseIndex15Threshold checks byte-size syntax and the shared lower bound.
func TestParseIndex15Threshold(t *testing.T) {
	for _, value := range []string{"", "0", "4095", "-8K", " -4K ", "invalid"} {
		if _, err := parseIndex15Threshold(value); err == nil {
			t.Fatalf("accepted invalid threshold %q", value)
		}
	}
	for _, value := range []string{"4096", "4K", " 4k "} {
		if got, err := parseIndex15Threshold(value); err != nil || got != kv.MinIndex15Threshold {
			t.Fatalf("threshold %q: got %d, %v", value, got, err)
		}
	}
}

// TestBuildIndexRejectsInvalidSeedIndex2 checks validation before any seeds are written.
func TestBuildIndexRejectsInvalidSeedIndex2(t *testing.T) {
	for _, partitions := range []int{1, 8, 4096} {
		outdir := filepath.Join(t.TempDir(), "index")
		opt := &IndexBuildingOptions{K: 32, Masks: 64, Chunks: 2, Partitions: partitions,
			RandSeed: 1, GenomeBatchSize: 3, SeedIndex2Threshold: kv.MinIndex15Threshold}
		if err := BuildIndex(outdir, []string{"missing.fa"}, opt); err == nil || !strings.Contains(err.Error(), "prefix") && !strings.Contains(err.Error(), "partitions") && !strings.Contains(err.Error(), "suffix") {
			t.Fatalf("unexpected validation result for %d partitions: %v", partitions, err)
		}
		if _, err := os.Stat(outdir); !os.IsNotExist(err) {
			t.Fatalf("created index before validating adaptive settings: %v", err)
		}
	}
}
