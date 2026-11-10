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
	"compress/gzip"
	"math/rand"
	"reflect"
	"strings"

	"bytes"
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"testing"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/seedposition"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/lexichash"
)

// A genome's completion log must precede the downstream owners that can recycle it.
func TestIndexDebugGenomeHandoff(t *testing.T) {
	dir := t.TempDir()
	db, logFile := filepath.Join(dir, "db"), filepath.Join(dir, "index.log")
	const genomes = 64
	args := []string{"index", "--quiet", "--debug", "--log", logFile,
		"-j", "4", "-m", "64", "-c", "1", "--partitions", "4", "-O", db}
	want := make(map[string]string, genomes)
	for i := range genomes {
		id := fmt.Sprintf("g%03d", i)
		want[id] = fmt.Sprint(i)
		args = append(args, writeTestInput(t, filepath.Join(dir, id+".fa"), ">s\n"+indexTestDNA(96+i)+"\n"))
	}
	runSeedIndexCommand(t, args...)
	data, err := os.ReadFile(logFile)
	if err != nil {
		t.Fatal(err)
	}
	pattern := regexp.MustCompile(`file #(\d+), (finished computing seeds|genome data saved|seeds collected): (\S+)`)
	finished := make(map[string]int, genomes)
	for _, match := range pattern.FindAllStringSubmatch(string(data), -1) {
		index, stage, id := match[1], match[2], match[3]
		if want[id] != index {
			t.Fatalf("debug log read a reused genome: index=%s id=%s stage=%s", index, id, stage)
		}
		if stage == "finished computing seeds" {
			finished[id]++
		} else if finished[id] != 1 {
			t.Fatalf("genome %s reached %q before its completion log", id, stage)
		}
	}
	for id := range want {
		if finished[id] != 1 {
			t.Errorf("genome %s has %d completion logs, want 1", id, finished[id])
		}
	}
	info, err := readIndexInfo(filepath.Join(db, FileInfo))
	if err != nil {
		t.Fatal(err)
	}
	if info.Genomes != genomes || info.InputGenomes != genomes {
		t.Fatalf("debug indexing lost genomes: stored=%d inputs=%d", info.Genomes, info.InputGenomes)
	}
}

// Several seed-free regions between adjacent seeds need only one boundary flag.
func TestIndexSeedPositionsAcrossUnseededRegions(t *testing.T) {
	for _, contigs := range []bool{false, true} {
		for _, gaps := range []int{2, 3} {
			for _, reverse := range []bool{false, true} {
				for _, noFill := range []bool{false, true} {
					t.Run(fmt.Sprintf("contigs=%t/gaps=%d/reverse=%t/noFill=%t", contigs, gaps, reverse, noFill), func(t *testing.T) {
						dir := t.TempDir()
						sequences := []string{indexTestDNA(1000)}
						for range gaps - 1 {
							sequences = append(sequences, strings.Repeat("A", 31))
						}
						sequences = append(sequences, indexTestDNA(1000))
						interval := 1000
						if !contigs {
							interval = 10
						}
						if reverse {
							for i, sequence := range sequences {
								sequences[i] = string((&seq.Seq{Alphabet: seq.DNA, Seq: []byte(sequence)}).RevCom().Seq)
							}
						}
						var input strings.Builder
						var ends []int
						position := 0
						for i, sequence := range sequences {
							if i > 0 {
								ends = append(ends, position+interval-1)
								position += interval
							}
							position += len(sequence)
							if contigs || i == 0 {
								fmt.Fprintf(&input, ">s%d\n", i)
							} else {
								input.WriteString(strings.Repeat("N", interval))
							}
							input.WriteString(sequence)
							input.WriteByte('\n')
						}
						file := writeTestInput(t, filepath.Join(dir, "g.fa"), input.String())
						db := filepath.Join(dir, "db")
						args := []string{"index", "--quiet", "-j", "1", "-m", "64", "-c", "1", "--partitions", "4", "--save-seed-pos", "-O", db, file}
						if noFill {
							args = append(args, "--no-desert-filling")
						}
						runSeedIndexCommand(t, args...)
						reader, err := seedposition.NewReader(filepath.Join(db, DirGenomes, batchDir(0), FileSeedPositions))
						if err != nil {
							t.Fatal(err)
						}
						var positions []uint32
						err = reader.SeedPositions(0, &positions)
						closeErr := reader.Close()
						if err != nil || closeErr != nil {
							t.Fatalf("reading positions: %v, %v", err, closeErr)
						}
						previous, flags := -1, 0
						for _, packed := range positions {
							position := int(packed >> 2)
							wantFlag := false
							for _, end := range ends {
								wantFlag = wantFlag || previous < end && end <= position
							}
							if (packed&1 != 0) != wantFlag {
								t.Errorf("seed at %d after %d: boundary flag=%d, want %t", position, previous, packed&1, wantFlag)
							}
							if wantFlag {
								flags++
							}
							previous = position
						}
						if flags != 1 || len(positions) < 4 {
							t.Fatalf("fixture did not cross all seed-free regions together: flags=%d, seeds=%d", flags, len(positions))
						}
						if contigs {
							output := runSeedIndexCommand(t, "utils", "seed-pos", "--quiet", "-j", "1", "-d", db, "-n", "g")
							rows := strings.Split(strings.TrimSuffix(string(output), "\n"), "\n")[1:]
							if len(rows) != len(positions) {
								t.Fatalf("seed-pos retained %d of %d positions", len(rows), len(positions))
							}
							for i, row := range rows {
								fields := strings.Split(row, "\t")
								packed := positions[i]
								if len(fields) != 6 || fields[2] != fmt.Sprint((packed>>2)+1) || fields[4] != string(lexichash.Strands[packed>>1&1]) {
									t.Fatalf("seed-pos changed position or strand at row %d: %q", i+1, row)
								}
							}
						}
					})
				}
			}
		}
	}
}

func TestIndexSeedChunkCountPreservesLayout(t *testing.T) {
	// Check every supported chunk count across divisible and uneven mask counts.
	for masks := 64; masks <= 1024; masks++ {
		for requested := 1; requested <= min(masks, 128); requested++ {
			chunks := indexSeedChunkCount(masks, requested)
			oldWidth := (masks + requested - 1) / requested
			width := (masks + chunks - 1) / chunks
			if chunks > requested || width != oldWidth || (chunks-1)*width >= masks || chunks*width < masks {
				t.Fatalf("masks=%d requested=%d: chunks=%d width=%d, legacy width=%d", masks, requested, chunks, width, oldWidth)
			}
			if (requested-1)*oldWidth < masks && chunks != requested {
				t.Fatalf("changed a valid legacy layout: masks=%d requested=%d", masks, requested)
			}
		}
	}
}

func TestIndexSeedChunks(t *testing.T) {
	for _, batch := range []string{"1", "2"} {
		t.Run("batch="+batch, func(t *testing.T) {
			dir := t.TempDir()
			input := writeTestInput(t, filepath.Join(dir, "g1.fa"), ">s\n"+indexTestDNA(1800)+"\n")
			other := writeTestInput(t, filepath.Join(dir, "g2.fa"), ">s\n"+indexTestDNA(1600)+"\n")
			var wantKmers, wantSearch []byte
			for _, tc := range []struct{ requested, actual int }{{32, 32}, {63, 32}, {3, 3}} {
				db := filepath.Join(dir, fmt.Sprintf("db%d", tc.requested))
				runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64",
					"-c", fmt.Sprint(tc.requested), "--partitions", "4", "-b", batch, "-O", db, input, other)
				info, err := readIndexInfo(filepath.Join(db, FileInfo))
				if err != nil || info.Chunks != tc.actual {
					t.Fatalf("incorrect chunk summary: %+v, %v", info, err)
				}
				width := (64 + tc.actual - 1) / tc.actual
				for c := 0; c < tc.actual; c++ {
					_, offset, masks, _, _, _, err := kv.ReadKVIndex(filepath.Join(db, DirSeeds, chunkFile(c)) + kv.KVIndexFileExt)
					if err != nil || offset != c*width || len(masks) != min(width, 64-c*width) {
						t.Fatalf("chunk %d: offset=%d masks=%d, %v", c, offset, len(masks), err)
					}
				}
				if _, err := os.Stat(filepath.Join(db, DirSeeds, chunkFile(tc.actual))); !os.IsNotExist(err) {
					t.Fatalf("unexpected extra chunk: %v", err)
				}
				// Both all-mask and selected-mask readers must agree across layouts.
				kmers := runSeedIndexCommand(t, "utils", "kmers", "--quiet", "-d", db, "-m", "0")
				last := runSeedIndexCommand(t, "utils", "kmers", "--quiet", "-d", db, "-m", "64")
				lines := bytes.SplitAfter(kmers, []byte("\n"))
				selected := bytes.Clone(lines[0])
				for _, line := range lines[1:] {
					if bytes.HasPrefix(line, []byte("64\t")) {
						selected = append(selected, line...)
					}
				}
				if !bytes.Equal(last, selected) || bytes.Count(last, []byte("\n")) < 2 {
					t.Fatal("selected mask differs from all-mask output")
				}
				if wantKmers == nil {
					wantKmers = kmers
				} else if !bytes.Equal(kmers, wantKmers) {
					t.Fatal("seed output changed across chunk layouts")
				}
				for _, load := range []bool{false, true} {
					args := []string{"search", "--quiet", "-j", "1", "-d", db, "-n", "0", "-N", "0", input}
					if load {
						args = append(args, "--load-whole-seeds")
					}
					out := runSeedIndexCommand(t, args...)
					if wantSearch == nil {
						wantSearch = out
						if bytes.Count(out, []byte("\n")) < 3 {
							t.Fatalf("missing search hits: %s", out)
						}
					} else if !bytes.Equal(out, wantSearch) {
						t.Fatalf("search output changed across chunk layouts or loading modes:\n%s\n%s", wantSearch, out)
					}
				}
			}
		})
	}
}

func TestIndexContigIntervalLengthLimit(t *testing.T) {
	for _, tc := range []struct {
		limit int
		want  []int
	}{{1800, []int{1800}}, {1799, []int{400, 400}}, {1000, []int{400, 400}}} {
		t.Run(fmt.Sprint(tc.limit), func(t *testing.T) {
			dir := t.TempDir()
			input := writeTestInput(t, filepath.Join(dir, "g.fa"), ">a\n"+indexTestDNA(400)+"\n>b\n"+indexTestDNA(400)+"\n")
			db := filepath.Join(dir, "db")
			runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-c", "2",
				"--partitions", "4", "-g", fmt.Sprint(tc.limit), "-O", db, input)
			info, err := readIndexInfo(filepath.Join(db, FileInfo))
			if err != nil || info.Genomes != len(tc.want) || info.InputGenomes != 1 {
				t.Fatalf("wrong chunk counts: %+v, %v", info, err)
			}
			reader, err := genome.NewReader(filepath.Join(db, DirGenomes, batchDir(0), FileGenomes))
			if err != nil {
				t.Fatal(err)
			}
			defer reader.Close()
			for i, length := range tc.want {
				g, err := reader.GenomeInfo(i)
				if err != nil {
					t.Fatal(err)
				}
				if g.Len != length || g.Len > tc.limit {
					t.Fatalf("chunk %d length=%d, want %d (limit %d)", i, g.Len, length, tc.limit)
				}
				genome.RecycleGenome(g)
			}
		})
	}
}

func TestIndexCustomMaskMinSequenceLength(t *testing.T) {
	dir := t.TempDir()
	lh, err := lexichash.NewWithSeed(16, 64, 1, 0)
	if err != nil {
		t.Fatal(err)
	}
	var masks strings.Builder
	for _, mask := range lh.Masks {
		fmt.Fprintf(&masks, "%s\n", lexichash.MustDecode(mask, 16))
	}
	maskFile := writeTestInput(t, filepath.Join(dir, "masks.txt"), masks.String())
	input := writeTestInput(t, filepath.Join(dir, "g.fa"), ">short\n"+indexTestDNA(25)+"\n")
	for _, tc := range []struct {
		name string
		args []string
		want int
	}{
		{"default", nil, 1},
		{"explicit_k", []string{"-k", "16"}, 1},
		{"explicit_min_below_cli_k", []string{"-l", "20"}, 1},
		{"explicit_min", []string{"-l", "26"}, 0},
	} {
		t.Run(tc.name, func(t *testing.T) {
			db := filepath.Join(dir, tc.name)
			args := []string{"index", "--quiet", "-j", "1", "-c", "2", "--partitions", "4",
				"-M", maskFile, "-O", db, input}
			runSeedIndexCommand(t, append(args, tc.args...)...)
			info, err := readIndexInfo(filepath.Join(db, FileInfo))
			if err != nil || info.Genomes != tc.want {
				t.Fatalf("genomes: %+v, %v; want %d", info, err, tc.want)
			}
		})
	}
}

func TestLoadIndexMasksUsesCustomMaskCount(t *testing.T) {
	lh, err := lexichash.NewWithSeed(31, 64, 1, 0)
	if err != nil {
		t.Fatal(err)
	}

	var data bytes.Buffer
	for _, mask := range lh.Masks {
		if _, err := fmt.Fprintf(&data, "%s\n", lexichash.MustDecode(mask, uint8(lh.K))); err != nil {
			t.Fatal(err)
		}
	}

	file := filepath.Join(t.TempDir(), "masks.txt")
	if err := os.WriteFile(file, data.Bytes(), 0o644); err != nil {
		t.Fatal(err)
	}

	opt := &IndexBuildingOptions{MaskFile: file, K: 21, Masks: 20_000, Chunks: 16}
	got, err := loadIndexMasks(opt)
	if err != nil {
		t.Fatal(err)
	}
	if got.K != 31 || opt.K != 31 {
		t.Fatalf("k-mer size: loaded=%d option=%d, want 31", got.K, opt.K)
	}
	if len(got.Masks) != 64 || opt.Masks != 64 {
		t.Fatalf("mask count: loaded=%d option=%d, want 64", len(got.Masks), opt.Masks)
	}
}

// indexTestDNA supplies reproducible, non-repetitive seeds for CLI regressions.
func indexTestDNA(n int) string {
	rng := rand.New(rand.NewSource(81))
	bases := make([]byte, n)
	for i := range bases {
		bases[i] = "ACGT"[rng.Intn(4)]
	}
	return string(bases)
}

func TestIndexRetainsEmptyBatches(t *testing.T) {
	for _, valid := range []int{0, 1, 2, -1} {
		t.Run(fmt.Sprintf("valid=%d", valid), func(t *testing.T) {
			dir := t.TempDir()
			query := writeTestInput(t, filepath.Join(dir, "query.fa"), ">q\n"+indexTestDNA(1500)+"\n")
			db := filepath.Join(dir, "db")
			args := []string{"index", "--quiet", "-j", "1", "-m", "64", "-c", "2", "-k", "16",
				"--partitions", "4", "-b", "1", "-O", db}
			for i := 0; i < 3; i++ {
				sequence := "ACGT"
				if i == valid {
					sequence = indexTestDNA(1500)
				}
				args = append(args, writeTestInput(t, filepath.Join(dir, fmt.Sprintf("g%d.fa", i)), ">s\n"+sequence+"\n"))
			}
			runSeedIndexCommand(t, args...)
			info, err := readIndexInfo(filepath.Join(db, FileInfo))
			if err != nil {
				t.Fatal(err)
			}
			want := 1
			if valid < 0 {
				want = 0
			}
			if info.GenomeBatches != 3 || info.InputGenomes != want || info.Genomes != want {
				t.Fatalf("wrong batch/genome counts: %+v", info)
			}
			for i := 0; i < 3; i++ {
				if _, err := os.Stat(filepath.Join(db, DirGenomes, batchDir(i), FileGenomes)); err != nil {
					t.Fatal(err)
				}
			}
			for _, load := range []bool{false, true} {
				args := []string{"search", "--quiet", "-j", "1", "-d", db, "-n", "0", "-N", "0", query}
				if load {
					args = append(args, "--load-whole-seeds")
				}
				if want == 0 {
					assertEmptyIndexSearchRejected(t, args...)
					continue
				}
				out := runSeedIndexCommand(t, args...)
				if (strings.Count(string(out), "\n") > 1) != (want > 0) {
					t.Fatalf("wrong search result: %s", out)
				}
			}
			for _, args := range [][]string{
				{"utils", "genome-details", "-i", "-e", "-d", db},
				{"genome", "search", "-m", "0", "-p", "12", "-d", db, query},
				{"genome", "pair", "-m", "0", "-p", "12", "-d", db},
			} {
				if want == 0 && args[0] == "genome" && args[1] == "search" {
					assertEmptyIndexSearchRejected(t, append(args, "--quiet", "-j", "1")...)
					continue
				}
				runSeedIndexCommand(t, append(args, "--quiet", "-j", "1")...)
			}
		})
	}
}

func TestIndexRejectsWholeOversizedGenome(t *testing.T) {
	for _, tc := range []struct {
		name    string
		lengths []int
		exclude bool
		want    int
		stdin   bool
	}{
		{"first", []int{1500, 700, 700}, false, 0, false},
		{"middle", []int{700, 1500, 700}, false, 0, false},
		{"last", []int{700, 700, 700, 1500}, false, 0, false},
		{"valid_split", []int{700, 700, 700}, false, 3, false},
		{"excluded_long_contig", []int{700, 700, 1500}, true, 2, false},
		{"stdin_valid_split", []int{700, 700, 700}, false, 3, true},
		{"stdin_gzip", []int{700, 700, 700}, false, 3, true},
		{"stdin_oversized", []int{700, 700, 700, 1500}, false, 0, true},
	} {
		t.Run(tc.name, func(t *testing.T) {
			dir := t.TempDir()
			var fasta strings.Builder
			for i, n := range tc.lengths {
				fmt.Fprintf(&fasta, ">c%d\n%s\n", i, indexTestDNA(n))
			}
			input := writeTestInput(t, filepath.Join(dir, "input.fa"), fasta.String())
			data := fasta.String()
			if tc.name == "valid_split" || tc.name == "stdin_gzip" {
				var compressed bytes.Buffer
				writer := gzip.NewWriter(&compressed)
				if _, err := writer.Write([]byte(fasta.String())); err != nil {
					t.Fatal(err)
				}
				if err := writer.Close(); err != nil {
					t.Fatal(err)
				}
				data = compressed.String()
				input = writeTestInput(t, input+".gz", data)
			}
			genomeID := "input"
			if tc.stdin {
				input, genomeID = "-", "-"
			}
			db, skipped := filepath.Join(dir, "db"), filepath.Join(dir, "skipped.tsv")
			args := []string{"index", "--quiet", "-j", "1", "-m", "64", "-c", "2", "-k", "16",
				"--partitions", "4", "-g", "1000", "--big-genomes", skipped, "-O", db, input}
			if tc.exclude {
				args = append(args, "--seq-name-filter", "^c2$")
			}
			_, stderr, err := runAuditCommand(t, data, args...)
			if err != nil {
				t.Fatalf("index: %v\n%s", err, stderr)
			}
			info, err := readIndexInfo(filepath.Join(db, FileInfo))
			if err != nil {
				t.Fatal(err)
			}
			if info.Genomes != tc.want {
				t.Fatalf("genomes=%d, want %d", info.Genomes, tc.want)
			}
			log, err := os.ReadFile(skipped)
			if err != nil {
				t.Fatal(err)
			}
			if strings.Contains(string(log), TOO_LARGE_GENOME) != (tc.want == 0) {
				t.Fatalf("inconsistent skip log: %s", log)
			}
			out := runSeedIndexCommand(t, "utils", "genomes", "--quiet", "-d", db)
			if strings.Contains(string(out), genomeID) != (tc.want > 0) {
				t.Fatalf("partial genome in mapping: %s", out)
			}
			query := writeTestInput(t, filepath.Join(dir, "query.fa"), ">q\n"+indexTestDNA(700)+"\n")
			searchArgs := []string{"search", "--quiet", "-j", "1", "-d", db, "-n", "0", "-N", "0", query}
			if tc.want == 0 {
				assertEmptyIndexSearchRejected(t, searchArgs...)
			} else {
				out = runSeedIndexCommand(t, searchArgs...)
				if !strings.Contains(string(out), "\t"+genomeID+"\t") {
					t.Fatalf("missing accepted genome in seeds: %s", out)
				}
			}
			leftovers, err := filepath.Glob(filepath.Join(db, ".genome-input-*"))
			if err != nil || len(leftovers) != 0 {
				t.Fatalf("staged input remains: %v, %v", leftovers, err)
			}
		})
	}
}

// Indexing may produce an empty index after filtering all input genomes;
// search must report that condition instead of accepting it as a usable index.
func assertEmptyIndexSearchRejected(t *testing.T, args ...string) {
	t.Helper()
	_, stderr, err := runAuditCommand(t, "", args...)
	if err == nil || !bytes.Contains(stderr, []byte("no sequences indexed in the index")) {
		t.Fatalf("empty index search: error=%v, stderr=%s", err, stderr)
	}
}

func TestWriteIndexInfoAtomic(t *testing.T) {
	dir := t.TempDir()
	file := filepath.Join(dir, FileInfo)
	want := &IndexInfo{K: 31, Genomes: 7, GenomeBatches: 3}
	if err := writeIndexInfo(file, want); err != nil {
		t.Fatal(err)
	}
	if err := os.Chmod(file, 0640); err != nil {
		t.Fatal(err)
	}
	want.Genomes++
	if err := writeIndexInfo(file, want); err != nil {
		t.Fatal(err)
	}
	got, err := readIndexInfo(file)
	if err != nil || !reflect.DeepEqual(got, want) {
		t.Fatalf("round trip: %+v, %v", got, err)
	}
	stat, err := os.Stat(file)
	if err != nil || stat.Mode().Perm() != 0640 {
		t.Fatalf("permissions changed: %v, %v", stat, err)
	}
	// A failed replacement must preserve its target and remove temporary files.
	if err := writeIndexInfo(dir, want); err == nil {
		t.Fatal("replaced a directory")
	}
	entries, err := os.ReadDir(dir)
	if err != nil || len(entries) != 1 || entries[0].Name() != FileInfo {
		t.Fatalf("temporary files remain: %v, %v", entries, err)
	}
}

func TestBuildIndexRejectsCurrentDirectory(t *testing.T) {
	cwd, err := os.Getwd()
	if err != nil {
		t.Fatal(err)
	}
	root := t.TempDir()
	alias := filepath.Join(root, "cwd")
	if err := os.Symlink(cwd, alias); err != nil {
		t.Fatal(err)
	}
	for _, path := range []string{".", cwd, alias} {
		if err := BuildIndex(path, nil, &IndexBuildingOptions{}); err == nil || !strings.Contains(err.Error(), "should not be current directory") {
			t.Fatalf("current output directory %q accepted: %v", path, err)
		}
	}
	db := filepath.Join(root, "db")
	if err := os.Symlink(cwd, db+ExtTmpDir); err != nil {
		t.Fatal(err)
	}
	opt := &IndexBuildingOptions{GenomeBatchSize: 1}
	if err := BuildIndex(db, []string{"a.fa", "b.fa"}, opt); err == nil || !strings.Contains(err.Error(), "should not be current directory") {
		t.Fatalf("current temporary directory accepted: %v", err)
	}
}

// Invalid construction budgets must fail before replacing the output directory.
func TestIndexRejectsInsufficientFileBudget(t *testing.T) {
	for _, positions := range []bool{false, true} {
		dir := t.TempDir()
		ref := writeTestInput(t, filepath.Join(dir, "ref.fa"), ">s\n"+indexTestDNA(1000)+"\n")
		db := filepath.Join(dir, "db")
		if err := os.MkdirAll(db, 0700); err != nil {
			t.Fatal(err)
		}
		sentinel := writeTestInput(t, filepath.Join(db, "keep"), "original\n")
		minimum := reservedOpenFiles + 4
		args := []string{"index", "--quiet", "--force", "-O", db, ref}
		if positions {
			minimum++
			args = append(args, "--save-seed-pos")
		}
		args = append(args, "--max-open-files", fmt.Sprint(minimum-1))
		_, stderr, err := runAuditCommand(t, "", args...)
		if err == nil || !strings.Contains(string(stderr), fmt.Sprintf("should be >= %d", minimum)) {
			t.Fatalf("invalid budget accepted: %v, %s", err, stderr)
		}
		got, err := os.ReadFile(sentinel)
		if err != nil || string(got) != "original\n" {
			t.Fatalf("invalid budget changed output directory: %q, %v", got, err)
		}
	}
}

func TestBuildIndexRejectsReportInsideOutput(t *testing.T) {
	for _, suffix := range []string{"", ExtTmpDir} {
		t.Run(suffix, func(t *testing.T) {
			db := filepath.Join(t.TempDir(), "db")
			if err := os.Mkdir(db+suffix, 0700); err != nil {
				t.Fatal(err)
			}
			report := writeTestInput(t, filepath.Join(db+suffix, "skipped.tsv"), "original\n")
			opt := &IndexBuildingOptions{GenomeBatchSize: 1, BigGenomeFile: report}
			if err := BuildIndex(db, []string{"a.fa", "b.fa"}, opt); err == nil || !strings.Contains(err.Error(), "big-genomes path") {
				t.Fatalf("unsafe report accepted: %v", err)
			}
			if got, err := os.ReadFile(report); err != nil || string(got) != "original\n" {
				t.Fatalf("report changed: %q, %v", got, err)
			}
		})
	}
}

// TestIndexBigGenomesReportSurvivesMerge covers an empty batch and final cleanup.
func TestIndexBigGenomesReportSurvivesMerge(t *testing.T) {
	root := t.TempDir()
	big := writeTestInput(t, filepath.Join(root, "big.fa"), ">big\n"+indexTestDNA(1020)+"\n")
	small := writeTestInput(t, filepath.Join(root, "small.fa"), ">small\n"+indexTestDNA(100)+"\n")
	db := filepath.Join(root, "db")
	report := filepath.Join(root, "skipped.tsv")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-c", "2", "--partitions", "4",
		"-b", "1", "-g", "1000", "-G", report, "-O", db, big, small)
	if got, err := os.ReadFile(report); err != nil || string(got) != big+"\ttoo_large_genome\n" {
		t.Fatalf("missing skipped genome: %q, %v", got, err)
	}
	info, err := readIndexInfo(filepath.Join(db, FileInfo))
	if err != nil || info.Genomes != 1 || info.GenomeBatches != 2 {
		t.Fatalf("incorrect merged index: %+v, %v", info, err)
	}
	if _, err := os.Stat(db + ExtTmpDir); !os.IsNotExist(err) {
		t.Fatalf("temporary batches remain: %v", err)
	}
}
