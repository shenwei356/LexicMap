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
	"testing"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
	"github.com/shenwei356/lexichash"
)

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
			out = runSeedIndexCommand(t, "search", "--quiet", "-j", "1", "-d", db, "-n", "0", "-N", "0", query)
			if strings.Contains(string(out), "\t"+genomeID+"\t") != (tc.want > 0) {
				t.Fatalf("partial genome in seeds: %s", out)
			}
			leftovers, err := filepath.Glob(filepath.Join(db, ".genome-input-*"))
			if err != nil || len(leftovers) != 0 {
				t.Fatalf("staged input remains: %v, %v", leftovers, err)
			}
		})
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
