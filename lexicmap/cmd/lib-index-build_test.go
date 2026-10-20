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

	"github.com/shenwei356/lexichash"
)

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
