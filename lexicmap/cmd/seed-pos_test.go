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
	"fmt"
	"os"
	"path/filepath"
	"slices"
	"strings"
	"testing"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/seedposition"
)

// writeSeedPosFixture stores real genome/position records without a seed index.
func writeSeedPosFixture(t *testing.T, batches int) string {
	t.Helper()
	dir := t.TempDir()
	sequence := []byte(indexTestDNA(128))
	var mapping bytes.Buffer
	for batch := 0; batch < batches; batch++ {
		name := fmt.Sprintf("g%d", batch)
		path := filepath.Join(dir, DirGenomes, batchDir(batch))
		if err := os.MkdirAll(path, 0700); err != nil {
			t.Fatal(err)
		}
		writer, err := genome.NewWriter(filepath.Join(path, FileGenomes), uint32(batch))
		if err != nil {
			t.Fatal(err)
		}
		id := []byte("s")
		g := &genome.Genome{ID: []byte(name), Seq: sequence, GenomeSize: len(sequence), Len: len(sequence), NumSeqs: 1, SeqSizes: []int{len(sequence)}, SeqIDs: []*[]byte{&id}}
		if err := writer.Write(g); err != nil {
			t.Fatal(err)
		}
		if err := writer.Close(); err != nil {
			t.Fatal(err)
		}
		positions, err := seedposition.NewWriter(filepath.Join(path, FileSeedPositions), uint32(batch))
		if err != nil {
			t.Fatal(err)
		}
		if err := positions.Write([]uint32{10 << 2, 30 << 2, 60 << 2}); err != nil {
			t.Fatal(err)
		}
		if err := positions.Close(); err != nil {
			t.Fatal(err)
		}
		binary.Write(&mapping, binary.BigEndian, uint16(len(name)))
		mapping.WriteString(name)
		binary.Write(&mapping, binary.BigEndian, uint64(batch)<<BITS_GENOME_IDX)
	}
	if err := os.WriteFile(filepath.Join(dir, FileGenomeIndex), mapping.Bytes(), 0600); err != nil {
		t.Fatal(err)
	}
	info := &IndexInfo{MainVersion: MainVersion, K: 31, Chunks: 2, GenomeBatches: batches, Genomes: batches, InputGenomes: batches}
	if err := writeIndexInfo(filepath.Join(dir, FileInfo), info); err != nil {
		t.Fatal(err)
	}
	return dir
}

func TestSeedPosReadsRequestedBatches(t *testing.T) {
	db := writeSeedPosFixture(t, 3)
	if err := os.Remove(filepath.Join(db, DirGenomes, batchDir(2), FileGenomes)); err != nil {
		t.Fatal(err)
	}
	output := runSeedIndexCommand(t, "utils", "seed-pos", "--quiet", "-d", db, "-n", "g0", "--max-open-files", "11")
	if strings.Count(string(output), "\n") != 4 {
		t.Fatalf("lost seed positions: %s", output)
	}
}

func TestSeedPosInvalidOffset(t *testing.T) {
	db := writeSeedPosFixture(t, 1)
	file := filepath.Join(db, DirGenomes, batchDir(0), FileSeedPositions) + seedposition.PositionsIndexFileExt
	data, err := os.ReadFile(file)
	if err != nil {
		t.Fatal(err)
	}
	binary.BigEndian.PutUint64(data[24:32], ^uint64(0)-7)
	if err := os.WriteFile(file, data, 0600); err != nil {
		t.Fatal(err)
	}
	_, stderr, err := runAuditCommand(t, "", "utils", "seed-pos", "--quiet", "-j", "1", "-d", db, "-n", "g0")
	if err == nil || bytes.Contains(stderr, []byte("panic:")) || !bytes.Contains(stderr, []byte(seedposition.ErrBrokenFile.Error())) {
		t.Fatalf("expected corrupt-index error without a panic: %v\n%s", err, stderr)
	}
}

// TestSeedPosZeroDistance keeps coincident seeds while emitting empty intervals.
func TestSeedPosZeroDistance(t *testing.T) {
	db := writeSeedPosFixture(t, 1)
	writer, err := seedposition.NewWriter(filepath.Join(db, DirGenomes, batchDir(0), FileSeedPositions), 0)
	if err != nil {
		t.Fatal(err)
	}
	// Both strands at the first base and at an internal base have zero distance.
	if err := writer.Write([]uint32{0, 2, 10 << 2, 10<<2 | 2, 30 << 2}); err != nil {
		t.Fatal(err)
	}
	if err := writer.Close(); err != nil {
		t.Fatal(err)
	}
	sequence := indexTestDNA(128)
	wantSeqs := []string{"", "", sequence[:10], "", sequence[10:30]}
	args := []string{"utils", "seed-pos", "--quiet", "-j", "1", "-d", db, "-n", "g0", "--max-open-files", "11"}
	plain := strings.Split(strings.TrimSuffix(string(runSeedIndexCommand(t, args...)), "\n"), "\n")
	verbose := strings.Split(strings.TrimSuffix(string(runSeedIndexCommand(t, append(args, "-v")...)), "\n"), "\n")
	if len(plain) != len(wantSeqs)+1 || len(verbose) != len(plain) {
		t.Fatalf("lost seed rows: plain=%q verbose=%q", plain, verbose)
	}
	for i, want := range wantSeqs {
		fields := strings.Split(verbose[i+1], "\t")
		if len(fields) != 8 || strings.Join(fields[:6], "\t") != plain[i+1] {
			t.Fatalf("row %d changed seed fields: %q", i, verbose[i+1])
		}
		if fields[5] != fmt.Sprint(len(want)) || fields[7] != want {
			t.Fatalf("row %d: distance=%s sequence=%q, want %d bases %q", i, fields[5], fields[7], len(want), want)
		}
		if want == "" && fields[6] != "0" {
			t.Fatalf("row %d: empty interval has len_aaa=%s", i, fields[6])
		}
	}
}

func TestSeedPosPlotsWithBoundedReaders(t *testing.T) {
	db := writeSeedPosFixture(t, 2)
	plots := filepath.Join(t.TempDir(), "plots")
	runSeedIndexCommand(t, "utils", "seed-pos", "--quiet", "-v", "-j", "8", "-d", db, "-a", "--max-open-files", "11", "--plot-dir", plots, "--plot-ext", ".svg", "--slid-window", "40", "--slid-step", "20")
	files, err := filepath.Glob(filepath.Join(plots, "*.svg"))
	if err != nil {
		t.Fatal(err)
	}
	if len(files) != 4 {
		t.Fatalf("expected two plots per genome, got %v", files)
	}
}

// Every complete window contributes once, including the last and empty windows.
func TestSeedPosWindowCounts(t *testing.T) {
	for _, tc := range []struct {
		name   string
		window string
		step   string
		locs   []uint32
		want   string
	}{
		{"whole-sequence", "128", "100", []uint32{10 << 2, 30 << 2, 60 << 2}, "min=3, median=3, max=3"},
		{"complete-kmers", "31", "31", []uint32{1 << 2, 31 << 2, 62 << 2}, "min=0, median=0, max=1"},
		{"short-window", "10", "10", []uint32{10 << 2, 30 << 2, 60 << 2}, "min=0, median=0, max=0"},
	} {
		t.Run(tc.name, func(t *testing.T) {
			db := writeSeedPosFixture(t, 1)
			writer, err := seedposition.NewWriter(filepath.Join(db, DirGenomes, batchDir(0), FileSeedPositions), 0)
			if err != nil {
				t.Fatal(err)
			}
			if err := writer.Write(tc.locs); err != nil {
				t.Fatal(err)
			}
			if err := writer.Close(); err != nil {
				t.Fatal(err)
			}
			plots := filepath.Join(t.TempDir(), "plots")
			runSeedIndexCommand(t, "utils", "seed-pos", "--quiet", "-d", db, "-n", "g0", "--plot-dir", plots,
				"--plot-ext", ".svg", "--slid-window", tc.window, "--slid-step", tc.step)
			data, err := os.ReadFile(filepath.Join(plots, "g0.seed_number.svg"))
			if err != nil {
				t.Fatal(err)
			}
			if !bytes.Contains(data, []byte(tc.want)) {
				t.Fatalf("seed-count plot does not contain %q", tc.want)
			}
		})
	}
}

// Exercise overlapping windows, gaps, coincident strands, and contig boundaries.
func TestSeedWindowCounts(t *testing.T) {
	locs := []uint32{0, 2, 10 << 2, 40<<2 | 1, 60 << 2, 80 << 2, 130<<2 | 1, 150 << 2, 150<<2 | 2, 169 << 2}
	regions := [][2]int{{0, 110}, {120, 199}}
	for window := uint32(1); window <= 112; window++ {
		for _, step := range []uint32{1, 7, 40, 128} {
			var want []float64
			// Independent exhaustive definition: inspect every seed for each window.
			for _, region := range regions {
				for start := region[0]; start+int(window)-1 <= region[1]; start += int(step) {
					n, gap := 0, false
					for _, loc := range locs {
						pos := int(loc >> 2)
						if pos >= start && pos+31 <= start+int(window) {
							n++
							gap = gap || loc&1 != 0
						}
					}
					if !gap {
						want = append(want, float64(n))
					}
				}
			}
			got := seedWindowCounts(locs, regions, 31, window, step, nil)
			if !slices.Equal(got, want) {
				t.Fatalf("window=%d step=%d: got %v, want %v", window, step, got, want)
			}
		}
	}
	// Reversing the sequence reverses an exactly tiled window grid's counts.
	locs = []uint32{0, 1 << 2, 31 << 2, 32 << 2, 63 << 2, 64 << 2, 95 << 2, 97 << 2}
	reverse := make([]uint32, len(locs))
	for i, loc := range locs {
		reverse[len(locs)-1-i] = (128-31-loc>>2)<<2 | 2
	}
	forwardCounts := seedWindowCounts(locs, [][2]int{{0, 127}}, 31, 32, 32, nil)
	reverseCounts := seedWindowCounts(reverse, [][2]int{{0, 127}}, 31, 32, 32, nil)
	slices.Reverse(reverseCounts)
	if !slices.Equal(forwardCounts, reverseCounts) {
		t.Fatalf("reverse complement changed window counts: %v != %v", forwardCounts, reverseCounts)
	}
}

// Chunks shorter than the counting window still need distinct distance plots.
func TestSeedPosChunkPlotsWithoutWindows(t *testing.T) {
	db := writeSeedPosFixture(t, 2)
	var mapping bytes.Buffer
	for batch := 0; batch < 2; batch++ {
		binary.Write(&mapping, binary.BigEndian, uint16(2))
		mapping.WriteString("g0")
		binary.Write(&mapping, binary.BigEndian, uint64(batch)<<BITS_GENOME_IDX)
	}
	if err := os.WriteFile(filepath.Join(db, FileGenomeIndex), mapping.Bytes(), 0600); err != nil {
		t.Fatal(err)
	}
	plots := filepath.Join(t.TempDir(), "plots")
	runSeedIndexCommand(t, "utils", "seed-pos", "--quiet", "-d", db, "-n", "g0", "--plot-dir", plots,
		"--plot-ext", ".svg", "--slid-window", "200")
	for _, name := range []string{"g0.svg", "g0_2.svg"} {
		if _, err := os.Stat(filepath.Join(plots, name)); err != nil {
			t.Fatal(err)
		}
	}
}

// Budget validation must precede replacement of a requested plot directory.
func TestSeedPosRejectsInsufficientFileBudget(t *testing.T) {
	db := writeSeedPosFixture(t, 1)
	plots := filepath.Join(t.TempDir(), "plots")
	if err := os.MkdirAll(plots, 0700); err != nil {
		t.Fatal(err)
	}
	sentinel := writeTestInput(t, filepath.Join(plots, "keep"), "original\n")
	_, stderr, err := runAuditCommand(t, "", "utils", "seed-pos", "--quiet", "-d", db, "-a", "--max-open-files", "10", "--plot-dir", plots, "--force")
	if err == nil || !strings.Contains(string(stderr), "should be >= 11") {
		t.Fatalf("invalid budget accepted: %v, %s", err, stderr)
	}
	got, err := os.ReadFile(sentinel)
	if err != nil || string(got) != "original\n" {
		t.Fatalf("invalid budget changed plot directory: %q, %v", got, err)
	}
}

func TestSeedPosPlotPathPreflight(t *testing.T) {
	for _, kind := range []string{"escaping-id", "escaping-extension", "result-collision", "plot-collision"} {
		t.Run(kind, func(t *testing.T) {
			db := writeSeedPosFixture(t, 2)
			dir := t.TempDir()
			plots := filepath.Join(dir, "plots")
			if err := os.Mkdir(plots, 0700); err != nil {
				t.Fatal(err)
			}
			sentinel := writeTestInput(t, filepath.Join(plots, "keep"), "original\n")
			logFile := writeTestInput(t, filepath.Join(dir, "run.log"), "original\n")
			output, ext := filepath.Join(dir, "result.tsv"), ".svg"
			switch kind {
			case "escaping-id":
				runSeedIndexCommand(t, "utils", "edit-genome-ids", "--quiet", "-d", db, "-p", "^g0$", "-r", "../outside")
			case "escaping-extension":
				ext = "/../../outside.svg"
			case "result-collision":
				output = filepath.Join(plots, "g0.svg")
			case "plot-collision":
				runSeedIndexCommand(t, "utils", "edit-genome-ids", "--quiet", "-d", db, "-p", "^g1$", "-r", "g0.seed_number")
			}
			_, stderr, err := runAuditCommand(t, "", "utils", "seed-pos", "--quiet", "-d", db, "-a", "--plot-dir", plots,
				"--plot-ext", ext, "--force", "-o", output, "--log", logFile)
			if err == nil || !strings.Contains(string(stderr), "plot") {
				t.Fatalf("plot path conflict accepted: %v, %s", err, stderr)
			}
			for _, path := range []string{sentinel, logFile} {
				if got, err := os.ReadFile(path); err != nil || string(got) != "original\n" {
					t.Fatalf("file changed before validation: %s: %q, %v", path, got, err)
				}
			}
			if _, err := os.Stat(filepath.Join(dir, "outside.svg")); !os.IsNotExist(err) {
				t.Fatalf("plot escaped directory: %v", err)
			}
		})
	}
}

func TestSeedPosResultBesidePlots(t *testing.T) {
	db := writeSeedPosFixture(t, 1)
	dir := t.TempDir()
	plots, target := filepath.Join(dir, "plots"), filepath.Join(dir, "old-plots")
	if err := os.Mkdir(target, 0700); err != nil {
		t.Fatal(err)
	}
	if err := os.Symlink(target, plots); err != nil {
		t.Fatal(err)
	}
	result := filepath.Join(plots, "seeds.tsv")
	runSeedIndexCommand(t, "utils", "seed-pos", "--quiet", "-d", db, "-a", "--plot-dir", plots, "--plot-ext", ".svg", "--slid-window", "40", "--slid-step", "20", "-o", result)
	data, err := os.ReadFile(result)
	if err != nil || bytes.Count(data, []byte("\n")) != 4 {
		t.Fatalf("result lost beside plots: %q, %v", data, err)
	}
	for _, name := range []string{"g0.svg", "g0.seed_number.svg"} {
		if _, err := os.Stat(filepath.Join(plots, name)); err != nil {
			t.Fatal(err)
		}
	}
}
