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
	"fmt"
	"math/rand"
	"os"
	"path/filepath"
	"reflect"
	"strconv"
	"strings"
	"testing"
)

func TestCombinationCount(t *testing.T) {
	if got, err := combinationCount(3, "genome IDs"); err != nil || got != 3 {
		t.Fatalf("count: got %d, err=%v, want 3", got, err)
	}
	if _, err := combinationCount(1, "genome IDs"); err == nil {
		t.Fatal("expected an error for fewer than two genome IDs")
	}
}

// TestComparePairThresholds exercises asymmetric AF/ANI through the CLI, with
// both input orders and both preparation paths. Filtering must use the pair.
func TestComparePairThresholds(t *testing.T) {
	dir := t.TempDir()
	rng := rand.New(rand.NewSource(84))
	randomSeq := func(n int) []byte {
		s := make([]byte, n)
		for i := range s {
			s[i] = "ACGT"[rng.Intn(4)]
		}
		return s
	}
	core := randomSeq(2040)
	write := func(name string, sequences ...[]byte) string {
		t.Helper()
		var data strings.Builder
		for i, s := range sequences {
			fmt.Fprintf(&data, ">s%d\n%s\n", i, s)
		}
		file := filepath.Join(dir, name+".fa")
		if err := os.WriteFile(file, []byte(data.String()), 0600); err != nil {
			t.Fatal(err)
		}
		return file
	}
	short := write("short", core)
	long := write("long", core, randomSeq(18360))
	contigs := [][]byte{core}
	for range 3 {
		mutated := append([]byte(nil), core...)
		for i := range mutated {
			if rng.Intn(10) == 0 {
				mutated[i] = "ACGT"[(strings.IndexByte("ACGT", mutated[i])+1)%4]
			}
		}
		contigs = append(contigs, mutated)
	}
	asymmetric := write("asymmetric", contigs...)
	for _, cache := range []string{"0", "1M"} {
		for _, reverse := range []bool{false, true} {
			for _, tc := range []struct {
				name, other, af, ani string
				ortho, wantMatch     bool
			}{
				{"ortho AF below", long, "50", "70", true, false},
				{"ortho AF equal", long, "10", "70", true, true},
				{"ordinary AF below", long, "50", "70", false, false},
				{"ordinary ANI either direction", asymmetric, "0", "99", false, true},
				{"ordinary ANI equal", asymmetric, "0", "100", false, true},
			} {
				t.Run(fmt.Sprintf("%s/cache=%s/reverse=%v", tc.name, cache, reverse), func(t *testing.T) {
					inputs := []string{short, tc.other}
					if reverse {
						inputs[0], inputs[1] = inputs[1], inputs[0]
					}
					args := []string{"genome", "compare", "--quiet", "-j", "2", "--max-genome-cache-memory", cache,
						"--min-af", tc.af, "--min-ani", tc.ani}
					if tc.ortho {
						args = append(args, "-O")
					}
					output := runSeedIndexCommand(t, append(args, inputs...)...)
					lines := strings.Split(strings.TrimSpace(string(output)), "\n")
					wantLines := 1
					if tc.wantMatch {
						wantLines++
					}
					if len(lines) != wantLines {
						t.Fatalf("got %d lines, want %d:\n%s", len(lines), wantLines, output)
					}
					if tc.other == asymmetric {
						fields := strings.Split(lines[1], "\t")
						a, _ := strconv.ParseFloat(fields[3], 64)
						b, _ := strconv.ParseFloat(fields[4], 64)
						if min(a, b) >= 99 || max(a, b) != 100 {
							t.Fatalf("fixture must have one ANI below 99 and one at 100: %v", fields[3:5])
						}
					}
				})
			}
		}
	}
}

func TestReadPairFilesSkipsHeaderInEachFile(t *testing.T) {
	dir := t.TempDir()
	file1 := filepath.Join(dir, "pairs1.tsv")
	file2 := filepath.Join(dir, "pairs2.tsv")
	if err := os.WriteFile(file1, []byte("genome1\tgenome2\nA\tB\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(file2, []byte("genome1\tgenome2\nC\tD\textra\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	got, err := readPairFiles([]string{file1, file2}, true)
	if err != nil {
		t.Fatal(err)
	}
	want := []string{"A", "B", "C", "D"}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("pairs: got %v, want %v", got, want)
	}
}

// TestCompareCaseInsensitive checks normalization before either cached or
// private fragment indexing, in both default and reciprocal-best-hit modes.
func TestCompareCaseInsensitive(t *testing.T) {
	dir := t.TempDir()
	rng := rand.New(rand.NewSource(719))
	sequence := make([]byte, 5100)
	for i := range sequence {
		sequence[i] = "ACGT"[rng.Intn(4)]
	}
	upper := string(sequence)
	files := make([]string, 0, 3)
	for i, s := range []string{upper, strings.ToLower(upper), upper[:2040] + strings.ToLower(upper[2040:3060]) + upper[3060:]} {
		file := filepath.Join(dir, fmt.Sprintf("case%d.fa", i))
		if err := os.WriteFile(file, []byte(">seq\n"+s+"\n"), 0600); err != nil {
			t.Fatal(err)
		}
		files = append(files, file)
	}
	for _, cache := range []string{"0", "1M"} {
		for _, ortho := range []bool{false, true} {
			t.Run(fmt.Sprintf("cache=%s/ortho=%v", cache, ortho), func(t *testing.T) {
				args := []string{"genome", "compare", "--quiet", "-j", "2", "--max-genome-cache-memory", cache}
				if ortho {
					args = append(args, "-O")
				}
				output := runSeedIndexCommand(t, append(args, files...)...)
				lines := strings.Split(strings.TrimSpace(string(output)), "\n")
				if len(lines) != 4 {
					t.Fatalf("expected three pairwise matches, got %s", output)
				}
				for _, line := range lines[1:] {
					fields := strings.Split(line, "\t")
					if strings.Join(fields[2:7], "\t") != "100.0000\t100.0000\t100.0000\t100.000\t100.000" {
						t.Fatalf("case changed ANI/AF: %s", line)
					}
				}
			})
		}
	}
}

// TestReadGenomeFromFileUppercase preserves ambiguity symbols and contig sizes.
func TestReadGenomeFromFileUppercase(t *testing.T) {
	file := filepath.Join(t.TempDir(), "mixed.fa")
	if err := os.WriteFile(file, []byte(">one\nacgtnryACGTNRY\n>two\nws\n"), 0600); err != nil {
		t.Fatal(err)
	}
	q, err := ReadGenomeFromFile(file, nil, false)
	if err != nil {
		t.Fatal(err)
	}
	defer RecycleGQuery(q)
	if q.genomeSize != 16 || len(q.seqs) != 2 || string(*q.seqs[0]) != "ACGTNRYACGTNRY" || string(*q.seqs[1]) != "WS" {
		t.Fatalf("unexpected normalization: size=%d seqs=%v", q.genomeSize, q.seqs)
	}
}
