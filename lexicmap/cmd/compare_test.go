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
	"compress/gzip"
	"encoding/json"
	"fmt"
	"math/rand"
	"os"
	"path/filepath"
	"reflect"
	"slices"
	"strconv"
	"strings"
	"testing"
)

// A zero-score terminal extension must not change ANI/AF or remove an AF-
// qualified pair when both input sequences are reverse-complemented.
func TestCompareZeroScoreExtension(t *testing.T) {
	data, err := os.ReadFile("testdata/chain-zero-score-extension.json")
	if err != nil {
		t.Fatal(err)
	}
	var fixture struct{ Query, Target string }
	if err := json.Unmarshal(data, &fixture); err != nil {
		t.Fatal(err)
	}
	const want = "89.4604\t89.1457\t89.1457\t96.614\t99.081\t1\t1004\t1\t979"
	for _, reverse := range []bool{false, true} {
		q, s := []byte(fixture.Query), []byte(fixture.Target)
		if reverse {
			RC(q)
			RC(s)
		}
		dir := t.TempDir()
		qfile := writeTestInput(t, filepath.Join(dir, "q.fa"), ">q\n"+string(q)+"\n")
		sfile := writeTestInput(t, filepath.Join(dir, "s.fa"), ">s\n"+string(s)+"\n")
		for _, cache := range []string{"0", "1M"} {
			for _, af := range []string{"0", "95"} {
				out := runSeedIndexCommand(t, "genome", "compare", "--quiet", "-j", "2", "--max-genome-cache-memory", cache, "--min-af", af, qfile, sfile)
				lines := strings.Split(strings.TrimSpace(string(out)), "\n")
				if len(lines) != 2 {
					t.Errorf("reverse=%t cache=%s min-af=%s: got %d hits, want one", reverse, cache, af, len(lines)-1)
					continue
				}
				if got := strings.Join(strings.Split(lines[1], "\t")[2:], "\t"); got != want {
					t.Errorf("reverse=%t cache=%s min-af=%s: got %s, want %s", reverse, cache, af, got, want)
				}
			}
		}
	}
}

// Variable anchor lengths must not drop a terminal region and fail the AF
// filter in just one orientation of the same sequence pair.
func TestCompareVariableAnchorCoverage(t *testing.T) {
	dir := t.TempDir()
	for _, reverse := range []bool{false, true} {
		q, s := []byte(variableAnchorQuery), []byte(variableAnchorTarget)
		if reverse {
			RC(q)
			RC(s)
		}
		qfile := writeTestInput(t, filepath.Join(dir, "q.fa"), ">q\n"+string(q)+"\n")
		sfile := writeTestInput(t, filepath.Join(dir, "s.fa"), ">s\n"+string(s)+"\n")
		for _, ortho := range []bool{false, true} {
			// Match the fragment size to this short fixture so OrthoANI's
			// shared-k-mer screen does not use a 1020-base count threshold.
			args := []string{"genome", "compare", "--quiet", "-j", "2", "--frag-size", "156", "--min-af", "95", "--max-genome-cache-memory", "0", qfile, sfile}
			if ortho {
				args = append(args, "-O")
			}
			output := runSeedIndexCommand(t, args...)
			lines := strings.Split(strings.TrimSpace(string(output)), "\n")
			if len(lines) != 2 {
				t.Errorf("reverse=%t ortho=%t: got %d hits, want one", reverse, ortho, len(lines)-1)
				continue
			}
			want := "91.0256\t91.0256\t91.0256\t100.000\t100.000\t1\t156\t1\t156"
			if got := strings.Join(strings.Split(lines[1], "\t")[2:], "\t"); got != want {
				t.Errorf("reverse=%t ortho=%t: got %s, want %s", reverse, ortho, got, want)
			}
		}
	}
}

// A 50-base tail is shorter than half a 101-base fragment and must merge,
// preserving the full 101-base shared region in both comparison directions.
func TestCompareOddFragmentBoundary(t *testing.T) {
	dir := t.TempDir()
	sequence := "CCCCGCACAACTCGTTGAGCTACGGGTATTACCCGCGTTTAGTGGCTAGCGACTCGTGGACTTGCTGTACTGTCTACGGGCGTCAACTTGATAATCCCAAAAAAGCTTGGAATTCTTAAAAAGATGGTGTACGTCCCTCATACTTCGTATT"
	want := "80.1587\t100.0000\t100.0000\t66.887\t100.000\t1\t151\t1\t101"
	for _, reverse := range []bool{false, true} {
		q, s := []byte(sequence), []byte(sequence[:101])
		if reverse {
			RC(q)
			RC(s)
		}
		qfile := writeTestInput(t, filepath.Join(dir, "q.fa"), ">q\n"+string(q)+"\n")
		sfile := writeTestInput(t, filepath.Join(dir, "s.fa"), ">s\n"+string(s)+"\n")
		for _, ortho := range []bool{false, true} {
			for _, cache := range []string{"0", "1M"} {
				args := []string{"genome", "compare", "--quiet", "-j", "2", "--frag-size", "101", "--min-af", "0", "--min-ani", "0", "--max-genome-cache-memory", cache, qfile, sfile}
				if ortho {
					args = append(args, "-O")
				}
				output := runSeedIndexCommand(t, args...)
				lines := strings.Split(strings.TrimSpace(string(output)), "\n")
				if len(lines) != 2 {
					t.Fatalf("reverse=%t ortho=%t cache=%s: expected one comparison", reverse, ortho, cache)
				}
				if got := strings.Join(strings.Split(lines[1], "\t")[2:], "\t"); got != want {
					t.Errorf("reverse=%t ortho=%t cache=%s: got %s, want %s", reverse, ortho, cache, got, want)
				}
			}
		}
	}
}

// Low-complexity filtering leaves a two-base tail outside the exact anchors.
// Extending that tail must retain full self-coverage and pass a 99% AF filter.
func TestCompareTwoBaseExtension(t *testing.T) {
	dir := t.TempDir()
	sequence := twoBaseExtensionSequence
	want := "100.0000\t100.0000\t100.0000\t100.000\t100.000\t1\t89\t1\t89"
	for _, reverse := range []bool{false, true} {
		seq := []byte(sequence)
		if reverse {
			RC(seq)
		}
		qfile := writeTestInput(t, filepath.Join(dir, "q.fa"), ">q\n"+string(seq)+"\n")
		sfile := writeTestInput(t, filepath.Join(dir, "s.fa"), ">s\n"+string(seq)+"\n")
		for _, cache := range []string{"0", "1M"} {
			for _, ortho := range []bool{false, true} {
				for _, af := range []string{"0", "99"} {
					args := []string{"genome", "compare", "--quiet", "-j", "2", "--max-genome-cache-memory", cache, "--min-af", af, qfile, sfile}
					if ortho {
						args = append(args, "-O")
					}
					output := runSeedIndexCommand(t, args...)
					lines := strings.Split(strings.TrimSpace(string(output)), "\n")
					if len(lines) != 2 {
						t.Errorf("reverse=%t cache=%s ortho=%t min-af=%s: expected one comparison, got %d", reverse, cache, ortho, af, len(lines)-1)
						continue
					}
					if got := strings.Join(strings.Split(lines[1], "\t")[2:], "\t"); got != want {
						t.Errorf("reverse=%t cache=%s ortho=%t min-af=%s: got %s, want %s", reverse, cache, ortho, af, got, want)
					}
				}
			}
		}
	}
}

// A 24-base overlap between 31-base anchors is below the 80% trimming cutoff.
// Keeping the terminal anchor preserves the full query and its AF-qualified hit.
func TestCompareTerminalAnchorOverlap(t *testing.T) {
	dir := t.TempDir()
	query := "TTTCGTCATGTGTGAGCCTAATGTGTGAGCCAAATCTCGGCGCTTGTGTGATTCTTGACATCGCGGCACAAACTTATTACGAGGTTACCAGATATAATGAGGCCTGAATACAGGTATACGGATGTGTGCCATCACGAGCCCCGTATTGTTACCCGTGCACGAACAGGTGGCCCTTTTATGCCGCCCAGACCAGTAAAGAAATGGGCTACATGTTGTGTAGAATCCCCCGAA"
	subject := query[:20] + query[7:] // Insert a 13-base tandem repeat.
	want := "97.2632\t94.6721\t94.6721\t100.000\t94.672\t1\t231\t1\t244"
	for _, reverse := range []bool{false, true} {
		q, s := []byte(query), []byte(subject)
		if reverse {
			RC(q)
			RC(s)
		}
		qfile := writeTestInput(t, filepath.Join(dir, "q.fa"), ">q\n"+string(q)+"\n")
		sfile := writeTestInput(t, filepath.Join(dir, "s.fa"), ">s\n"+string(s)+"\n")
		for _, cache := range []string{"0", "1M"} {
			for _, af := range []string{"0", "94"} {
				output := runSeedIndexCommand(t, "genome", "compare", "--quiet", "-j", "2", "--max-genome-cache-memory", cache, "--min-af", af, qfile, sfile)
				lines := strings.Split(strings.TrimSpace(string(output)), "\n")
				if len(lines) != 2 {
					t.Errorf("reverse=%t cache=%s min-af=%s: expected one comparison, got %d", reverse, cache, af, len(lines)-1)
					continue
				}
				if got := strings.Join(strings.Split(lines[1], "\t")[2:], "\t"); got != want {
					t.Errorf("reverse=%t cache=%s min-af=%s: got %s, want %s", reverse, cache, af, got, want)
				}
			}
		}
	}
}

// A repeated terminal anchor must not change AF filtering when both inputs
// are reverse-complemented, with or without the prepared-genome cache.
func TestCompareTerminalRepeatOrientation(t *testing.T) {
	dir := t.TempDir()
	query := "GGCCCAGTCCAGATCCTCGGAAGTCCCAACTAACGAATAAGTAGATCCTTCTAAATAGTAGTATACGAATGTAATCAACATATTACCGCACCCAACTTTGATTGGGTCATAAACAAACATCATGTTTCTGGAGGCGGTTTGTGAGACCGACAAACGCGATAGAGTAATCTTATAAATCAGGCAAATACTCGGTCCTTTCTACGGGTCGCGGCCGCGTGGTCATAGGCAGTT"
	subject := "GGCCCAGTCCAGATCCTCGGAAGTCCCAACTAACGATTAAGTTGATCCGTCTAAGTAGTATTATACCAATGTCATCAATAGGCCCAGTCCAGATCCTCGGAAGTCCCAACTAACGAACATAATGTTCCTGGACGCGGTGTGTGACACCGATAAACGGGATAGCGTAATTTTATAGATCAGCCAAATGCTCGGGCCTTTGTACGGGTCGCGGCCGCGTGGTCATAGGCAGTT"
	var want string
	for _, cache := range []string{"0", "1M"} {
		for _, reverse := range []bool{false, true} {
			q, s := []byte(query), []byte(subject)
			if reverse {
				RC(q)
				RC(s)
			}
			qfile := writeTestInput(t, filepath.Join(dir, "q.fa"), ">q\n"+string(q)+"\n")
			sfile := writeTestInput(t, filepath.Join(dir, "s.fa"), ">s\n"+string(s)+"\n")
			for _, af := range []string{"0", "50"} {
				output := runSeedIndexCommand(t, "genome", "compare", "--quiet", "-j", "2", "--max-genome-cache-memory", cache, "--min-af", af, qfile, sfile)
				lines := strings.Split(strings.TrimSpace(string(output)), "\n")
				if af == "50" {
					if len(lines) != 1 {
						t.Errorf("cache=%s reverse=%t: comparison below AF threshold was retained", cache, reverse)
					}
					continue
				}
				if len(lines) != 2 {
					t.Fatalf("cache=%s reverse=%t: expected one comparison, got %d", cache, reverse, len(lines)-1)
				}
				got := strings.Join(strings.Split(lines[1], "\t")[2:], "\t")
				if want == "" {
					want = got
				} else if got != want {
					t.Errorf("cache=%s reverse=%t: metrics %q differ from %q", cache, reverse, got, want)
				}
			}
		}
	}
}

// Extended HSPs must survive -l in both comparison modes and both orientations.
func TestCompareExtendedAlignmentLength(t *testing.T) {
	dir := t.TempDir()
	for _, reverse := range []bool{false, true} {
		q, s := terminalMismatchSequences()
		if reverse {
			RC(q)
			RC(s)
		}
		qfile := writeTestInput(t, filepath.Join(dir, "q.fa"), ">q\n"+string(q)+"\n")
		sfile := writeTestInput(t, filepath.Join(dir, "s.fa"), ">s\n"+string(s)+"\n")
		for _, ortho := range []bool{false, true} {
			for _, minimum := range []int{50, 980, 1000, 1001} {
				args := []string{"genome", "compare", "--quiet", "-j", "2", "-l", strconv.Itoa(minimum), qfile, sfile}
				if ortho {
					args = append(args, "--OrthoANI", "--frag-size", "1000")
				}
				out := runSeedIndexCommand(t, args...)
				lines := strings.Split(strings.TrimSpace(string(out)), "\n")
				want := 2
				if minimum > 1000 {
					want = 1
				}
				if len(lines) != want {
					t.Fatalf("reverse=%t/ortho=%t/minimum=%d: got %d lines, want %d", reverse, ortho, minimum, len(lines), want)
				}
			}
		}
	}
}

func TestCombinationCount(t *testing.T) {
	if got, err := combinationCount(3, "genome IDs"); err != nil || got != 3 {
		t.Fatalf("count: got %d, err=%v, want 3", got, err)
	}
	if _, err := combinationCount(1, "genome IDs"); err == nil {
		t.Fatal("expected an error for fewer than two genome IDs")
	}
}

// TestCompareMinAlignmentLength checks ordinary and reciprocal-best-hit
// comparisons after extending a short anchor in 100-base inputs. The final
// alignment spans 84 bases after mismatching ends are trimmed.
func TestCompareMinAlignmentLength(t *testing.T) {
	dir := t.TempDir()
	query := []byte("TGGTCCGCGTAGTACAGAGACTGCCTTGAAATTTCGTGATTATGCAGCGCTACTTTGGAGTCCCGACTATGACTGGCTTTCTCTTTGGCGGAGGCCTACC")
	subject := []byte("TGGTCCGGGTAGTACATAGACTGCCTTGAAATTTGGTGATTATTCAGCGCTAGTTTGGAGTGCCGACTATTACTGGCTTACTCTTTGGGGGAGGCCTCCC")
	for _, reverse := range []bool{false, true} {
		q, s := slices.Clone(query), slices.Clone(subject)
		if reverse {
			RC(q)
			RC(s)
		}
		qfile := writeTestInput(t, filepath.Join(dir, fmt.Sprintf("q-%t.fa", reverse)), ">q\n"+string(q)+"\n")
		sfile := writeTestInput(t, filepath.Join(dir, fmt.Sprintf("s-%t.fa", reverse)), ">s\n"+string(s)+"\n")
		for _, ortho := range []bool{false, true} {
			for _, minimum := range []int{0, 84, 85} {
				t.Run(fmt.Sprintf("reverse=%t/ortho=%t/minimum=%d", reverse, ortho, minimum), func(t *testing.T) {
					args := []string{"genome", "compare", "--quiet", "-j", "2"}
					if minimum != 0 { // also exercise the default threshold
						args = append(args, "-l", strconv.Itoa(minimum))
					}
					if ortho {
						args = append(args, "-O", "--frag-size", "100")
					}
					output := runSeedIndexCommand(t, append(args, qfile, sfile)...)
					lines := strings.Split(strings.TrimSpace(string(output)), "\n")
					want := 2
					if minimum > 84 {
						want = 1
					}
					if len(lines) != want {
						t.Fatalf("got %d lines, want %d:\n%s", len(lines), want, output)
					}
				})
			}
		}
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

// TestCompareStdinReplay covers both pair roles, cache eviction, and concurrent pairs.
func TestCompareStdinReplay(t *testing.T) {
	dir := t.TempDir()
	data := ">seq\n" + indexTestDNA(2040) + "\n"
	files := make([]string, 3)
	for i := range files {
		files[i] = writeTestInput(t, filepath.Join(dir, fmt.Sprintf("g%d.fa", i)), data)
	}
	normalize := func(output []byte, input int) string {
		lines := strings.Split(strings.TrimSpace(string(output)), "\n")
		if len(lines) != 4 {
			t.Fatalf("expected three pairs: %s", output)
		}
		for i := 1; i < len(lines); i++ {
			fields := strings.Split(lines[i], "\t")
			for j := 0; j < 2; j++ {
				if fields[j] == fmt.Sprintf("g%d", input) {
					fields[j] = "-"
				}
			}
			lines[i] = strings.Join(fields, "\t")
		}
		slices.Sort(lines[1:])
		return strings.Join(lines, "\n")
	}
	for _, ortho := range []bool{false, true} {
		args := []string{"genome", "compare", "--quiet", "-j", "2"}
		if ortho {
			args = append(args, "-O")
		}
		baseline := runSeedIndexCommand(t, append(slices.Clone(args), files...)...)
		for _, cache := range []string{"0", "1", "4G"} {
			for _, threads := range []string{"2", "4"} {
				for position := 0; position < 3; position++ {
					inputs := slices.Clone(files)
					inputs[position] = "-"
					command := append(slices.Clone(args), "--max-genome-cache-memory", cache, "-j", threads)
					stdin := data
					if position == 2 {
						var compressed bytes.Buffer
						w := gzip.NewWriter(&compressed)
						if _, err := w.Write([]byte(data)); err != nil {
							t.Fatal(err)
						}
						if err := w.Close(); err != nil {
							t.Fatal(err)
						}
						stdin = compressed.String()
					}
					output, stderr, err := runAuditCommand(t, stdin, append(command, inputs...)...)
					if err != nil {
						t.Fatalf("cache=%s threads=%s stdin=%d: %v %s", cache, threads, position, err, stderr)
					}
					if normalize(output, position) != normalize(baseline, position) {
						t.Fatalf("stdin changed results: %s vs %s", output, baseline)
					}
				}
			}
		}
	}
}

func TestCompareSeedlessGenomes(t *testing.T) {
	dir := t.TempDir()
	long := ">seq\n" + indexTestDNA(2040) + "\n"
	a := writeTestInput(t, filepath.Join(dir, "a.fa"), long)
	b := writeTestInput(t, filepath.Join(dir, "b.fa"), long)
	for _, ortho := range []bool{false, true} {
		for _, cache := range []string{"0", "1", "4G"} {
			args := []string{"genome", "compare", "--quiet", "-j", "4", "--max-genome-cache-memory", cache}
			if ortho {
				args = append(args, "-O")
			}
			baseline := runSeedIndexCommand(t, append(slices.Clone(args), a, b)...)
			for i, data := range []string{"", ">short\nACGT\n", strings.Repeat(">short\nACGT\n", 8)} {
				short := writeTestInput(t, filepath.Join(dir, fmt.Sprintf("short%d.fa", i)), data)
				output := runSeedIndexCommand(t, append(slices.Clone(args), short, a, b)...)
				if string(output) != string(baseline) {
					t.Fatalf("seedless genome changed valid pair: %s vs %s", output, baseline)
				}
			}
		}
	}
}

func TestReadPairsSkipsShortRows(t *testing.T) {
	for _, header := range []bool{false, true} {
		data := "short\na\tb\nshort\nc\td\textra\n"
		if header {
			data = "query\tsubject\n" + data
		}
		file := writeTestInput(t, filepath.Join(t.TempDir(), "pairs.tsv"), data)
		got, err := readPairs(file, header)
		if err != nil || !reflect.DeepEqual(got, []string{"a", "b", "c", "d"}) {
			t.Fatalf("header=%v: %v, %v", header, got, err)
		}
	}
}
