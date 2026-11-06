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
	"encoding/csv"
	"errors"
	"fmt"
	"io"
	"math/rand"
	"os"
	"os/exec"
	"path/filepath"
	"reflect"
	"slices"
	"strconv"
	"strings"
	"testing"
)

// mergeTestInput exposes close failures without relying on filesystem behavior.
type mergeTestInput struct {
	io.Reader
	closeErr error
	closes   int
}

func (r *mergeTestInput) Close() error { r.closes++; return r.closeErr }

// Direct reads must work without a temporary directory whenever the configured
// budget can hold every input, including equality after reserved descriptors.
func TestMergeSearchResultsFileBudget(t *testing.T) {
	dir := t.TempDir()
	files := make([]string, 200)
	for i := range files {
		row := fmt.Sprintf("q\t4\t1\tg%d\ts\t100\t1\t1\t100\t4\t100\t0\t1\t4\t1\t4\t+\t4\t0\t%d\n", i, 1000-i)
		files[i] = writeTestInput(t, filepath.Join(dir, fmt.Sprintf("%d.tsv", i)), auditSearchHeader+"\n"+row)
	}
	// A regular file cannot contain a spool file. Direct reading must ignore it.
	unusableTemp := writeTestInput(t, filepath.Join(dir, "not-a-directory"), "preserve\n")
	for _, name := range []string{"TMPDIR", "TMP", "TEMP"} {
		t.Setenv(name, unusableTemp)
	}
	var expected []byte
	for _, budget := range []int{0, len(files) + reservedOpenFiles, len(files) + reservedOpenFiles - 1} {
		args := []string{"utils", "merge-search-results", "--quiet"}
		if budget != 0 {
			args = append(args, "--max-open-files", strconv.Itoa(budget))
		}
		out, stderr, err := runAuditCommand(t, "", append(args, files...)...)
		if budget == len(files)+reservedOpenFiles-1 {
			if err == nil || !strings.Contains(string(stderr), "failed to create temporary merge input") {
				t.Fatalf("insufficient direct-reader budget did not use temporary storage: %v: %s", err, stderr)
			}
			continue
		}
		if err != nil || bytes.Count(out, []byte("\n")) != len(files)+1 {
			t.Fatalf("budget=%d: direct reading failed: %v: %s", budget, err, stderr)
		}
		if expected == nil {
			expected = out
		} else if !bytes.Equal(out, expected) {
			t.Fatal("file budget changed merged output")
		}
	}
}

// Reject an unusable budget before opening either output or log files.
func TestMergeSearchResultsInvalidFileBudget(t *testing.T) {
	dir := t.TempDir()
	output := writeTestInput(t, filepath.Join(dir, "output.tsv"), "preserve\n")
	logFile := writeTestInput(t, filepath.Join(dir, "run.log"), "preserve\n")
	for _, budget := range []int{0, reservedOpenFiles} {
		_, stderr, err := runAuditCommand(t, "", "utils", "merge-search-results", "--quiet",
			"--max-open-files", strconv.Itoa(budget), "--log", logFile, "-o", output)
		if err == nil || !strings.Contains(string(stderr), "--max-open-files") {
			t.Fatalf("invalid file budget accepted: %v: %s", err, stderr)
		}
		for _, file := range []string{output, logFile} {
			if got, err := os.ReadFile(file); err != nil || string(got) != "preserve\n" {
				t.Fatalf("file changed before budget validation: %s: %q: %v", file, got, err)
			}
		}
	}
}

// EOF, validation failure, and early cancellation each close the source once.
func TestSearchResultReaderClosesOnce(t *testing.T) {
	row := "q\t4\t2\ta\ts\t100\t1\t1\t100\t4\t100\t0\t1\t4\t1\t4\t+\t4\t0\t10\n"
	closeErr := errors.New("close failed")
	for _, mode := range []string{"eof", "parse-error", "early-close"} {
		body := auditSearchHeader + "\n" + row
		if mode == "parse-error" {
			body += "invalid\n"
		}
		if mode == "early-close" {
			body += strings.Replace(row, "\ta\t", "\tb\t", 1)
		}
		source := &mergeTestInput{Reader: strings.NewReader(body), closeErr: closeErr}
		reader := newSearchResultReader("input.tsv", "", 65536, source)
		result := reader.Next()
		if mode == "early-close" {
			if result == nil {
				t.Fatal("missing first group")
			}
			RecycleSearchResultOfAGenome(result)
			if err := reader.Close(); !errors.Is(err, closeErr) {
				t.Fatalf("lost close error: %v", err)
			}
		} else {
			if result != nil {
				RecycleSearchResultOfAGenome(result)
				t.Fatal("returned a group after input failure")
			}
			if !errors.Is(reader.Err(), closeErr) {
				t.Fatalf("lost close error: %v", reader.Err())
			}
			if mode == "parse-error" && !strings.Contains(reader.Err().Error(), "columns") {
				t.Fatal("lost parse error")
			}
		}
		if err := reader.Close(); err != nil || source.closes != 1 || reader.Next() != nil {
			t.Fatalf("%s: close count=%d, error=%v", mode, source.closes, err)
		}
	}
}

func TestMergeSearchResultsCommand(t *testing.T) {
	if os.Getenv("LEXICMAP_MERGE_TEST_HELPER") != "1" {
		return
	}
	for i, arg := range os.Args {
		if arg == "--" {
			RootCmd.SetArgs(os.Args[i+1:])
			if err := RootCmd.Execute(); err != nil {
				os.Exit(1)
			}
			// Behave like the CLI; the test runner's PASS line is not TSV output.
			os.Exit(0)
		}
	}
	t.Fatal("missing command arguments")
}

// Mixed formats must preserve the selected width without borrowing optional
// alignment data from an earlier record, including within one subject group.
func TestMergeSearchResultsMixedFormats(t *testing.T) {
	dir := t.TempDir()
	const extraHeader = "\tcigar\tqseq\tsseq\talign"
	const extra = "\t4M\tACGT\tACGT\t||||"
	row := func(genome string, score int) string {
		return fmt.Sprintf("q\t4\t2\t%s\ts\t100\t1\t1\t100\t4\t100\t0\t1\t4\t1\t4\t+\t4\t0\t%d", genome, score)
	}
	for _, fullFirst := range []bool{false, true} {
		for _, together := range []bool{false, true} {
			for _, sameGenome := range []bool{false, true} {
				if sameGenome && !together {
					continue // distinct genome IDs are required across indexes
				}
				for _, topN := range []string{"0", "2"} {
					t.Run(fmt.Sprintf("fullFirst=%t/together=%t/sameGenome=%t/topN=%s", fullFirst, together, sameGenome, topN), func(t *testing.T) {
						secondGenome := "b"
						if sameGenome {
							secondGenome = "a"
						}
						first, second := row("a", 100), row(secondGenome, 90)
						firstHeader, secondHeader := auditSearchHeader, auditSearchHeader
						if fullFirst {
							first += extra
							firstHeader += extraHeader
						} else {
							second += extra
							secondHeader += extraHeader
						}
						files := []string{filepath.Join(dir, "first.tsv")}
						body := firstHeader + "\n" + first + "\n"
						if together {
							body += second + "\n"
						} else {
							files = append(files, writeTestInput(t, filepath.Join(dir, "second.tsv"), secondHeader+"\n"+second+"\n"))
						}
						writeTestInput(t, files[0], body)
						output := runSeedIndexCommand(t, append([]string{"utils", "merge-search-results", "--quiet", "-n", topN}, files...)...)
						// Keep trailing tabs: missing optional fields are four empty cells.
						lines := strings.Split(strings.TrimSuffix(string(output), "\n"), "\n")
						if len(lines) != 3 {
							t.Fatalf("got %d lines, want header plus two records", len(lines))
						}
						width := 20
						if fullFirst {
							width = 24
						}
						for i, line := range lines {
							fields := strings.Split(line, "\t")
							if len(fields) != width {
								t.Fatalf("line %d has %d columns, want %d", i+1, len(fields), width)
							}
							if i == 2 && fullFirst && strings.Join(fields[20:], "") != "" {
								t.Fatalf("record without alignment data inherited %q", fields[20:])
							}
							if i == 1 && fullFirst && "\t"+strings.Join(fields[20:], "\t") != extra {
								t.Fatal("alignment data changed")
							}
						}
					})
				}
			}
		}
	}
}

func TestMergeSearchResultsTopN(t *testing.T) {
	header := "query\tqlen\thits\tsgenome\tsseqid\tqcovGnm\tcls\thsp\tqcovHSP\talenHSP\tpident\tgaps\tqstart\tqend\tsstart\tsend\tsstr\tslen\tevalue\tbitscore"
	row := func(genome, hsp string, bitscore, pident int) string {
		return fmt.Sprintf("q\t100\t3\t%s\tseq%s\t100\t1\t%s\t100\t100\t%d\t0\t1\t100\t1\t100\t+\t100\t1e-10\t%d\n", genome, hsp, hsp, pident, bitscore)
	}
	// Genome a's best HSP is not its first row. Genomes b and c have equal
	// alignment scores despite different bitscores and identities.
	shardA := row("a", "1", 10, 100) + row("a", "2", 100, 100) + row("d", "1", 100, 75) + row("f", "1", 60, 100)
	shardB := row("b", "1", 100, 90) + row("c", "1", 90, 100) + row("e", "1", 70, 100)
	for _, extra := range []bool{false, true} {
		t.Run(fmt.Sprintf("extra=%v", extra), func(t *testing.T) {
			dir := t.TempDir()
			format := func(body string) string {
				h := header
				if extra {
					h += "\tcigar\tqseq\tsseq\talign"
					body = strings.ReplaceAll(body, "\n", "\t100M\tACGT\tACGT\t||||\n")
				}
				return h + "\n" + body
			}
			paths := []string{filepath.Join(dir, "a.tsv"), filepath.Join(dir, "b.tsv.gz"), filepath.Join(dir, "empty.tsv")}
			inputs := []string{format(shardA), format(shardB), format("")}
			for i, path := range paths {
				data := []byte(inputs[i])
				if strings.HasSuffix(path, ".gz") {
					var compressed bytes.Buffer
					gz := gzip.NewWriter(&compressed)
					if _, err := gz.Write(data); err != nil {
						t.Fatal(err)
					}
					if err := gz.Close(); err != nil {
						t.Fatal(err)
					}
					data = compressed.Bytes()
				}
				if err := os.WriteFile(path, data, 0600); err != nil {
					t.Fatal(err)
				}
			}
			for _, tt := range []struct {
				name  string
				topN  int
				files []string
				query string
				stdin string
				want  []string
			}{
				{name: "unlimited", files: paths[:2], want: []string{"a", "b", "c", "d", "e", "f"}},
				{name: "top one", topN: 1, files: paths[:2], want: []string{"a"}},
				{name: "cutoff ties", topN: 2, files: paths[:2], want: []string{"a", "b", "c"}},
				{name: "tie group ends at N", topN: 3, files: paths[:2], want: []string{"a", "b", "c"}},
				{name: "unique cutoff", topN: 4, files: paths[:2], want: []string{"a", "b", "c", "d"}},
				{name: "below limit", topN: 20, files: paths[:2], want: []string{"a", "b", "c", "d", "e", "f"}},
				{name: "single file ties", topN: 1, files: paths[1:2], want: []string{"b", "c"}},
				{name: "stdin ties", topN: 1, files: []string{"-"}, stdin: format(shardB), want: []string{"b", "c"}},
				{name: "empty shard", topN: 1, files: []string{paths[2], paths[0]}, want: []string{"a"}},
				{name: "selected query", topN: 2, query: "q", files: paths[:2], want: []string{"a", "b", "c"}},
				{name: "missing query", topN: 2, query: "missing", files: paths[:2], want: []string{}},
			} {
				t.Run(tt.name, func(t *testing.T) {
					args := []string{"utils", "merge-search-results", "--quiet", "-b", "64K", "-n", strconv.Itoa(tt.topN), "-q", tt.query}
					args = append(args, tt.files...)
					command := exec.Command(os.Args[0], append([]string{"-test.run=^TestMergeSearchResultsCommand$", "--"}, args...)...)
					// Filtering must work even when no temporary directory is available.
					command.Env = append(os.Environ(), "LEXICMAP_MERGE_TEST_HELPER=1", "TMPDIR="+filepath.Join(dir, "unavailable"))
					command.Stdin = strings.NewReader(tt.stdin)
					var output, stderr bytes.Buffer
					command.Stdout, command.Stderr = &output, &stderr
					if err := command.Run(); err != nil {
						t.Fatalf("merge failed: %v\n%s", err, stderr.String())
					}
					parser := csv.NewReader(&output)
					parser.Comma = '\t'
					records, err := parser.ReadAll()
					if err != nil {
						t.Fatal(err)
					}
					seen := make(map[string]int)
					previousScore := float64(1e10)
					for _, record := range records[min(1, len(records)):] {
						if record[2] != strconv.Itoa(len(tt.want)) {
							t.Fatalf("hits=%s, want %d", record[2], len(tt.want))
						}
						seen[record[3]]++
						if seen[record[3]] == 1 {
							score := map[string]float64{"a": 10000, "b": 9000, "c": 9000, "d": 7500, "e": 7000, "f": 6000}[record[3]]
							if score > previousScore {
								t.Fatal("genomes are not sorted by descending alignment score")
							}
							previousScore = score
						}
						original := append([]string(nil), record...)
						original[2] = "3"
						if !strings.Contains(inputs[0]+inputs[1], strings.Join(original, "\t")+"\n") {
							t.Fatalf("alignment row changed: %v", record)
						}
						if extra && strings.Join(record[20:], "\t") != "100M\tACGT\tACGT\t||||" {
							t.Fatal("alignment output fields changed")
						}
					}
					got := make([]string, 0, len(seen))
					for genome := range seen {
						got = append(got, genome)
					}
					slices.Sort(got)
					if !reflect.DeepEqual(got, tt.want) {
						t.Fatalf("got genomes %v, want %v", got, tt.want)
					}
					if seen["a"] != 0 && seen["a"] != 2 {
						t.Fatal("selected genome lost an HSP")
					}
				})
			}
		})
	}
}

// Negative scores remain ordered; zero is not an implicit alignment in a genome.
func TestMergeSearchResultsNegativeScores(t *testing.T) {
	for _, separate := range []bool{false, true} {
		for _, topN := range []int{1, 2, 3} {
			t.Run(fmt.Sprintf("separate=%t/topN=%d", separate, topN), func(t *testing.T) {
				dir := t.TempDir()
				row := func(genome string, score int) string {
					return fmt.Sprintf("q\t4\t4\t%s\ts\t100\t1\t1\t100\t4\t75.000\t0\t1\t4\t1\t4\t+\t4\t1e10\t%d\n", genome, score)
				}
				// A later HSP raises a's negative maximum; b ties that maximum.
				groups := []string{row("a", -20) + row("a", -10), row("c", -30), row("b", -10), row("z", 0)}
				var files []string
				if separate {
					for i, group := range groups {
						files = append(files, writeTestInput(t, filepath.Join(dir, fmt.Sprintf("%d.tsv", i)), auditSearchHeader+"\n"+strings.ReplaceAll(group, "\t4\t4\t", "\t4\t1\t")))
					}
				} else {
					files = []string{writeTestInput(t, filepath.Join(dir, "all.tsv"), auditSearchHeader+"\n"+strings.Join(groups, ""))}
				}
				out := runSeedIndexCommand(t, append([]string{"utils", "merge-search-results", "--quiet", "-n", strconv.Itoa(topN)}, files...)...)
				lines := strings.Split(strings.TrimSpace(string(out)), "\n")
				want := map[string]int{"z": 1}
				if topN > 1 {
					want["a"], want["b"] = 2, 1
				}
				got := make(map[string]int)
				for i, line := range lines[1:] {
					fields := strings.Split(line, "\t")
					if len(fields) != 20 || fields[2] != strconv.Itoa(len(want)) {
						t.Fatalf("invalid merged row: %q", line)
					}
					if i == 0 && fields[3] != "z" {
						t.Fatalf("zero score must precede negative scores: %q", line)
					}
					got[fields[3]]++
				}
				if !reflect.DeepEqual(got, want) {
					t.Fatalf("retained HSP counts = %v, want %v", got, want)
				}
			})
		}
	}
}

// Search orders full-precision scores, while TSV identity has three decimals.
// These real score combinations tie before formatting but reverse after it.
func TestMergeSearchResultsTopNRoundedIdentity(t *testing.T) {
	for _, fullWinner := range []bool{false, true} {
		for _, split := range []bool{false, true} {
			t.Run(fmt.Sprintf("fullWinner=%t/split=%t", fullWinner, split), func(t *testing.T) {
				dir := t.TempDir()
				row := func(genome string, aligned, matched, bitscore int, full bool) string {
					c := &AlignmentResult{AlignedLength: aligned,
						AlignedFraction: float64(aligned), PIdent: float64(matched) / float64(aligned) * 100,
						BitScore: bitscore, QEnd: aligned - 1, TEnd: aligned - 1}
					sd := &SimilarityDetail{SeqID: []byte("s"), SeqLen: 100}
					text := string(appendSearchResultRow(nil, []byte("q"), 100, 2, []byte(genome), "",
						&SearchResult{AlignedFraction: float64(aligned)}, sd, c, 1, 1, searchOutputOptions{}, 0))
					if full {
						text += fmt.Sprintf("\t%dM%dX\t%s\t%s\t%s", matched, aligned-matched,
							strings.Repeat("A", aligned), strings.Repeat("A", matched)+strings.Repeat("T", aligned-matched),
							strings.Repeat("|", matched)+strings.Repeat(" ", aligned-matched))
					}
					return text + "\n"
				}
				// Both exact scores are 7800; exported scores are 7799.975 and 7799.988.
				first := row("a", 95, 78, 95, !fullWinner)
				winner := row("b", 56, 52, 84, fullWinner)
				files := []string{filepath.Join(dir, "a.tsv")}
				body := auditSearchHeader + "\n" + first
				if split {
					// Hide the winner behind a lower-scoring head in the other input too.
					files = append(files, writeTestInput(t, filepath.Join(dir, "b.tsv"), auditSearchHeader+"\n"+row("c", 95, 78, 95, false)+winner))
				} else {
					body += winner
				}
				writeTestInput(t, files[0], body)
				out := runSeedIndexCommand(t, append([]string{"utils", "merge-search-results", "--quiet", "-n", "1"}, files...)...)
				lines := strings.Split(strings.TrimSuffix(string(out), "\n"), "\n")
				if len(lines) != 2 {
					t.Fatalf("expected only the higher exported score, got %d records", len(lines)-1)
				}
				fields := strings.Split(lines[1], "\t")
				if fields[3] != "b" || fields[2] != "1" {
					t.Fatalf("wrong Top-1 result: genome=%s hits=%s", fields[3], fields[2])
				}
				width := 20
				if fullWinner {
					width = 24
				}
				if len(fields) != width || len(strings.Split(lines[0], "\t")) != width {
					t.Fatalf("discarded result determined output width: %d, want %d", len(fields), width)
				}
			})
		}
	}
}

// Decimal TSV scores can tie even when float multiplication differs by one ULP.
func TestMergeSearchResultsDecimalCutoffTies(t *testing.T) {
	for _, identity := range []string{"75.900", "75.9000", "7.59e1", "75.9000001"} {
		for _, split := range []bool{false, true} {
			t.Run(identity+"/"+fmt.Sprint(split), func(t *testing.T) {
				dir := t.TempDir()
				row := func(genome, pident, bitscore string) string {
					return fmt.Sprintf("q\t2000\t2\t%s\ts\t100\t1\t1\t100\t2000\t%s\t0\t1\t2000\t1\t2000\t+\t2000\t0\t%s\n", genome, pident, bitscore)
				}
				first, second := row("a", "76.175", "1104"), row("b", identity, "1108")
				files := []string{filepath.Join(dir, "a.tsv")}
				body := auditSearchHeader + "\n" + first
				if split {
					files = append(files, writeTestInput(t, filepath.Join(dir, "b.tsv"), auditSearchHeader+"\n"+second))
				} else {
					body += second
				}
				writeTestInput(t, files[0], body)
				out := runSeedIndexCommand(t, append([]string{"utils", "merge-search-results", "--quiet", "-n", "1"}, files...)...)
				lines := strings.Split(strings.TrimSuffix(string(out), "\n"), "\n")
				want := 2
				if identity == "75.9000001" {
					want = 1
				}
				if len(lines) != want+1 {
					t.Fatalf("identity=%s: got %d records, want %d", identity, len(lines)-1, want)
				}
				for _, line := range lines[1:] {
					fields := strings.Split(line, "\t")
					if fields[2] != strconv.Itoa(want) || want == 1 && fields[3] != "b" {
						t.Fatalf("wrong retained hit count: %s", line)
					}
				}
			})
		}
	}
}

// A full stable sort is the reference for every arrival order and cutoff tie.
func TestTopSearchResultsMatchesFullSort(t *testing.T) {
	rng := rand.New(rand.NewSource(42))
	for _, n := range []int{1, 2, 5, 20, 100} {
		for repeat := 0; repeat < 20; repeat++ {
			scores := make([]int, 80)
			want := make([]int, len(scores))
			for i := range scores {
				scores[i], want[i] = rng.Intn(12), i
			}
			slices.SortStableFunc(want, func(i, j int) int { return scores[j] - scores[i] })
			cut := min(n, len(want))
			for cut < len(want) && scores[want[cut]] == scores[want[cut-1]] {
				cut++
			}
			want = want[:cut]
			var entries []*SearchResultOfAGenome
			selected := topSearchResults{SearchResultsHeap: SearchResultsHeap{entries: &entries}, limit: n}
			for i, score := range scores {
				selected.Add(&SearchResultOfAGenome{Score: float64(score), order: i, Sgenome: strconv.Itoa(i)})
			}
			got := make([]int, 0, len(want))
			for _, r := range selected.Results() {
				id, err := strconv.Atoi(r.Sgenome)
				if err != nil {
					t.Fatal(err)
				}
				got = append(got, id)
				RecycleSearchResultOfAGenome(r)
			}
			if !slices.Equal(got, want) {
				t.Fatalf("N=%d repeat=%d: got %v, want %v", n, repeat, got, want)
			}
		}
	}
}

// TestMergeSearchResultsEmptyInputs covers both copying and parsed merging.
func TestMergeSearchResultsEmptyInputs(t *testing.T) {
	dir := t.TempDir()
	header := writeTestInput(t, filepath.Join(dir, "header.tsv"), auditSearchHeader+"\n")
	empty := writeTestInput(t, filepath.Join(dir, "empty.tsv"), "")
	var compressed bytes.Buffer
	gz := gzip.NewWriter(&compressed)
	if err := gz.Close(); err != nil {
		t.Fatal(err)
	}
	emptyGzip := writeTestInput(t, filepath.Join(dir, "empty.tsv.gz"), compressed.String())
	for _, file := range []string{empty, emptyGzip, "-"} {
		for _, mode := range []struct {
			name string
			args []string
		}{
			{name: "copy", args: []string{file}},
			{name: "merge", args: []string{header, file}},
			{name: "query", args: []string{"-q", "q", file}},
			{name: "top-n", args: []string{"-n", "1", file}},
		} {
			t.Run(filepath.Base(file)+"/"+mode.name, func(t *testing.T) {
				args := append([]string{"utils", "merge-search-results"}, mode.args...)
				_, stderr, err := runAuditCommand(t, "", args...)
				want := fmt.Sprintf("empty search result file %q", file)
				if err != nil || !strings.Contains(string(stderr), want) || !strings.Contains(string(stderr), "[WARN]") {
					t.Fatalf("expected warning naming skipped empty input: %v\n%s", err, stderr)
				}
			})
		}
	}
	// Report every empty input while preserving valid hits and file indexing.
	row := "q\t100\t1\tg\ts\t100\t1\t1\t100\t100\t100\t0\t1\t100\t1\t100\t+\t100\t0\t100\n"
	valid := writeTestInput(t, filepath.Join(dir, "valid.tsv"), auditSearchHeader+"\n"+row)
	for _, flags := range [][]string{nil, {"-q", "q"}, {"-n", "1"}} {
		args := append([]string{"utils", "merge-search-results"}, flags...)
		args = append(args, empty, valid, emptyGzip, "-")
		out, stderr, err := runAuditCommand(t, "", args...)
		if err != nil || string(out) != auditSearchHeader+"\n"+row {
			t.Fatalf("merge with empty inputs changed valid results: %v\n%s\n%s", err, out, stderr)
		}
		for _, file := range []string{empty, emptyGzip, "-"} {
			if strings.Count(string(stderr), fmt.Sprintf("empty search result file %q", file)) != 1 {
				t.Fatalf("expected one warning for %q:\n%s", file, stderr)
			}
		}
	}
	out, stderr, err := runAuditCommand(t, "", "utils", "merge-search-results", empty, emptyGzip)
	if err != nil || len(out) != 0 || strings.Count(string(stderr), "[WARN]") != 2 {
		t.Fatalf("all-empty inputs did not produce two warnings and no results: %v\n%s\n%s", err, out, stderr)
	}
	_, stderr, err = runAuditCommand(t, "", "utils", "merge-search-results", "--quiet", empty, emptyGzip)
	if err != nil || len(stderr) != 0 {
		t.Fatalf("quiet mode did not suppress warnings: %v\n%s", err, stderr)
	}
	// A completed search without matches still has a header and remains valid.
	for _, files := range [][]string{{header}, {header, header}} {
		_, stderr, err := runAuditCommand(t, "", append([]string{"utils", "merge-search-results", "--quiet"}, files...)...)
		if err != nil {
			t.Fatalf("header-only input rejected: %v\n%s", err, stderr)
		}
	}
	// A truncated compressed stream is a read error, not a clean empty input.
	broken := writeTestInput(t, filepath.Join(dir, "broken.tsv.gz"), compressed.String()[:5])
	_, stderr, err = runAuditCommand(t, "", "utils", "merge-search-results", "--quiet", broken)
	if err == nil || !strings.Contains(string(stderr), broken) ||
		!strings.Contains(string(stderr), "unexpected EOF") || strings.Contains(string(stderr), "empty search result") {
		t.Fatalf("compressed input error lost or misclassified: %v\n%s", err, stderr)
	}
}

func TestMergeSearchResultsLongAlignmentRows(t *testing.T) {
	dir := t.TempDir()
	header := auditSearchHeader + "\tcigar\tqseq\tsseq\talign\n"
	sequence := strings.Repeat("ACGT", 1<<18) // each row exceeds 3 MiB
	rows, files := make([]string, 2), make([]string, 2)
	for i := range files {
		rows[i] = fmt.Sprintf("q\t1048576\t1\tg%d\ts\t100\t1\t1\t100\t1048576\t100\t0\t1\t1048576\t1\t1048576\t+\t1048576\t0\t%d\t1048576M\t%s\t%s\t%s\n",
			i, 100-i, sequence, sequence, strings.Repeat("|", len(sequence)))
		files[i] = writeTestInput(t, filepath.Join(dir, fmt.Sprintf("%d.tsv", i)), header+rows[i])
	}
	args := []string{"utils", "merge-search-results", "--quiet"}
	got := runSeedIndexCommand(t, append(args, files...)...)
	want := header
	for _, row := range rows {
		want += strings.Replace(row, "\t1\tg", "\t2\tg", 1)
	}
	if string(got) != want {
		t.Fatal("long alignment rows changed")
	}
	_, stderr, err := runAuditCommand(t, "", append(append(args, "-b", "1M"), files...)...)
	if err == nil || !strings.Contains(string(stderr), "token too long") {
		t.Fatalf("maximum row length was not enforced: %v %s", err, stderr)
	}
}

func TestMergeSearchResultsRejectsMixedQueries(t *testing.T) {
	for _, sameSubject := range []bool{false, true} {
		dir := t.TempDir()
		subject2 := "other"
		if sameSubject {
			subject2 = "subject"
		}
		row := func(query, subject string) string {
			return fmt.Sprintf("%s\t100\t1\t%s\ts\t100\t1\t1\t100\t100\t100\t0\t1\t100\t1\t100\t+\t100\t0\t100\n", query, subject)
		}
		file := writeTestInput(t, filepath.Join(dir, "mixed.tsv"), auditSearchHeader+"\n"+row("query1", "subject")+row("query2", subject2))
		for _, topN := range []string{"0", "1"} {
			_, stderr, err := runAuditCommand(t, "", "utils", "merge-search-results", "--quiet", "-n", topN, file)
			if err == nil || !strings.Contains(string(stderr), "queries") {
				t.Errorf("mixed queries accepted with -n %s: %v %s", topN, err, stderr)
			}
		}
		for _, query := range []string{"query1", "query2"} {
			output := runSeedIndexCommand(t, "utils", "merge-search-results", "--quiet", "-n", "1", "-q", query, file)
			lines := strings.Split(strings.TrimSpace(string(output)), "\n")
			if len(lines) != 2 || !strings.HasPrefix(lines[1], query+"\t") {
				t.Fatalf("query filter changed: %s", output)
			}
		}
	}
}

// A reused query ID does not make different query lengths interchangeable.
// Validate even discarded Top-N rows, but leave unselected queries alone.
func TestMergeSearchResultsQueryLengths(t *testing.T) {
	dir := t.TempDir()
	row := func(query, length, subject string, score int) string {
		return fmt.Sprintf("%s\t%s\t1\t%s\ts\t100\t1\t1\t100\t4\t100\t0\t1\t4\t1\t4\t+\t100\t0\t%d\n", query, length, subject, score)
	}
	for _, layout := range []string{"same subject", "different subjects", "different files"} {
		secondSubject := "a"
		if layout != "same subject" {
			secondSubject = "b"
		}
		first, second := row("q", "4", "a", 100), row("q", "9", secondSubject, 90)
		files := []string{filepath.Join(dir, "first.tsv")}
		body := auditSearchHeader + "\n" + first
		if layout == "different files" {
			files = append(files, writeTestInput(t, filepath.Join(dir, "second.tsv"), auditSearchHeader+"\n"+second))
		} else {
			body += second
		}
		writeTestInput(t, files[0], body)
		for _, topN := range []string{"0", "1"} {
			for _, query := range []string{"", "q"} {
				_, stderr, err := runAuditCommand(t, "", append([]string{"utils", "merge-search-results", "--quiet", "-n", topN, "-q", query}, files...)...)
				if err == nil || !strings.Contains(string(stderr), "inconsistent query lengths") {
					t.Errorf("%s topN=%s query=%q: got %v %s", layout, topN, query, err, stderr)
				}
			}
		}
	}
	for _, length := range []string{"0", "bad", "-1", "18446744073709551616"} {
		file := writeTestInput(t, filepath.Join(dir, "invalid.tsv"), auditSearchHeader+"\n"+row("q", length, "a", 100))
		_, stderr, err := runAuditCommand(t, "", "utils", "merge-search-results", "--quiet", file)
		if err == nil || !strings.Contains(string(stderr), "invalid query length") {
			t.Errorf("length=%q: got %v %s", length, err, stderr)
		}
	}
	file := writeTestInput(t, filepath.Join(dir, "selected.tsv"), auditSearchHeader+"\n"+row("other", "bad", "x", 101)+row("q", "4", "a", 100)+row("q", "004", "b", 90))
	output := runSeedIndexCommand(t, "utils", "merge-search-results", "--quiet", "-q", "q", file)
	if len(strings.Split(strings.TrimSpace(string(output)), "\n")) != 3 {
		t.Fatal("selected query with numerically equal lengths lost records")
	}
}

// Invalid numeric fields must be rejected before grouping or Top-N filtering.
func TestMergeSearchResultsNumericFields(t *testing.T) {
	dir := t.TempDir()
	fields := strings.Split("q\t100\t1\ta\ts\t100\t1\t1\t100\t100\t100\t0\t1\t100\t1\t100\t+\t100\t0\t100", "\t")
	valid := strings.Join(fields, "\t") + "\n"
	file := filepath.Join(dir, "input.tsv")
	for _, field := range []struct {
		name   string
		column int
		values []string
	}{
		{"pident", 10, []string{"", "bad", "NaN", "Inf", "-Inf", "1e309", "-1", "100.1"}},
		{"hits", 2, []string{"", "bad", "0", "-1", "1.5", "18446744073709551616"}},
		{"bitscore", 19, []string{"", "bad", "1.5", "18446744073709551616"}},
	} {
		for _, value := range field.values {
			for _, sameGenome := range []bool{false, true} {
				for _, topN := range []string{"0", "1"} {
					t.Run(fmt.Sprintf("%s=%s/sameGenome=%t/topN=%s", field.name, value, sameGenome, topN), func(t *testing.T) {
						bad := slices.Clone(fields)
						if !sameGenome {
							bad[3] = "b"
						}
						bad[field.column] = value
						writeTestInput(t, file, auditSearchHeader+"\n\n"+valid+strings.Join(bad, "\t")+"\n")
						_, stderr, err := runAuditCommand(t, "", "utils", "merge-search-results", "--quiet", "-n", topN, file)
						if err == nil || !strings.Contains(string(stderr), file+":4") || !strings.Contains(string(stderr), field.name) || strings.Contains(string(stderr), "panic:") {
							t.Fatalf("invalid %s was not diagnosed at its source: %v\n%s", field.name, err, stderr)
						}
					})
				}
			}
		}
	}
	// Keep valid boundary values and numeric spellings unchanged in output.
	for _, identity := range []string{"0", "100", "1e2", "75.9000001"} {
		row := slices.Clone(fields)
		row[2], row[10], row[19] = "01", identity, "-10"
		text := strings.Join(row, "\t") + "\n"
		writeTestInput(t, file, auditSearchHeader+"\n"+text)
		got := runSeedIndexCommand(t, "utils", "merge-search-results", "--quiet", "-n", "1", file)
		row[2] = "1" // The merger reports the retained genome count.
		if string(got) != auditSearchHeader+"\n"+strings.Join(row, "\t")+"\n" {
			t.Fatalf("valid numeric values changed for identity %s", identity)
		}
	}
	// An unselected query's numeric fields remain outside the merge operation.
	other := slices.Clone(fields)
	other[0], other[2], other[10], other[19] = "other", "bad", "NaN", "bad"
	writeTestInput(t, file, auditSearchHeader+"\n"+strings.Join(other, "\t")+"\n"+valid)
	got := runSeedIndexCommand(t, "utils", "merge-search-results", "--quiet", "-q", "q", file)
	if string(got) != auditSearchHeader+"\n"+valid {
		t.Fatal("numeric validation changed query selection")
	}
}

// Valid single inputs retain their fields, ordering, and empty-result header.
func TestMergeSearchResultsSingleInput(t *testing.T) {
	for _, extra := range []bool{false, true} {
		for _, empty := range []bool{false, true} {
			t.Run(fmt.Sprintf("extra=%t/empty=%t", extra, empty), func(t *testing.T) {
				header, fields := auditSearchHeader, ""
				if extra {
					header += "\tcigar\tqseq\tsseq\talign"
					fields = "\t4M\tACGT\tACGT\t||||"
				}
				data := header + "\n"
				if !empty {
					for _, subject := range []string{"a", "b"} {
						data += "q\t100\t2\t" + subject + "\ts\t4\t1\t1\t4\t4\t100\t0\t1\t4\t1\t4\t+\t100\t1e-10\t100" + fields + "\n"
					}
				}
				file := writeTestInput(t, filepath.Join(t.TempDir(), "input.tsv"), data)
				for _, stdin := range []bool{false, true} {
					input, stream := file, ""
					if stdin {
						input, stream = "-", data
					}
					out, stderr, err := runAuditCommand(t, stream, "utils", "merge-search-results", "--quiet", input)
					if err != nil || string(out) != data {
						t.Fatalf("single input changed (stdin=%t): %v\n%s\n%s", stdin, err, out, stderr)
					}
				}
			})
		}
	}
}
