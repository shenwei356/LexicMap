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
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"reflect"
	"slices"
	"strconv"
	"strings"
	"testing"
)

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
		_, stderr, err := runAuditCommand(t, "", "utils", "merge-search-results", "--quiet", "-n", "1", file)
		if err == nil || !strings.Contains(string(stderr), "queries") {
			t.Fatalf("mixed queries accepted: %v %s", err, stderr)
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
