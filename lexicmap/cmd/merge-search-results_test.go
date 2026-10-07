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
			return
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
