//go:build linux

package cmd

import (
	"bytes"
	"compress/gzip"
	"context"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"slices"
	"strconv"
	"strings"
	"testing"
	"time"
)

// Many input shards must retain every field and cutoff tie under a small FD
// limit, including compressed shards and mixed alignment-column formats.
func TestMergeSearchResultsManyShards(t *testing.T) {
	dir := t.TempDir()
	temp := t.TempDir()
	t.Setenv("TMPDIR", temp)
	const header = auditSearchHeader + "\tcigar\tqseq\tsseq\talign\n"
	var combined strings.Builder
	combined.WriteString(header)
	files := make([]string, 200)
	for i := range files {
		score := 1000 - i/2
		if i == 0 {
			score = 2000 // The first retained row selects 24-column output.
		}
		row := fmt.Sprintf("q\t4\t1\tg%d\ts\t100\t1\t1\t100\t4\t100\t0\t1\t4\t1\t4\t+\t4\t0\t%d", i, score)
		if i%2 == 0 {
			row += "\t4M\tACGT\tACGT\t||||"
		}
		row += "\n"
		combined.WriteString(strings.Replace(row, "\t1\tg", fmt.Sprintf("\t%d\tg", len(files)), 1))
		if i == len(files)-1 {
			row = strings.TrimSuffix(row, "\n")
		}
		data := []byte(header + row)
		file := filepath.Join(dir, fmt.Sprintf("shard-%03d.tsv", i))
		if i%3 == 0 {
			var compressed bytes.Buffer
			writer := gzip.NewWriter(&compressed)
			if _, err := writer.Write(data); err != nil {
				t.Fatal(err)
			}
			if err := writer.Close(); err != nil {
				t.Fatal(err)
			}
			data, file = compressed.Bytes(), file+".gz"
		}
		if err := os.WriteFile(file, data, 0600); err != nil {
			t.Fatal(err)
		}
		files[i] = file
	}
	oracle := writeTestInput(t, filepath.Join(dir, "combined.tsv"), combined.String())
	for _, topN := range []string{"0", "3"} {
		args := []string{"utils", "merge-search-results", "--quiet", "-n", topN}
		want := runSeedIndexCommand(t, append(args, oracle)...)
		for _, budget := range []int{reservedOpenFiles + 1, 64, reservedOpenFiles + len(files)} {
			output, logfile := filepath.Join(dir, "output.tsv"), filepath.Join(dir, "run.log")
			limited := append(slices.Clone(args), "--max-open-files", strconv.Itoa(budget), "-o", output, "--log", logfile)
			_, stderr, err := runResourceLimitedCommand(t, "open-files="+strconv.Itoa(budget), append(limited, files...)...)
			if err != nil {
				t.Fatalf("topN=%s, budget=%d: %v: %s", topN, budget, err, stderr)
			}
			got, err := os.ReadFile(output)
			if err != nil {
				t.Fatal(err)
			}
			// Equal-score genomes already have nondeterministic cross-file order.
			actual := strings.Split(strings.TrimSuffix(string(got), "\n"), "\n")
			expected := strings.Split(strings.TrimSuffix(string(want), "\n"), "\n")
			slices.Sort(actual[1:])
			slices.Sort(expected[1:])
			if !slices.Equal(actual, expected) {
				t.Fatalf("topN=%s, budget=%d: merged records differ from the single-input oracle", topN, budget)
			}
		}
	}
	// A late parse failure must remove any temporary input storage.
	bad := writeTestInput(t, filepath.Join(dir, "bad.tsv"), header+"invalid\n")
	args := append([]string{"utils", "merge-search-results", "--quiet", "--max-open-files", "64"}, files...)
	_, stderr, err := runResourceLimitedCommand(t, "open-files", append(args, bad)...)
	if err == nil || !strings.Contains(string(stderr), "columns") {
		t.Fatalf("parse error lost: %v: %s", err, stderr)
	}
	// Spool-write and final-output failures must also release the shared file.
	_, stderr, err = runResourceLimitedCommand(t, "file-size", args...)
	if err == nil || !strings.Contains(string(stderr), "file too large") {
		t.Fatalf("temporary write error lost: %v: %s", err, stderr)
	}
	_, stderr, err = runResourceLimitedCommand(t, "open-files", append(args, "-o", "/dev/full")...)
	if err == nil || !strings.Contains(string(stderr), "no space left on device") {
		t.Fatalf("output error lost: %v: %s", err, stderr)
	}
	entries, err := os.ReadDir(temp)
	if err != nil || len(entries) != 0 {
		t.Fatalf("temporary inputs retained: %v, %v", entries, err)
	}
}

// A closed stdout pipe must report EPIPE and remove staged inputs, including
// when the first write happens during a deferred buffer flush or gzip close.
func TestMergeSearchResultsBrokenPipeCleansSpool(t *testing.T) {
	for _, rows := range []int{1, 1500} {
		for _, compressed := range []bool{false, true} {
			t.Run(fmt.Sprintf("rows=%d/gzip=%t", rows, compressed), func(t *testing.T) {
				dir, temp := t.TempDir(), t.TempDir()
				args := []string{"utils", "merge-search-results", "--quiet", "--max-open-files", "10"}
				if compressed {
					args = append(args, "--force-gzip-output")
				}
				for file := range 2 {
					var input strings.Builder
					input.WriteString(auditSearchHeader + "\n")
					for row := range rows {
						fmt.Fprintf(&input, "q\t4\t%d\tg%d_%d\ts\t100\t1\t1\t100\t4\t100\t0\t1\t4\t1\t4\t+\t4\t0\t100\n", rows, file, row)
					}
					args = append(args, writeTestInput(t, filepath.Join(dir, fmt.Sprintf("%d.tsv", file)), input.String()))
				}
				read, write, err := os.Pipe()
				if err != nil {
					t.Fatal(err)
				}
				if err = read.Close(); err != nil {
					t.Fatal(err)
				}
				defer write.Close()
				ctx, cancel := context.WithTimeout(context.Background(), 30*time.Second)
				defer cancel()
				command := exec.CommandContext(ctx, os.Args[0], append([]string{"-test.run=^TestSeedIndexCommand$", "--"}, args...)...)
				command.Env = append(os.Environ(), "LEXICMAP_SEED_INDEX_TEST_HELPER=1", "TMPDIR="+temp)
				command.Stdout = write
				var stderr bytes.Buffer
				command.Stderr = &stderr
				if err := command.Run(); err == nil || !strings.Contains(stderr.String(), "broken pipe") {
					t.Fatalf("missing output error: %v: %s", err, stderr.String())
				}
				entries, err := os.ReadDir(temp)
				if err != nil || len(entries) != 0 {
					t.Fatalf("temporary inputs retained: %v, %v", entries, err)
				}
			})
		}
	}
}
