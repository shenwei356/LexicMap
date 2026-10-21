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
	"os"
	"path/filepath"
	"reflect"

	"bytes"
	"strings"
	"testing"
	"time"

	"github.com/shenwei356/go-logging"
	"github.com/spf13/cobra"
)

func TestInputListPreservesExplicitStdin(t *testing.T) {
	dir := t.TempDir()
	for _, tc := range []struct {
		name string
		args []string
		list string
		want []string
	}{
		{"implicit", nil, "listed.fa\n", []string{"listed.fa"}},
		{"explicit", []string{"-"}, "listed.fa\n", []string{"-", "listed.fa"}},
		{"mixed", []string{"first.fa", "-"}, "listed.fa\n", []string{"first.fa", "-", "listed.fa"}},
		{"listed_stdin", nil, "-\nlisted.fa\n", []string{"-", "listed.fa"}},
		{"empty_list", []string{"-"}, "", []string{"-"}},
	} {
		t.Run(tc.name, func(t *testing.T) {
			list := writeTestInput(t, filepath.Join(dir, tc.name), tc.list)
			cmd := &cobra.Command{}
			cmd.Flags().String("infile-list", list, "")
			got := getFileListFromArgsAndFile(cmd, tc.args, false, "infile-list", false)
			if !reflect.DeepEqual(got, tc.want) {
				t.Fatalf("inputs=%v, want %v", got, tc.want)
			}
		})
	}

	sequence := indexTestDNA(1500)
	input := writeTestInput(t, filepath.Join(dir, "g.fa"), ">g\n"+sequence+"\n")
	query := writeTestInput(t, filepath.Join(dir, "query.fa"), ">file_query\n"+sequence+"\n")
	list := writeTestInput(t, filepath.Join(dir, "queries.txt"), query+"\n")
	db := filepath.Join(dir, "db")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-c", "2", "--partitions", "4", "-O", db, input)
	out, stderr, err := runAuditCommand(t, ">stdin_query\n"+sequence+"\n", "search", "--quiet", "-j", "1",
		"-d", db, "-n", "0", "-N", "0", "-X", list, "-")
	if err != nil || !bytes.Contains(out, []byte("stdin_query\t")) || !bytes.Contains(out, []byte("file_query\t")) {
		t.Fatalf("missing mixed input: %v\n%s\n%s", err, stderr, out)
	}
}

// TestLogStartupInfoDate checks that the local date follows the command line.
func TestLogStartupInfoDate(t *testing.T) {
	var output bytes.Buffer
	previousLog := log
	log = logging.MustGetLogger("startup-test")
	backend := logging.NewBackendFormatter(logging.NewLogBackend(&output, "", 0),
		logging.MustStringFormatter("%{message}"))
	log.SetBackend(logging.AddModuleLevel(backend))
	t.Cleanup(func() { log = previousLog })

	before := time.Now().Format("2006-01-02")
	logStartupInfo()
	after := time.Now().Format("2006-01-02")
	lines := strings.Split(output.String(), "\n")
	for i, line := range lines {
		if !strings.HasPrefix(line, " CMD: ") {
			continue
		}
		if i+1 >= len(lines) {
			t.Fatal("missing date after CMD")
		}
		// Accept either date if logging crosses local midnight.
		if date := lines[i+1]; date != "DATE: "+before && date != "DATE: "+after {
			t.Fatalf("unexpected line after CMD: %q", date)
		}
		return
	}
	t.Fatal("missing CMD in startup log")
}

// writeTestInput writes a small input without depending on local databases.
func writeTestInput(t *testing.T, file, data string) string {
	t.Helper()
	if err := os.WriteFile(file, []byte(data), 0600); err != nil {
		t.Fatal(err)
	}
	return file
}

func TestCommandsRejectRepeatedStdin(t *testing.T) {
	dir := t.TempDir()
	list := writeTestInput(t, filepath.Join(dir, "files.txt"), "-\n")
	for _, command := range [][]string{
		{"index", "-O", filepath.Join(dir, "db")}, {"search", "-d", dir}, {"genome", "search", "-d", dir},
		{"genome", "compare"}, {"utils", "2blast"}, {"utils", "2sam"}, {"utils", "merge-search-results"},
	} {
		for _, inputs := range [][]string{{"-", "-"}, {"-X", list, "-"}, {"-X", "-", "-"}} {
			// Logging must not truncate an existing file before reporting invalid stdin.
			logFile := writeTestInput(t, filepath.Join(dir, "command.log"), "preserve\n")
			args := append(append([]string{}, command...), "--quiet", "--log", logFile)
			args = append(args, inputs...)
			out, stderr, err := runAuditCommand(t, "", args...)
			if err == nil || !bytes.Contains(stderr, []byte("stdin (-) may only be used once")) || len(out) != 0 {
				t.Fatalf("%v: %v, %s, %s", args, err, out, stderr)
			}
			if got, err := os.ReadFile(logFile); err != nil || string(got) != "preserve\n" {
				t.Fatalf("log changed: %q, %v", got, err)
			}
		}
	}
	for _, args := range [][]string{
		{"utils", "merge-search-results", "-", "-"},  // stdout-only path
		{"utils", "merge-search-results", "-X", "-"}, // stdin list itself contains stdin
		{"utils", "2blast", "--kv-file-genome", "-", "-"},
		{"genome", "compare", "-d", dir, "-P", "-", "-P", "-"},
	} {
		_, stderr, err := runAuditCommand(t, "-\n", args...)
		if err == nil || !bytes.Contains(stderr, []byte("stdin (-) may only be used once")) {
			t.Fatalf("%v: %v, %s", args, err, stderr)
		}
	}

	output := writeTestInput(t, filepath.Join(dir, "output.tsv"), "preserve\n")
	_, stderr, err := runAuditCommand(t, "", "utils", "merge-search-results", "-", "-", "-o", output)
	if err == nil || !bytes.Contains(stderr, []byte("stdin (-) may only be used once")) {
		t.Fatalf("output preflight: %v, %s", err, stderr)
	}
	if got, err := os.ReadFile(output); err != nil || string(got) != "preserve\n" {
		t.Fatalf("output changed: %q, %v", got, err)
	}
	// Even an empty stdin file list would otherwise fall back to reading stdin again.
	_, stderr, err = runAuditCommand(t, "", "utils", "merge-search-results", "-X", "-")
	if err == nil || !bytes.Contains(stderr, []byte("stdin (-) may only be used once")) {
		t.Fatalf("empty stdin list: %v, %s", err, stderr)
	}
}

func TestStdinPreflightDistinguishesGenomeIDs(t *testing.T) {
	cmd := &cobra.Command{Use: "compare"}
	cmd.Flags().String("index", filepath.Join(t.TempDir(), "db"), "")
	if err := checkCommandFiles(cmd, []string{"-", "-"}); err != nil {
		t.Fatalf("genome IDs treated as stdin: %v", err)
	}
	// File-list entries in index mode are IDs as well.
	list := writeTestInput(t, filepath.Join(t.TempDir(), "ids.txt"), "-\n-\n")
	cmd.Flags().String("infile-list", list, "")
	if err := checkCommandFiles(cmd, nil); err != nil {
		t.Fatalf("listed genome IDs treated as stdin: %v", err)
	}
}
