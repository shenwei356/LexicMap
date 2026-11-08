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
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// Index consumers reject invalid summary fields before truncating result files.
// A real index also provides successful controls for every command.
func TestIndexCommandsRejectInvalidSummary(t *testing.T) {
	dir := t.TempDir()
	query := writeTestInput(t, filepath.Join(dir, "g.fa"), ">s\n"+indexTestDNA(2040)+"\n")
	db := filepath.Join(dir, "db")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-c", "2", "--partitions", "4", "--save-seed-pos", "-O", db, query)
	summary := filepath.Join(db, FileInfo)
	original, err := os.ReadFile(summary)
	if err != nil {
		t.Fatal(err)
	}
	commands := []struct {
		name           string
		args           []string
		seeds, genomes bool
	}{
		{"search", []string{"search", "-p", "15", "-P", "15", query}, true, true},
		{"search-memory", []string{"search", "-p", "15", "-P", "15", "-w", query}, true, true},
		{"genome-search", []string{"genome", "search", "-m", "0", "-p", "15", query}, true, true},
		{"compare", []string{"genome", "compare", "g", "g"}, false, true},
		{"kmers", []string{"utils", "kmers", "-m", "1"}, true, false},
		{"genome-seqs", []string{"utils", "genome-seqs", "-n", "g"}, false, true},
		{"subseq", []string{"utils", "subseq", "-n", "g", "-s", "s", "-r", "1:31"}, false, true},
		{"seed-pos", []string{"utils", "seed-pos", "-n", "g"}, false, true},
		{"genome-details", []string{"utils", "genome-details"}, false, true},
	}
	for _, command := range commands {
		t.Run(command.name, func(t *testing.T) {
			args := append(append([]string{}, command.args...), "--quiet", "-j", "1", "-d", db)
			if output := runSeedIndexCommand(t, args...); len(output) == 0 {
				t.Fatal("valid index produced no output")
			}
			cases := []struct{ old, replacement string }{{"main-version = 3", "main-version = 0"}}
			if command.seeds {
				cases = append(cases, []struct{ old, replacement string }{
					{"chunks = 2", "chunks = 0"}, {"chunks = 2", "chunks = 129"},
					{"max-K = 31", "max-K = 0"}, {"max-K = 31", "max-K = 32"},
					{"masks = 64", "masks = 63"},
				}...)
			}
			if command.genomes {
				cases = append(cases, []struct{ old, replacement string }{
					{"genome-batches = 1", "genome-batches = 0"},
					{"genome-batches = 1", "genome-batches = -1"},
					{"genome-batches = 1", "genome-batches = 131073"},
					{"\ngenomes = 1", "\ngenomes = -1"},
					{"input-bases = 2040", "input-bases = -1"},
				}...)
			}
			for _, tc := range cases {
				t.Run(tc.replacement, func(t *testing.T) {
					data := strings.Replace(string(original), tc.old, tc.replacement, 1)
					if data == string(original) {
						t.Fatalf("missing field %q", tc.old)
					}
					if err := os.WriteFile(summary, []byte(data), 0600); err != nil {
						t.Fatal(err)
					}
					t.Cleanup(func() {
						if err := os.WriteFile(summary, original, 0600); err != nil {
							t.Error(err)
						}
					})
					output := writeTestInput(t, filepath.Join(t.TempDir(), "out.txt"), "keep existing output\n")
					_, stderr, err := runAuditCommand(t, "", append(append([]string{}, args...), "-o", output)...)
					if err == nil || bytes.Contains(stderr, []byte("panic:")) {
						t.Errorf("expected summary rejection without panic: %v\n%s", err, stderr)
					}
					got, err := os.ReadFile(output)
					if err != nil || string(got) != "keep existing output\n" {
						t.Errorf("invalid summary changed output: %v (%d bytes)", err, len(got))
					}
				})
			}
		})
	}
}
