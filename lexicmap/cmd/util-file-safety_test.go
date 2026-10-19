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
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestCommandInputOutputCollisions(t *testing.T) {
	for _, alias := range []string{"same", "symlink", "hardlink"} {
		for _, kind := range []string{"merge-output", "merge-log", "sam", "blast", "search", "genome-search", "compare", "subseq", "pair-file", "mask-file", "kv-file", "taxid-file", "genome2taxid", "file-list", "list-entry", "index-log", "index-dir", "index-metadata", "index-data"} {
			t.Run(alias+"/"+kind, func(t *testing.T) {
				dir := t.TempDir()
				input := filepath.Join(dir, "input.fa")
				data := []byte("input must survive\n")
				if err := os.WriteFile(input, data, 0600); err != nil {
					t.Fatal(err)
				}
				output := input
				if alias != "same" {
					output = filepath.Join(dir, "alias")
					var err error
					if alias == "symlink" {
						err = os.Symlink(input, output)
					} else {
						err = os.Link(input, output)
					}
					if err != nil {
						t.Fatal(err)
					}
				}
				var args []string
				switch kind {
				case "merge-output":
					args = []string{"utils", "merge-search-results", input, "-o", output}
				case "merge-log":
					args = []string{"utils", "merge-search-results", input, "--log", output}
				case "sam":
					args = []string{"utils", "2sam", input, "-o", output}
				case "blast":
					args = []string{"utils", "2blast", input, "-o", output}
				case "search":
					args = []string{"search", "-d", filepath.Join(dir, "missing"), input, "-o", output}
				case "genome-search":
					args = []string{"genome", "search", "-d", filepath.Join(dir, "missing"), input, "--log", output}
				case "compare":
					args = []string{"genome", "compare", input, "--log", output}
				case "subseq":
					args = []string{"utils", "subseq", "-d", dir, "-f", input, "-o", output}
				case "pair-file":
					args = []string{"genome", "compare", "-d", dir, "-P", input, "--log", output}
				case "mask-file":
					args = []string{"index", "-M", input, "-O", filepath.Join(dir, "index"), "--log", output}
				case "kv-file":
					args = []string{"utils", "2blast", "--kv-file-genome", input, "--log", output}
				case "taxid-file":
					args = []string{"search", "-d", dir, "--taxid-file", input, "--log", output}
				case "genome2taxid":
					args = []string{"genome", "search", "-d", dir, "-G", input, "--log", output}
				case "file-list":
					args = []string{"utils", "merge-search-results", "-X", input, "--log", output}
				case "list-entry":
					list := filepath.Join(dir, "files.txt")
					if err := os.WriteFile(list, []byte(input+"\n"), 0600); err != nil {
						t.Fatal(err)
					}
					args = []string{"utils", "merge-search-results", "-X", list, "-o", output}
				case "index-log":
					args = []string{"index", input, "-O", filepath.Join(dir, "index"), "--log", output}
				case "index-dir":
					args = []string{"index", "-I", dir, "-O", filepath.Join(dir, "index"), "--log", output}
				case "index-metadata":
					if err := os.Link(input, filepath.Join(dir, FileGenomeIndex)); err != nil {
						t.Fatal(err)
					}
					args = []string{"utils", "genomes", "-d", dir, "-o", output}
				case "index-data":
					seeds := filepath.Join(dir, DirSeeds)
					if err := os.Mkdir(seeds, 0700); err != nil {
						t.Fatal(err)
					}
					if err := os.Link(input, filepath.Join(seeds, "chunk_000.bin")); err != nil {
						t.Fatal(err)
					}
					args = []string{"search", "-d", dir, "-o", output}
				}
				args = append(args, "--quiet")
				_, stderr, err := runAuditCommand(t, "", args...)
				if err == nil || !bytes.Contains(stderr, []byte("conflicts with input")) {
					t.Fatalf("collision not rejected: %v, %s", err, stderr)
				}
				got, err := os.ReadFile(input)
				if err != nil || !bytes.Equal(got, data) {
					t.Fatalf("input changed: %q, %v", got, err)
				}
			})
		}
	}
}

func TestInputListStdinConsumedOnce(t *testing.T) {
	dir := t.TempDir()
	input := filepath.Join(dir, "input.tsv")
	output := filepath.Join(dir, "output.tsv")
	data := []byte("one file copied unchanged\n")
	if err := os.WriteFile(input, data, 0600); err != nil {
		t.Fatal(err)
	}
	_, stderr, err := runAuditCommand(t, input+"\n", "utils", "merge-search-results", "--quiet", "-X", "-", "-o", output)
	if err != nil {
		t.Fatalf("stdin list failed: %v, %s", err, stderr)
	}
	got, err := os.ReadFile(output)
	if err != nil || !bytes.Equal(got, data) {
		t.Fatalf("stdin list was lost: %q, %v", got, err)
	}
}

func TestIndexOutputCannotContainInputs(t *testing.T) {
	dir := t.TempDir()
	index := filepath.Join(dir, "index")
	if err := os.Mkdir(index, 0700); err != nil {
		t.Fatal(err)
	}
	input := filepath.Join(index, "input.fa")
	data := []byte(">seq\nACGT\n")
	if err := os.WriteFile(input, data, 0600); err != nil {
		t.Fatal(err)
	}
	_, stderr, err := runAuditCommand(t, "", "index", "--quiet", "--force", "-O", index, input)
	if err == nil || !bytes.Contains(stderr, []byte("conflicts with input")) {
		t.Fatalf("destructive output directory accepted: %v, %s", err, stderr)
	}
	got, err := os.ReadFile(input)
	if err != nil || !bytes.Equal(got, data) {
		t.Fatalf("input removed: %q, %v", got, err)
	}
}

// The derived merge directory is destructive even when -O does not exist yet.
func TestIndexTemporaryDirectoryCannotContainInputs(t *testing.T) {
	for _, mode := range []string{"positional", "list-entry", "list-file", "directory", "mask", "symlink", "parent-symlink"} {
		t.Run(mode, func(t *testing.T) {
			dir := t.TempDir()
			index := filepath.Join(dir, "db.lmi")
			tmp := index + ExtTmpDir
			if err := os.Mkdir(tmp, 0700); err != nil {
				t.Fatal(err)
			}
			input := filepath.Join(tmp, "input.fa")
			if mode == "list-file" {
				input = filepath.Join(dir, "input.fa") // Only the list is inside the merge directory.
			}
			data := []byte(">seq\nACGTACGT\n")
			if err := os.WriteFile(input, data, 0600); err != nil {
				t.Fatal(err)
			}
			var inputs []string
			switch mode {
			case "positional":
				inputs = []string{input}
			case "list-entry", "list-file":
				list := filepath.Join(dir, "files.txt")
				if mode == "list-file" {
					list = filepath.Join(tmp, "files.txt")
				}
				if err := os.WriteFile(list, []byte(input+"\n"), 0600); err != nil {
					t.Fatal(err)
				}
				inputs = []string{"-X", list}
			case "directory":
				inputs = []string{"-I", tmp}
			case "mask":
				inputs = []string{"-M", input}
			case "symlink":
				link := filepath.Join(dir, "input.fa")
				if err := os.Symlink(input, link); err != nil {
					t.Fatal(err)
				}
				inputs = []string{link}
			case "parent-symlink":
				link := filepath.Join(t.TempDir(), "parent")
				if err := os.Symlink(dir, link); err != nil {
					t.Fatal(err)
				}
				index = filepath.Join(link, "db.lmi")
				inputs = []string{input}
			}
			logFile := filepath.Join(dir, "run.log")
			args := []string{"index", "--quiet", "-j", "1", "-m", "64", "--partitions", "4", "-O", index, "--log", logFile}
			_, stderr, err := runAuditCommand(t, "", append(args, inputs...)...)
			if err == nil || !bytes.Contains(stderr, []byte("conflicts with input")) {
				t.Fatalf("temporary directory conflict accepted: %v, %s", err, stderr)
			}
			got, err := os.ReadFile(input)
			if err != nil || !bytes.Equal(got, data) {
				t.Fatalf("input changed: %q, %v", got, err)
			}
			if _, err := os.Stat(logFile); !os.IsNotExist(err) {
				t.Fatalf("log opened before validation: %v", err)
			}
		})
	}
}

func TestSingleBatchIndexPreservesUnusedTemporaryDirectory(t *testing.T) {
	dir := t.TempDir()
	index := filepath.Join(dir, "db.lmi")
	if err := os.Mkdir(index+ExtTmpDir, 0700); err != nil {
		t.Fatal(err)
	}
	sentinel := filepath.Join(index+ExtTmpDir, "keep.txt")
	if err := os.WriteFile(sentinel, []byte("keep"), 0600); err != nil {
		t.Fatal(err)
	}
	input := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(input, []byte(">seq\n"+strings.Repeat("ACGT", 100)), 0600); err != nil {
		t.Fatal(err)
	}
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-k", "16", "--partitions", "4", "-O", index, input)
	got, err := os.ReadFile(sentinel)
	if err != nil || string(got) != "keep" {
		t.Fatalf("unused temporary directory changed: %q, %v", got, err)
	}
}

func TestIndexedGenomeIDsAreNotFilePaths(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "id")
	if err := os.WriteFile(path, []byte("not a FASTA input"), 0600); err != nil {
		t.Fatal(err)
	}
	// These are genome IDs; only the missing database should be rejected.
	_, stderr, err := runAuditCommand(t, "", "genome", "compare", "--quiet", "-d", filepath.Join(dir, "missing"), path, "other", "-o", path)
	if err == nil || bytes.Contains(stderr, []byte("conflicts with input")) {
		t.Fatalf("genome ID treated as file: %v, %s", err, stderr)
	}
}

func TestAuxiliaryOutputsCannotShareLog(t *testing.T) {
	file := filepath.Join(t.TempDir(), "shared")
	_, stderr, err := runAuditCommand(t, "", "index", "--quiet", "-O", file+".index", "--big-genomes", file, "--log", file)
	if err == nil || !strings.Contains(string(stderr), "output paths should not be the same") {
		t.Fatalf("outputs not checked: %v, %s", err, stderr)
	}
	if _, err = os.Stat(file); !os.IsNotExist(err) {
		t.Fatal(fmt.Sprintf("output opened before checking: %v", err))
	}
}
