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
	"os/exec"
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

// A symlink can be removed with its directory even when its target is outside.
func TestIndexPreservesInputSymlinkInOutput(t *testing.T) {
	for _, suffix := range []string{"", ExtTmpDir} {
		for _, parentAlias := range []bool{false, true} {
			t.Run(fmt.Sprintf("suffix=%s/alias=%t", suffix, parentAlias), func(t *testing.T) {
				root := t.TempDir()
				db := filepath.Join(root, "db")
				replaced := db + suffix
				if err := os.Mkdir(replaced, 0700); err != nil {
					t.Fatal(err)
				}
				input := writeTestInput(t, filepath.Join(root, "ref.fa"), ">s\n"+indexTestDNA(100)+"\n")
				link := filepath.Join(replaced, "input.fa")
				if err := os.Symlink(input, link); err != nil {
					t.Fatal(err)
				}
				passed := link
				if parentAlias {
					alias := filepath.Join(root, "alias")
					if err := os.Symlink(replaced, alias); err != nil {
						t.Fatal(err)
					}
					passed = filepath.Join(alias, "input.fa")
				}
				sentinel := writeTestInput(t, filepath.Join(replaced, "keep"), "original\n")
				logFile := writeTestInput(t, filepath.Join(root, "run.log"), "original\n")
				_, stderr, err := runAuditCommand(t, "", "index", "--quiet", "--force", "-j", "1", "-m", "64",
					"-c", "2", "--partitions", "4", "-b", "1", "-O", db, "--log", logFile, passed, input)
				if err == nil || !bytes.Contains(stderr, []byte("conflicts with input")) {
					t.Errorf("input symlink conflict accepted: %v, %s", err, stderr)
				}
				if target, err := os.Readlink(link); err != nil || target != input {
					t.Errorf("input symlink changed: %q, %v", target, err)
				}
				for _, file := range []string{sentinel, logFile} {
					if got, err := os.ReadFile(file); err != nil || string(got) != "original\n" {
						t.Errorf("file changed: %s: %q, %v", file, got, err)
					}
				}
				if got, err := os.ReadFile(input); err != nil || string(got) != ">s\n"+indexTestDNA(100)+"\n" {
					t.Errorf("source changed: %v", err)
				}
			})
		}
	}
}

// Nested aliases must not hide an input or auxiliary output entry from cleanup.
func TestIndexPreservesNestedSymlinks(t *testing.T) {
	for _, suffix := range []string{"", ExtTmpDir} {
		for _, mode := range []string{"input", "link-chain", "relative-link-chain", "log", "big-genomes"} {
			t.Run(suffix+"/"+mode, func(t *testing.T) {
				root := t.TempDir()
				real, outside := filepath.Join(root, "real"), filepath.Join(root, "outside")
				replaced := filepath.Join(real, "db") + suffix
				for _, dir := range []string{replaced, outside} {
					if err := os.MkdirAll(dir, 0700); err != nil {
						t.Fatal(err)
					}
				}
				alias := filepath.Join(root, "alias")
				link := filepath.Join(replaced, "inputs")
				for name, target := range map[string]string{alias: real, link: outside} {
					if err := os.Symlink(target, name); err != nil {
						t.Fatal(err)
					}
				}
				data := ">s\n" + indexTestDNA(100) + "\n"
				input := writeTestInput(t, filepath.Join(outside, "q.fa"), data)
				sentinel := writeTestInput(t, filepath.Join(replaced, "keep"), "original\n")
				logFile := writeTestInput(t, filepath.Join(root, "run.log"), "original\n")
				passed := filepath.Join(alias, "db"+suffix, "inputs", "q.fa")
				if mode == "link-chain" || mode == "relative-link-chain" {
					entry := filepath.Join(root, "entry")
					target := filepath.Join(alias, "db"+suffix, "inputs")
					if mode == "relative-link-chain" {
						target = filepath.Join("alias", "db"+suffix, "inputs")
					}
					if err := os.Symlink(target, entry); err != nil {
						t.Fatal(err)
					}
					passed = filepath.Join(entry, "q.fa")
				}
				args := []string{"index", "--quiet", "--force", "-j", "1", "-m", "64", "-c", "2", "--partitions", "4", "-O", filepath.Join(alias, "db")}
				wantError := "conflicts with input"
				if mode == "log" || mode == "big-genomes" {
					// Use a different output spelling so lexical containment is insufficient.
					args[len(args)-1] = filepath.Join(real, "db")
					// An absent auxiliary file must also retain every parent link entry.
					args = append(args, "--"+mode, filepath.Join(filepath.Dir(passed), "new.log"), input)
					wantError = "must be outside output directory"
				} else {
					args = append(args, "--log", logFile, passed)
				}
				_, stderr, err := runAuditCommand(t, "", args...)
				if err == nil || !bytes.Contains(stderr, []byte(wantError)) {
					t.Errorf("conflict accepted: %v, %s", err, stderr)
				}
				for file, want := range map[string]string{sentinel: "original\n", logFile: "original\n", input: data} {
					if got, err := os.ReadFile(file); err != nil || string(got) != want {
						t.Errorf("file changed: %s: %v", file, err)
					}
				}
				if target, err := os.Readlink(link); err != nil || target != outside {
					t.Errorf("input entry changed: %q, %v", target, err)
				}
				if _, err := os.Stat(filepath.Join(outside, "new.log")); !os.IsNotExist(err) {
					t.Errorf("auxiliary file created before validation: %v", err)
				}
			})
		}
	}
}

func TestCommandOutputCurrentDirectory(t *testing.T) {
	inputs := t.TempDir()
	writeTestInput(t, filepath.Join(inputs, "ref.fa"), ">ref\n"+strings.Repeat("ACGT", 300)+"\n")
	db := writeSeedPosFixture(t, 1)

	for _, commandKind := range []string{"index", "plot", "remerge"} {
		for _, pathKind := range []string{"dot", "dot-slash", "absolute", "normalized-relative", "symlink"} {
			for _, force := range []bool{false, true} {
				if commandKind == "remerge" && force {
					continue // remerge has no --force flag.
				}
				forceName := "without-force"
				if force {
					forceName = "with-force"
				}
				t.Run(commandKind+"/"+pathKind+"/"+forceName, func(t *testing.T) {
					root := t.TempDir()
					cwd := filepath.Join(root, "work")
					if err := os.Mkdir(cwd, 0700); err != nil {
						t.Fatal(err)
					}
					sentinel := writeTestInput(t, filepath.Join(cwd, "sentinel"), "preserve")
					var outDir string
					switch pathKind {
					case "dot":
						outDir = "."
					case "dot-slash":
						outDir = "." + string(filepath.Separator)
					case "absolute":
						outDir = cwd
					case "normalized-relative":
						outDir = "nested" + string(filepath.Separator) + ".."
					case "symlink":
						outDir = filepath.Join(root, "cwd-alias")
						if err := os.Symlink(cwd, outDir); err != nil {
							t.Fatal(err)
						}
					}

					// Keep the input outside cwd so input-alias checks cannot mask output validation.
					var args []string
					if commandKind == "index" {
						args = []string{"index", "--quiet", "-j", "1", "-m", "64", "-k", "16", "--partitions", "4", "-I", inputs, "-O", outDir}
					} else if commandKind == "remerge" {
						args = []string{"utils", "remerge", "--quiet", "-d", outDir}
					} else {
						args = []string{"utils", "seed-pos", "--quiet", "-d", db, "-a", "--plot-dir", outDir, "--plot-ext", ".svg"}
					}
					// Directory validation must run before an existing log can be truncated.
					args = append(args, "--log", sentinel)
					if force {
						args = append(args, "--force")
					}
					command := exec.Command(os.Args[0], append([]string{"-test.run=^TestSeedIndexCommand$", "--"}, args...)...)
					command.Dir = cwd
					command.Env = append(os.Environ(), "LEXICMAP_SEED_INDEX_TEST_HELPER=1")
					command.Stdin = strings.NewReader("")
					var stdout, stderr bytes.Buffer
					command.Stdout, command.Stderr = &stdout, &stderr
					err := command.Run()
					if err == nil || !bytes.Contains(stderr.Bytes(), []byte("should not be current directory")) {
						t.Fatalf("current directory accepted: %v\n%s\n%s", err, stderr.Bytes(), stdout.Bytes())
					}
					if got, err := os.ReadFile(sentinel); err != nil || string(got) != "preserve" {
						t.Fatalf("sentinel changed: %q, %v", got, err)
					}
					entries, err := os.ReadDir(cwd)
					if err != nil || len(entries) != 1 || entries[0].Name() != "sentinel" {
						t.Fatalf("unexpected output in cwd: %v, %v", entries, err)
					}
				})
			}
		}
	}
}

func TestCommandOutputDirectoryExistingBehavior(t *testing.T) {
	root := t.TempDir()
	input := writeTestInput(t, filepath.Join(root, "ref.fa"), ">ref\n"+strings.Repeat("ACGT", 300)+"\n")
	baseArgs := []string{"index", "--quiet", "-j", "1", "-m", "64", "-k", "16", "--partitions", "4"}

	empty := filepath.Join(root, "empty")
	if err := os.Mkdir(empty, 0700); err != nil {
		t.Fatal(err)
	}
	args := append(append([]string{}, baseArgs...), "-O", empty, input)
	if _, stderr, err := runAuditCommand(t, "", args...); err != nil {
		t.Fatalf("empty output directory rejected: %v, %s", err, stderr)
	}
	if _, err := os.Stat(filepath.Join(empty, FileInfo)); err != nil {
		t.Fatalf("index was not written to empty output directory: %v", err)
	}

	nonempty := filepath.Join(root, "nonempty")
	if err := os.Mkdir(nonempty, 0700); err != nil {
		t.Fatal(err)
	}
	sentinel := writeTestInput(t, filepath.Join(nonempty, "sentinel"), "preserve without force")
	args = append(append([]string{}, baseArgs...), "-O", nonempty, input)
	if _, stderr, err := runAuditCommand(t, "", args...); err == nil || !bytes.Contains(stderr, []byte("not empty")) {
		t.Fatalf("nonempty output directory accepted without --force: %v, %s", err, stderr)
	}
	if got, err := os.ReadFile(sentinel); err != nil || string(got) != "preserve without force" {
		t.Fatalf("nonempty output directory changed without --force: %q, %v", got, err)
	}

	args = append(append([]string{}, baseArgs...), "-O", nonempty, input, "--force")
	if _, stderr, err := runAuditCommand(t, "", args...); err != nil {
		t.Fatalf("nonempty output directory rejected with --force: %v, %s", err, stderr)
	}
	if _, err := os.Stat(sentinel); !os.IsNotExist(err) {
		t.Fatalf("--force did not replace old output directory: %v", err)
	}
	if _, err := os.Stat(filepath.Join(nonempty, FileInfo)); err != nil {
		t.Fatalf("index was not written to forced output directory: %v", err)
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

func TestRemergeProtectsTemporaryBatches(t *testing.T) {
	for _, name := range []string{FileMasks, filepath.Join(DirGenomes, batchDir(0), FileGenomes), filepath.Join(DirSeeds, "chunk_000.bin")} {
		for _, alias := range []string{"same", "symlink", "hardlink"} {
			t.Run(name+"/"+alias, func(t *testing.T) {
				dir := t.TempDir()
				db := filepath.Join(dir, "db")
				input := filepath.Join(db+ExtTmpDir, "batch_0000", name)
				if err := os.MkdirAll(filepath.Dir(input), 0700); err != nil {
					t.Fatal(err)
				}
				original := []byte("recovery input must survive\n")
				if err := os.WriteFile(input, original, 0600); err != nil {
					t.Fatal(err)
				}
				output := input
				if alias != "same" {
					output = filepath.Join(dir, "run.log")
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
				_, stderr, err := runAuditCommand(t, "", "utils", "remerge", "--quiet", "-d", db, "--log", output)
				if err == nil || !bytes.Contains(stderr, []byte("conflicts with input")) {
					t.Fatalf("recovery alias accepted: %v %s", err, stderr)
				}
				got, err := os.ReadFile(input)
				if err != nil || !bytes.Equal(got, original) {
					t.Fatalf("recovery input changed: %v %q", err, got)
				}
			})
		}
	}
}

func TestSeedPosPlotDirectoryPreservesIndex(t *testing.T) {
	for _, location := range []string{"same", "parent", "symlink", "parent-symlink"} {
		t.Run(location, func(t *testing.T) {
			dir := t.TempDir()
			parent := filepath.Join(dir, "parent")
			db := filepath.Join(parent, "db")
			if err := os.MkdirAll(db, 0700); err != nil {
				t.Fatal(err)
			}
			original := []byte("input index must survive\n")
			input := writeTestInput(t, filepath.Join(db, FileInfo), string(original))
			output := db
			if location == "parent" {
				output = parent
			}
			if strings.Contains(location, "symlink") {
				output = filepath.Join(dir, "alias")
				target := db
				if location == "parent-symlink" {
					target = parent
				}
				if err := os.Symlink(target, output); err != nil {
					t.Fatal(err)
				}
			}
			_, stderr, err := runAuditCommand(t, "", "utils", "seed-pos", "--quiet", "-d", db, "-a", "--plot-dir", output, "--force")
			if err == nil || !bytes.Contains(stderr, []byte("conflicts with input")) {
				t.Fatalf("destructive plot directory accepted: %v %s", err, stderr)
			}
			got, err := os.ReadFile(input)
			if err != nil || !bytes.Equal(got, original) {
				t.Fatalf("input index changed: %v %q", err, got)
			}
		})
	}
}

// Temporary directories are output targets too, including remerge's input batches.
func TestCommandTemporaryDirectoryCurrentDirectory(t *testing.T) {
	for _, kind := range []string{"index", "remerge"} {
		for _, alias := range []bool{false, true} {
			t.Run(fmt.Sprintf("%s/symlink=%t", kind, alias), func(t *testing.T) {
				root := t.TempDir()
				db := filepath.Join(root, "db")
				cwd := db + ExtTmpDir
				if alias {
					cwd = filepath.Join(root, "work")
				}
				if err := os.Mkdir(cwd, 0700); err != nil {
					t.Fatal(err)
				}
				if alias {
					if err := os.Symlink(cwd, db+ExtTmpDir); err != nil {
						t.Fatal(err)
					}
				}
				sentinel := writeTestInput(t, filepath.Join(cwd, "sentinel"), "preserve")
				args := []string{"utils", "remerge", "-d", db}
				if kind == "index" {
					ref := writeTestInput(t, filepath.Join(root, "ref.fa"), ">ref\n"+strings.Repeat("ACGT", 300)+"\n")
					args = []string{"index", "-O", db, "--force", "-b", "1", ref, ref}
				}
				args = append(args, "--quiet", "--log", sentinel)
				command := exec.Command(os.Args[0], append([]string{"-test.run=^TestSeedIndexCommand$", "--"}, args...)...)
				command.Dir = cwd
				command.Env = append(os.Environ(), "LEXICMAP_SEED_INDEX_TEST_HELPER=1")
				output, err := command.CombinedOutput()
				if err == nil || !bytes.Contains(output, []byte("should not be current directory")) {
					t.Fatalf("current temporary directory accepted: %v\n%s", err, output)
				}
				if got, err := os.ReadFile(sentinel); err != nil || string(got) != "preserve" {
					t.Fatalf("sentinel changed: %q, %v", got, err)
				}
				if _, err := os.Stat(db); !os.IsNotExist(err) {
					t.Fatalf("output index created before validation: %v", err)
				}
			})
		}
	}
}

func TestLogOutsideReplacedDirectories(t *testing.T) {
	for _, kind := range []string{"index", "index-tmp", "plot", "remerge", "remerge-tmp"} {
		for _, alias := range []bool{false, true} {
			for _, existingLog := range []bool{false, true} {
				for _, force := range []bool{false, true} {
					if force && strings.HasPrefix(kind, "remerge") {
						continue
					}
					t.Run(fmt.Sprintf("%s/alias=%v/existing=%v/force=%v", kind, alias, existingLog, force), func(t *testing.T) {
						root := t.TempDir()
						outDir := filepath.Join(root, "output")
						replaced := outDir
						if strings.HasSuffix(kind, "-tmp") {
							replaced += ExtTmpDir
						}
						if err := os.Mkdir(replaced, 0700); err != nil {
							t.Fatal(err)
						}
						sentinel := writeTestInput(t, filepath.Join(replaced, "sentinel"), "preserve\n")
						logFile := filepath.Join(replaced, "run.log")
						if existingLog {
							writeTestInput(t, logFile, "original log\n")
						}
						passedLog := logFile
						if alias {
							link := filepath.Join(root, "alias")
							if err := os.Symlink(replaced, link); err != nil {
								t.Fatal(err)
							}
							passedLog = filepath.Join(link, "run.log")
						}
						var args []string
						switch kind {
						case "index", "index-tmp":
							args = []string{"index", "-O", outDir}
						case "plot":
							args = []string{"utils", "seed-pos", "-d", filepath.Join(root, "db"), "--plot-dir", outDir}
						default:
							args = []string{"utils", "remerge", "-d", outDir}
						}
						args = append(args, "--log", passedLog, "--quiet")
						if force {
							args = append(args, "--force")
						}
						_, stderr, err := runAuditCommand(t, "", args...)
						if err == nil || !bytes.Contains(stderr, []byte("must be outside output directory")) {
							t.Fatalf("%v: %v, %s", args, err, stderr)
						}
						if got, err := os.ReadFile(sentinel); err != nil || string(got) != "preserve\n" {
							t.Fatalf("directory changed: %q, %v", got, err)
						}
						got, err := os.ReadFile(logFile)
						if existingLog {
							if err != nil || string(got) != "original log\n" {
								t.Fatalf("log changed: %q, %v", got, err)
							}
						} else if !os.IsNotExist(err) {
							t.Fatalf("log created: %v", err)
						}
					})
				}
			}
		}
	}
}

// Reports must survive directory replacement, including aliases in either direction.
func TestBigGenomesOutsideReplacedDirectories(t *testing.T) {
	for _, suffix := range []string{"", ExtTmpDir} {
		for _, form := range []string{"existing", "missing", "parent-alias", "leaf-symlink", "aliased-leaf-symlink"} {
			t.Run(suffix+"/"+form, func(t *testing.T) {
				root := t.TempDir()
				db := filepath.Join(root, "db")
				replaced := db + suffix
				report := filepath.Join(replaced, "skipped.tsv")
				var sentinel string
				if form != "missing" {
					if err := os.Mkdir(replaced, 0700); err != nil {
						t.Fatal(err)
					}
					sentinel = writeTestInput(t, filepath.Join(replaced, "keep"), "original\n")
				}
				passedReport := report
				if form == "parent-alias" || form == "aliased-leaf-symlink" {
					alias := filepath.Join(root, "alias")
					if err := os.Symlink(replaced, alias); err != nil {
						t.Fatal(err)
					}
					passedReport = filepath.Join(alias, "skipped.tsv")
				}
				switch form {
				case "leaf-symlink", "aliased-leaf-symlink":
					target := writeTestInput(t, filepath.Join(root, "external.tsv"), "original report\n")
					if err := os.Symlink(target, report); err != nil {
						t.Fatal(err)
					}
				case "existing":
					writeTestInput(t, report, "original report\n")
				}
				logFile := writeTestInput(t, filepath.Join(root, "run.log"), "original log\n")
				input := writeTestInput(t, filepath.Join(root, "ref.fa"), ">s\n"+indexTestDNA(100)+"\n")
				big := writeTestInput(t, filepath.Join(root, "big.fa"), ">big\n"+indexTestDNA(1020)+"\n")
				_, stderr, err := runAuditCommand(t, "", "index", "--quiet", "--force", "-j", "1", "-m", "64",
					"-c", "2", "--partitions", "4", "-b", "1", "-g", "1000", "-O", db, "-G", passedReport, "--log", logFile, big, input)
				if err == nil || !bytes.Contains(stderr, []byte("big-genomes path")) || !bytes.Contains(stderr, []byte("must be outside output directory")) {
					t.Fatalf("unsafe report accepted: %v, %s", err, stderr)
				}
				if got, err := os.ReadFile(logFile); err != nil || string(got) != "original log\n" {
					t.Fatalf("log changed: %q, %v", got, err)
				}
				if sentinel != "" {
					if got, err := os.ReadFile(sentinel); err != nil || string(got) != "original\n" {
						t.Fatalf("directory changed: %q, %v", got, err)
					}
				} else if _, err := os.Stat(replaced); !os.IsNotExist(err) {
					t.Fatalf("directory created: %v", err)
				}
				got, err := os.ReadFile(report)
				if form == "existing" || strings.HasSuffix(form, "leaf-symlink") {
					if err != nil || string(got) != "original report\n" {
						t.Fatalf("report changed: %q, %v", got, err)
					}
				} else if !os.IsNotExist(err) {
					t.Fatalf("report created: %v", err)
				}
			})
		}
	}
}
