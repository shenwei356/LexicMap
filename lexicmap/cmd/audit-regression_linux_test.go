//go:build linux

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
	"math/rand"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

func TestCommandOutputWriteFailures(t *testing.T) {
	dir := t.TempDir()
	input := filepath.Join(dir, "input.tsv")
	if err := os.WriteFile(input, []byte("buffered result\n"), 0600); err != nil {
		t.Fatal(err)
	}
	for _, compressed := range []bool{false, true} {
		output := "/dev/full"
		if compressed {
			output = filepath.Join(dir, "full.gz")
			if err := os.Symlink("/dev/full", output); err != nil {
				t.Fatal(err)
			}
		}
		_, stderr, err := runAuditCommand(t, "", "utils", "merge-search-results", "--quiet", input, "-o", output)
		if err == nil || !bytes.Contains(stderr, []byte("no space left")) {
			t.Fatalf("lost output error: %v, %s", err, stderr)
		}
	}
	// Verify the real search output path, including the header-only case.
	rng := rand.New(rand.NewSource(91))
	sequence := make([]byte, 600)
	for i := range sequence {
		sequence[i] = "ACGT"[rng.Intn(4)]
	}
	ref := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(ref, append([]byte(">seq\n"), sequence...), 0600); err != nil {
		t.Fatal(err)
	}
	db := filepath.Join(dir, "db")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-c", "2", "-m", "64", "-k", "16", "--partitions", "4", "-O", db, ref)
	for _, query := range []string{">match\n" + string(sequence), ">short\nACGT\n"} {
		_, stderr, err := runAuditCommand(t, query, "search", "--quiet", "-j", "1", "-d", db, "-p", "12", "-P", "12", "-n", "0", "-N", "0", "-o", "/dev/full")
		if err == nil || !bytes.Contains(stderr, []byte("no space left")) {
			t.Fatalf("lost search output error: %v, %s", err, stderr)
		}
	}
}

func TestRedirectedStdinOutputCollision(t *testing.T) {
	file := filepath.Join(t.TempDir(), "input.tsv")
	data := []byte("original input\n")
	if err := os.WriteFile(file, data, 0600); err != nil {
		t.Fatal(err)
	}
	input, err := os.Open(file)
	if err != nil {
		t.Fatal(err)
	}
	defer input.Close()
	command := exec.Command(os.Args[0], "-test.run=^TestSeedIndexCommand$", "--", "utils", "merge-search-results", "--quiet", "-o", file)
	command.Env = append(os.Environ(), "LEXICMAP_SEED_INDEX_TEST_HELPER=1")
	command.Stdin = input
	out, err := command.CombinedOutput()
	if err == nil || !strings.Contains(string(out), "conflicts with stdin") {
		t.Fatalf("stdin alias accepted: %v, %s", err, out)
	}
	got, err := os.ReadFile(file)
	if err != nil || !bytes.Equal(got, data) {
		t.Fatalf("stdin overwritten: %q, %v", got, err)
	}
}
