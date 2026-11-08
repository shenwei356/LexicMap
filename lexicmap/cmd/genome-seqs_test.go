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

// Unknown genome names must not create or truncate extraction outputs.
func TestGenomeSeqsMissingNamePreservesOutput(t *testing.T) {
	dir := t.TempDir()
	sequence := indexTestDNA(500)
	input := writeTestInput(t, filepath.Join(dir, "g.fa"), ">s\n"+sequence+"\n")
	db := filepath.Join(dir, "db")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-c", "2", "--partitions", "4", "-O", db, input)
	for _, suffix := range []string{".fa", ".fa.gz"} {
		for _, existing := range []bool{false, true} {
			t.Run(fmt.Sprintf("%s/existing=%t", suffix, existing), func(t *testing.T) {
				output := filepath.Join(t.TempDir(), "output"+suffix)
				const original = "keep existing output\n"
				if existing {
					writeTestInput(t, output, original)
				}
				_, stderr, err := runAuditCommand(t, "", "utils", "genome-seqs", "--quiet", "-d", db, "-n", "missing", "-o", output)
				if err == nil || !bytes.Contains(stderr, []byte("reference name not found: missing")) {
					t.Fatalf("expected missing-name error: %v\n%s", err, stderr)
				}
				got, err := os.ReadFile(output)
				if existing {
					if err != nil || string(got) != original {
						t.Fatalf("missing name changed existing output: %q, %v", got, err)
					}
				} else if !os.IsNotExist(err) {
					t.Fatalf("missing name created output: %v", err)
				}
			})
		}
	}
	got := runSeedIndexCommand(t, "utils", "genome-seqs", "--quiet", "-d", db, "-n", "g", "-w", "0")
	if string(got) != ">s\n"+sequence+"\n" {
		t.Fatal("valid genome extraction changed")
	}
}

// writeChunkExtractionFixture uses a real split index shared by extraction tests.
func writeChunkExtractionFixture(t *testing.T) (string, string) {
	t.Helper()
	dir := t.TempDir()
	sequence := indexTestDNA(1000)
	var fasta strings.Builder
	for i := 0; i < 96; i++ {
		fmt.Fprintf(&fasta, ">s%d\n%s\n", i, sequence)
	}
	input := writeTestInput(t, filepath.Join(dir, "many.fa"), fasta.String())
	db := filepath.Join(dir, "db")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-c", "2", "--partitions", "4", "-g", "1500", "-O", db, input)
	return db, sequence
}
