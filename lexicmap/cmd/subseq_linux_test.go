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
	"slices"
	"strings"
	"testing"
)

func TestSubseqFileBudget(t *testing.T) {
	db, query := buildMultiBatchSearchFixture(t)
	dir := t.TempDir()
	results := runSeedIndexCommand(t, "search", "--quiet", "-j", "1", "-d", db, "-p", "12", "-P", "12", "-n", "0", "-N", "0", query)
	input := writeTestInput(t, filepath.Join(dir, "hits.tsv"), string(results))
	args := []string{"utils", "subseq", "--quiet", "-j", "32", "-d", db, "-f", input, "-w", "0"}
	want := runSeedIndexCommand(t, args...)
	if got, records := bytes.Count(want, []byte("\n"))/2, bytes.Count(results, []byte("\n"))-1; got != records || got < 2 {
		t.Fatalf("expected one extracted sequence per input alignment: got %d, want %d", got, records)
	}
	for _, budget := range []int{10, 11, 24} {
		out, logFile := filepath.Join(dir, "out.fa"), filepath.Join(dir, "subseq.log")
		limited := append(slices.Clone(args), "--max-open-files", fmt.Sprint(budget), "-o", out, "--log", logFile)
		_, stderr, err := runResourceLimitedCommand(t, "open-files="+fmt.Sprint(budget), limited...)
		if err != nil {
			log, _ := os.ReadFile(logFile)
			t.Fatalf("budget=%d: %v, %s\n%s", budget, err, stderr, log)
		}
		got, err := os.ReadFile(out)
		if err != nil || !bytes.Equal(got, want) {
			t.Fatalf("budget=%d: extraction changed: %v", budget, err)
		}
	}
	// Direct region extraction also validates the process allowance.
	out, stderr, err := runResourceLimitedCommand(t, "open-files=16", "utils", "subseq", "--quiet", "-d", db,
		"-n", "ref", "-r", "1:10", "--max-open-files", "9")
	if err == nil || len(out) != 0 || !strings.Contains(string(stderr), "should be >= 10") {
		t.Fatalf("direct extraction accepted invalid budget: %v, %s, %s", err, out, stderr)
	}
}

func TestSubseqClosesUnmatchedChunkReaders(t *testing.T) {
	db, sequence := writeChunkExtractionFixture(t)
	t.Setenv("GOGC", "off")
	output, stderr, err := runResourceLimitedCommand(t, "open-files", "utils", "subseq", "--quiet", "-d", db, "-n", "many", "-s", "s95", "-r", "1:10", "-w", "0")
	if err != nil {
		t.Fatalf("last chunk extraction failed: %v %s", err, stderr)
	}
	if string(output) != ">s95:1-10:+\n"+sequence[:10]+"\n" {
		t.Fatalf("wrong sequence: %s", output)
	}
}
