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
	"fmt"
	"os"
	"path/filepath"
	"slices"
	"strings"
	"testing"
)

// A near-full batch pool must leave room for the process and command output.
func TestCompareResourceFileBudget(t *testing.T) {
	db := writeDetailsFixture(t, 22)
	args := []string{"genome", "compare", "--quiet", "-j", "4", "-d", db, "g0", "g1"}
	want := runSeedIndexCommand(t, args...)
	for _, budget := range []string{"31", "32"} {
		got, stderr, err := runResourceLimitedCommand(t, "open-files="+budget, append(args, "--max-open-files", budget)...)
		if err != nil || !bytes.Equal(got, want) {
			t.Fatalf("budget %s: error=%v stdout=%s stderr=%s", budget, err, got, stderr)
		}
	}
	for _, budget := range []string{"30", "24", "10"} {
		got, stderr, err := runResourceLimitedCommand(t, "open-files=32", append(args, "--max-open-files", budget)...)
		if err == nil || len(got) != 0 || !strings.Contains(string(stderr), "should be >= 31") {
			t.Fatalf("budget %s: missing pool was accepted: %v stdout=%s stderr=%s", budget, err, got, stderr)
		}
	}
}

func TestCompareRawInputFileBudget(t *testing.T) {
	dir := t.TempDir()
	files := make([]string, 24)
	for i := range files {
		files[i] = writeTestInput(t, filepath.Join(dir, fmt.Sprintf("g%d.fa", i)), ">s\n"+strings.Repeat("ACGT", 50000)+"\n")
	}
	baseline := append([]string{"genome", "compare", "--quiet", "-j", "2", "--max-genome-cache-memory", "0"}, files...)
	want := strings.Split(string(runSeedIndexCommand(t, baseline...)), "\n")
	slices.Sort(want)
	for _, cache := range []string{"0", "64M"} {
		outfile := filepath.Join(dir, "compare-"+cache+".tsv")
		args := append([]string{"genome", "compare", "--quiet", "-j", "64", "--max-open-files", "10", "--max-genome-cache-memory", cache, "--log", outfile + ".log", "-o", outfile}, files...)
		_, stderr, err := runResourceLimitedCommand(t, "open-files=10", args...)
		if err != nil {
			t.Fatalf("cache=%s: error=%v stderr=%s", cache, err, stderr)
		}
		out, err := os.ReadFile(outfile)
		got := strings.Split(string(out), "\n")
		slices.Sort(got)
		if err != nil || !slices.Equal(got, want) {
			t.Fatalf("cache=%s: error=%v stderr=%s output differs=%v", cache, err, stderr, !slices.Equal(got, want))
		}
	}
}
