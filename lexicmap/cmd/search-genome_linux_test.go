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
	"fmt"
	"os"
	"path/filepath"
	"slices"
	"strings"
	"testing"
)

// Multiple genome query files must coexist with every persistent index pool,
// command output, and logging under the actual OS limit.
func TestGenomeSearchQueryFileBudget(t *testing.T) {
	db, ref := buildMultiBatchSearchFixture(t)
	data, err := os.ReadFile(ref)
	if err != nil {
		t.Fatal(err)
	}
	dir := t.TempDir()
	files := make([]string, 24)
	for i := range files {
		files[i] = writeTestInput(t, filepath.Join(dir, fmt.Sprintf("q%02d.fa", i)), string(data))
	}
	args := append([]string{"genome", "search", "--quiet", "-j", "32", "-J", "32", "-d", db, "-m", "0", "-S", "-p", "12"}, files...)
	normalize := func(data []byte) []string {
		lines := strings.Split(strings.TrimRight(string(data), "\r\n"), "\n")
		if len(lines) != 49 {
			t.Fatalf("expected two hits per query: %d rows", len(lines))
		}
		slices.Sort(lines[1:]) // concurrent query completion order is incidental
		return lines
	}
	want := normalize(runSeedIndexCommand(t, args...))
	for _, memory := range []bool{false, true} {
		for _, budget := range []int{15, 24} {
			if memory {
				budget -= 4
			}
			out := filepath.Join(dir, "out.tsv")
			limited := append(slices.Clone(args), "--max-open-files", fmt.Sprint(budget), "--log", filepath.Join(dir, "search.log"), "-o", out)
			if memory {
				limited = append(limited, "-w")
			}
			_, stderr, err := runResourceLimitedCommand(t, "open-files="+fmt.Sprint(budget), limited...)
			if err != nil {
				log, _ := os.ReadFile(filepath.Join(dir, "search.log"))
				t.Fatalf("budget=%d memory=%v: %v, %s\n%s", budget, memory, err, stderr, log)
			}
			data, err := os.ReadFile(out)
			if err != nil {
				t.Fatal(err)
			}
			if !slices.Equal(normalize(data), want) {
				t.Fatalf("budget=%d memory=%v: output changed", budget, memory)
			}
		}
	}
}
