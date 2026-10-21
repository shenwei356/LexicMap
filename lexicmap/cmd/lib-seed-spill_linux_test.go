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

// Multiple queries share readers and spill files under the actual process limit.
func TestSeedSpillSearchResourceFileBudget(t *testing.T) {
	db, ref := buildMultiBatchSearchFixture(t)
	data, err := os.ReadFile(ref)
	if err != nil {
		t.Fatal(err)
	}
	var queries strings.Builder
	for i := range 4 {
		queries.WriteString(strings.Replace(string(data), ">seq", fmt.Sprintf(">q%d", i), 1))
		queries.WriteByte('\n')
	}
	query := writeTestInput(t, filepath.Join(t.TempDir(), "queries.fa"), queries.String())
	args := []string{"search", "--quiet", "-j", "4", "-J", "4", "-d", db, "-p", "12", "-P", "12", "-n", "0", "-N", "0", query}
	// Only query completion order is incidental; retain all columns and duplicates.
	normalize := func(data []byte) string {
		lines := strings.Split(strings.TrimSpace(string(data)), "\n")
		if len(lines) < 5 {
			t.Fatalf("fixture lost query hits: %s", data)
		}
		slices.SortStableFunc(lines[1:], func(a, b string) int {
			queryA, _, _ := strings.Cut(a, "\t")
			queryB, _, _ := strings.Cut(b, "\t")
			return strings.Compare(queryA, queryB)
		})
		return strings.Join(lines, "\n")
	}
	want := normalize(runSeedIndexCommand(t, args...))
	for _, tc := range []struct {
		budget string
		memory bool
	}{
		{"24", false}, {"18", false}, {"14", true},
	} {
		limited := append(slices.Clone(args), "--max-open-files", tc.budget, "--max-seed-memory", "96")
		if tc.memory {
			limited = append(limited, "-w")
		}
		got, stderr, err := runResourceLimitedCommand(t, "open-files="+tc.budget, limited...)
		if err != nil {
			t.Fatalf("%+v: %v, %s", tc, err, stderr)
		}
		if normalize(got) != want {
			t.Fatalf("%+v: file-budgeted spill changed search results", tc)
		}
	}
	// Reject an impossible reservation before starting collection or producing TSV.
	out, stderr, err := runResourceLimitedCommand(t, "open-files=24", append(slices.Clone(args), "--max-open-files", "17", "--max-seed-memory", "96")...)
	if err == nil || !strings.Contains(string(stderr), "should be >= 18") || len(out) != 0 {
		t.Fatalf("invalid spill budget: error=%v stdout=%s stderr=%s", err, out, stderr)
	}
}
