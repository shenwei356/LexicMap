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
	"slices"
	"strings"
	"testing"
)

func TestSeedPosFileBudget(t *testing.T) {
	db := writeSeedPosFixture(t, 64)
	normalize := func(output []byte) string {
		lines := strings.Split(strings.TrimSpace(string(output)), "\n")
		if len(lines) != 193 {
			t.Fatalf("lost positions: %d lines", len(lines))
		}
		slices.Sort(lines[1:])
		return strings.Join(lines, "\n")
	}
	for _, verbose := range []bool{false, true} {
		args := []string{"utils", "seed-pos", "--quiet", "-d", db, "-a"}
		if verbose {
			args = append(args, "-v")
		}
		baseline := normalize(runSeedIndexCommand(t, append(slices.Clone(args), "-j", "1")...))
		for _, threads := range []string{"1", "16"} {
			dir := t.TempDir()
			out, logFile := filepath.Join(dir, "positions.tsv"), filepath.Join(dir, "seed-pos.log")
			_, stderr, err := runResourceLimitedCommand(t, "open-files=11", append(slices.Clone(args), "-j", threads, "--max-open-files", "11", "-o", out, "--log", logFile)...)
			if err != nil {
				log, _ := os.ReadFile(logFile)
				t.Fatalf("bounded seed-pos failed: %v %s\n%s", err, stderr, log)
			}
			output, err := os.ReadFile(out)
			if err != nil {
				t.Fatal(err)
			}
			if normalize(output) != baseline {
				t.Fatalf("file budget changed seed positions")
			}
		}
	}
}
