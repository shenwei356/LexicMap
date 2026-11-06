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
	"path/filepath"
	"strings"
	"testing"
)

// Coverage filtering must retain the full hit when both input sequences are
// reverse complemented, including when the chain has unequal anchor lengths.
func TestSearchVariableAnchorCoverage(t *testing.T) {
	for _, reverse := range []bool{false, true} {
		t.Run(fmt.Sprintf("reverse=%t", reverse), func(t *testing.T) {
			dir := t.TempDir()
			q, s := []byte(variableAnchorQuery), []byte(variableAnchorTarget)
			if reverse {
				RC(q)
				RC(s)
			}
			qfile := writeTestInput(t, filepath.Join(dir, "q.fa"), ">q\n"+string(q)+"\n")
			sfile := writeTestInput(t, filepath.Join(dir, "s.fa"), ">s\n"+string(s)+"\n")
			db := filepath.Join(dir, "db")
			runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-k", "21", "-c", "1", "--partitions", "16", "-O", db, sfile)
			output := runSeedIndexCommand(t, "search", "--quiet", "-j", "2", "-d", db, "-n", "0", "-N", "0", "-Q", "95", qfile)
			lines := strings.Split(strings.TrimSpace(string(output)), "\n")
			if len(lines) != 2 {
				t.Fatalf("got %d hits, want one", len(lines)-1)
			}
			fields := strings.Split(lines[1], "\t")
			for _, i := range []int{5, 8} {
				if fields[i] != "100.000" {
					t.Errorf("coverage column %d = %s, want 100.000", i, fields[i])
				}
			}
			if got := strings.Join(fields[12:16], "\t"); got != "1\t156\t1\t156" {
				t.Errorf("coordinates=%s, want full-length hit", got)
			}
		})
	}
}
