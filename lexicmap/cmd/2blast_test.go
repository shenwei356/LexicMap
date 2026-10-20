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

func TestBLASTGroupingBoundaries(t *testing.T) {
	dir := t.TempDir()
	row := func(query, genome string) string {
		return fmt.Sprintf("%s\t4\t2\t%s\ts\t100\t1\t1\t100\t4\t100\t0\t1\t4\t1\t4\t+\t4\t0\t8\t4M\tACGT\tACGT\t||||\n", query, genome)
	}
	header := auditSearchHeader + "\tcigar\tqseq\tsseq\talign\n"
	// The old sentinel must also be accepted as an ordinary query ID.
	a := writeTestInput(t, filepath.Join(dir, "a.tsv"), header+row("shenwei356", "g1")+row("q2", "g1")+row("q2", "g2"))
	b := writeTestInput(t, filepath.Join(dir, "b.tsv"), header+row("q2", "g2"))
	out := string(runSeedIndexCommand(t, "utils", "2blast", "--quiet", a, b))
	if strings.Count(out, "Query = ") != 3 || strings.Count(out, "[Subject genome #") != 4 || strings.Count(out, ">s ") != 4 {
		t.Fatalf("missing group headings:\n%s", out)
	}
	if strings.Count(out, "[Subject genome #1/") != 3 || strings.Count(out, "[Subject genome #2/") != 1 {
		t.Fatalf("wrong group numbering:\n%s", out)
	}
}
