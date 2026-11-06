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

func TestBLASTAlignmentLengths(t *testing.T) {
	for _, tc := range []struct {
		name         string
		q, s, middle int
		valid        bool
	}{
		{"short", 4, 4, 4, true},
		{"line boundary", 60, 60, 60, true},
		{"multiple lines", 61, 61, 61, true},
		{"short subject", 61, 60, 61, false},
		{"long subject", 60, 61, 60, false},
		{"short annotation", 61, 61, 60, false},
		{"long annotation", 60, 60, 61, false},
		{"empty", 0, 0, 0, false},
	} {
		t.Run(tc.name, func(t *testing.T) {
			row := fmt.Sprintf("q\t61\t1\tg\ts\t100\t1\t1\t100\t%d\t100\t0\t1\t%d\t1\t%d\t+\t61\t0\t100\t%dM\t%s\t%s\t%s\n", tc.q, tc.q, tc.s, tc.q, strings.Repeat("A", tc.q), strings.Repeat("A", tc.s), strings.Repeat("|", tc.middle))
			input := auditSearchHeader + "\tcigar\tqseq\tsseq\talign\n\n" + row
			out, stderr, err := runAuditCommand(t, input, "utils", "2blast", "--quiet", "-")
			if tc.valid {
				if err != nil || strings.Count(string(out), "Query  ") != (tc.q+59)/60 {
					t.Fatalf("valid alignment failed: %v\n%s\n%s", err, stderr, out)
				}
			} else if err == nil || strings.Contains(string(stderr), "panic:") || !strings.Contains(string(stderr), "invalid alignment lengths") || !strings.Contains(string(stderr), "line 3") {
				t.Fatalf("expected alignment-length error at line 3, got %v\n%s", err, stderr)
			}
		})
	}
}

func TestBLASTGappedCoordinates(t *testing.T) {
	for _, strand := range []string{"+", "-"} {
		t.Run(strand, func(t *testing.T) {
			// Gaps cross the 60-column display boundary; each strand must end at
			// the corresponding inclusive coordinate despite different spans.
			q, s := strings.Repeat("A", 59)+"--AC", strings.Repeat("A", 59)+"GT-C"
			row := fmt.Sprintf("q\t100\t1\tg\ts\t61\t1\t1\t61\t63\t95.2\t3\t3\t63\t10\t71\t%s\t100\t0\t100\t59M2D1I1M\t%s\t%s\t%s\n", strand, q, s, strings.Repeat("|", 59)+"   |")
			out := runSeedIndexCommand(t, "utils", "2blast", "--quiet", writeTestInput(t, filepath.Join(t.TempDir(), "input.tsv"), auditSearchHeader+"\tcigar\tqseq\tsseq\talign\n"+row))
			var query, subject string
			for _, line := range strings.Split(string(out), "\n") {
				if strings.HasPrefix(line, "Query  ") {
					query = strings.Join(strings.Fields(line), " ")
				} else if strings.HasPrefix(line, "Sbjct  ") {
					subject = strings.Join(strings.Fields(line), " ")
				}
			}
			wantSubject := "Sbjct 70 T-C 71"
			if strand == "-" {
				wantSubject = "Sbjct 11 T-C 10"
			}
			if query != "Query 62 -AC 63" || subject != wantSubject {
				t.Fatalf("wrong terminal coordinates: %s; %s", query, subject)
			}
		})
	}
}

func TestBLASTInvalidCoordinates(t *testing.T) {
	base := strings.Split("q\t4\t1\tg\ts\t100\t1\t1\t100\t4\t100\t0\t1\t4\t1\t4\t+\t4\t0\t8\t4M\tACGT\tACGT\t||||", "\t")
	for _, tc := range []struct {
		name, value, message string
		column               int
	}{
		{"bad query start", "bad", "invalid query coordinates", 12},
		{"zero query start", "0", "invalid query coordinates", 12},
		{"reversed query", "5", "invalid query coordinates", 12},
		{"query overrun", "5", "invalid query coordinates", 13},
		{"query length overflow", "18446744073709551616", "invalid query coordinates", 1},
		{"bad subject start", "bad", "invalid subject coordinates", 14},
		{"zero subject start", "0", "invalid subject coordinates", 14},
		{"reversed subject", "5", "invalid subject coordinates", 14},
		{"subject overrun", "5", "invalid subject coordinates", 15},
		{"subject length overflow", "18446744073709551616", "invalid subject coordinates", 17},
		{"invalid strand", "?", "invalid strand", 16},
		{"query span", "3", "alignment spans", 13},
		{"subject span", "3", "alignment spans", 15},
	} {
		t.Run(tc.name, func(t *testing.T) {
			fields := append([]string(nil), base...)
			fields[tc.column] = tc.value
			input := auditSearchHeader + "\tcigar\tqseq\tsseq\talign\n\n" + strings.Join(fields, "\t") + "\n"
			_, stderr, err := runAuditCommand(t, input, "utils", "2blast", "--quiet", "-")
			if err == nil || strings.Contains(string(stderr), "panic:") || !strings.Contains(string(stderr), tc.message) || !strings.Contains(string(stderr), "line 3") {
				t.Fatalf("expected %s at line 3: %v\n%s", tc.message, err, stderr)
			}
		})
	}
}
