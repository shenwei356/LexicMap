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

// TestSAMAlignmentLengths distinguishes the three lengths and preserves tag order.
func TestSAMAlignmentLengths(t *testing.T) {
	file := filepath.Join(t.TempDir(), "span.tsv")
	header := auditSearchHeader + "\tcigar\tqseq\tsseq\talign\n"
	row := "q\t12\t1\tg\ts\t58.3\t1\t1\t58.3\t9\t55.6\t3\t3\t9\t10\t17\t+\t100\t0\t100\t2M2D1X1I3M\tAC--GTACG\tACTTA-ACG\t||    |||\n"
	writeTestInput(t, file, header+row)
	out := runSeedIndexCommand(t, "utils", "2sam", "--quiet", file)
	records := 0
	for _, line := range strings.Split(strings.TrimSpace(string(out)), "\n") {
		if strings.HasPrefix(line, "@") {
			continue
		}
		fields := strings.Split(line, "\t")
		records++
		if len(fields) != 16 || fields[8] != "0" || fields[5] != "2S2M2D1X1I3M3S" || strings.Join(fields[13:], "\t") != "ZA:i:9\tZQ:i:7\tZR:i:8" {
			t.Fatalf("expected alignment length 9, query span 7, reference span 8: %s", line)
		}
	}
	if records != 1 {
		t.Fatalf("expected one SAM record: %s", out)
	}
	for _, coords := range []string{"0\t9", "3\t2", "bad\t9", "3\tbad", "3\t13"} {
		writeTestInput(t, file, header+strings.Replace(row, "\t3\t9\t10\t17\t", "\t"+coords+"\t10\t17\t", 1))
		_, stderr, err := runAuditCommand(t, "", "utils", "2sam", "--quiet", file)
		if err == nil || !strings.Contains(string(stderr), "invalid query coordinates") {
			t.Fatalf("accepted invalid query coordinates %q: %v %s", coords, err, stderr)
		}
	}
	for _, coords := range []string{"0\t17", "10\t9", "bad\t17", "10\tbad"} {
		writeTestInput(t, file, header+strings.Replace(row, "\t10\t17\t", "\t"+coords+"\t", 1))
		_, stderr, err := runAuditCommand(t, "", "utils", "2sam", "--quiet", file)
		if err == nil || !strings.Contains(string(stderr), "invalid subject coordinates") {
			t.Fatalf("accepted invalid coordinates %q: %v %s", coords, err, stderr)
		}
	}
}

func TestSAMQualityDecoratedLength(t *testing.T) {
	dir := t.TempDir()
	var body strings.Builder
	body.WriteString(auditSearchHeader + "\tcigar\tqseq\tssseq\talign\n")
	for i, score := range []int{100, 90} {
		fmt.Fprintf(&body, "q\t100\t2\tg%d\ts\t100\t1\t1\t100\t100\t100\t0\t1\t100\t1\t100\t+\t100\t0\t%d\t100M\t%s\t%s\t%s\n",
			i, score, strings.Repeat("A", 100), strings.Repeat("A", 100), strings.Repeat("|", 100))
	}
	plain := writeTestInput(t, filepath.Join(dir, "plain.tsv"), body.String())
	quality := writeTestInput(t, filepath.Join(dir, "quality.tsv"),
		strings.ReplaceAll(body.String(), "\t100\t100\t0\t", "\t100:40.0\t100\t0\t"))
	want := runSeedIndexCommand(t, "utils", "2sam", "--quiet", plain)
	got := runSeedIndexCommand(t, "utils", "2sam", "--quiet", quality)
	if !bytes.Equal(got, want) {
		t.Fatalf("quality suffix changed SAM:\n%s\n%s", got, want)
	}
	if bytes.Count(got, []byte("\tZA:i:100\tZQ:i:100\tZR:i:100\n")) != 2 {
		t.Fatalf("quality suffix changed alignment length tags: %s", got)
	}
	for _, invalid := range []string{"bad", "0", "-1", "bad:40.0"} {
		writeTestInput(t, quality, strings.ReplaceAll(body.String(), "\t100\t100\t0\t", "\t"+invalid+"\t100\t0\t"))
		_, stderr, err := runAuditCommand(t, "", "utils", "2sam", "--quiet", quality)
		if err == nil || !strings.Contains(string(stderr), "invalid alenHSP") {
			t.Fatalf("accepted invalid length %q: %v %s", invalid, err, stderr)
		}
	}
}

// Raw scores cannot be recovered from the rounded bitscore in search output.
func TestSAMRawAlignmentScore(t *testing.T) {
	for _, tc := range []struct {
		cigar                                       string
		qlen, tlen, aligned, score, edits, bitscore int
	}{
		{"1000M", 1000, 1000, 1000, 2000, 0, 1804},
		{"500M1X499M", 1000, 1000, 1000, 1995, 1, 1799},
		{"500M1I499M", 1000, 999, 1000, 1991, 1, 1795},
		{"500M2D500M", 1000, 1002, 1002, 1991, 2, 1795},
		{"500M1I1D499M", 1000, 1000, 1001, 1984, 2, 1790},
	} {
		for _, strand := range []string{"+", "-"} {
			t.Run(tc.cigar+strand, func(t *testing.T) {
				file := filepath.Join(t.TempDir(), "input.tsv")
				row := fmt.Sprintf("q\t%d\t1\tg\ts\t100\t1\t1\t100\t%d\t100\t0\t2\t%d\t1\t%d\t%s\t%d\t0\t%d\t%s\t%s\t*\t*\n",
					tc.qlen+2, tc.aligned, tc.qlen+1, tc.tlen, strand, tc.tlen, tc.bitscore, tc.cigar, strings.Repeat("A", tc.qlen))
				writeTestInput(t, file, auditSearchHeader+"\tcigar\tqseq\tsseq\talign\n"+row)
				out := runSeedIndexCommand(t, "utils", "2sam", "--quiet", file)
				want := fmt.Sprintf("\tNM:i:%d\tAS:i:%d\t", tc.edits, tc.score)
				if !strings.Contains(string(out), want) {
					t.Fatalf("missing %s in SAM output", want)
				}
			})
		}
	}
}

func TestSAMMappingQualityAtQueryBoundaryAndEOF(t *testing.T) {
	file := filepath.Join(t.TempDir(), "input.tsv")
	var input strings.Builder
	input.WriteString(auditSearchHeader + "\tcigar\tqseq\tsseq\talign\n")
	for _, query := range []string{"first", "last"} {
		for i, matches := range []int{100, 62, 96} {
			// Two HSPs from the best genome, then the second-best genome.
			genome := "g1"
			if i == 2 {
				genome = "g2"
			}
			// Raw scores are 200, 10, and 180, so the strongest competitor gives MAPQ 4.
			cigar := fmt.Sprintf("%dM", matches)
			if matches < 100 {
				cigar += fmt.Sprintf("%dX", 100-matches)
			}
			bitscore := []int{181, 10, 163}[i]
			fmt.Fprintf(&input, "%s\t100\t2\t%s\ts\t100\t1\t%d\t100\t100\t%d\t0\t1\t100\t1\t100\t+\t100\t0\t%d\t%s\t%s\t%s\t%s\n", query, genome, i+1, matches, bitscore, cigar, strings.Repeat("A", 100), strings.Repeat("A", matches)+strings.Repeat("T", 100-matches), strings.Repeat("|", matches)+strings.Repeat(" ", 100-matches))
		}
	}
	if err := os.WriteFile(file, []byte(input.String()), 0600); err != nil {
		t.Fatal(err)
	}
	out := runSeedIndexCommand(t, "utils", "2sam", "--quiet", file)
	primary := 0
	for _, line := range strings.Split(strings.TrimSpace(string(out)), "\n") {
		if strings.HasPrefix(line, "@") {
			continue
		}
		r := strings.Split(line, "\t")
		if r[1] == "0" {
			primary++
			if r[4] != "4" {
				t.Fatalf("wrong MAPQ: %s", line)
			}
		} else if r[1] != "256" {
			t.Fatalf("wrong secondary flag: %s", line)
		}
	}
	if primary != 2 {
		t.Fatalf("primary alignments: %d", primary)
	}
}
