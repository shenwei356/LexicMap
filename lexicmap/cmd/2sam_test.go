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

// Invalid MAPQ inputs must not silently become zero or inflate confidence.
func TestSAMNumericFields(t *testing.T) {
	header := auditSearchHeader + "\tcigar\tqseq\tsseq\talign\n"
	base := strings.Split("q\t100\t2\tg\ts\t100\t1\t1\t100\t100\t100\t0\t1\t100\t1\t100\t+\t100\t0\t200\t100M\t"+
		strings.Repeat("A", 100)+"\t"+strings.Repeat("A", 100)+"\t"+strings.Repeat("|", 100), "\t")
	other := append([]string(nil), base...)
	other[7], other[10], other[19], other[20] = "2", "99", "195", "50M1X49M"
	other[22] = strings.Repeat("A", 50) + "T" + strings.Repeat("A", 49)
	other[23] = strings.Repeat("|", 50) + " " + strings.Repeat("|", 49)
	for _, field := range []struct {
		name   string
		column int
		bad    []string
	}{
		{"qcovHSP", 8, []string{"", "bad", "NaN", "Inf", "1e309", "-1", "100.1"}},
		{"pident", 10, []string{"", "bad", "NaN", "Inf", "1e309", "-1", "100.1"}},
		{"gaps", 11, []string{"", "bad", "-1", "101", "1.5", "999999999999999999999"}},
	} {
		for _, value := range field.bad {
			for row := 0; row < 2; row++ {
				t.Run(fmt.Sprintf("%s=%s/row=%d", field.name, value, row), func(t *testing.T) {
					rows := [][]string{append([]string(nil), base...), append([]string(nil), other...)}
					rows[row][field.column] = value
					file := writeTestInput(t, filepath.Join(t.TempDir(), "input.tsv"), header+"\n"+
						strings.Join(rows[0], "\t")+"\n"+strings.Join(rows[1], "\t")+"\n")
					_, stderr, err := runAuditCommand(t, "", "utils", "2sam", "--quiet", file)
					if err == nil || !strings.Contains(string(stderr), "invalid "+field.name) ||
						!strings.Contains(string(stderr), fmt.Sprintf("input.tsv: line %d", row+3)) {
						t.Fatalf("invalid field was not rejected with its location: %v\n%s", err, stderr)
					}
				})
			}
		}
	}
	file := writeTestInput(t, filepath.Join(t.TempDir(), "valid.tsv"), header+
		strings.Join(base, "\t")+"\n"+strings.Join(other, "\t")+"\n")
	want := runSeedIndexCommand(t, "utils", "2sam", "--quiet", file)
	if !bytes.Contains(want, []byte("q\t0\ts\t1\t1\t100M\t")) {
		t.Fatalf("valid primary MAPQ changed: %s", want)
	}
	// Keep valid numeric spellings and --show-avg-qual support unchanged.
	base[8], base[9], base[10], base[11] = "1e2", "100:40.0", "100.000", "00"
	writeTestInput(t, file, header+strings.Join(base, "\t")+"\n"+strings.Join(other, "\t")+"\n")
	got := runSeedIndexCommand(t, "utils", "2sam", "--quiet", file)
	if !bytes.Equal(got, want) {
		t.Fatalf("valid numeric spelling changed SAM output:\n%s\n%s", got, want)
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

// One SAM reference name cannot describe subjects of different lengths.
func TestSAMReferenceLengths(t *testing.T) {
	header := auditSearchHeader + "\tcigar\tqseq\tsseq\talign\n"
	row := func(genome, length, start, end string) string {
		return fmt.Sprintf("q\t4\t2\t%s\ts\t100\t1\t1\t100\t4\t100\t0\t1\t4\t%s\t%s\t+\t%s\t0\t8\t4M\tACGT\tACGT\t||||\n", genome, start, end, length)
	}
	for _, separate := range []bool{false, true} {
		for _, reverse := range []bool{false, true} {
			for _, concat := range []bool{false, true} {
				t.Run(fmt.Sprintf("separate=%t/reverse=%t/concat=%t", separate, reverse, concat), func(t *testing.T) {
					dir := t.TempDir()
					rows := []string{row("g1", "4", "1", "4"), row("g2", "10", "7", "10")}
					if reverse {
						rows[0], rows[1] = rows[1], rows[0]
					}
					args := []string{"utils", "2sam", "--quiet"}
					if concat {
						args = append(args, "-c")
					}
					if separate {
						for i, record := range rows {
							args = append(args, writeTestInput(t, filepath.Join(dir, fmt.Sprintf("input%d.tsv", i)), header+record))
						}
					} else {
						args = append(args, writeTestInput(t, filepath.Join(dir, "input.tsv"), header+strings.Join(rows, "")))
					}
					out, stderr, err := runAuditCommand(t, "", args...)
					if !concat {
						line := "line 3"
						if separate {
							line = "line 2"
							if !strings.Contains(string(stderr), "input1.tsv") {
								t.Errorf("missing conflicting file: %s", stderr)
							}
						}
						if err == nil || len(out) != 0 || !strings.Contains(string(stderr), "inconsistent lengths for SAM reference") || !strings.Contains(string(stderr), line) {
							t.Fatalf("expected a reference-length error before the SAM header: %v\n%s\n%s", err, stderr, out)
						}
					} else if err != nil || !bytes.Contains(out, []byte("@SQ\tSN:g1~s\tLN:4\n")) || !bytes.Contains(out, []byte("@SQ\tSN:g2~s\tLN:10\n")) {
						t.Fatalf("distinct references were not preserved: %v\n%s\n%s", err, stderr, out)
					}
				})
			}
		}
	}
	// Repeated references are valid; compare numeric lengths, not their spelling.
	dir := t.TempDir()
	first := writeTestInput(t, filepath.Join(dir, "first.tsv"), header+row("g1", "4", "1", "4"))
	second := writeTestInput(t, filepath.Join(dir, "second.tsv"), header+row("g2", "004", "1", "4"))
	out := runSeedIndexCommand(t, "utils", "2sam", "--quiet", first, second)
	if bytes.Count(out, []byte("@SQ\t")) != 1 || !bytes.Contains(out, []byte("@SQ\tSN:s\tLN:4\n")) || bytes.Count(out, []byte("\nq\t")) != 2 {
		t.Fatalf("repeated reference changed: %s", out)
	}
}

func TestSAMNames(t *testing.T) {
	header := auditSearchHeader + "\tcigar\tqseq\tsseq\talign\n"
	base := strings.Split("q\t4\t1\tg\ts\t100\t1\t1\t100\t4\t100\t0\t1\t4\t1\t4\t+\t4\t0\t8\t4M\tACGT\tACGT\t||||", "\t")
	for _, tc := range []struct {
		name, value, message string
		column               int
		flags                []string
	}{
		{"empty query", "", "query", 0, nil},
		{"query starts with at", "@q", "query", 0, nil},
		{"query contains at", "q@1", "query", 0, nil},
		{"query space", "q 1", "query", 0, nil},
		{"long query", strings.Repeat("q", 255), "query", 0, nil},
		{"non-ASCII query", "qé", "query", 0, nil},
		{"query control", "q\x7f", "query", 0, nil},
		{"empty reference", "", "reference", 4, nil},
		{"reference star", "*s", "reference", 4, nil},
		{"reference equals", "=s", "reference", 4, nil},
		{"reference space", "s 1", "reference", 4, nil},
		{"reference comma", "s,1", "reference", 4, nil},
		{"reference backslash", "s\\1", "reference", 4, nil},
		{"reference quote", "s\"1", "reference", 4, nil},
		{"reference apostrophe", "s'1", "reference", 4, nil},
		{"reference backtick", "s`1", "reference", 4, nil},
		{"reference parentheses", "s(1)", "reference", 4, nil},
		{"reference brackets", "s[1]", "reference", 4, nil},
		{"reference braces", "s{1}", "reference", 4, nil},
		{"reference angles", "s<1>", "reference", 4, nil},
		{"non-ASCII reference", "sé", "reference", 4, nil},
		{"reference control", "s\x7f", "reference", 4, nil},
		{"concatenated genome", "g one", "reference", 3, []string{"-c"}},
		{"separator", "g", "reference", 3, []string{"-c", "-s", " "}},
	} {
		t.Run(tc.name, func(t *testing.T) {
			fields := append([]string(nil), base...)
			fields[tc.column] = tc.value
			file := writeTestInput(t, filepath.Join(t.TempDir(), "input.tsv"), header+"\n"+strings.Join(fields, "\t")+"\n")
			args := append([]string{"utils", "2sam", "--quiet", file}, tc.flags...)
			out, stderr, err := runAuditCommand(t, "", args...)
			if err == nil || len(out) != 0 || !strings.Contains(string(stderr), "invalid SAM "+tc.message+" name") || !strings.Contains(string(stderr), "line 3") {
				t.Fatalf("invalid name was not rejected before the header: %v\n%s\n%s", err, stderr, out)
			}
		})
	}
	// Query names allow punctuation forbidden in references, and @ is allowed
	// in reference names. Genome names matter only when -c includes them.
	for _, names := range [][2]string{{strings.Repeat("q", 254), "@ref:1|*=~"}, {"q[](){}<>,'\"`\\", "ref/1"}} {
		fields := append([]string(nil), base...)
		fields[0], fields[3], fields[4] = names[0], "g with spaces", names[1]
		file := writeTestInput(t, filepath.Join(t.TempDir(), "valid.tsv"), header+strings.Join(fields, "\t")+"\n")
		out := runSeedIndexCommand(t, "utils", "2sam", "--quiet", file)
		if !bytes.Contains(out, []byte("@SQ\tSN:"+names[1]+"\tLN:4\n")) || !bytes.Contains(out, []byte("\n"+names[0]+"\t0\t"+names[1]+"\t")) {
			t.Fatalf("valid names changed: %s", out)
		}
	}
}

// File boundaries and header-only inputs do not end a continuing query group.
func TestSAMQueryGroupsAcrossFiles(t *testing.T) {
	header := auditSearchHeader + "\tcigar\tqseq\tsseq\talign\n"
	row := func(query, cigar, identity, target, alignment string) string {
		return fmt.Sprintf("%s\t4\t1\tg\ts\t100\t1\t1\t100\t4\t%s\t0\t1\t4\t1\t4\t+\t4\t0\t8\t%s\tACGT\t%s\t%s\n", query, identity, cigar, target, alignment)
	}
	first := row("first", "4M", "100", "ACGT", "||||") + row("shared", "4M", "100", "ACGT", "||||")
	second := row("shared", "3M1X", "75", "ACGA", "||| ") + row("last", "4M", "100", "ACGT", "||||")
	dir := t.TempDir()
	combined := writeTestInput(t, filepath.Join(dir, "combined.tsv"), header+first+second)
	a := writeTestInput(t, filepath.Join(dir, "first.tsv"), header+first)
	empty := writeTestInput(t, filepath.Join(dir, "empty.tsv"), header)
	b := writeTestInput(t, filepath.Join(dir, "second.tsv"), header+second)
	want := runSeedIndexCommand(t, "utils", "2sam", "--quiet", combined)
	for _, files := range [][]string{{a, b}, {a, empty, b}} {
		got := runSeedIndexCommand(t, append([]string{"utils", "2sam", "--quiet"}, files...)...)
		if !bytes.Equal(got, want) {
			t.Fatalf("file boundaries changed query flags, MAPQ, or sequence:\n%s\nwant:\n%s", got, want)
		}
	}
	if !bytes.Contains(want, []byte("\nshared\t0\ts\t1\t25\t")) || !bytes.Contains(want, []byte("\nshared\t256\ts\t1\t0\t")) {
		t.Fatalf("wrong shared-query flags or MAPQ: %s", want)
	}
}

func TestSAMQueryGroupValidation(t *testing.T) {
	header := auditSearchHeader + "\tcigar\tqseq\tsseq\talign\n"
	row := "q\t4\t1\tg\ts\t100\t1\t1\t100\t4\t100\t0\t1\t4\t1\t4\t+\t4\t0\t8\t4M\tACGT\tACGT\t||||\n"
	for _, separate := range []bool{false, true} {
		for _, inconsistentLength := range []bool{false, true} {
			t.Run(fmt.Sprintf("separate=%t/length=%t", separate, inconsistentLength), func(t *testing.T) {
				dir := t.TempDir()
				first, second := row+strings.Replace(row, "q\t", "other\t", 1), row
				message := "noncontiguous records for SAM query"
				if inconsistentLength {
					first, second = row, strings.Replace(row, "q\t4\t", "q\t5\t", 1)
					message = "inconsistent lengths for SAM query"
				}
				args := []string{"utils", "2sam", "--quiet"}
				if separate {
					args = append(args, writeTestInput(t, filepath.Join(dir, "first.tsv"), header+first), writeTestInput(t, filepath.Join(dir, "second.tsv"), header+second))
				} else {
					args = append(args, writeTestInput(t, filepath.Join(dir, "input.tsv"), header+first+second))
				}
				out, stderr, err := runAuditCommand(t, "", args...)
				if err == nil || len(out) != 0 || !strings.Contains(string(stderr), message) {
					t.Fatalf("expected %s before output: %v\n%s\n%s", message, err, stderr, out)
				}
			})
		}
	}
	// Equivalent numeric lengths do not split or invalidate a query group.
	file := writeTestInput(t, filepath.Join(t.TempDir(), "valid.tsv"), header+row+strings.Replace(row, "q\t4\t", "q\t004\t", 1))
	out := runSeedIndexCommand(t, "utils", "2sam", "--quiet", file)
	if bytes.Count(out, []byte("\nq\t0\t")) != 1 || bytes.Count(out, []byte("\nq\t256\t")) != 1 {
		t.Fatalf("numeric query-length equivalence changed grouping: %s", out)
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

// A nonpositive primary score cannot give positive confidence over alternatives.
func TestSAMNonpositivePrimaryScore(t *testing.T) {
	type alignment struct {
		cigar                             string
		identity, gaps, subjectLen, score int
	}
	negative := alignment{"30M1I18X1I19X31M", 61, 2, 98, -3}
	weaker := alignment{"30M1I18X1I20X30M", 60, 2, 98, -8}
	zero := alignment{"60M40X", 60, 0, 100, 0}
	positive := alignment{"61M39X", 61, 0, 100, 5}
	for _, tc := range []struct {
		name           string
		primary, other alignment
	}{
		{"negative versus zero", negative, zero},
		{"negative versus positive", negative, positive},
		{"negative versus better negative", weaker, negative},
		{"negative versus worse negative", negative, weaker},
		{"zero versus negative", zero, negative},
		{"zero versus zero", zero, zero},
	} {
		t.Run(tc.name, func(t *testing.T) {
			var input strings.Builder
			input.WriteString(auditSearchHeader + "\tcigar\tqseq\tsseq\talign\n")
			for i, a := range []alignment{tc.primary, tc.other} {
				fmt.Fprintf(&input, "q\t100\t2\tg%d\ts%d\t100\t1\t1\t100\t100\t%d\t%d\t1\t100\t1\t%d\t+\t%d\t1e10\t0\t%s\t%s\t*\t*\n",
					i, i, a.identity, a.gaps, a.subjectLen, a.subjectLen, a.cigar, strings.Repeat("A", 100))
			}
			file := writeTestInput(t, filepath.Join(t.TempDir(), "input.tsv"), input.String())
			out := runSeedIndexCommand(t, "utils", "2sam", "--quiet", file)
			records := 0
			for _, line := range strings.Split(strings.TrimSpace(string(out)), "\n") {
				if strings.HasPrefix(line, "@") {
					continue
				}
				fields := strings.Split(line, "\t")
				if records == 0 && (fields[1] != "0" || fields[4] != "0" || fields[12] != fmt.Sprintf("AS:i:%d", tc.primary.score)) {
					t.Fatalf("nonpositive primary score gained confidence: %s", line)
				}
				if records == 1 && fields[1] != "256" {
					t.Fatalf("secondary flag changed: %s", line)
				}
				records++
			}
			if records != 2 {
				t.Fatalf("expected two alignments, got %d", records)
			}
		})
	}
}

// TestSAMInvalidAlignments rejects records that cannot describe a valid SAM read.
func TestSAMInvalidAlignments(t *testing.T) {
	base := strings.Split("q\t4\t1\tg\ts\t100\t1\t1\t100\t4\t100\t0\t1\t4\t1\t4\t+\t4\t0\t8\t4M\tACGT\tACGT\t||||", "\t")
	for _, tc := range []struct {
		name, value, message string
		column               int
	}{
		{"missing cigar", "", "invalid CIGAR", 20},
		{"trailing garbage", "4Mgarbage", "invalid CIGAR", 20},
		{"leading garbage", "garbage4M", "invalid CIGAR", 20},
		{"missing count", "M", "invalid CIGAR", 20},
		{"zero count", "0I4M", "invalid CIGAR", 20},
		{"count overflow", "18446744073709551616M", "invalid CIGAR", 20},
		{"unknown operation", "4Z", "invalid CIGAR", 20},
		{"clipped input", "1S4M", "invalid CIGAR", 20},
		{"query span", "3M1D", "CIGAR spans", 20},
		{"subject span", "3M1I", "CIGAR spans", 20},
		{"alignment length", "5", "CIGAR spans", 9},
		{"missing query sequence", "", "query sequence length", 21},
		{"short query sequence", "ACG", "query sequence length", 21},
		{"long query sequence", "ACGTA", "query sequence length", 21},
		{"invalid strand", "?", "invalid strand", 16},
		{"invalid query length", "bad", "invalid query coordinates", 1},
		{"query length overflow", "18446744073709551616", "invalid query coordinates", 1},
		{"invalid subject length", "bad", "invalid subject coordinates", 17},
		{"subject overrun", "3", "invalid subject coordinates", 17},
	} {
		t.Run(tc.name, func(t *testing.T) {
			fields := append([]string(nil), base...)
			fields[tc.column] = tc.value
			file := writeTestInput(t, filepath.Join(t.TempDir(), "input.tsv"), auditSearchHeader+"\tcigar\tqseq\tsseq\talign\n\n"+strings.Join(fields, "\t")+"\n")
			_, stderr, err := runAuditCommand(t, "", "utils", "2sam", "--quiet", file)
			if err == nil || strings.Contains(string(stderr), "panic:") || !strings.Contains(string(stderr), tc.message) || !strings.Contains(string(stderr), "line 3") {
				t.Fatalf("expected %s at line 3: %v\n%s", tc.message, err, stderr)
			}
		})
	}
}

func TestSAMMergedMissingAlignment(t *testing.T) {
	dir := t.TempDir()
	row := "q\t4\t1\tg\ts\t100\t1\t1\t100\t4\t100\t0\t1\t4\t1\t4\t+\t4\t0\t8"
	full := writeTestInput(t, filepath.Join(dir, "full.tsv"), auditSearchHeader+"\tcigar\tqseq\tsseq\talign\n"+row+"\t4M\tACGT\tACGT\t||||\n")
	// The first retained hit selects the format; avoid an unordered score tie.
	shortRow := strings.TrimSuffix(strings.Replace(row, "\tg\t", "\tg2\t", 1), "\t8") + "\t7"
	short := writeTestInput(t, filepath.Join(dir, "short.tsv"), auditSearchHeader+"\n"+shortRow+"\n")
	merged := runSeedIndexCommand(t, "utils", "merge-search-results", "--quiet", full, short)
	lines := strings.Split(strings.TrimSuffix(string(merged), "\n"), "\n")
	if len(lines) != 3 {
		t.Fatalf("expected a header and two merged records: %s", merged)
	}
	for _, line := range lines {
		if len(strings.Split(line, "\t")) != 24 {
			t.Fatalf("expected merged alignment columns: %s", merged)
		}
	}
	file := writeTestInput(t, filepath.Join(dir, "merged.tsv"), string(merged))
	out, stderr, err := runAuditCommand(t, "", "utils", "2sam", "--quiet", file)
	if err == nil || !strings.Contains(string(stderr), "invalid CIGAR") {
		t.Fatalf("accepted missing alignment from mixed-format results: %v\n%s\n%s", err, stderr, out)
	}
}
