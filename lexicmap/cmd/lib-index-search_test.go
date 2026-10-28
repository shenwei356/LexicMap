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
	"encoding/csv"
	"fmt"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
)

// Spare seed/spill/query capacity must never displace the last reader of a batch.
func TestSearchFilePlan(t *testing.T) {
	for _, memory := range []bool{false, true} {
		for _, seedSpill := range []bool{false, true} {
			for _, payload := range []bool{false, true} {
				opt := DefaultIndexSearchingOptions
				opt.NumCPUs, opt.MaxSeedSearchingConcurrency = 16, 16
				opt.MaxQueryFileConcurrency, opt.MaxQueryConcurrency = 16, 16
				opt.InMemorySearch, opt.OutputSeq = memory, payload
				opt.MaxSeedMemory, opt.MaxAlignResultMemory = 0, 0
				minimum := reservedOpenFiles + 4 // four genome batches
				if !memory {
					minimum += 2 // two seed chunks
				}
				if seedSpill {
					opt.MaxSeedMemory = 96
					minimum += 3
				}
				if payload {
					opt.MaxAlignResultMemory = 1
					minimum++
				}
				for budget := minimum - 1; budget <= minimum+128; budget++ {
					opt.MaxOpenFiles = budget
					p, err := planSearchFiles(&opt, 2, 4)
					if budget < minimum {
						if err == nil || !strings.Contains(err.Error(), fmt.Sprintf("should be >= %d", minimum)) {
							t.Fatalf("accepted budget %d below minimum %d: %+v, %v", budget, minimum, p, err)
						}
						continue
					}
					if err != nil || p.genomeReaders < 1 || p.seedReaders < 1 || p.queryReaders < 1 {
						t.Fatalf("missing reader pool: budget=%d plan=%+v err=%v", budget, p, err)
					}
					used := reservedOpenFiles + 4*p.genomeReaders + p.queryReaders - 1 + p.seedSpillFiles + p.alignmentFiles
					if !memory {
						used += 2 * p.seedReaders
					}
					if used > budget || (seedSpill && p.seedSpillFiles < 3) || (payload && p.alignmentFiles < 1) {
						t.Fatalf("invalid reservation: budget=%d used=%d plan=%+v", budget, used, p)
					}
					if budget == minimum && (p.genomeReaders != 1 || p.queryReaders != 1 || (!memory && p.seedReaders != 1)) {
						t.Fatalf("minimum budget did not preserve one reader per pool: %+v", p)
					}
				}
			}
		}
	}
}

func TestSearchRejectsInsufficientPoolsBeforeRecount(t *testing.T) {
	db, query := buildMultiBatchSearchFixture(t)
	file := filepath.Join(db, FileInfo)
	info, err := readIndexInfo(file)
	if err != nil {
		t.Fatal(err)
	}
	info.InputBases = 0 // old indexes recount this field during initialization
	if err := writeIndexInfo(file, info); err != nil {
		t.Fatal(err)
	}
	for _, command := range [][]string{{"search"}, {"genome", "search", "-m", "0"}} {
		for _, memory := range []bool{false, true} {
			minimum := reservedOpenFiles + info.GenomeBatches
			if !memory {
				minimum += info.Chunks
			}
			args := append(append([]string{}, command...), "--quiet", "-d", db, "--max-open-files", strconv.Itoa(minimum-1), query)
			if memory {
				args = append(args, "-w")
			}
			out, stderr, err := runAuditCommand(t, "", args...)
			if err == nil || len(out) != 0 || !strings.Contains(string(stderr), fmt.Sprintf("should be >= %d", minimum)) {
				t.Fatalf("invalid budget accepted: %v, %s, %s", err, out, stderr)
			}
			got, err := readIndexInfo(file)
			if err != nil || got.InputBases != 0 {
				t.Fatalf("invalid budget recounted the index: %+v, %v", got, err)
			}
		}
	}
}

// writeSearchContigFixture gives each contig different random bases and length.
func writeSearchContigFixture(t *testing.T) (string, string, []string) {
	t.Helper()
	dir := t.TempDir()
	dna := indexTestDNA(4800)
	contigs := []string{dna[:1000], dna[1000:3200], dna[3200:]}
	var fasta strings.Builder
	for i, s := range contigs {
		fmt.Fprintf(&fasta, ">%c\n%s\n", 'a'+i, s)
	}
	ref := writeTestInput(t, filepath.Join(dir, "multi.fa"), fasta.String())
	return dir, ref, contigs
}

// TestSearchCrossContigCoordinates checks both strands after a contig transition.
// The large extension deliberately places several contigs in one alignment window.
func TestSearchCrossContigCoordinates(t *testing.T) {
	dir, ref, contigs := writeSearchContigFixture(t)
	db := filepath.Join(dir, "db")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-c", "2", "--partitions", "4", "-O", db, ref)
	queries := make(map[string]string)
	expected := make(map[string]int)
	var fasta strings.Builder
	for _, tc := range []struct {
		id, sequence string
		hits         int
	}{
		{"abc", strings.Join(contigs, ""), 3}, {"ab", contigs[0] + contigs[1], 2}, {"b", contigs[1], 1},
	} {
		for _, reverse := range []bool{false, true} {
			name, sequence := tc.id, []byte(tc.sequence)
			if reverse {
				name += "_rc"
				RC(sequence)
			}
			queries[name], expected[name] = string(sequence), tc.hits
			fmt.Fprintf(&fasta, ">%s\n%s\n", name, sequence)
		}
	}
	query := writeTestInput(t, filepath.Join(dir, "query.fa"), fasta.String())
	output := runSeedIndexCommand(t, "search", "--quiet", "-j", "1", "-d", db, "-n", "0", "-N", "0", "-Q", "0", "-q", "0", "-a", "--align-ext-len", "10000", query)
	reader := csv.NewReader(strings.NewReader(string(output)))
	reader.Comma = '\t'
	rows, err := reader.ReadAll()
	if err != nil {
		t.Fatal(err)
	}
	seen := make(map[string]int)
	for _, row := range rows[1:] {
		q, ok := queries[row[0]]
		if !ok {
			t.Fatalf("unknown query: %q", row[0])
		}
		if len(row[4]) != 1 || row[4][0] < 'a' || row[4][0] > 'c' {
			t.Fatalf("wrong contig: %v", row[:20])
		}
		subject := contigs[int(row[4][0]-'a')]
		number := func(i int) int {
			n, e := strconv.Atoi(row[i])
			if e != nil {
				t.Fatal(e)
			}
			return n
		}
		qs, qe, ss, se := number(12), number(13), number(14), number(15)
		if qs < 1 || qe > len(q) || qs > qe || ss != 1 || se != len(subject) || number(17) != len(subject) {
			t.Fatalf("wrong coordinates: %v", row[:21])
		}
		target := []byte(subject[ss-1 : se])
		if row[16] == "-" {
			RC(target)
		}
		if strings.ReplaceAll(row[21], "-", "") != q[qs-1:qe] || strings.ReplaceAll(row[22], "-", "") != string(target) {
			t.Fatalf("alignment text does not match coordinates: %v", row[:21])
		}
		if row[20] != fmt.Sprintf("%dM", len(subject)) || qe-qs+1 != len(subject) {
			t.Fatalf("exact contig lost bases: %v", row[:21])
		}
		seen[row[0]]++
	}
	for query, n := range expected {
		if seen[query] != n {
			t.Fatalf("%s: got %d HSPs, want %d", query, seen[query], n)
		}
	}
}

// TestSearchRepeatRefinement checks forward/reverse contig coordinates, CIGAR
// text, and coverage fallback through the ordinary search command.
func TestSearchRepeatRefinement(t *testing.T) {
	dir := t.TempDir()
	rc := []byte(repeatWindowTarget)
	RC(rc)
	ref := writeTestInput(t, filepath.Join(dir, "repeat.fa"), ">forward\n"+repeatWindowTarget+"\n>reverse\n"+string(rc)+"\n")
	rcq := []byte(repeatWindowQuery)
	RC(rcq)
	query := writeTestInput(t, filepath.Join(dir, "query.fa"), ">query\n"+repeatWindowQuery+"\n>rc-query\n"+string(rcq)+"\n")
	db := filepath.Join(dir, "db")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-c", "2", "--partitions", "4", "-O", db, ref)
	for _, coverage := range []struct{ hsp, genome, gap string }{{"0", "0", "20"}, {"0", "0", "100"}, {"90", "0", "100"}, {"0", "90", "100"}} {
		output := runSeedIndexCommand(t, "search", "--quiet", "-j", "1", "-d", db, "-n", "0", "-N", "0", "-Q", coverage.genome, "-q", coverage.hsp, "-a", "--align-max-gap", coverage.gap, query)
		reader := csv.NewReader(strings.NewReader(string(output)))
		reader.Comma = '\t'
		rows, err := reader.ReadAll()
		if err != nil {
			t.Fatal(err)
		}
		seen := map[string]bool{}
		for _, row := range rows[1:] {
			number := func(i int) int {
				n, e := strconv.Atoi(row[i])
				if e != nil {
					t.Fatal(e)
				}
				return n
			}
			qs, qe, ss, se := number(12), number(13), number(14), number(15)
			querySeq := repeatWindowQuery
			wantBegin, wantEnd, wantCIGAR := 46, 399, "5M3X17M2X7M2X2M3X313M"
			if row[0] == "rc-query" {
				querySeq = string(rcq)
				wantBegin, wantEnd, wantCIGAR = 19, 372, "313M3X2M2X7M2X17M3X5M"
			}
			target := repeatWindowTarget
			if row[4] == "reverse" {
				target = string(rc)
			}
			subject := []byte(target[ss-1 : se])
			if row[16] == "-" {
				RC(subject)
			}
			if strings.ReplaceAll(row[21], "-", "") != querySeq[qs-1:qe] || strings.ReplaceAll(row[22], "-", "") != string(subject) {
				t.Fatalf("inconsistent sequence/coordinates: %v", row)
			}
			if ss != 1 || se != 354 {
				continue
			}
			if coverage.hsp == "0" && coverage.genome == "0" && qs == wantBegin && qe == wantEnd && number(9) == 354 && number(11) == 0 {
				if row[20] != wantCIGAR {
					t.Fatalf("wrong refined CIGAR: %s", row[20])
				}
				seen[row[0]+"/"+row[4]] = true
			}
			if coverage.hsp == "90" || coverage.genome == "90" {
				// The existing pseudo-alignment already clips the reversed query
				// below 90% coverage; the new refinement must not change filtering.
				if row[0] == "rc-query" {
					t.Fatal("reverse query unexpectedly passed the existing coverage filter")
				}
				if qs == 1 && qe == 399 && number(9) == 408 && number(11) == 63 {
					seen[row[0]+"/"+row[4]] = true
				}
			}
		}
		wantHits := 4
		if coverage.hsp == "90" || coverage.genome == "90" {
			wantHits = 2
		}
		if len(seen) != wantHits {
			t.Fatalf("coverage=%+v: missing expected results: %s", coverage, output)
		}
	}
}
