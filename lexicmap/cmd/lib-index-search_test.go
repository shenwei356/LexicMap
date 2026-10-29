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
	"math/rand"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
)

// TestSearchMinAlignmentLength checks the final HSP length after extending a
// single exact anchor, including the inclusive boundary and both strands.
func TestSearchMinAlignmentLength(t *testing.T) {
	dir := t.TempDir()
	rng := rand.New(rand.NewSource(193))
	var refs, queries strings.Builder
	wantLengths := make(map[string]int)
	for _, length := range []int{31, 49, 50, 51} {
		q := make([]byte, length)
		for i := range q {
			q[i] = "ACGT"[rng.Intn(4)]
		}
		subject := append([]byte(nil), q...)
		if length > 31 {
			// Only bases 18..34 form a seed of at least 11 bp. Extension
			// must recover the surrounding matches before applying -l.
			for _, pos := range []int{7, 16, 34, 43} {
				subject[pos] = "ACGT"[(strings.IndexByte("ACGT", subject[pos])+1)%4]
			}
		}
		fmt.Fprintf(&refs, ">s%d\n%s\n", length, subject)
		for _, reverse := range []bool{false, true} {
			id := fmt.Sprintf("q%d_rc%t", length, reverse)
			sequence := append([]byte(nil), q...)
			if reverse {
				RC(sequence)
			}
			fmt.Fprintf(&queries, ">%s\n%s\n", id, sequence)
			wantLengths[id] = length
		}
	}
	ref := writeTestInput(t, filepath.Join(dir, "ref.fa"), refs.String())
	query := writeTestInput(t, filepath.Join(dir, "query.fa"), queries.String())
	db := filepath.Join(dir, "db")
	// A 17-mer index makes the short common seed available to coarse search.
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-k", "17", "-m", "1024", "-c", "2", "--partitions", "4", "-O", db, ref)
	for _, minimum := range []int{30, 50, 51} {
		for _, all := range []bool{false, true} {
			t.Run(fmt.Sprintf("minimum=%d/all=%t", minimum, all), func(t *testing.T) {
				args := []string{"search", "--quiet", "-j", "1", "-d", db, "-n", "0", "-N", "0"}
				if minimum != 50 { // exercise the default as well as explicit -l
					args = append(args, "-l", strconv.Itoa(minimum))
				}
				if all {
					args = append(args, "-a")
				}
				output := runSeedIndexCommand(t, append(args, query)...)
				reader := csv.NewReader(strings.NewReader(string(output)))
				reader.Comma = '\t'
				rows, err := reader.ReadAll()
				if err != nil {
					t.Fatal(err)
				}
				seen := make(map[string]int)
				for _, row := range rows[1:] {
					length, ok := wantLengths[row[0]]
					if !ok || row[4] != fmt.Sprintf("s%d", length) || row[9] != strconv.Itoa(length) {
						t.Fatalf("unexpected alignment: %v", row[:20])
					}
					if length < minimum {
						t.Errorf("HSP shorter than -l %d: %v", minimum, row[:20])
					}
					strand := "+"
					if strings.HasSuffix(row[0], "rctrue") {
						strand = "-"
					}
					if row[16] != strand || row[12] != "1" || row[13] != strconv.Itoa(length) || row[14] != "1" || row[15] != strconv.Itoa(length) {
						t.Fatalf("wrong coordinates or strand: %v", row[:20])
					}
					seen[row[0]]++
				}
				for id, length := range wantLengths {
					want := 0
					if length >= minimum {
						want = 1
					}
					if seen[id] != want {
						t.Errorf("%s: got %d hits, want %d", id, seen[id], want)
					}
				}
			})
		}
	}
}

// TestSearchMinAlignmentLengthAcrossContigs exercises the separate completion
// path for the preceding contig, including an HSP shortened by refinement.
func TestSearchMinAlignmentLengthAcrossContigs(t *testing.T) {
	dir := t.TempDir()
	long := indexTestDNA(1000)
	ref := writeTestInput(t, filepath.Join(dir, "ref.fa"), ">repeat\n"+repeatWindowTarget+"\n>long\n"+long+"\n")
	q := repeatWindowQuery + strings.Repeat("N", 100) + long
	query := writeTestInput(t, filepath.Join(dir, "query.fa"), ">forward\n"+q+"\n")
	db := filepath.Join(dir, "db")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-c", "2", "--partitions", "4", "-O", db, ref)
	for _, minimum := range []int{350, 360} {
		for _, all := range []bool{false, true} {
			t.Run(fmt.Sprintf("minimum=%d/all=%t", minimum, all), func(t *testing.T) {
				args := []string{"search", "--quiet", "-j", "1", "-d", db, "-n", "0", "-N", "0",
					"--align-ext-len", "10000", "--align-max-gap", "100", "-l", strconv.Itoa(minimum)}
				if all {
					args = append(args, "-a")
				}
				output := runSeedIndexCommand(t, append(args, query)...)
				reader := csv.NewReader(strings.NewReader(string(output)))
				reader.Comma = '\t'
				rows, err := reader.ReadAll()
				if err != nil {
					t.Fatal(err)
				}
				seen := make(map[string]bool)
				for _, row := range rows[1:] {
					length, err := strconv.Atoi(row[9])
					if err != nil || length < minimum {
						t.Fatalf("HSP shorter than -l %d: %v", minimum, row[:20])
					}
					seen[row[4]] = true
				}
				if !seen["long"] {
					t.Error("lost the long alignment")
				}
				if minimum == 350 && !seen["repeat"] {
					t.Error("missing refinement control")
				}
			})
		}
	}
}

// TestSearchSparseAnchors retains a 34-point chain of two 17-bp anchors even
// when -l and -i would previously have raised its score threshold above 34.
func TestSearchSparseAnchors(t *testing.T) {
	dir := t.TempDir()
	q := []byte(indexTestDNA(100))
	target := append([]byte(nil), q...)
	for _, pos := range []int{10, 28, 39, 50, 61, 79, 90} {
		target[pos] = "ACGT"[(strings.IndexByte("ACGT", target[pos])+1)%4]
	}
	rc := append([]byte(nil), q...)
	RC(rc)
	ref := writeTestInput(t, filepath.Join(dir, "ref.fa"), ">target\n"+string(target)+"\n")
	query := writeTestInput(t, filepath.Join(dir, "query.fa"), ">forward\n"+string(q)+"\n>reverse\n"+string(rc)+"\n")
	db := filepath.Join(dir, "db")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-k", "17", "-m", "1024", "-c", "2", "--partitions", "4", "-O", db, ref)
	for _, limits := range [][2]int{{50, 70}, {60, 90}} {
		t.Run(fmt.Sprintf("length=%d/identity=%d", limits[0], limits[1]), func(t *testing.T) {
			output := runSeedIndexCommand(t, "search", "--quiet", "-j", "1", "-d", db, "-n", "0", "-N", "0",
				"-l", strconv.Itoa(limits[0]), "-i", strconv.Itoa(limits[1]), query)
			reader := csv.NewReader(strings.NewReader(string(output)))
			reader.Comma = '\t'
			rows, err := reader.ReadAll()
			if err != nil {
				t.Fatal(err)
			}
			seen := make(map[string]int)
			for _, row := range rows[1:] {
				strand := "+"
				if row[0] == "reverse" {
					strand = "-"
				}
				if row[4] != "target" || row[9] != "100" || row[10] != "93.000" || row[12] != "1" || row[13] != "100" ||
					row[14] != "1" || row[15] != "100" || row[16] != strand {
					t.Fatalf("unexpected alignment: %v", row[:20])
				}
				seen[row[0]]++
			}
			if len(rows) != 3 || seen["forward"] != 1 || seen["reverse"] != 1 {
				t.Fatalf("lost sparse-anchor alignments: %v", rows)
			}
		})
	}
}

// TestSetSeqCompareOptionsChainingScore covers sampled-anchor callers that
// bypass Compare, and proves per-call overrides do not leak into Chain.
func TestSetSeqCompareOptionsChainingScore(t *testing.T) {
	opt := DefaultSeqComparatorOptions
	opt.MinPrefix = 13
	opt.Chaining2Options = DefaultChaining2Options
	before := opt
	idx := &Index{}
	idx.SetSeqCompareOptions(&opt)
	chainer := idx.poolChainers2.Get().(*Chainer2)
	defer idx.poolChainers2.Put(chainer)
	// Two 13-bp anchors with a one-base gap score 25: pass 22, fail 26.
	subs := []SubstrPair{{Len: 13}, {QBegin: 56, TBegin: 57, Len: 13}}
	for _, override := range []bool{false, true, false} {
		var chains *[]*Chain2Result
		if override {
			chains, _, _, _, _, _, _, _ = chainer.chainWithMinScore(&subs, 22)
		} else {
			chains, _, _, _, _, _, _, _ = chainer.Chain(&subs)
		}
		got := chains != nil
		if chains != nil {
			RecycleChaining2Result(chains)
		}
		if got != override {
			t.Errorf("override=%t: retained=%t, want %t", override, got, override)
		}
	}
	if opt != before {
		t.Fatal("score derivation mutated caller options")
	}
}

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
