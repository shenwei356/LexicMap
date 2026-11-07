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
	"encoding/json"
	"errors"
	"fmt"
	"math/rand"
	"os"
	"path/filepath"
	"slices"
	"strconv"
	"strings"
	"testing"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
	"github.com/shenwei356/wfa"
)

// Prefix matches can all be shorter than k. Their longest anchors must keep
// the reference's ordering and duplicate strand choice in both orientations.
func TestClearSubstrPairsObservedMaximum(t *testing.T) {
	rng := rand.New(rand.NewSource(20261107))
	for _, k := range []int{13, 31, 32, 255} {
		for trial := 0; trial < 120; trial++ {
			longest := trial % k // include zero-length and k-1 cases
			input := make([]SubstrPair, 100+trial)
			for i := range input {
				length := longest
				if i%3 == 0 {
					length = rng.Intn(longest + 1)
				}
				input[i] = SubstrPair{QBegin: int32(rng.Intn(20)), TBegin: int32(rng.Intn(100)), Len: uint8(length), QRC: i%2 == 0, TRC: i%3 == 0}
				if i > 0 && i%5 == 0 {
					input[i] = input[i-1]
					input[i].QRC = !input[i].QRC
				}
			}
			checkClearSubstrPairsReference(t, input, k)
			for i := range input {
				input[i].QBegin = 1000 - input[i].QBegin - int32(input[i].Len)
				input[i].TBegin = 1000 - input[i].TBegin - int32(input[i].Len)
				input[i].QRC = !input[i].QRC
				input[i].TRC = !input[i].TRC
			}
			checkClearSubstrPairsReference(t, input, k)
		}
	}
	// Repeated matches at the same query start retain every distinct target.
	input := make([]SubstrPair, 4000)
	for i := range input {
		input[i] = SubstrPair{TBegin: int32(40 * i), Len: 30}
	}
	want := slices.Clone(input)
	ClearSubstrPairs(&input, 31)
	if !slices.Equal(input, want) {
		t.Fatal("distinct targets were removed or reordered")
	}
}

// Measure the complete cleanup call, including sorting and pooled markers.
func BenchmarkClearSubstrPairsObservedMaximum(b *testing.B) {
	for _, length := range []uint8{30, 31} {
		input := make([]SubstrPair, 4000)
		for i := range input {
			input[i] = SubstrPair{TBegin: int32(40 * i), Len: length}
		}
		b.Run(fmt.Sprintf("length=%d", length), func(b *testing.B) {
			work := make([]SubstrPair, len(input))
			b.ReportAllocs()
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				work = work[:len(input)]
				copy(work, input)
				ClearSubstrPairs(&work, 31)
			}
		})
	}
}

// Closing several batches can report errors concurrently. Empty readers inject
// os.ErrInvalid without closing a descriptor owned by another live reader.
func TestIndexCloseReaderErrors(t *testing.T) {
	file := filepath.Join(t.TempDir(), "genomes.bin")
	w, err := genome.NewWriter(file, 0)
	if err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	directory, err := genome.ReadIndex(file)
	if err != nil {
		t.Fatal(err)
	}
	idx := &Index{
		opt:            &IndexSearchingOptions{},
		hasGenomeRdrs:  true,
		poolGenomeRdrs: make([]chan *genome.Reader, 16),
	}
	idx.genomeReaderIndexes.Store(0, &sharedGenomeReaderIndex{index: directory})
	for i := range idx.poolGenomeRdrs {
		pool := make(chan *genome.Reader, 2)
		idx.poolGenomeRdrs[i] = pool
		pool <- &genome.Reader{}
		rdr, err := directory.NewReader()
		if err != nil {
			t.Fatal(err)
		}
		pool <- rdr
	}
	if err := idx.Close(); !errors.Is(err, os.ErrInvalid) {
		t.Fatalf("Close returned %v, want the injected reader error", err)
	}
	for batch, pool := range idx.poolGenomeRdrs {
		if _, ok := <-pool; ok {
			t.Fatalf("batch %d still contains an unclosed reader", batch)
		}
	}
	idx.genomeReaderIndexes.Range(func(_, _ any) bool {
		t.Error("Close retained shared genome metadata")
		return false
	})
}

// Permissive E-values can retain negative bitscores, which must still sort by
// the best actual HSP score within a cluster and within a genome.
func TestSearchNegativeSimilarityScores(t *testing.T) {
	dir := t.TempDir()
	rng := rand.New(rand.NewSource(93))
	query := make([]byte, 1000)
	for i := range query {
		query[i] = "ACGT"[rng.Intn(4)]
	}
	var refs []string
	var combined strings.Builder
	for j, shape := range [][2]int{{44, 17}, {44, 18}, {44, 19}, {38, 16}, {38, 17}, {38, 18}} {
		subject := append([]byte(nil), query...)
		for i, base := range subject {
			if i%shape[0] >= shape[1] {
				subject[i] = "ACGT"[(strings.IndexByte("ACGT", base)+1+rng.Intn(3))%4]
			}
		}
		fasta := fmt.Sprintf(">s%d\n%s\n", j, subject)
		combined.WriteString(fasta)
		refs = append(refs, writeTestInput(t, filepath.Join(dir, fmt.Sprintf("g%d.fa", j)), fasta))
	}
	for _, sameGenome := range []bool{false, true} {
		db := filepath.Join(dir, fmt.Sprintf("db-%t", sameGenome))
		inputs := refs
		if sameGenome {
			inputs = []string{writeTestInput(t, filepath.Join(dir, "combined.fa"), combined.String())}
		}
		runSeedIndexCommand(t, append([]string{"index", "--quiet", "-j", "1", "-k", "16", "-m", "256", "-c", "1", "--partitions", "4", "-O", db}, inputs...)...)
		for _, reverse := range []bool{false, true} {
			t.Run(fmt.Sprintf("sameGenome=%t/reverse=%t", sameGenome, reverse), func(t *testing.T) {
				q := append([]byte(nil), query...)
				if reverse {
					RC(q)
				}
				file := writeTestInput(t, filepath.Join(dir, "q.fa"), ">q\n"+string(q)+"\n")
				out := runSeedIndexCommand(t, "search", "--quiet", "-j", "1", "-d", db, "-p", "11", "-P", "11", "-n", "0", "-N", "0", "-i", "60", "-e", "1e100", "--align-max-gap", "100", "-a", file)
				reader := csv.NewReader(strings.NewReader(string(out)))
				reader.Comma = '\t'
				rows, err := reader.ReadAll()
				if err != nil {
					t.Fatal(err)
				}
				var keys []string
				scores := make(map[string]float64)
				for _, row := range rows[1:] {
					key := row[3]
					if sameGenome {
						key += ":" + row[6] // HSP cluster within the genome
					}
					bitscore, err := strconv.Atoi(row[19])
					if err != nil {
						t.Fatal(err)
					}
					// Use unrounded identity from the alignment, like search itself.
					score := float64(bitscore) * float64(strings.Count(row[23], "|")) / float64(len(row[23])) * 100
					previous, ok := scores[key]
					if !ok {
						keys = append(keys, key)
					}
					if !ok || score > previous {
						scores[key] = score
					}
				}
				negative := 0
				for i, key := range keys {
					if scores[key] < 0 {
						negative++
					}
					if i > 0 && scores[keys[i-1]] < scores[key] {
						t.Errorf("score order increased: %s=%g, %s=%g", keys[i-1], scores[keys[i-1]], key, scores[key])
					}
				}
				if negative < 2 {
					t.Fatalf("fixture needs at least two negative-score groups, got %d", negative)
				}
			})
		}
	}
}

// A borderline DUST score must not remove an exact hit on either strand.
func TestSearchDustOrientation(t *testing.T) {
	dir := t.TempDir()
	sequence := []byte("AAAATCAAAATAAACAAAACCATTAAACAAA")
	var previous string
	for _, reverseTarget := range []bool{false, true} {
		ref := writeTestInput(t, filepath.Join(dir, "g.fa"), ">s\n"+string(sequence)+"\n")
		db := filepath.Join(dir, fmt.Sprintf("db-%t", reverseTarget))
		runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-k", "31", "-m", "64", "-c", "2", "--partitions", "4", "-O", db, ref)
		for _, reverseQuery := range []bool{false, true} {
			query := []byte("AAAATCAAAATAAACAAAACCATTAAACAAA")
			if reverseQuery {
				RC(query)
			}
			input := writeTestInput(t, filepath.Join(dir, "query.fa"), ">q\n"+string(query)+"\n")
			output := runSeedIndexCommand(t, "search", "--quiet", "-j", "1", "-d", db, "-p", "11", "-P", "11", "-n", "0", "-N", "0", "-l", "20", input)
			reader := csv.NewReader(strings.NewReader(string(output)))
			reader.Comma = '\t'
			rows, err := reader.ReadAll()
			if err != nil || len(rows) != 2 {
				t.Fatalf("targetRC=%t queryRC=%t: got %d rows, err=%v; want one hit", reverseTarget, reverseQuery, len(rows)-1, err)
			}
			row := rows[1]
			strand := "+"
			if reverseTarget != reverseQuery {
				strand = "-"
			}
			if row[9] != "31" || row[10] != "100.000" || row[12] != "1" || row[13] != "31" || row[14] != "1" || row[15] != "31" || row[16] != strand {
				t.Fatalf("wrong exact-hit coordinates or strand: %v", row)
			}
			row[16] = "" // only the relative strand changes
			result := strings.Join(row, "\t")
			if previous != "" && result != previous {
				t.Fatalf("orientation changed exact-hit fields: %s != %s", result, previous)
			}
			previous = result
		}
		RC(sequence)
	}
}

// A two-base extension must reach -l and preserve complete coordinates/CIGAR
// on both strands, even when the terminal k-mers were filtered as low complexity.
func TestSearchTwoBaseExtension(t *testing.T) {
	dir := t.TempDir()
	sequence := twoBaseExtensionSequence
	ref := writeTestInput(t, filepath.Join(dir, "s.fa"), ">s\n"+sequence+"\n")
	db := filepath.Join(dir, "db")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-k", "16", "-m", "64", "-c", "1", "--partitions", "4", "-O", db, ref)
	for _, reverse := range []bool{false, true} {
		q := []byte(sequence)
		strand := "+"
		if reverse {
			RC(q)
			strand = "-"
		}
		query := writeTestInput(t, filepath.Join(dir, "q.fa"), ">q\n"+string(q)+"\n")
		output := runSeedIndexCommand(t, "search", "--quiet", "-j", "1", "-d", db, "-p", "11", "-P", "11", "-n", "0", "-N", "0", "-l", "89", "-a", query)
		reader := csv.NewReader(strings.NewReader(string(output)))
		reader.Comma = '\t'
		rows, err := reader.ReadAll()
		if err != nil || len(rows) != 2 {
			t.Errorf("reverse=%t: got %d rows, err=%v; want one hit", reverse, len(rows)-1, err)
			continue
		}
		row := rows[1]
		if row[5] != "100.000" || row[9] != "89" || row[10] != "100.000" || row[12] != "1" || row[13] != "89" || row[14] != "1" || row[15] != "89" || row[16] != strand || row[20] != "89M" {
			t.Errorf("reverse=%t: wrong full-length alignment: %v", reverse, row[:21])
		}
	}
}

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

// terminalMismatchSequences leaves a 979-base exact core in a 1000-base HSP.
func terminalMismatchSequences() ([]byte, []byte) {
	q := []byte(indexTestDNA(1000))
	s := append([]byte(nil), q...)
	for _, pos := range []int{10, 990} {
		s[pos] = "ACGT"[(strings.IndexByte("ACGT", s[pos])+1)%4]
	}
	return q, s
}

// TestSearchExtendedAlignmentLength applies -l after extension on both strands.
func TestSearchExtendedAlignmentLength(t *testing.T) {
	dir := t.TempDir()
	q, s := terminalMismatchSequences()
	ref := writeTestInput(t, filepath.Join(dir, "ref.fa"), ">s\n"+string(s)+"\n")
	forward := string(q)
	RC(q)
	query := writeTestInput(t, filepath.Join(dir, "query.fa"), ">forward\n"+forward+"\n>reverse\n"+string(q)+"\n")
	db := filepath.Join(dir, "db")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-c", "2", "--partitions", "4", "-O", db, ref)
	for _, minimum := range []int{50, 980, 1000, 1001} {
		for _, all := range []bool{false, true} {
			args := []string{"search", "--quiet", "-j", "1", "-d", db, "-n", "0", "-N", "0", "-l", strconv.Itoa(minimum), query}
			if all {
				args = append(args, "-a")
			}
			output := runSeedIndexCommand(t, args...)
			lines := strings.Split(strings.TrimSpace(string(output)), "\n")
			want := 3
			if minimum > 1000 {
				want = 1
			}
			if len(lines) != want {
				t.Fatalf("minimum=%d/all=%t: got %d lines, want %d", minimum, all, len(lines), want)
			}
			for _, line := range lines[1:] {
				row := strings.Split(line, "\t")
				strand := "+"
				if row[0] == "reverse" {
					strand = "-"
				}
				if strings.Join(row[9:17], ",") != "1000,99.800,0,1,1000,1,1000,"+strand {
					t.Fatalf("unexpected HSP: %v", row[:20])
				}
			}
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
	// 85% lies between the refined 354/417 and unrefined 372/417 query spans.
	for _, coverage := range []struct{ hsp, genome, gap string }{{"0", "0", "20"}, {"0", "0", "100"}, {"85", "0", "100"}, {"0", "85", "100"}} {
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
			for field, minimum := range map[int]string{5: coverage.genome, 8: coverage.hsp} {
				value, err := strconv.ParseFloat(row[field], 64)
				if err != nil || (minimum == "85" && value < 85) {
					t.Fatalf("coverage filter violated: %v", row[:20])
				}
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
			if coverage.hsp == "85" || coverage.genome == "85" {
				// Both query orientations retain the longer window when the
				// refined 354-base alignment would fail the coverage threshold.
				wantBegin, wantEnd, wantCIGAR = 28, 399, "11M3X9M9I12M1X11M9I307M"
				if row[0] == "rc-query" {
					wantBegin, wantEnd, wantCIGAR = 19, 390, "313M9I5M1X18M3X9M9I5M"
				}
				if qs == wantBegin && qe == wantEnd && number(9) == 372 && number(11) == 18 && row[20] == wantCIGAR {
					seen[row[0]+"/"+row[4]] = true
				}
			}
		}
		if len(seen) != 4 {
			t.Fatalf("coverage=%+v: missing expected results: %s", coverage, output)
		}
	}
}

// Short queries and targets in a small-k index must not abort other queries.
func TestSearchShortSequences(t *testing.T) {
	for _, targetLen := range []int{30, 100} {
		t.Run(strconv.Itoa(targetLen), func(t *testing.T) {
			dir := t.TempDir()
			sequence := indexTestDNA(100)
			ref := writeTestInput(t, filepath.Join(dir, "ref.fa"), ">target\n"+sequence[:targetLen]+"\n")
			var queries strings.Builder
			for _, length := range []int{17, 30, 31, 100} {
				fmt.Fprintf(&queries, ">q%d\n%s\n", length, sequence[:length])
			}
			query := writeTestInput(t, filepath.Join(dir, "query.fa"), queries.String())
			db := filepath.Join(dir, "db")
			runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-k", "17", "-m", "64", "-c", "2", "--partitions", "4", "-O", db, ref)
			out := runSeedIndexCommand(t, "search", "--quiet", "-j", "2", "-d", db, "-n", "0", "-N", "0", "-l", "17", query)
			lines := strings.Split(strings.TrimSpace(string(out)), "\n")
			want := 1
			if targetLen == 100 {
				want = 3 // 31 and 100 bp queries remain searchable.
			}
			if len(lines) != want || strings.Contains(string(out), "\nq17\t") || strings.Contains(string(out), "\nq30\t") {
				t.Fatalf("unexpected results: %s", out)
			}
		})
	}
}

// Distinct seed chains can converge after extension; retain each final HSP once.
func TestSearchFinalAlignmentDeduplication(t *testing.T) {
	dir := t.TempDir()
	data, err := os.ReadFile("testdata/final-alignment-duplicates.json")
	if err != nil {
		t.Fatal(err)
	}
	var fixture struct{ Unit, Other, Left, Right string }
	if err := json.Unmarshal(data, &fixture); err != nil {
		t.Fatal(err)
	}
	target := strings.Repeat(fixture.Unit, 3)
	forward := fixture.Left + target + fixture.Right
	reverse := []byte(forward)
	RC(reverse)
	// A second synthetic contig competes for masks, leaving sparse seed chains.
	ref := writeTestInput(t, filepath.Join(dir, "ref.fa"), ">target\n"+target+"\n>other\n"+fixture.Other+"\n")
	query := writeTestInput(t, filepath.Join(dir, "query.fa"), ">forward\n"+forward+"\n>reverse\n"+string(reverse)+"\n")
	db := filepath.Join(dir, "db")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "1024", "-c", "2", "--partitions", "4", "-O", db, ref)
	for _, all := range []bool{false, true} {
		args := []string{"search", "--quiet", "-j", "1", "-d", db, "-n", "0", "-N", "0", query}
		if all {
			args = append(args, "-a")
		}
		out := runSeedIndexCommand(t, args...)
		seen := make(map[string]bool)
		full := make(map[string]int)
		for _, line := range strings.Split(strings.TrimSpace(string(out)), "\n")[1:] {
			row := strings.Split(line, "\t")
			key := strings.Join(append(append([]string{}, row[:6]...), row[8:]...), "\t")
			if seen[key] {
				t.Fatalf("duplicate final HSP (-a=%t): %v", all, row[:20])
			}
			seen[key] = true
			if row[9] == "1500" && row[14] == "1" && row[15] == "1500" {
				full[row[0]]++
			}
		}
		if full["forward"] != 1 || full["reverse"] != 1 {
			t.Fatalf("missing full-length controls (-a=%t): %v", all, full)
		}
	}
}

// Exact CIGAR comparison keeps alternatives and distinct loci, and owns its keys
// after the caller reuses both the scratch buffer and the WFA operation slice.
func TestDuplicateFinalAlignmentIdentity(t *testing.T) {
	seen := make(map[finalAlignmentKey]struct{})
	var scratch []byte
	c := &Chain2Result{QEnd: 99, TEnd: 99, AlignedLength: 100, MatchedBases: 98, Gaps: 2, Score: 180, BitScore: 100}
	cigar := &wfa.AlignmentResult{Ops: []uint64{OpM<<32 | 50, OpI<<32 | 1, OpD<<32 | 1, OpM<<32 | 49}}
	if duplicateFinalAlignment(seen, &scratch, c, cigar, 0, false) || !duplicateFinalAlignment(seen, &scratch, c, cigar, 0, false) {
		t.Fatal("exact duplicate was not removed")
	}
	if duplicateFinalAlignment(seen, &scratch, c, cigar, 1, false) || duplicateFinalAlignment(seen, &scratch, c, cigar, 0, true) {
		t.Fatal("distinct contig or strand was removed")
	}
	c.TBegin++
	c.TEnd++
	if duplicateFinalAlignment(seen, &scratch, c, cigar, 0, false) {
		t.Fatal("distinct locus was removed")
	}
	c.TBegin--
	c.TEnd--
	cigar.Ops[0] = OpM<<32 | 49
	cigar.Ops[3] = OpM<<32 | 50
	if duplicateFinalAlignment(seen, &scratch, c, cigar, 0, false) {
		t.Fatal("different CIGAR was removed")
	}
	cigar.Ops[0] = OpM<<32 | 50
	cigar.Ops[3] = OpM<<32 | 49
	if !duplicateFinalAlignment(seen, &scratch, c, cigar, 0, false) {
		t.Fatal("operation buffer reuse changed a retained key")
	}
}

// Reusing deduplication scratch must not hide a later reference's identical HSP.
func TestFinalAlignmentScratchReset(t *testing.T) {
	s := poolFinalAlignmentScratch.Get().(*finalAlignmentScratch)
	defer func() { s.reset(); poolFinalAlignmentScratch.Put(s) }()
	if len(s.keys) != 0 || len(s.ops) != 0 {
		t.Fatal("borrowed dirty scratch")
	}
	c := &Chain2Result{QEnd: 99, TEnd: 99, AlignedLength: 100, MatchedBases: 100}
	cigar := &wfa.AlignmentResult{Ops: []uint64{OpM<<32 | 100}}
	if duplicateFinalAlignment(s.keys, &s.ops, c, cigar, 0, false) || !duplicateFinalAlignment(s.keys, &s.ops, c, cigar, 0, false) {
		t.Fatal("deduplication failed before reset")
	}
	capacity := cap(s.ops)
	s.reset()
	if len(s.keys) != 0 || len(s.ops) != 0 || cap(s.ops) != capacity {
		t.Fatal("reset did not clear keys and retain scratch capacity")
	}
	if duplicateFinalAlignment(s.keys, &s.ops, c, cigar, 0, false) {
		t.Fatal("reset retained an HSP from the preceding reference")
	}
}

// Workers must exclusively own their scratch until it is cleared and returned.
func TestFinalAlignmentScratchParallel(t *testing.T) {
	for worker := 0; worker < 8; worker++ {
		t.Run(strconv.Itoa(worker), func(t *testing.T) {
			t.Parallel()
			c := &Chain2Result{QEnd: 99, TEnd: 99, AlignedLength: 100, MatchedBases: 100}
			cigar := &wfa.AlignmentResult{Ops: []uint64{OpM<<32 | 100}}
			for i := 0; i < 100; i++ {
				s := poolFinalAlignmentScratch.Get().(*finalAlignmentScratch)
				dirty := len(s.keys) != 0 || len(s.ops) != 0
				first := duplicateFinalAlignment(s.keys, &s.ops, c, cigar, 0, false)
				second := duplicateFinalAlignment(s.keys, &s.ops, c, cigar, 0, false)
				s.reset()
				poolFinalAlignmentScratch.Put(s)
				if dirty || first || !second {
					t.Fatal("deduplication state crossed worker lifetimes")
				}
			}
		})
	}
}
