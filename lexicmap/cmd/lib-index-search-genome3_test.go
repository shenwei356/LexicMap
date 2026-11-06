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
	"errors"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"sync"
	"testing"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
	"github.com/shenwei356/wfa"
)

// Contained variable-length anchors must not select different alignment windows
// after reversing both sequences, which could change final identity filtering.
func TestAlignChainVariableAnchorContainment(t *testing.T) {
	idx := compareTestIndex(t)
	idx.seqCompareOption.MaxGap = 20
	cpr := NewSeqComparator(idx.seqCompareOption, idx.poolChainers2)
	algn := wfa.New(wfa.DefaultPenalties, &wfa.Options{GlobalAlignment: true})
	defer wfa.RecycleAligner(algn)
	query := []byte("ACTTCAGAACATCTGCGTTTACGGAATGCATCTTAGAAAGAATGGGTGTTTCACTACCAGTGATTCGGCCATCCGCGGCCTCAGGGAAG")
	target := []byte("ACTTCAGAACATCACTTCAGAACATCTGCGTTTACGGATGCGTTTACGGCATGCATCTCAGAAAGGATGGGTGTTCCACTACCAGTGATCCAGCAACCCGCGGCCTCCGGGCAG")
	for _, reverse := range []bool{false, true} {
		if err := cpr.Index(query); err != nil {
			t.Fatal(err)
		}
		score := scoreAndEvalue(2, -3, 5, 2, len(target), 0.625, 0.41)
		chain := &Chain2Result{QEnd: len(query) - 1, TEnd: len(target) - 1}
		matched, aligned, gaps, _, _, ok := alignChain(query, nil, target, chain,
			&subjectSketch{seqLen: len(target)}, algn, cpr, 0, 50, 70, 30, idx, &score)
		cpr.RecycleIndex()
		if !ok || [3]int{matched, aligned, gaps} != [3]int{75, 94, 12} {
			t.Errorf("reverse=%t: accepted=%t, matches/aligned/gaps=%d/%d/%d; want true, 75/94/12",
				reverse, ok, matched, aligned, gaps)
		}
		RC(query)
		RC(target)
	}
}

// The coverage bound must allow terminal extension, equality, and later chains.
func TestChainsCanReachQueryCoverage(t *testing.T) {
	short := &Chain2Result{QBegin: 100, QEnd: 139, TBegin: 100, TEnd: 139}
	long := &Chain2Result{QBegin: 100, QEnd: 399, TBegin: 100, TEnd: 399}
	// Evaluate the percentage in the same order as the final coverage filter.
	maxPaired, queryLen := short.QEnd-short.QBegin+101, 1000
	boundary := float64(maxPaired) / float64(queryLen) * 100
	for _, tc := range []struct {
		name   string
		chains []*Chain2Result
		qlen   int
		tlen   int
		ext    int
		qcov   float64
		want   bool
	}{
		{"short", []*Chain2Result{short}, 1000, 1000, 50, 30, false},
		{"extended-equality", []*Chain2Result{short}, 1000, 1000, 50, boundary, true},
		{"above-bound", []*Chain2Result{short}, 1000, 1000, 50, math.Nextafter(boundary, 100), false},
		{"query-edge", []*Chain2Result{{QBegin: 0, QEnd: 39, TBegin: 100, TEnd: 139}}, 1000, 1000, 50, 10, false},
		{"target-edge", []*Chain2Result{short}, 1000, 140, 50, 10, false},
		{"later-chain", []*Chain2Result{short, long}, 1000, 1000, 50, 30, true},
		{"disabled", []*Chain2Result{short}, 1000, 1000, 0, 0, true},
		{"unlimited-extension", []*Chain2Result{short}, 1000, 1000, math.MaxInt, 100, true},
	} {
		t.Run(tc.name, func(t *testing.T) {
			if got := chainsCanReachQueryCoverage(tc.chains, tc.qlen, tc.tlen, tc.ext, tc.qcov); got != tc.want {
				t.Fatalf("can reach coverage=%t, want %t", got, tc.want)
			}
		})
	}
}

// Disabling the early bound and applying coverage to the completed alignment
// is an independent oracle for rejection, statistics, and reverse complements.
func TestAlignChainCoverageMatchesFinalFilter(t *testing.T) {
	idx := compareTestIndex(t)
	idx.opt.MaxEvalue = 10
	q := []byte(indexTestDNA(500))
	for _, length := range []int{40, 150, 300, 500} {
		for _, gapped := range []bool{false, true} {
			for _, reverse := range []bool{false, true} {
				query := append([]byte(nil), q...)
				target := append([]byte(nil), q[(500-length)/2:(500+length)/2]...)
				if gapped {
					mid := len(target) / 2
					inserted := make([]byte, len(target)+7)
					copy(inserted, target[:mid])
					copy(inserted[mid:], "GATTACA")
					copy(inserted[mid+7:], target[mid:])
					target = inserted
				}
				if reverse {
					RC(query)
					RC(target)
				}
				cpr := NewSeqComparator(idx.seqCompareOption, idx.poolChainers2)
				if err := cpr.Index(query); err != nil {
					t.Fatal(err)
				}
				algn := wfa.New(wfa.DefaultPenalties, &wfa.Options{GlobalAlignment: true})
				score := scoreAndEvalue(2, -3, 5, 2, len(target), 0.625, 0.41)
				chain := &Chain2Result{QBegin: 240, QEnd: 259, TBegin: len(target)/2 - 10, TEnd: len(target)/2 + 9}
				for _, ext := range []int{0, 2, 50} {
					wm, wa, wg, wp, wpa, wok := alignChain(query, nil, target, chain, &subjectSketch{seqLen: len(target)},
						algn, cpr, 500, ext, 70, 0, idx, &score)
					coverage := float64(wa-wg) / float64(len(query)) * 100
					for _, minimum := range []float64{0, 30, 100, coverage, math.Nextafter(coverage, 101)} {
						m, a, g, pi, pia, ok := alignChain(query, nil, target, chain, &subjectSketch{seqLen: len(target)},
							algn, cpr, 500, ext, 70, minimum, idx, &score)
						want := wok && coverage >= minimum
						if ok != want || ok && (m != wm || a != wa || g != wg || pi != wp || pia != wpa) {
							t.Fatalf("length=%d gapped=%t reverse=%t ext=%d coverage=%g: got %t/%d/%d/%d/%g/%g; want %t/%d/%d/%d/%g/%g",
								length, gapped, reverse, ext, minimum, ok, m, a, g, pi, pia, want, wm, wa, wg, wp, wpa)
						}
					}
				}
				wfa.RecycleAligner(algn)
				cpr.RecycleIndex()
			}
		}
	}
}

// OrthoANI must retain the final two matching bases before applying AF filters.
func TestGenomeSearchTwoBaseExtension(t *testing.T) {
	dir := t.TempDir()
	sequence := indexTestDNA(1000) + twoBaseExtensionSequence
	ref := writeTestInput(t, filepath.Join(dir, "s.fa"), ">s\n"+sequence+"\n")
	db := filepath.Join(dir, "db")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-k", "16", "-m", "64", "-c", "1", "--partitions", "4", "-O", db, ref)
	for _, reverse := range []bool{false, true} {
		q := []byte(sequence)
		if reverse {
			RC(q)
		}
		query := writeTestInput(t, filepath.Join(dir, "q.fa"), ">q\n"+string(q)+"\n")
		out := runSeedIndexCommand(t, "genome", "search", "--quiet", "-j", "1", "-d", db, "-p", "11", "-m", "0", "-n", "0", "-O", "--min-af", "99.9", query)
		lines := strings.Split(strings.TrimSpace(string(out)), "\n")
		if len(lines) != 2 {
			t.Errorf("reverse=%t: got %d hits, want one", reverse, len(lines)-1)
			continue
		}
		row := strings.Split(lines[1], "\t")
		if row[2] != "100.0000" || row[3] != "100.000" || row[4] != "100.000" {
			t.Errorf("reverse=%t: incomplete self-alignment: %v", reverse, row)
		}
	}
}

// Genome search must apply -l after extending its sampled and dense chains.
func TestGenomeSearchExtendedAlignmentLength(t *testing.T) {
	dir := t.TempDir()
	q, s := terminalMismatchSequences()
	ref := writeTestInput(t, filepath.Join(dir, "ref.fa"), ">s\n"+string(s)+"\n")
	db := filepath.Join(dir, "db")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-c", "2", "--partitions", "4", "-O", db, ref)
	for _, reverse := range []bool{false, true} {
		if reverse {
			RC(q)
		}
		query := writeTestInput(t, filepath.Join(dir, "query.fa"), ">q\n"+string(q)+"\n")
		for _, ortho := range []bool{false, true} {
			for _, minimum := range []int{50, 980, 1000, 1001} {
				args := []string{"genome", "search", "--quiet", "-j", "2", "-d", db, "-m", "0", "--frag-size", "1000", "-l", strconv.Itoa(minimum), query}
				if ortho {
					args = append(args, "--OrthoANI")
				}
				out := runSeedIndexCommand(t, args...)
				lines := strings.Split(strings.TrimSpace(string(out)), "\n")
				want := 2
				if minimum > 1000 {
					want = 1
				}
				if len(lines) != want {
					t.Fatalf("reverse=%t/ortho=%t/minimum=%d: got %d lines, want %d", reverse, ortho, minimum, len(lines), want)
				}
			}
		}
	}
}

// A chain's inclusive end must include its last base in the comparator window.
func TestAlignChainQueryWindowEnd(t *testing.T) {
	for _, reverse := range []bool{false, true} {
		t.Run(fmt.Sprintf("reverse=%t", reverse), func(t *testing.T) {
			idx := compareTestIndex(t)
			idx.opt.MaxEvalue = 10
			q := []byte(indexTestDNA(31))
			if reverse {
				RC(q)
			}
			cpr := NewSeqComparator(idx.seqCompareOption, &sync.Pool{New: func() interface{} {
				return NewChainer2(&idx.seqCompareOption.Chaining2Options)
			}})
			if err := cpr.Index(q); err != nil {
				t.Fatal(err)
			}
			defer cpr.RecycleIndex()
			algn := wfa.New(wfa.DefaultPenalties, &wfa.Options{GlobalAlignment: true})
			defer wfa.RecycleAligner(algn)
			chain := &Chain2Result{QBegin: 0, QEnd: 30, TBegin: 0, TEnd: 30}
			score := scoreAndEvalue(2, -3, 5, 2, len(q), 0.625, 0.41)
			matched, aligned, gaps, _, _, ok := alignChain(q, nil, q, chain, &subjectSketch{seqLen: len(q)},
				algn, cpr, 0, 0, 100, 100, idx, &score)
			if !ok || matched != len(q) || aligned != len(q) || gaps != 0 {
				t.Fatalf("terminal match lost: ok=%t matched/aligned/gaps=%d/%d/%d", ok, matched, aligned, gaps)
			}
		})
	}
}

// TestAlignVotedWindowMinLength checks the repetitive-window fallback,
// including the inclusive boundary, reverse complements, and pooled reuse.
func TestAlignVotedWindowMinLength(t *testing.T) {
	for _, minimum := range []int{30, 50} {
		for _, length := range []int{minimum - 1, minimum, minimum + 1} {
			for _, reverse := range []bool{false, true} {
				t.Run(fmt.Sprintf("minimum=%d/length=%d/reverse=%t", minimum, length, reverse), func(t *testing.T) {
					idx := compareTestIndex(t)
					idx.seqCompareOption.MinAlignLen = minimum
					// A permissive E-value isolates the final length filter.
					idx.opt.MaxEvalue = 10
					q := []byte(indexTestDNA(length))
					if reverse {
						RC(q)
					}
					target := append([]byte(nil), q...)
					algn := wfa.New(wfa.DefaultPenalties, &wfa.Options{GlobalAlignment: true})
					defer wfa.RecycleAligner(algn)
					sketch := &subjectSketch{seqLen: length}
					score := scoreAndEvalue(2, -3, 5, 2, length, 0.625, 0.41)
					subs := []SubstrPair{{QBegin: 0, TBegin: 0, Len: 11}, {QBegin: 4, TBegin: 4, Len: 11}, {QBegin: 8, TBegin: 8, Len: 11}}
					for trial := 0; trial < 2; trial++ {
						matched, aligned, gaps, _, _, ok := alignVotedWindow(q, nil, &subs, sketch, target, algn, 50, 3, 90, 70, 30, idx, &score)
						if ok != (length >= minimum) {
							t.Errorf("trial=%d: accepted=%t, aligned=%d, minimum=%d", trial, ok, aligned, minimum)
						}
						if ok && (matched != length || aligned != length || gaps != 0) {
							t.Errorf("trial=%d: got matched/aligned/gaps=%d/%d/%d, want %d/%d/0", trial, matched, aligned, gaps, length, length)
						}
					}
				})
			}
		}
	}
}

func TestAddSampledKmerPosition(t *testing.T) {
	kmers := make(map[uint64]uint32)
	repeated := make(map[uint64]repeatedKmerList)
	positions := make([]repeatedKmerPosition, 0)

	addSampledKmerPosition(&kmers, &repeated, &positions, 42, 10)
	addSampledKmerPosition(&kmers, &repeated, &positions, 42, 20)
	addSampledKmerPosition(&kmers, &repeated, &positions, 42, 30)

	if got := kmers[42]; got != 11 {
		t.Fatalf("first encoded position: got %d, want 11", got)
	}

	index := repeated[42].head
	want := []uint32{20, 30}
	for i, expected := range want {
		if index == 0 {
			t.Fatalf("position %d is missing", i)
		}
		position := positions[index-1]
		if position.position != expected {
			t.Fatalf("position %d: got %d, want %d", i, position.position, expected)
		}
		index = position.next
	}
	if index != 0 {
		t.Fatalf("unexpected additional position at index %d", index)
	}
}

func TestSampledKmerMapCapacity(t *testing.T) {
	if got := sampledKmerMapCapacity(176, 4); got != 44 {
		t.Fatalf("capacity: got %d, want 44", got)
	}
}

func TestSubjectContigBoundsForwardAndReverse(t *testing.T) {
	sketch := &subjectSketch{
		seqLen:       260,
		forwardLen:   100,
		rcStart:      160,
		contigBounds: [][2]int{{0, 40}, {60, 100}},
	}

	if start, end := subjectContigBounds(sketch, 20); start != 0 || end != 40 {
		t.Fatalf("forward bounds: got [%d,%d), want [0,40)", start, end)
	}
	if start, end := subjectContigBounds(sketch, 180); start != 160 || end != 200 {
		t.Fatalf("reverse bounds: got [%d,%d), want [160,200)", start, end)
	}
	if start, end := subjectContigBounds(sketch, 240); start != 220 || end != 260 {
		t.Fatalf("reverse bounds: got [%d,%d), want [220,260)", start, end)
	}

	linear := func(position int) (int, int) {
		if position >= sketch.rcStart {
			for i := len(sketch.contigBounds) - 1; i >= 0; i-- {
				start := sketch.rcStart + sketch.forwardLen - sketch.contigBounds[i][1]
				end := sketch.rcStart + sketch.forwardLen - sketch.contigBounds[i][0]
				if position >= start && position < end {
					return start, end
				}
			}
		} else {
			for _, bounds := range sketch.contigBounds {
				if position >= bounds[0] && position < bounds[1] {
					return bounds[0], bounds[1]
				}
			}
		}
		return 0, sketch.seqLen
	}
	for position := 0; position < sketch.seqLen; position++ {
		wantStart, wantEnd := linear(position)
		start, end := subjectContigBounds(sketch, position)
		if start != wantStart || end != wantEnd {
			t.Fatalf("position %d: got [%d,%d), want [%d,%d)", position, start, end, wantStart, wantEnd)
		}
	}
}

func TestReadGenomeRejectsAllChunksWhenCombinedSizeExceedsLimit(t *testing.T) {
	dbDir := t.TempDir()
	dir := filepath.Join(dbDir, DirGenomes, batchDir(0))
	if err := os.MkdirAll(dir, 0o755); err != nil {
		t.Fatal(err)
	}
	wtr, err := genome.NewWriter(filepath.Join(dir, FileGenomes), 0)
	if err != nil {
		t.Fatal(err)
	}
	for i, sequence := range []string{"ACGTACGT", "TGCATGCA"} {
		g := genome.PoolGenome.Get().(*genome.Genome)
		g.Reset()
		g.ID = append(g.ID, byte('0'+i))
		g.Seq = append(g.Seq, sequence...)
		g.GenomeSize = len(g.Seq)
		g.Len = len(g.Seq)
		g.NumSeqs = 1
		g.SeqSizes = append(g.SeqSizes, len(g.Seq))
		seqID := []byte("seq")
		g.SeqIDs = append(g.SeqIDs, &seqID)
		if err = wtr.Write(g); err != nil {
			genome.RecycleGenome(g)
			t.Fatal(err)
		}
		genome.RecycleGenome(g)
	}
	if err = wtr.Close(); err != nil {
		t.Fatal(err)
	}

	idx := &Index{
		path:           dbDir,
		opt:            &IndexSearchingOptions{MaxSubjectGenomeSize: 10},
		info:           &IndexInfo{GenomeBatches: 1},
		openFileTokens: make(chan int, 1),
	}
	chunks := []uint64{0, 1}
	q, err := idx.ReadGenome(&chunks, "chunked")
	if q != nil || !errors.Is(err, errGenomeTooLarge) {
		t.Fatalf("oversized genome: query=%v err=%v", q, err)
	}

	firstChunk := []uint64{0}
	q, err = idx.ReadGenome(&firstChunk, "chunked")
	if err != nil {
		t.Fatal(err)
	}
	if q.genomeSize != 8 || len(q.seqs) != 1 {
		t.Fatalf("reused query contains partial oversized data: size=%d contigs=%d", q.genomeSize, len(q.seqs))
	}
	RecycleGQuery(q)
}
