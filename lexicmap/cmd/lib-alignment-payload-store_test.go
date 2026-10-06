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
	"testing"
)

func newPayloadTestResult(qseq, tseq string) *SearchResult {
	alignments := []AlignmentResult{{
		CIGAR: []byte("2M"), QSeq: []byte(qseq),
		TSeq: []byte(tseq), Alignment: []byte("|."),
	}}
	sd := &SimilarityDetail{Alignments: alignments}
	details := []*SimilarityDetail{sd}
	return &SearchResult{SimilarityDetails: &details}
}

func TestAlignmentPayloadStoreStaysInMemoryBelowGlobalBudget(t *testing.T) {
	budget := newAlignmentPayloadBudget(1024)
	budget.tempDir = t.TempDir()
	store := newAlignmentPayloadStore(budget)
	r := newPayloadTestResult("AC", "AG")
	previous := []*SearchResult{}
	wantBytes := searchResultPayloadBytes(r)

	if err := store.retain(r, &previous); err != nil {
		t.Fatal(err)
	}
	a := &(*r.SimilarityDetails)[0].Alignments[0]
	if store.spilled || store.file != nil || store.refs != nil {
		t.Fatal("payload below the budget was spilled")
	}
	if string(a.CIGAR) != "2M" || string(a.QSeq) != "AC" ||
		string(a.TSeq) != "AG" || string(a.Alignment) != "|." {
		t.Fatal("in-memory payload changed")
	}
	if budget.usedBytes() != wantBytes {
		t.Fatalf("budget usage: got %d, want %d", budget.usedBytes(), wantBytes)
	}
	if err := store.close(); err != nil {
		t.Fatal(err)
	}
	if budget.usedBytes() != 0 {
		t.Fatalf("closed store retained %d budget bytes", budget.usedBytes())
	}
}

func TestAlignmentPayloadStoreSpillsEarlierAndCurrentResults(t *testing.T) {
	r1 := newPayloadTestResult("AC", "AG")
	r2 := newPayloadTestResult("GT", "GC")
	oneResultBytes := searchResultPayloadBytes(r1)
	budget := newAlignmentPayloadBudget(oneResultBytes)
	budget.tempDir = t.TempDir()
	store := newAlignmentPayloadStore(budget)
	previous := []*SearchResult{}

	if err := store.retain(r1, &previous); err != nil {
		t.Fatal(err)
	}
	previous = append(previous, r1)
	if err := store.retain(r2, &previous); err != nil {
		t.Fatal(err)
	}
	if !store.spilled || store.file == nil {
		t.Fatal("budget overflow did not switch the query to spill mode")
	}
	if budget.usedBytes() != 0 {
		t.Fatalf("spilled query retained %d budget bytes", budget.usedBytes())
	}

	spillName := store.file.Name()
	for _, r := range []*SearchResult{r1, r2} {
		a := &(*r.SimilarityDetails)[0].Alignments[0]
		if _, ok := store.refs[a]; !ok || a.CIGAR != nil || a.QSeq != nil || a.TSeq != nil || a.Alignment != nil {
			t.Fatal("spilled alignment retained its in-memory payload")
		}
	}
	if err := store.finalize(); err != nil {
		t.Fatal(err)
	}

	var scratch []byte
	for i, test := range []struct {
		r          *SearchResult
		qseq, tseq string
	}{{r1, "AC", "AG"}, {r2, "GT", "GC"}} {
		a := &(*test.r.SimilarityDetails)[0].Alignments[0]
		cigar, qseq, tseq, alignment, nextScratch, err := store.payload(a, scratch)
		if err != nil {
			t.Fatal(err)
		}
		scratch = nextScratch
		if string(cigar) != "2M" || string(qseq) != test.qseq ||
			string(tseq) != test.tseq || string(alignment) != "|." {
			t.Fatalf("payload %d changed after spill: %q %q %q %q", i, cigar, qseq, tseq, alignment)
		}
	}

	if err := store.close(); err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(spillName); !os.IsNotExist(err) {
		t.Fatalf("spill file was not removed: %s", spillName)
	}
}

func TestAlignmentPayloadBudgetIsSharedAcrossQueries(t *testing.T) {
	r1 := newPayloadTestResult("AC", "AG")
	r2 := newPayloadTestResult("GT", "GC")
	oneResultBytes := searchResultPayloadBytes(r1)
	budget := newAlignmentPayloadBudget(oneResultBytes)
	budget.tempDir = t.TempDir()
	store1 := newAlignmentPayloadStore(budget)
	store2 := newAlignmentPayloadStore(budget)
	empty := []*SearchResult{}

	if err := store1.retain(r1, &empty); err != nil {
		t.Fatal(err)
	}
	if err := store2.retain(r2, &empty); err != nil {
		t.Fatal(err)
	}
	if store1.spilled || !store2.spilled {
		t.Fatalf("global budget was not shared: store1 spilled=%v, store2 spilled=%v", store1.spilled, store2.spilled)
	}
	if budget.usedBytes() != oneResultBytes {
		t.Fatalf("global usage: got %d, want %d", budget.usedBytes(), oneResultBytes)
	}
	if err := store1.close(); err != nil {
		t.Fatal(err)
	}
	if err := store2.close(); err != nil {
		t.Fatal(err)
	}
}

func TestAlignmentPayloadStoreReconcilesFilteredResults(t *testing.T) {
	r1 := newPayloadTestResult("AC", "AG")
	r2 := newPayloadTestResult("GT", "GC")
	oneResultBytes := searchResultPayloadBytes(r1)
	budget := newAlignmentPayloadBudget(oneResultBytes * 2)
	store := newAlignmentPayloadStore(budget)
	previous := []*SearchResult{}
	if err := store.retain(r1, &previous); err != nil {
		t.Fatal(err)
	}
	previous = append(previous, r1)
	if err := store.retain(r2, &previous); err != nil {
		t.Fatal(err)
	}

	live := []*SearchResult{r1}
	store.reconcile(&live)
	if budget.usedBytes() != oneResultBytes {
		t.Fatalf("filtered payload remained charged: got %d, want %d", budget.usedBytes(), oneResultBytes)
	}
	if err := store.close(); err != nil {
		t.Fatal(err)
	}
}

func TestAlignmentPayloadStoreDropsFilteredSpillReferences(t *testing.T) {
	r1 := newPayloadTestResult("AC", "AG")
	r2 := newPayloadTestResult("GT", "GC")
	budget := newAlignmentPayloadBudget(searchResultPayloadBytes(r1))
	budget.tempDir = t.TempDir()
	store := newAlignmentPayloadStore(budget)
	previous := []*SearchResult{}
	if err := store.retain(r1, &previous); err != nil {
		t.Fatal(err)
	}
	previous = append(previous, r1)
	if err := store.retain(r2, &previous); err != nil {
		t.Fatal(err)
	}

	a1 := &(*r1.SimilarityDetails)[0].Alignments[0]
	a2 := &(*r2.SimilarityDetails)[0].Alignments[0]
	live := []*SearchResult{r1}
	store.reconcile(&live)
	if len(store.refs) != 1 {
		t.Fatalf("got %d spill references after filtering, want 1", len(store.refs))
	}
	if _, ok := store.refs[a1]; !ok {
		t.Fatal("live spill reference was removed")
	}
	if _, ok := store.refs[a2]; ok {
		t.Fatal("filtered spill reference was retained")
	}
	if err := store.finalize(); err != nil {
		t.Fatal(err)
	}
	_, qseq, _, _, _, err := store.payload(a1, nil)
	if err != nil || string(qseq) != "AC" {
		t.Fatalf("live payload was not readable after reconciliation: %q, %v", qseq, err)
	}
	if err := store.close(); err != nil {
		t.Fatal(err)
	}
}
