// Copyright © 2026 Wei Shen <shenwei356@gmail.com>

package cmd

import (
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/shenwei356/wfa"
)

func TestSeqs2FragmentsWithQual(t *testing.T) {
	seq1, seq2 := []byte("AAAACCCCGG"), []byte("TTTTAAA")
	qual1, qual2 := []byte("1234567890"), []byte("ABCDEFG")
	seqs := []*[]byte{&seq1, &seq2}
	quals := []*[]byte{&qual1, &qual2}

	frags, qualFrags, n, err := seqs2fragmentsWithQual(&seqs, &quals, 4, 3)
	if err != nil {
		t.Fatal(err)
	}
	defer recycleFragments(frags)
	defer recycleFragments(qualFrags)

	if n != 15 {
		t.Fatalf("fragment bases: got %d, want 15", n)
	}
	wantSeqs := []string{"AAAA", "CCCC", "TTTT", "AAA"}
	wantQuals := []string{"1234", "5678", "ABCD", "EFG"}
	if len(*frags) != len(wantSeqs) || len(*qualFrags) != len(wantQuals) {
		t.Fatalf("fragment counts: got %d sequences and %d qualities", len(*frags), len(*qualFrags))
	}
	for i := range wantSeqs {
		if got := string((*frags)[i]); got != wantSeqs[i] {
			t.Errorf("sequence fragment %d: got %q, want %q", i, got, wantSeqs[i])
		}
		if got := string((*qualFrags)[i]); got != wantQuals[i] {
			t.Errorf("quality fragment %d: got %q, want %q", i, got, wantQuals[i])
		}
	}
}

func TestSeqs2FragmentsWithQualRejectsMismatchedLengths(t *testing.T) {
	sequence := []byte("ACGT")
	quality := []byte("III")
	seqs := []*[]byte{&sequence}
	quals := []*[]byte{&quality}

	if _, _, _, err := seqs2fragmentsWithQual(&seqs, &quals, 4, 1); err == nil {
		t.Fatal("expected an error for mismatched sequence and quality lengths")
	}
}

func TestAdjustPIdentByQual(t *testing.T) {
	const (
		trueIdentity = 0.95
		errorRate    = 0.01 // Q20
	)
	observed := (trueIdentity*(1-errorRate) + (1-trueIdentity)*errorRate/3) * 100
	qual := []byte(strings.Repeat("5", 100)) // ASCII 53 is Phred+33 Q20

	if got := adjustPIdentByQual(observed, qual); math.Abs(got-trueIdentity*100) > 1e-12 {
		t.Fatalf("adjusted identity: got %.12f, want %.12f", got, trueIdentity*100)
	}
	if got := adjustPIdentByQual(100, qual); got != 100 {
		t.Fatalf("upper clamp: got %g, want 100", got)
	}
	if got := adjustPIdentByQual(87.5, nil); got != 87.5 {
		t.Fatalf("empty quality fallback: got %g, want 87.5", got)
	}
	if got := adjustPIdentByQual(87.5, []byte("!")); got != 87.5 {
		t.Fatalf("non-identifiable Q0 fallback: got %g, want 87.5", got)
	}
}

func TestAdjustedPIdentForAlignmentUsesAlignedQueryRegion(t *testing.T) {
	qual := []byte("!!!!555IIII")
	cigar := &wfa.AlignmentResult{QBegin: 2, QEnd: 4}

	got := adjustedPIdentForAlignment(90, qual, 3, cigar)
	want := adjustPIdentByQual(90, qual[4:7])
	if math.Abs(got-want) > 1e-12 {
		t.Fatalf("adjusted identity: got %.12f, want %.12f", got, want)
	}
	if got := adjustedPIdentForAlignment(90, qual, -10, cigar); got != 90 {
		t.Fatalf("invalid coordinate fallback: got %g, want 90", got)
	}
}

func TestGenomeReaderQualityLoading(t *testing.T) {
	dir := t.TempDir()
	fastq := filepath.Join(dir, "query.fastq")
	if err := os.WriteFile(fastq, []byte("@r1\nACGT\n+\nI5+!\n"), 0o600); err != nil {
		t.Fatal(err)
	}

	reader := NewGenomeReader(3, nil)
	query, err := reader.ReadWithQual(fastq, false, false)
	if err != nil {
		t.Fatal(err)
	}
	defer RecycleGQuery(query)
	if len(query.seqs) != 1 || len(query.quals) != 1 {
		t.Fatalf("records: got %d sequences and %d qualities", len(query.seqs), len(query.quals))
	}
	if got := string(*query.quals[0]); got != "I5+!" {
		t.Fatalf("quality: got %q, want %q", got, "I5+!")
	}
}

func TestGenomeReaderQualityLoadingRequiresFASTQ(t *testing.T) {
	file := filepath.Join(t.TempDir(), "query.fasta")
	if err := os.WriteFile(file, []byte(">r1\nACGT\n"), 0o600); err != nil {
		t.Fatal(err)
	}

	reader := NewGenomeReader(3, nil)
	if query, err := reader.ReadWithQual(file, false, false); err == nil {
		if query != nil {
			RecycleGQuery(query)
		}
		t.Fatal("expected FASTA input to be rejected when quality adjustment is enabled")
	}
}

func TestGSearchResultResetQualityAdjustedFields(t *testing.T) {
	r := &GSearchResult{PidentsAdjustedSum: 93.5, ANIAdjusted: 0.935}
	r.Reset()
	if r.PidentsAdjustedSum != 0 || r.ANIAdjusted != 0 {
		t.Fatalf("quality-adjusted fields were not reset: %+v", r)
	}
}
