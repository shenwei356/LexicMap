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
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/shenwei356/wfa"
)

// Tail merging preserves every base and uses the exact half-fragment boundary,
// including odd fragment lengths and seed lengths larger than a full fragment.
func TestSeqs2FragmentsTailBoundaries(t *testing.T) {
	for _, tc := range []struct {
		length, fragment, seed int
		lengths                []int
	}{
		{149, 100, 13, []int{149}},
		{150, 100, 13, []int{100, 50}},
		{151, 100, 13, []int{100, 51}},
		{150, 101, 13, []int{150}},
		{151, 101, 13, []int{151}},
		{152, 101, 13, []int{101, 51}},
		{14, 4, 5, []int{4, 4, 6}},
	} {
		t.Run(fmt.Sprintf("length=%d/fragment=%d/seed=%d", tc.length, tc.fragment, tc.seed), func(t *testing.T) {
			sequence := []byte(indexTestDNA(tc.length))
			quality := make([]byte, tc.length)
			for i := range quality {
				quality[i] = byte(33 + i%94)
			}
			seqs, quals := []*[]byte{&sequence}, []*[]byte{&quality}
			plain, n := seqs2fragments(&seqs, tc.fragment, tc.seed)
			defer recycleFragments(plain)
			withQual, qual, qn, err := seqs2fragmentsWithQual(&seqs, &quals, tc.fragment, tc.seed)
			if err != nil {
				t.Fatal(err)
			}
			defer recycleFragments(withQual)
			defer recycleFragments(qual)
			if n != tc.length || qn != tc.length {
				t.Fatalf("base count: %d/%d, want %d", n, qn, tc.length)
			}
			for i, fragments := range []*[][]byte{plain, withQual, qual} {
				if len(*fragments) != len(tc.lengths) {
					t.Errorf("result %d: got %d fragments, want %v", i, len(*fragments), tc.lengths)
					continue
				}
				for j, fragment := range *fragments {
					if len(fragment) != tc.lengths[j] {
						t.Errorf("result %d fragment %d: length %d, want %d", i, j, len(fragment), tc.lengths[j])
					}
				}
				original := sequence
				if i == 2 {
					original = quality
				}
				if !bytes.Equal(bytes.Join(*fragments, nil), original) {
					t.Errorf("result %d did not preserve all original bytes in order", i)
				}
			}
		})
	}
}

func TestSeqs2FragmentsInvalidLength(t *testing.T) {
	sequence, quality := []byte("ACGT"), []byte("IIII")
	seqs, quals := []*[]byte{&sequence}, []*[]byte{&quality}
	for _, length := range []int{0, -1} {
		frags, n := seqs2fragments(&seqs, length, 1)
		if frags != nil || n != 0 {
			t.Fatalf("invalid length %d produced fragments", length)
		}
		for _, values := range []*[]*[]byte{nil, &quals} {
			frags, qual, n, err := seqs2fragmentsWithQual(&seqs, values, length, 1)
			if err == nil || frags != nil || qual != nil || n != 0 {
				t.Fatalf("invalid length %d: fragments=%v, qualities=%v, bases=%d, err=%v", length, frags, qual, n, err)
			}
		}
	}
}

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

	if n != 17 {
		t.Fatalf("fragment bases: got %d, want 17", n)
	}
	wantSeqs := []string{"AAAA", "CCCCGG", "TTTT", "AAA"}
	wantQuals := []string{"1234", "567890", "ABCD", "EFG"}
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

// TestSeqs2FragmentsRetainsShortEnds checks that a seedless tail joins its
// preceding fragment while a short whole contig stays separate.
func TestSeqs2FragmentsRetainsShortEnds(t *testing.T) {
	seq1, seq2 := []byte("AAAAC"), []byte("TT")
	seqs := []*[]byte{&seq1, &seq2}
	frags, n := seqs2fragments(&seqs, 4, 3)
	defer recycleFragments(frags)
	if n != 7 {
		t.Fatalf("fragment bases: got %d, want 7", n)
	}
	want := []string{"AAAAC", "TT"}
	if len(*frags) != len(want) {
		t.Fatalf("fragment count: got %d, want %d", len(*frags), len(want))
	}
	for i, seq := range want {
		if got := string((*frags)[i]); got != seq {
			t.Errorf("fragment %d: got %q, want %q", i, got, seq)
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
