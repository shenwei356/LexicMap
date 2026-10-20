package cmd

import (
	"fmt"
	"testing"
)

// TestExtendRightFirstAnchor retains the sole matching terminal 2-mer.
func TestExtendRightFirstAnchor(t *testing.T) {
	if q, s := _extendRight([]byte("ACA"), []byte("ACC")); q != 2 || s != 2 {
		t.Fatalf("extension = (%d, %d), want (2, 2)", q, s)
	}
}

// TestExtendMatchFirstAnchor checks terminal matches within sequence and contig bounds on both strands.
func TestExtendMatchFirstAnchor(t *testing.T) {
	for _, tc := range []struct {
		name                       string
		seq1, seq2                 string
		start1, end1, start2, end2 int
		leftLimit, rightLimit      int
		want                       string
		ext                        [4]int // Query left/right, subject left/right extension lengths.
	}{
		{
			name: "right",
			seq1: "GGACA", seq2: "GGACC",
			start1: 0, end1: 2, start2: 0, end2: 2,
			rightLimit: 50,
			want:       "GGAC", ext: [4]int{0, 2, 0, 2},
		},
		{
			name: "left",
			seq1: "ACAGG", seq2: "CCAGG",
			start1: 3, end1: 5, start2: 3, end2: 5,
			leftLimit: 50,
			want:      "CAGG", ext: [4]int{2, 0, 2, 0},
		},
		{
			name: "both clipped to sequence",
			seq1: "ACAGGACA", seq2: "CCAGGACC",
			start1: 3, end1: 5, start2: 3, end2: 5,
			leftLimit: 50, rightLimit: 50,
			want: "CAGGAC", ext: [4]int{2, 2, 2, 2},
		},
		{
			name: "right contig boundary",
			seq1: "GGACA", seq2: "GGACC",
			start1: 0, end1: 2, start2: 0, end2: 2,
			rightLimit: 2,
			want:       "GG",
		},
		{
			name: "left contig boundary",
			seq1: "ACAGG", seq2: "CCAGG",
			start1: 3, end1: 5, start2: 3, end2: 5,
			leftLimit: 2,
			want:      "GG",
		},
		{
			name: "two bases remaining",
			seq1: "GGAC", seq2: "GGAC",
			start1: 0, end1: 2, start2: 0, end2: 2,
			leftLimit: 50, rightLimit: 50,
			want: "GG",
		},
		{
			name: "different query and subject coordinates",
			seq1: "TTGGACA", seq2: "GGACC",
			start1: 2, end1: 4, start2: 0, end2: 2,
			rightLimit: 50,
			want:       "GGAC", ext: [4]int{0, 2, 0, 2},
		},
	} {
		for _, rc := range []bool{false, true} {
			t.Run(fmt.Sprintf("%s/rc=%t", tc.name, rc), func(t *testing.T) {
				tBegin, maxExtLen := tc.leftLimit, tc.rightLimit
				if rc { // The subject's contig limits swap in reverse-complement orientation.
					tBegin, maxExtLen = maxExtLen, tBegin
				}
				q, s, qLeft, qRight, sLeft, sRight, err := extendMatch(
					[]byte(tc.seq1), []byte(tc.seq2), tc.start1, tc.end1, tc.start2, tc.end2,
					50, tBegin, maxExtLen, rc,
				)
				if err != nil {
					t.Fatal(err)
				}
				if got := [4]int{qLeft, qRight, sLeft, sRight}; got != tc.ext {
					t.Fatalf("extension = %v, want %v", got, tc.ext)
				}
				if string(q) != tc.want || string(s) != tc.want {
					t.Fatalf("extended sequences = (%q, %q), want (%q, %q)", q, s, tc.want, tc.want)
				}
			})
		}
	}
}
