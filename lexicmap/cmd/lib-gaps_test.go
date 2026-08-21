// Copyright © 2023-2026 Wei Shen <shenwei356@gmail.com>

package cmd

import (
	"bytes"
	"math/rand"
	"reflect"
	"regexp"
	"testing"
)

func TestFindGapRegions(t *testing.T) {
	reference := regexp.MustCompile(`[Nn]{5,}`)
	tests := [][]byte{
		nil,
		[]byte("ACGT"),
		[]byte("NNNN"),
		[]byte("NNNNN"),
		[]byte("nnNNn"),
		[]byte("NNNNNA"),
		[]byte("ANNNNN"),
		[]byte("ANNNNNANNNNANNNNNN"),
		[]byte("NNNNNnnnnnNNNNN"),
	}

	rng := rand.New(rand.NewSource(1))
	alphabet := []byte("ACGTNnX")
	for length := 0; length <= 4096; length += 17 {
		seq := make([]byte, length)
		for i := range seq {
			seq[i] = alphabet[rng.Intn(len(alphabet))]
		}
		tests = append(tests, seq)
	}

	for _, seq := range tests {
		want0 := reference.FindAllIndex(seq, -1)
		var want [][2]int
		for _, gap := range want0 {
			want = append(want, [2]int{gap[0], gap[1]})
		}
		got0 := findGapRegions(seq)
		var got [][2]int
		if got0 != nil {
			for _, gap := range *got0 {
				start, end := unpackGapRegion(gap)
				got = append(got, [2]int{start, end})
			}
			recycleGapRegions(got0)
		}
		if !reflect.DeepEqual(got, want) {
			t.Fatalf("sequence %q: got %v, want %v", seq, got, want)
		}
	}
}

func benchmarkGapSequence() []byte {
	seq := bytes.Repeat([]byte("ACGT"), 1<<18)
	for i := 10000; i+6 <= len(seq); i += 10000 {
		copy(seq[i:i+6], "NNNNNN")
	}
	return seq
}

func BenchmarkFindGapRegionsPacked(b *testing.B) {
	seq := benchmarkGapSequence()
	b.SetBytes(int64(len(seq)))
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		gaps := findGapRegions(seq)
		recycleGapRegions(gaps)
	}
}

func BenchmarkFindGapRegionsRegexp(b *testing.B) {
	seq := benchmarkGapSequence()
	reference := regexp.MustCompile(`[Nn]{5,}`)
	b.SetBytes(int64(len(seq)))
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = reference.FindAllSubmatchIndex(seq, -1)
	}
}
