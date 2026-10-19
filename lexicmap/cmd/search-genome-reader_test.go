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
	"math/rand"
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"

	"github.com/shenwei356/lexichash"
)

// writeGenomeReaderInput keeps the same contigs in FASTA and Phred+33 FASTQ fixtures.
func writeGenomeReaderInput(t *testing.T, file string, contigs []string, fastq bool) {
	t.Helper()
	var data strings.Builder
	for i, sequence := range contigs {
		if fastq {
			fmt.Fprintf(&data, "@seq%d\n%s\n+\n%s\n", i, sequence, strings.Repeat("I", len(sequence)))
		} else {
			fmt.Fprintf(&data, ">seq%d\n%s\n", i, sequence)
		}
	}
	if err := os.WriteFile(file, []byte(data.String()), 0600); err != nil {
		t.Fatal(err)
	}
}

// Spacers and real N runs can overlap, but must form one ordered set of intervals.
func TestGenomeReaderContigBoundaries(t *testing.T) {
	contigs := []string{"ACGTNNNNN", "NNNNNGTAC", "TGCATGCA"}
	for _, k := range []int{3, 31} {
		for _, fastq := range []bool{false, true} {
			for _, convert := range []bool{false, true} {
				t.Run(fmt.Sprintf("k=%d/fastq=%v/convert=%v", k, fastq, convert), func(t *testing.T) {
					file := filepath.Join(t.TempDir(), "query.fa")
					writeGenomeReaderInput(t, file, contigs, fastq)
					reader := NewGenomeReader(k, nil)
					var q *GQuery
					var err error
					if fastq {
						q, err = reader.ReadWithQual(file, convert, false)
					} else {
						q, err = reader.Read(file, convert, false)
					}
					if err != nil {
						t.Fatal(err)
					}
					defer RecycleGQuery(q)
					wantSeqs := append([]string(nil), contigs...)
					wantRegions := []int{4, k + 13, k + 18, 2*k + 17}
					if convert {
						for i := range wantSeqs {
							wantSeqs[i] = strings.ReplaceAll(wantSeqs[i], "N", "A")
						}
						wantRegions = []int{9, k + 8, k + 18, 2*k + 17}
					}
					want := strings.Join(wantSeqs, strings.Repeat("N", k))
					if string(q.bigSeq) != want || !reflect.DeepEqual(q.skipRegions, wantRegions) {
						t.Fatalf("screening input: %s, %v; want %s, %v", q.bigSeq, q.skipRegions, want, wantRegions)
					}
					if q.genomeSize != 26 || len(q.seqs) != len(contigs) {
						t.Fatalf("spacers changed genome size/contigs: %d/%d", q.genomeSize, len(q.seqs))
					}
					for i, sequence := range q.seqs {
						if string(*sequence) != wantSeqs[i] {
							t.Fatalf("alignment contig %d changed: %s", i, *sequence)
						}
						if fastq && string(*q.quals[i]) != strings.Repeat("I", len(contigs[i])) {
							t.Fatalf("qualities changed for contig %d", i)
						}
					}
				})
			}
		}
	}
}

func TestGenomeReaderScreeningIgnoresContigOrder(t *testing.T) {
	rng := rand.New(rand.NewSource(147))
	contigs := make([]string, 40)
	for i := range contigs {
		sequence := make([]byte, 250)
		for j := range sequence {
			sequence[j] = "ACGT"[rng.Intn(4)]
		}
		contigs[i] = string(sequence)
	}
	lh, err := lexichash.NewWithSeed(31, 64, 1, 0)
	if err != nil {
		t.Fatal(err)
	}
	if err := lh.IndexMasks(3); err != nil {
		t.Fatal(err)
	}
	if err := lh.IndexMasksWithDistinctPrefixes(4); err != nil {
		t.Fatal(err)
	}
	selection := make([]int, len(lh.Masks))
	for i := range selection {
		selection[i] = i
	}
	sketcher, err := lh.NewSketcher(selection)
	if err != nil {
		t.Fatal(err)
	}
	var want []uint64
	for _, reverse := range []bool{false, true} {
		if reverse {
			for i, j := 0, len(contigs)-1; i < j; i, j = i+1, j-1 {
				contigs[i], contigs[j] = contigs[j], contigs[i]
			}
		}
		file := filepath.Join(t.TempDir(), "query.fa")
		writeGenomeReaderInput(t, file, contigs, false)
		q, err := NewGenomeReader(31, nil).Read(file, true, false)
		if err != nil {
			t.Fatal(err)
		}
		got, err := sketcher.Mask(q.bigSeq, q.skipRegions, true)
		RecycleGQuery(q)
		if err != nil {
			t.Fatal(err)
		}
		if want == nil {
			want = append([]uint64(nil), got...)
		} else if !reflect.DeepEqual(got, want) {
			t.Fatalf("reordering contigs changed screening seeds:\ngot %v\nwant %v", got, want)
		}
	}
}

// Neither the current sketcher nor legacy masking may invent seeds across records.
func TestGenomeScreenRejectsJunctionSeeds(t *testing.T) {
	dir := t.TempDir()
	rng := rand.New(rand.NewSource(7))
	sequence := make([]byte, 5000)
	for i := range sequence {
		sequence[i] = "ACGT"[rng.Intn(4)]
	}
	ref := filepath.Join(dir, "ref.fa")
	writeGenomeReaderInput(t, ref, []string{string(sequence)}, false)
	index := filepath.Join(dir, "db.lmi")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "--partitions", "4", "-O", index, ref)
	var contigs []string
	for i := 0; i < len(sequence); i += 16 {
		contigs = append(contigs, string(sequence[i:min(i+16, len(sequence))]))
	}
	info, err := readIndexInfo(filepath.Join(index, FileInfo))
	if err != nil {
		t.Fatal(err)
	}
	for _, version := range []uint8{5, 4} {
		info.MinorVersion = version
		if err := writeIndexInfo(filepath.Join(index, FileInfo), info); err != nil {
			t.Fatal(err)
		}
		for _, fastq := range []bool{false, true} {
			query := filepath.Join(dir, "query.fa")
			writeGenomeReaderInput(t, query, contigs, fastq)
			out := runSeedIndexCommand(t, "genome", "search", "--quiet", "-j", "1", "-d", index, "-m", "0", "-n", "0", "-S", query)
			if bytes.Count(out, []byte("\n")) != 1 {
				t.Fatalf("version %d, FASTQ %v: seed hits from contigs shorter than k: %s", version, fastq, out)
			}
		}
	}
}
