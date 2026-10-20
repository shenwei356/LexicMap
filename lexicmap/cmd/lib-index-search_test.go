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
