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
	"strings"
	"testing"
)

// TestSearchLazyGenomeMeta exercises default search with full alignment output,
// multiple contigs, reverse-strand matches, and an unmatched query.
func TestSearchLazyGenomeMeta(t *testing.T) {
	dir := t.TempDir()
	rng := rand.New(rand.NewSource(17))
	var queries bytes.Buffer
	expectedQueries := make(map[string]string)
	refs := make([]string, 2)
	for record, contigs := range []int{1, 4} {
		var fasta bytes.Buffer
		for contig := 0; contig < contigs; contig++ {
			sequence := make([]byte, 800)
			for i := range sequence {
				sequence[i] = "ACGT"[rng.Intn(4)]
			}
			fmt.Fprintf(&fasta, ">contig-%d-%d\n%s\n", record, contig, sequence)
			if contig == contigs-1 {
				query := bytes.Clone(sequence[200:600])
				if record == 1 {
					RC(query)
				}
				expectedQueries[fmt.Sprintf("query-%d", record)] = string(query)
				fmt.Fprintf(&queries, ">query-%d\n%s\n", record, query)
			}
		}
		refs[record] = filepath.Join(dir, fmt.Sprintf("ref%d.fa", record))
		if err := os.WriteFile(refs[record], fasta.Bytes(), 0600); err != nil {
			t.Fatal(err)
		}
	}
	fmt.Fprintf(&queries, ">unmatched\n%s\n", bytes.Repeat([]byte("A"), 400))
	queryFile := filepath.Join(dir, "query.fa")
	if err := os.WriteFile(queryFile, queries.Bytes(), 0600); err != nil {
		t.Fatal(err)
	}
	db := filepath.Join(dir, "test.lmi")
	args := []string{"index", "--quiet", "-j", "1", "-c", "2", "-m", "64", "-k", "16", "--partitions", "4", "-O", db}
	runSeedIndexCommand(t, append(args, refs...)...)
	args = []string{"search", "--quiet", "-j", "1", "-J", "1", "-d", db, "-p", "12", "-P", "12", "-n", "0", "-N", "0", "-a", queryFile}
	output := runSeedIndexCommand(t, args...)
	rows := strings.Split(strings.TrimSpace(string(output)), "\n")
	if len(rows) != 3 { // Header plus the two exact matches; no hit for unmatched.
		t.Fatalf("unexpected fixture matches:\n%s", output)
	}
	for i, row := range rows[1:] {
		fields := strings.Split(row, "\t")
		queryID := fmt.Sprintf("query-%d", i)
		contigID, strand := "contig-0-0", "+"
		if i == 1 {
			contigID, strand = "contig-1-3", "-"
		}
		// Metadata must resolve the original contig and coordinates on both strands.
		if len(fields) != 24 || fields[0] != queryID || fields[3] != fmt.Sprintf("ref%d", i) || fields[4] != contigID ||
			strings.Join(fields[12:18], "\t") != "1\t400\t201\t600\t"+strand+"\t800" ||
			fields[21] != expectedQueries[queryID] || fields[22] != expectedQueries[queryID] {
			t.Fatalf("unexpected alignment for %s:\n%s", queryID, row)
		}
	}
}
