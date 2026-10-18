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
	"testing"
)

// TestSearchLazyGenomeMeta exercises the public flag with full alignment output,
// multiple contigs, reverse-strand matches, and an unmatched query.
func TestSearchLazyGenomeMeta(t *testing.T) {
	flag := mapCmd.Flags().Lookup("lazy-genome-meta")
	if flag == nil || flag.DefValue != "false" || DefaultIndexSearchingOptions.LazyGenomeMeta {
		t.Fatal("lazy metadata must be opt-in at both CLI and library levels")
	}
	dir := t.TempDir()
	rng := rand.New(rand.NewSource(17))
	var queries bytes.Buffer
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
	eager := runSeedIndexCommand(t, args...)
	if !bytes.Contains(eager, []byte("query-0\t")) || !bytes.Contains(eager, []byte("query-1\t")) || bytes.Contains(eager, []byte("unmatched\t")) {
		t.Fatalf("unexpected fixture matches:\n%s", eager)
	}
	for _, option := range []string{"--lazy-genome-meta=false", "--lazy-genome-meta"} {
		output := runSeedIndexCommand(t, append(append([]string{}, args...), option)...)
		if !bytes.Equal(eager, output) {
			t.Fatalf("%s changed search output\neager:\n%s\nactual:\n%s", option, eager, output)
		}
	}
}
