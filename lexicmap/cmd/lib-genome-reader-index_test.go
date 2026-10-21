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
	"math/rand"
	"os"
	"path/filepath"
	"strings"
	"sync"
	"testing"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
)

func writeSharedReaderFixture(t *testing.T, root string) string {
	t.Helper()
	file := filepath.Join(root, DirGenomes, batchDir(0), FileGenomes)
	if err := os.MkdirAll(filepath.Dir(file), 0755); err != nil {
		t.Fatal(err)
	}
	w, err := genome.NewWriter(file, 0)
	if err != nil {
		t.Fatal(err)
	}
	g := genome.PoolGenome.Get().(*genome.Genome)
	g.Reset()
	g.ID = append(g.ID, "shared"...)
	g.Seq = append(g.Seq, "ACGTACGT"...)
	g.GenomeSize, g.Len, g.NumSeqs = 8, 8, 1
	g.SeqSizes = append(g.SeqSizes, 8)
	id := []byte("contig")
	g.SeqIDs = append(g.SeqIDs, &id)
	err = w.Write(g)
	genome.RecycleGenome(g)
	if err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	return file
}

func TestSharedGenomeReaderFallbackAndRetry(t *testing.T) {
	root := t.TempDir()
	idx := &Index{path: root, opt: &IndexSearchingOptions{}, info: &IndexInfo{GenomeBatches: 1}, openFileTokens: make(chan int, 2)}
	if _, err := idx.acquireGenomeReader(0); err == nil || len(idx.openFileTokens) != 0 {
		t.Fatal("failed metadata load must release its descriptor token")
	}
	file := writeSharedReaderFixture(t, root)
	r, err := idx.acquireGenomeReader(0)
	if err != nil {
		t.Fatal(err)
	}
	directory := &r.Index[0]
	if err := idx.releaseGenomeReader(0, r); err != nil {
		t.Fatal(err)
	}
	// The fallback must reopen only data, even after a previously failed load.
	if err := os.Rename(file+genome.GenomeIndexFileExt, file+".saved-index"); err != nil {
		t.Fatal(err)
	}
	var wg sync.WaitGroup
	for range 8 {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for range 10 {
				r, err := idx.acquireGenomeReader(0)
				if err != nil {
					t.Error(err)
					return
				}
				if &r.Index[0] != directory {
					t.Error("fallback reader did not reuse its batch directory")
				}
				g, err := r.SubSeq(0, 1, 6)
				if err != nil {
					t.Error(err)
				} else {
					if !bytes.Equal(g.Seq, []byte("CGTACG")) {
						t.Errorf("unexpected sequence: %s", g.Seq)
					}
					genome.RecycleGenome(g)
				}
				if err := idx.releaseGenomeReader(0, r); err != nil {
					t.Error(err)
				}
			}
		}()
	}
	wg.Wait()
	if err := idx.Close(); err != nil {
		t.Fatal(err)
	}
	idx.genomeReaderIndexes.Range(func(_, _ any) bool {
		t.Error("closed index retains a batch directory")
		return false
	})
	if len(idx.openFileTokens) != 0 {
		t.Fatal("fallback leaked descriptor tokens")
	}
}

func TestSharedGenomeReaderDirectoryIsIndexLocal(t *testing.T) {
	a, b := &Index{path: t.TempDir()}, &Index{path: t.TempDir()}
	writeSharedReaderFixture(t, a.path)
	writeSharedReaderFixture(t, b.path)
	ra, err := a.newGenomeReader(0)
	if err != nil {
		t.Fatal(err)
	}
	defer ra.Close()
	rb, err := b.newGenomeReader(0)
	if err != nil {
		t.Fatal(err)
	}
	defer rb.Close()
	if &ra.Index[0] == &rb.Index[0] {
		t.Fatal("different indexes share a batch directory")
	}
}

func TestSearchSeedFileBudget(t *testing.T) {
	for _, tc := range []struct {
		budget, chunks, request, want int
		memory                        bool
	}{
		{32, 4, 16, 5, false}, {13, 4, 16, 1, false}, {12, 4, 16, 0, false},
		{1024, 4, 0, 2, false}, {9, 16, 16, 16, true}, {8, 16, 16, 0, true},
	} {
		got, err := searchSeedConcurrency(tc.budget, tc.chunks, tc.request, tc.memory)
		if got != tc.want || (err != nil) != (tc.want == 0) {
			t.Fatalf("%+v: got %d, %v", tc, got, err)
		}
	}
}

func TestGenomeComparatorFileBudget(t *testing.T) {
	db := writeDetailsFixture(t, 2)
	for _, budget := range []int{8, 9, 10, 16} {
		opt := DefaultIndexSearchingOptions
		opt.MaxOpenFiles, opt.NumCPUs = budget, 4
		idx, err := NewGenomeComparator(db, &opt)
		if budget == 8 {
			if err == nil {
				idx.Close()
				t.Fatal("accepted a budget without room for a genome reader")
			}
			continue
		}
		if err != nil {
			t.Fatal(err)
		}
		if cap(idx.openFileTokens) != budget-searchReservedFiles || idx.hasGenomeRdrs != (budget >= 10) {
			t.Fatalf("budget %d: incorrect process allowance or reader pool", budget)
		}
		ids := []uint64{0, 1 << BITS_GENOME_IDX}
		q, err := idx.ReadGenome(&ids, "fixture")
		if err != nil {
			t.Fatal(err)
		}
		if q == nil || q.genomeSize == 0 {
			t.Fatal("bounded readers lost sequences")
		}
		RecycleGQuery(q)
		if err := idx.Close(); err != nil {
			t.Fatal(err)
		}
	}
}

// buildMultiBatchSearchFixture builds four chunks from deterministic temporary data.
func buildMultiBatchSearchFixture(t *testing.T) (string, string) {
	t.Helper()
	dir := t.TempDir()
	query := filepath.Join(dir, "ref.fa")
	seq := make([]byte, 16000)
	rng := rand.New(rand.NewSource(13))
	for i := range seq {
		seq[i] = "ACGT"[rng.Intn(4)]
	}
	if err := os.WriteFile(query, append([]byte(">seq\n"), seq...), 0600); err != nil {
		t.Fatal(err)
	}
	reference := filepath.Join(dir, "ref2.fa")
	for i := 0; i < len(seq); i += 37 {
		seq[i] = "ACGT"[(strings.IndexByte("ACGT", seq[i])+1)%4]
	}
	if err := os.WriteFile(reference, append([]byte(">seq2\n"), seq...), 0600); err != nil {
		t.Fatal(err)
	}
	db := filepath.Join(dir, "db.lmi")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-c", "4", "-m", "64", "-k", "16", "--partitions", "4", "-b", "1", "-O", db, query, reference)
	return db, query
}
