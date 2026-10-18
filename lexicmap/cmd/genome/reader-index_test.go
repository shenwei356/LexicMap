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

package genome

import (
	"bytes"
	"fmt"
	"os"
	"path/filepath"
	"sync"
	"testing"
)

func writeReaderIndexFixture(t testing.TB, records int) string {
	t.Helper()
	file := filepath.Join(t.TempDir(), "genomes.bin")
	w, err := NewWriter(file, 7)
	if err != nil {
		t.Fatal(err)
	}
	for i := range records {
		g := PoolGenome.Get().(*Genome)
		g.Reset()
		g.ID = append(g.ID, fmt.Sprintf("genome-%d", i)...)
		g.Seq = append(g.Seq, "ACGTACGTTGCATGCA"...)
		g.GenomeSize, g.Len, g.NumSeqs = len(g.Seq), len(g.Seq), 2
		g.SeqSizes = append(g.SeqSizes, 8, 8)
		for _, name := range []string{"first", "second"} {
			id := []byte(name)
			g.SeqIDs = append(g.SeqIDs, &id)
		}
		err = w.Write(g)
		RecycleGenome(g)
		if err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	return file
}

func TestReaderIndexIndependentReaders(t *testing.T) {
	file := writeReaderIndexFixture(t, 3)
	index, err := ReadIndex(file)
	if err != nil {
		t.Fatal(err)
	}
	if !filepath.IsAbs(index.file) {
		t.Fatal("shared directory must bind readers to an absolute data path")
	}
	a, err := index.NewReader()
	if err != nil {
		t.Fatal(err)
	}
	b, err := index.NewReader()
	if err != nil {
		t.Fatal(err)
	}
	defer b.Close()
	if &a.Index[0] != &b.Index[0] || a.fhData == b.fhData || a.bufReader == b.bufReader {
		t.Fatal("readers must share only the immutable directory")
	}
	if err := a.Close(); err != nil {
		t.Fatal(err)
	}
	if a.Index != nil {
		t.Fatal("closed pooled reader retains the shared directory")
	}
	// New readers use the loaded directory even when the companion is unavailable.
	if err := os.Rename(file+GenomeIndexFileExt, file+".saved-index"); err != nil {
		t.Fatal(err)
	}
	c, err := index.NewReader()
	if err != nil {
		t.Fatal(err)
	}
	defer c.Close()
	for _, r := range []*Reader{b, c} {
		g, err := r.SubSeq3(2, 3, 12, nil)
		if err != nil {
			t.Fatal(err)
		}
		if !bytes.Equal(g.Seq, []byte("TACGTTGCAT")) || g.NumSeqs != 2 || len(g.SeqSizes) != 2 {
			t.Fatalf("wrong shared-reader extraction: sequence=%s contigs=%v", g.Seq, g.SeqSizes)
		}
		RecycleGenome(g)
	}
	if _, err := NewReader(file); err == nil {
		t.Fatal("standalone reader unexpectedly skipped loading its own directory")
	}
}

func TestReaderIndexConcurrentExtraction(t *testing.T) {
	file := writeReaderIndexFixture(t, 4)
	index, err := ReadIndex(file)
	if err != nil {
		t.Fatal(err)
	}
	want := []byte("ACGTACGTTGCATGCA")
	var wg sync.WaitGroup
	for range 8 {
		wg.Add(1)
		go func() {
			defer wg.Done()
			r, err := index.NewReader()
			if err != nil {
				t.Error(err)
				return
			}
			defer r.Close()
			for i := range 64 {
				start := i % len(want)
				g, err := r.SubSeq(i%4, start, len(want)-1)
				if err != nil {
					t.Error(err)
					return
				}
				if !bytes.Equal(g.Seq, want[start:]) {
					t.Errorf("start=%d sequence=%s", start, g.Seq)
				}
				RecycleGenome(g)
			}
		}()
	}
	wg.Wait()
}

func TestReaderIndexErrorsAndDataOpenRetry(t *testing.T) {
	file := writeReaderIndexFixture(t, 1)
	for _, path := range []string{"", file + GenomeIndexFileExt} {
		if _, err := ReadIndex(path); err == nil {
			t.Fatalf("accepted invalid genome data path %q", path)
		}
	}
	original, err := os.ReadFile(file + GenomeIndexFileExt)
	if err != nil {
		t.Fatal(err)
	}
	badMagic, badVersion, badCount := bytes.Clone(original), bytes.Clone(original), bytes.Clone(original)
	badMagic[0] ^= 255
	badVersion[8]++
	be.PutUint32(badCount[20:24], ^uint32(0))
	for _, data := range [][]byte{nil, original[:7], original[:15], original[:23], original[:len(original)-1], badMagic, badVersion, badCount} {
		if err := os.WriteFile(file+GenomeIndexFileExt, data, 0600); err != nil {
			t.Fatal(err)
		}
		if _, err := ReadIndex(file); err == nil {
			t.Fatalf("accepted malformed directory of length %d", len(data))
		}
	}
	if err := os.WriteFile(file+GenomeIndexFileExt, original, 0600); err != nil {
		t.Fatal(err)
	}
	index, err := ReadIndex(file)
	if err != nil {
		t.Fatal(err)
	}
	if err := os.Rename(file, file+".saved-data"); err != nil {
		t.Fatal(err)
	}
	if _, err := index.NewReader(); err == nil {
		t.Fatal("accepted missing data file")
	}
	if err := os.Rename(file+".saved-data", file); err != nil {
		t.Fatal(err)
	}
	r, err := index.NewReader()
	if err != nil {
		t.Fatal(err)
	}
	r.Close()
}

func BenchmarkReaderIndexOpen(b *testing.B) {
	file := writeReaderIndexFixture(b, 10000)
	index, err := ReadIndex(file)
	if err != nil {
		b.Fatal(err)
	}
	for _, shared := range []bool{false, true} {
		b.Run(fmt.Sprintf("shared=%t", shared), func(b *testing.B) {
			b.ReportAllocs()
			for range b.N {
				var r *Reader
				var err error
				if shared {
					r, err = index.NewReader()
				} else {
					r, err = NewReader(file)
				}
				if err != nil {
					b.Fatal(err)
				}
				if err := r.Close(); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}
