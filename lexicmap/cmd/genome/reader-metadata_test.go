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
	"reflect"
	"testing"
)

// writeDeferredFixture includes an incomplete final packed byte and variable IDs.
func writeDeferredFixture(t testing.TB, contigs int) (string, [][]byte) {
	t.Helper()
	file := filepath.Join(t.TempDir(), "genomes.bin")
	w, err := NewWriter(file, 0)
	if err != nil {
		t.Fatal(err)
	}
	sequences := make([][]byte, 3)
	for record := range sequences {
		g := PoolGenome.Get().(*Genome)
		g.Reset()
		g.ID = append(g.ID, fmt.Sprintf("record-%d", record)...)
		for i := 0; i < contigs; i++ {
			length := 5 + (record+i)%4
			g.SeqSizes = append(g.SeqSizes, length)
			id := []byte(fmt.Sprintf("contig-%d-%d-%s", record, i, bytes.Repeat([]byte("x"), i%130)))
			g.SeqIDs = append(g.SeqIDs, &id)
			for j := 0; j < length; j++ {
				g.Seq = append(g.Seq, "ACGT"[(j+i+record)%4])
			}
			g.GenomeSize += length
			if i+1 < contigs {
				g.Seq = append(g.Seq, "AAA"...)
			}
		}
		g.Len, g.NumSeqs = len(g.Seq), contigs
		sequences[record] = bytes.Clone(g.Seq)
		if err := w.Write(g); err != nil {
			t.Fatal(err)
		}
		RecycleGenome(g)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	return file, sequences
}

func TestDeferredMetadataMatchesEager(t *testing.T) {
	for _, contigs := range []int{1, 4, 128} {
		file, sequences := writeDeferredFixture(t, contigs)
		r, err := NewReader(file)
		if err != nil {
			t.Fatal(err)
		}
		for idx, sequence := range sequences {
			// Short reads exercise all four packed-base offsets; later reads reuse g.
			for _, bounds := range [][2]int{{0, 0}, {1, 2}, {2, 4}, {3, 6}, {-3, 3}, {0, len(sequence) + 9}, {len(sequence) - 2, len(sequence) + 2}} {
				eager, err := r.SubSeq3(idx, bounds[0], bounds[1], nil)
				if err != nil {
					t.Fatal(err)
				}
				lazy, err := r.SubSeqWithoutMetadata(idx, bounds[0], bounds[1], nil)
				if err != nil {
					t.Fatal(err)
				}
				if lazy.metadataLoaded || len(lazy.SeqSizes) != 0 || len(lazy.SeqIDs) != 0 || lazy.NumSeqs != 0 {
					t.Fatal("sequence-only extraction parsed metadata")
				}
				if !bytes.Equal(lazy.Seq, eager.Seq) || lazy.SeqOffSet != eager.SeqOffSet || lazy.Len != eager.Len {
					t.Fatalf("contigs=%d record=%d bounds=%v: extraction differs", contigs, idx, bounds)
				}
				begin, end := max(bounds[0], 0), min(bounds[1], len(sequence)-1)
				if !bytes.Equal(lazy.Seq, sequence[begin:end+1]) {
					t.Fatal("incorrect decoded bases")
				}
				lazy, err = r.LoadMetadata(idx, lazy)
				if err != nil {
					t.Fatal(err)
				}
				if lazy.GenomeSize != eager.GenomeSize || lazy.NumSeqs != eager.NumSeqs || lazy.Len != eager.Len || !reflect.DeepEqual(lazy.SeqSizes, eager.SeqSizes) || !reflect.DeepEqual(lazy.SeqIDs, eager.SeqIDs) {
					t.Fatal("hydrated metadata differs")
				}
				id := lazy.SeqIDs[0]
				// Loading metadata again must not duplicate IDs or overwrite a modified strand.
				lazy.Seq[0] = 'T'
				saved := bytes.Clone(lazy.Seq)
				lazy, err = r.LoadMetadata(idx, lazy)
				if err != nil {
					t.Fatal(err)
				}
				if lazy.SeqIDs[0] != id || !bytes.Equal(lazy.Seq, saved) {
					t.Fatal("reloading mutated sequence or metadata")
				}
				lazy, err = r.SubSeqWithoutMetadata(idx, 0, len(sequence)-1, lazy)
				if err != nil {
					t.Fatal(err)
				}
				if !bytes.Equal(lazy.Seq, sequence) || lazy.SeqIDs[0] != id {
					t.Fatal("hydrated reader reuse differs")
				}
				RecycleGenome(eager)
				RecycleGenome(lazy)
			}
		}
		if err := r.Close(); err != nil {
			t.Fatal(err)
		}
	}
}

func TestDeferredMetadataRejectAndLaterLoad(t *testing.T) {
	file, sequences := writeDeferredFixture(t, 128)
	r, err := NewReader(file)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	var g *Genome
	for attempt := 0; attempt < 3; attempt++ {
		// Several rejected chain regions reuse sequence state but never parse IDs.
		g, err = r.SubSeqWithoutMetadata(1, attempt, attempt+3, g)
		if err != nil {
			t.Fatal(err)
		}
		if g.metadataLoaded || len(g.SeqIDs) != 0 {
			t.Fatal("metadata loaded before a surviving chain")
		}
	}
	// The first surviving region can be on the reverse strand.
	for i, j := 0, len(g.Seq)-1; i < j; i, j = i+1, j-1 {
		g.Seq[i], g.Seq[j] = g.Seq[j], g.Seq[i]
	}
	saved := bytes.Clone(g.Seq)
	g, err = r.LoadMetadata(1, g)
	if err != nil {
		t.Fatal(err)
	}
	if !bytes.Equal(g.Seq, saved) || !g.metadataLoaded || g.NumSeqs != 128 {
		t.Fatal("metadata load lost strand or contigs")
	}
	RecycleGenome(g)
	// Recycling a hydrated object must not make the next record appear loaded.
	g, err = r.SubSeqWithoutMetadata(2, 0, len(sequences[2])-1, nil)
	if err != nil {
		t.Fatal(err)
	}
	if g.metadataLoaded || len(g.SeqIDs) != 0 || !bytes.Equal(g.Seq, sequences[2]) {
		t.Fatal("pooled metadata leaked across records")
	}
	RecycleGenome(g)
}

// TestDeferredMetadataInvalidStartRecycles covers transferred eager and lazy state.
func TestDeferredMetadataInvalidStartRecycles(t *testing.T) {
	file, sequences := writeDeferredFixture(t, 2)
	r, err := NewReader(file)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	for _, read := range []func(int, int, int, *Genome) (*Genome, error){r.SubSeq3, r.SubSeqWithoutMetadata} {
		for idx, sequence := range sequences {
			for _, start := range []int{len(sequence), len(sequence) + 5} {
				g, err := read(idx, 0, 0, nil)
				if err != nil {
					t.Fatal(err)
				}
				g, err = r.LoadMetadata(idx, g)
				if err != nil {
					t.Fatal(err)
				}
				// Keep the ID slots to verify error-path recycling without pool ordering.
				ids := g.SeqIDs
				g, err = read(idx, start, start+1, g)
				if err == nil || g != nil {
					t.Fatal("invalid start must return an error and nil genome")
				}
				for _, id := range ids {
					if id != nil {
						t.Fatal("invalid start retained a borrowed contig ID")
					}
				}
				g, err = read(idx, 0, len(sequence)+10, nil)
				if err != nil {
					t.Fatal(err)
				}
				if !bytes.Equal(g.Seq, sequence) {
					t.Fatal("read after an invalid start returned stale sequence")
				}
				RecycleGenome(g)
			}
		}
	}
}

func TestDeferredMetadataFallbackAndErrors(t *testing.T) {
	t.Run("trailing-data", func(t *testing.T) {
		file, sequences := writeDeferredFixture(t, 2)
		f, err := os.OpenFile(file, os.O_APPEND|os.O_WRONLY, 0)
		if err != nil {
			t.Fatal(err)
		}
		if _, err = f.Write([]byte("trailing bytes")); err != nil {
			t.Fatal(err)
		}
		f.Close()
		r, err := NewReader(file)
		if err != nil {
			t.Fatal(err)
		}
		defer r.Close()
		g, err := r.SubSeqWithoutMetadata(2, 0, len(sequences[2])-1, nil)
		if err != nil {
			t.Fatal(err)
		}
		// Without a header read, EOF-derived offsets are checked on metadata load.
		g, err = r.LoadMetadata(2, g)
		if err != ErrBrokenFile || g != nil {
			t.Fatal("trailing data must fail offset validation and recycle g")
		}
	})
	t.Run("bad-contig-id", func(t *testing.T) {
		file, _ := writeDeferredFixture(t, 2)
		r, err := NewReader(file)
		if err != nil {
			t.Fatal(err)
		}
		defer r.Close()
		f, err := os.OpenFile(file, os.O_RDWR, 0)
		if err != nil {
			t.Fatal(err)
		}
		// Corrupt the first contig ID length without touching the packed sequence.
		off := int64(r.Index[0]) + 2 + int64(len("record-0")) + 12 + 4
		if _, err = f.WriteAt([]byte{255, 255}, off); err != nil {
			t.Fatal(err)
		}
		f.Close()
		g, err := r.SubSeqWithoutMetadata(0, 0, 3, nil)
		if err != nil {
			t.Fatal(err)
		}
		if len(g.SeqIDs) != 0 {
			t.Fatal("rejected candidate parsed malformed metadata")
		}
		g, err = r.LoadMetadata(0, g)
		if err == nil || g != nil {
			t.Fatal("malformed metadata must fail when requested and recycle g")
		}
	})
	t.Run("truncated-sequence", func(t *testing.T) {
		file, _ := writeDeferredFixture(t, 2)
		info, err := os.Stat(file)
		if err != nil {
			t.Fatal(err)
		}
		if err = os.Truncate(file, info.Size()-2); err != nil {
			t.Fatal(err)
		}
		r, err := NewReader(file)
		if err != nil {
			t.Fatal(err)
		}
		defer r.Close()
		g, err := r.SubSeqWithoutMetadata(2, 0, 1000, nil)
		if err != nil {
			t.Fatal(err)
		}
		// A truncated final record shifts the inferred start, even if extraction
		// reaches EOF successfully. Metadata must reject that shifted offset.
		g, err = r.LoadMetadata(2, g)
		if err != ErrBrokenFile || g != nil {
			t.Fatal("truncated sequence must fail offset validation and recycle g")
		}
		for _, idx := range []int{-1, 3} {
			if _, err = r.SubSeqWithoutMetadata(idx, 0, 3, nil); err == nil {
				t.Fatal("invalid record accepted")
			}
		}
	})
	t.Run("invalid-record-bounds", func(t *testing.T) {
		file, sequences := writeDeferredFixture(t, 2)
		r, err := NewReader(file)
		if err != nil {
			t.Fatal(err)
		}
		defer r.Close()
		// An impossible next-record offset must select the original parser.
		r.Index[2] = r.Index[0]
		g, err := r.SubSeqWithoutMetadata(0, 0, len(sequences[0])-1, nil)
		if err != nil {
			t.Fatal(err)
		}
		defer RecycleGenome(g)
		if !g.metadataLoaded || !bytes.Equal(g.Seq, sequences[0]) {
			t.Fatal("invalid bounds did not use eager fallback")
		}
	})
}

func BenchmarkDeferredMetadata(b *testing.B) {
	for _, contigs := range []int{1, 4096} {
		file, _ := writeDeferredFixture(b, contigs)
		r, err := NewReader(file)
		if err != nil {
			b.Fatal(err)
		}
		defer r.Close()
		for _, mode := range []string{"eager", "reject", "accept"} {
			b.Run(fmt.Sprintf("contigs=%d/%s", contigs, mode), func(b *testing.B) {
				b.ReportAllocs()
				for range b.N {
					var g *Genome
					var err error
					if mode == "eager" {
						g, err = r.SubSeq3(1, 0, 4, nil)
					} else {
						g, err = r.SubSeqWithoutMetadata(1, 0, 4, nil)
					}
					if err != nil {
						b.Fatal(err)
					}
					if mode == "accept" {
						g, err = r.LoadMetadata(1, g)
						if err != nil {
							b.Fatal(err)
						}
					}
					if len(g.Seq) != 5 {
						b.Fatal("short extraction")
					}
					RecycleGenome(g)
				}
			})
		}
	}
}
