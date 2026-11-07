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
	"encoding/binary"
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
	"github.com/shenwei356/lexichash"
)

// Invalid summaries and chunk headers must not panic or truncate existing output.
func TestKmersRejectsInvalidLayoutBeforeOutput(t *testing.T) {
	dir := t.TempDir()
	input := writeTestInput(t, filepath.Join(dir, "g.fa"), ">s\n"+indexTestDNA(500)+"\n")
	db := filepath.Join(dir, "db")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-c", "4", "--partitions", "4", "-O", db, input)
	summary := filepath.Join(db, FileInfo)
	original, err := os.ReadFile(summary)
	if err != nil {
		t.Fatal(err)
	}
	for _, tc := range []struct{ name, old, replacement string }{
		{"zero chunks", "chunks = 4", "chunks = 0"},
		{"negative chunks", "chunks = 4", "chunks = -1"},
		{"too many chunks", "chunks = 4", "chunks = 65"},
		{"missing chunks", "chunks = 4", ""},
		{"wrong chunks", "chunks = 4", "chunks = 1"},
		{"wrong masks", "masks = 64", "masks = 63"},
		{"wrong k", "max-K = 31", "max-K = 32"},
	} {
		t.Run(tc.name, func(t *testing.T) {
			data := strings.Replace(string(original), tc.old, tc.replacement, 1)
			if data == string(original) {
				t.Fatal("fixture field not found")
			}
			if err := os.WriteFile(summary, []byte(data), 0600); err != nil {
				t.Fatal(err)
			}
			t.Cleanup(func() {
				if err := os.WriteFile(summary, original, 0600); err != nil {
					t.Error(err)
				}
			})
			assertKmersLayoutError(t, db)
		})
	}
	// Corrupt a later chunk so preflight must cover more than the first chunk.
	seedFile := filepath.Join(db, DirSeeds, chunkFile(3))
	for _, suffix := range []string{"", kv.KVIndexFileExt} {
		for _, field := range []struct {
			name   string
			offset int
			value  uint64
		}{
			{"first mask", 16, 0}, {"mask count", 24, 15},
		} {
			t.Run(suffix+field.name, func(t *testing.T) {
				file := seedFile + suffix
				data, err := os.ReadFile(file)
				if err != nil {
					t.Fatal(err)
				}
				original := bytes.Clone(data)
				binary.BigEndian.PutUint64(data[field.offset:field.offset+8], field.value)
				if err := os.WriteFile(file, data, 0600); err != nil {
					t.Fatal(err)
				}
				t.Cleanup(func() {
					if err := os.WriteFile(file, original, 0600); err != nil {
						t.Error(err)
					}
				})
				assertKmersLayoutError(t, db)
			})
		}
	}
}

// assertKmersLayoutError exercises the public CLI and checks output preservation.
func assertKmersLayoutError(t *testing.T, db string) {
	t.Helper()
	output := writeTestInput(t, filepath.Join(t.TempDir(), "out.tsv"), "keep existing output\n")
	_, stderr, err := runAuditCommand(t, "", "utils", "kmers", "--quiet", "-d", db, "-m", "0", "-o", output)
	if err == nil || bytes.Contains(stderr, []byte("panic:")) {
		t.Errorf("expected a normal layout error: %v\n%s", err, stderr)
	}
	got, err := os.ReadFile(output)
	if err != nil || string(got) != "keep existing output\n" {
		t.Errorf("invalid layout changed output (%d bytes): %v", len(got), err)
	}
}

// Export preserves every field, occurrence, and row order for both index formats
// and position widths, including empty masks and negative/reversed postings.
func TestKmersIndexFormats(t *testing.T) {
	const k = 31
	const header = "mask\tkmer\tprefix\tnumber\tref\tpos\tstrand\treversed\n"
	for _, narrow := range []bool{false, true} {
		for _, threshold := range []uint64{0, kv.MinIndex15Threshold} {
			t.Run(fmt.Sprintf("narrow=%t/threshold=%d", narrow, threshold), func(t *testing.T) {
				dir := t.TempDir()
				lh, err := lexichash.NewWithSeed(k, 64, 1, 0)
				if err != nil {
					t.Fatal(err)
				}
				if _, err := lh.WriteToFile(filepath.Join(dir, FileMasks)); err != nil {
					t.Fatal(err)
				}
				info := &IndexInfo{MainVersion: MainVersion, K: k, Masks: 64, Chunks: 2, Genomes: 1, GenomeBatches: 1}
				if err := writeIndexInfo(filepath.Join(dir, FileInfo), info); err != nil {
					t.Fatal(err)
				}
				genomeCode := uint64(0)
				if !narrow {
					genomeCode = uint64(512) << BITS_GENOME_IDX
				}
				mapping := make([]byte, 11)
				be.PutUint16(mapping[:2], 1)
				mapping[2] = 'g'
				be.PutUint64(mapping[3:], genomeCode)
				if err := os.WriteFile(filepath.Join(dir, FileGenomeIndex), mapping, 0600); err != nil {
					t.Fatal(err)
				}
				if err := os.Mkdir(filepath.Join(dir, DirSeeds), 0700); err != nil {
					t.Fatal(err)
				}
				var rows [64][]string
				for chunk := 0; chunk < 2; chunk++ {
					file := filepath.Join(dir, DirSeeds, chunkFile(chunk))
					w, err := kv.NewWriterWithIndex15(k, chunk*32, 32, file, 3, 8, narrow, threshold)
					if err != nil {
						t.Fatal(err)
					}
					for mask := chunk * 32; mask < (chunk+1)*32; mask++ {
						data := make(map[uint64]*[]uint64)
						if mask == 0 || mask == 32 {
							prefix := lh.Masks[mask] >> 40 << 40
							for _, suffix := range []uint64{0, 1, 4} {
								key := prefix | suffix
								values := make([]uint64, 700)
								for i := range values {
									pos, flag := i%257, i%4
									values[i] = genomeCode<<BITS_NONE_IDX | uint64(pos)<<2 | uint64(flag)
									strand := byte('+')
									if flag&2 != 0 {
										strand = '-'
									}
									reversed := "no"
									if flag&1 != 0 {
										reversed = "yes"
									}
									rows[mask] = append(rows[mask], fmt.Sprintf("%d\t%s\t%d\t%d\tg\t%d\t%c\t%s\n",
										mask+1, lexichash.MustDecoder()(key, k), util.MustKmerLongestPrefix(key, lh.Masks[mask], k, k), len(values), pos+1, strand, reversed))
								}
								data[key] = &values
							}
						}
						if err := w.WriteDataOfAMask(data); err != nil {
							t.Fatal(err)
						}
					}
					if err := w.Close(); err != nil {
						t.Fatal(err)
					}
					if threshold > 0 && w.Index15Stats().IndexedBlocks == 0 {
						t.Fatal("fixture did not exercise adaptive blocks")
					}
				}
				for _, mask := range []int{0, 33, 64} {
					for _, forward := range []bool{false, true} {
						args := []string{"utils", "kmers", "--quiet", "-j", "1", "-d", dir, "-m", strconv.Itoa(mask)}
						if forward {
							args = append(args, "-f")
						}
						var want strings.Builder
						want.WriteString(header)
						for i, records := range rows {
							if mask != 0 && i != mask-1 {
								continue
							}
							for _, row := range records {
								if !forward || strings.HasSuffix(row, "\tno\n") {
									want.WriteString(row)
								}
							}
						}
						got := runSeedIndexCommand(t, args...)
						if string(got) != want.String() {
							t.Fatalf("mask=%d forward=%t: output changed (%d vs %d bytes)", mask, forward, len(got), want.Len())
						}
					}
				}
			})
		}
	}
}
