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
	"strings"
	"testing"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/seedposition"
)

// writeSeedPosFixture stores real genome/position records without a seed index.
func writeSeedPosFixture(t *testing.T, batches int) string {
	t.Helper()
	dir := t.TempDir()
	sequence := []byte(indexTestDNA(128))
	var mapping bytes.Buffer
	for batch := 0; batch < batches; batch++ {
		name := fmt.Sprintf("g%d", batch)
		path := filepath.Join(dir, DirGenomes, batchDir(batch))
		if err := os.MkdirAll(path, 0700); err != nil {
			t.Fatal(err)
		}
		writer, err := genome.NewWriter(filepath.Join(path, FileGenomes), uint32(batch))
		if err != nil {
			t.Fatal(err)
		}
		id := []byte("s")
		g := &genome.Genome{ID: []byte(name), Seq: sequence, GenomeSize: len(sequence), Len: len(sequence), NumSeqs: 1, SeqSizes: []int{len(sequence)}, SeqIDs: []*[]byte{&id}}
		if err := writer.Write(g); err != nil {
			t.Fatal(err)
		}
		if err := writer.Close(); err != nil {
			t.Fatal(err)
		}
		positions, err := seedposition.NewWriter(filepath.Join(path, FileSeedPositions), uint32(batch))
		if err != nil {
			t.Fatal(err)
		}
		if err := positions.Write([]uint32{10 << 2, 30 << 2, 60 << 2}); err != nil {
			t.Fatal(err)
		}
		if err := positions.Close(); err != nil {
			t.Fatal(err)
		}
		binary.Write(&mapping, binary.BigEndian, uint16(len(name)))
		mapping.WriteString(name)
		binary.Write(&mapping, binary.BigEndian, uint64(batch)<<BITS_GENOME_IDX)
	}
	if err := os.WriteFile(filepath.Join(dir, FileGenomeIndex), mapping.Bytes(), 0600); err != nil {
		t.Fatal(err)
	}
	info := &IndexInfo{MainVersion: MainVersion, K: 31, Chunks: 2, GenomeBatches: batches, Genomes: batches, InputGenomes: batches}
	if err := writeIndexInfo(filepath.Join(dir, FileInfo), info); err != nil {
		t.Fatal(err)
	}
	return dir
}

func TestSeedPosReadsRequestedBatches(t *testing.T) {
	db := writeSeedPosFixture(t, 3)
	if err := os.Remove(filepath.Join(db, DirGenomes, batchDir(2), FileGenomes)); err != nil {
		t.Fatal(err)
	}
	output := runSeedIndexCommand(t, "utils", "seed-pos", "--quiet", "-d", db, "-n", "g0", "--max-open-files", "4")
	if strings.Count(string(output), "\n") != 4 {
		t.Fatalf("lost seed positions: %s", output)
	}
}

func TestSeedPosPlotsWithBoundedReaders(t *testing.T) {
	db := writeSeedPosFixture(t, 2)
	plots := filepath.Join(t.TempDir(), "plots")
	runSeedIndexCommand(t, "utils", "seed-pos", "--quiet", "-v", "-j", "8", "-d", db, "-a", "--max-open-files", "4", "--plot-dir", plots, "--plot-ext", ".svg", "--slid-window", "40", "--slid-step", "20")
	files, err := filepath.Glob(filepath.Join(plots, "*.svg"))
	if err != nil {
		t.Fatal(err)
	}
	if len(files) != 4 {
		t.Fatalf("expected two plots per genome, got %v", files)
	}
}
