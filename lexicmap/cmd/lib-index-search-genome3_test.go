// Copyright © 2026 Wei Shen <shenwei356@gmail.com>

package cmd

import (
	"errors"
	"os"
	"path/filepath"
	"testing"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
)

func TestAddSampledKmerPosition(t *testing.T) {
	kmers := make(map[uint64]uint32)
	repeated := make(map[uint64]uint64)
	positions := make([]repeatedKmerPosition, 0)

	addSampledKmerPosition(&kmers, &repeated, &positions, 42, 10)
	addSampledKmerPosition(&kmers, &repeated, &positions, 42, 20)
	addSampledKmerPosition(&kmers, &repeated, &positions, 42, 30)

	if got := kmers[42]; got != 11 {
		t.Fatalf("first encoded position: got %d, want 11", got)
	}

	index := uint32(repeated[42] >> 32)
	want := []uint32{20, 30}
	for i, expected := range want {
		if index == 0 {
			t.Fatalf("position %d is missing", i)
		}
		position := positions[index-1]
		if position.position != expected {
			t.Fatalf("position %d: got %d, want %d", i, position.position, expected)
		}
		index = position.next
	}
	if index != 0 {
		t.Fatalf("unexpected additional position at index %d", index)
	}
}

func TestSampledKmerMapCapacity(t *testing.T) {
	if got := sampledKmerMapCapacity(176, 4); got != 44 {
		t.Fatalf("capacity: got %d, want 44", got)
	}
}

func TestSubjectContigBoundsForwardAndReverse(t *testing.T) {
	sketch := &subjectSketch{
		seqLen:       260,
		forwardLen:   100,
		rcStart:      160,
		contigBounds: [][2]int{{0, 40}, {60, 100}},
	}

	if start, end := subjectContigBounds(sketch, 20); start != 0 || end != 40 {
		t.Fatalf("forward bounds: got [%d,%d), want [0,40)", start, end)
	}
	if start, end := subjectContigBounds(sketch, 180); start != 160 || end != 200 {
		t.Fatalf("reverse bounds: got [%d,%d), want [160,200)", start, end)
	}
	if start, end := subjectContigBounds(sketch, 240); start != 220 || end != 260 {
		t.Fatalf("reverse bounds: got [%d,%d), want [220,260)", start, end)
	}

	linear := func(position int) (int, int) {
		if position >= sketch.rcStart {
			for i := len(sketch.contigBounds) - 1; i >= 0; i-- {
				start := sketch.rcStart + sketch.forwardLen - sketch.contigBounds[i][1]
				end := sketch.rcStart + sketch.forwardLen - sketch.contigBounds[i][0]
				if position >= start && position < end {
					return start, end
				}
			}
		} else {
			for _, bounds := range sketch.contigBounds {
				if position >= bounds[0] && position < bounds[1] {
					return bounds[0], bounds[1]
				}
			}
		}
		return 0, sketch.seqLen
	}
	for position := 0; position < sketch.seqLen; position++ {
		wantStart, wantEnd := linear(position)
		start, end := subjectContigBounds(sketch, position)
		if start != wantStart || end != wantEnd {
			t.Fatalf("position %d: got [%d,%d), want [%d,%d)", position, start, end, wantStart, wantEnd)
		}
	}
}

func TestReadGenomeRejectsAllChunksWhenCombinedSizeExceedsLimit(t *testing.T) {
	dbDir := t.TempDir()
	dir := filepath.Join(dbDir, DirGenomes, batchDir(0))
	if err := os.MkdirAll(dir, 0o755); err != nil {
		t.Fatal(err)
	}
	wtr, err := genome.NewWriter(filepath.Join(dir, FileGenomes), 0)
	if err != nil {
		t.Fatal(err)
	}
	for i, sequence := range []string{"ACGTACGT", "TGCATGCA"} {
		g := genome.PoolGenome.Get().(*genome.Genome)
		g.Reset()
		g.ID = append(g.ID, byte('0'+i))
		g.Seq = append(g.Seq, sequence...)
		g.GenomeSize = len(g.Seq)
		g.Len = len(g.Seq)
		g.NumSeqs = 1
		g.SeqSizes = append(g.SeqSizes, len(g.Seq))
		seqID := []byte("seq")
		g.SeqIDs = append(g.SeqIDs, &seqID)
		if err = wtr.Write(g); err != nil {
			genome.RecycleGenome(g)
			t.Fatal(err)
		}
		genome.RecycleGenome(g)
	}
	if err = wtr.Close(); err != nil {
		t.Fatal(err)
	}

	idx := &Index{
		path:           dbDir,
		opt:            &IndexSearchingOptions{MaxSubjectGenomeSize: 10},
		info:           &IndexInfo{GenomeBatches: 1},
		openFileTokens: make(chan int, 1),
	}
	chunks := []uint64{0, 1}
	q, err := idx.ReadGenome(&chunks, "chunked")
	if q != nil || !errors.Is(err, errGenomeTooLarge) {
		t.Fatalf("oversized genome: query=%v err=%v", q, err)
	}

	firstChunk := []uint64{0}
	q, err = idx.ReadGenome(&firstChunk, "chunked")
	if err != nil {
		t.Fatal(err)
	}
	if q.genomeSize != 8 || len(q.seqs) != 1 {
		t.Fatalf("reused query contains partial oversized data: size=%d contigs=%d", q.genomeSize, len(q.seqs))
	}
	RecycleGQuery(q)
}
