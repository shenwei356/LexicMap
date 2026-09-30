// Copyright © 2026 Wei Shen <shenwei356@gmail.com>

package cmd

import "testing"

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
}
