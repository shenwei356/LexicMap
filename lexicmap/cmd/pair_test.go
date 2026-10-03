package cmd

import (
	"math"
	"reflect"
	"testing"
)

func TestDensePairIndexRoundTrip(t *testing.T) {
	for nGenomes := uint32(2); nGenomes <= 100; nGenomes++ {
		seen := make(map[int]struct{}, nGenomes*(nGenomes-1)/2)
		for g1 := uint32(0); g1 < nGenomes; g1++ {
			for g2 := g1 + 1; g2 < nGenomes; g2++ {
				idx := densePairIndex(g1, g2, nGenomes)
				if _, ok := seen[idx]; ok {
					t.Fatalf("duplicate index %d for (%d, %d), n=%d", idx, g1, g2, nGenomes)
				}
				seen[idx] = struct{}{}
				got1, got2 := densePairGenomes(uint32(idx), nGenomes)
				if got1 != g1 || got2 != g2 {
					t.Fatalf("round trip index %d: got (%d, %d), want (%d, %d)", idx, got1, got2, g1, g2)
				}
			}
		}
		want := int(nGenomes * (nGenomes - 1) / 2)
		if len(seen) != want {
			t.Fatalf("n=%d: got %d indices, want %d", nGenomes, len(seen), want)
		}
	}
}

func TestDenseGenomeLimits(t *testing.T) {
	pairs := func(n uint64) uint64 { return n * (n - 1) / 2 }
	if got := pairs(defaultMaxDenseGenomes); got != 199_990_000 {
		t.Fatalf("default dense genome limit: got %d pairs, want 199990000", got)
	}
	if pairs(maxSupportedDenseGenomes) > math.MaxUint32 || pairs(maxSupportedDenseGenomes+1) <= math.MaxUint32 {
		t.Fatal("maxSupportedDenseGenomes does not match the uint32 pair-index limit")
	}
}

func TestProcessKmerWithWindowDenseMatchesMap(t *testing.T) {
	const nGenomes = uint32(4)
	run := func(useDense bool) map[uint64]uint8 {
		counts := make(map[uint64]uint8)
		pairStarts := make([]uint32, nGenomes)
		for g := range nGenomes {
			pairStarts[g] = uint32(densePairRowStart(g, nGenomes))
		}
		var countsPtr *map[uint64]uint8
		var dense *DenseMaskCounts
		if useDense {
			dense = newDenseMaskCounts(int(nGenomes * (nGenomes - 1) / 2))
		} else {
			countsPtr = &counts
		}
		window := &KmerWindow{}
		inputs := []struct {
			code    uint64
			genomes []uint32
		}{
			{0x100, []uint32{0, 2}},
			{0x101, []uint32{1, 3}},
			{0x102, []uint32{0, 3}},
		}
		for _, input := range inputs {
			genomes := append([]uint32(nil), input.genomes...)
			processKmerWithWindow(input.code, &genomes, window, countsPtr, dense, pairStarts,
				math.MaxUint64, 0, 20)
		}
		for _, record := range window.records {
			record.genomes = record.genomes[:0]
			poolKmerRecord.Put(record)
		}
		if useDense {
			for g1 := uint32(0); g1 < nGenomes; g1++ {
				for g2 := g1 + 1; g2 < nGenomes; g2++ {
					prefix := dense.prefixes[densePairIndex(g1, g2, nGenomes)]
					if prefix > 0 {
						counts[uint64(g1)<<32|uint64(g2)] = prefix
					}
				}
			}
		}
		return counts
	}

	want := run(false)
	got := run(true)
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("dense counts: got %v, want %v", got, want)
	}
}

func TestProcessKmerWithWindowResetsAtPrefixBucketBoundary(t *testing.T) {
	window := &KmerWindow{}
	counts := make(map[uint64]uint8)
	genomes1 := []uint32{1}
	genomes2 := []uint32{2}

	processKmerWithWindow(0xff, &genomes1, window, &counts, nil, nil, 0x100, 0, 28)
	processKmerWithWindow(0x100, &genomes2, window, &counts, nil, nil, 0x100, 0, 28)

	if len(counts) != 0 {
		t.Fatalf("cross-bucket k-mers were compared: %v", counts)
	}
	if window.head != 0 || len(window.records) != 1 || window.records[0].code != 0x100 {
		t.Fatalf("window was not reset at bucket boundary: head=%d records=%v", window.head, window.records)
	}
	window.records[0].genomes = window.records[0].genomes[:0]
	poolKmerRecord.Put(window.records[0])
}

func TestMergeDenseMaskCountsResetsAndAccumulates(t *testing.T) {
	const nGenomes = uint32(4)
	nPairs := int(nGenomes * (nGenomes - 1) / 2)
	stats := newDensePairStats(nPairs)
	counts := newDenseMaskCounts(nPairs)
	pairStarts := make([]uint32, nGenomes)
	for g := range nGenomes {
		pairStarts[g] = uint32(densePairRowStart(g, nGenomes))
	}

	setDensePairPrefix(counts, pairStarts, 0, 1, 21)
	setDensePairPrefix(counts, pairStarts, 0, 1, 25)
	setDensePairPrefix(counts, pairStarts, 2, 3, 22)
	mergeDenseMaskCounts(stats, counts, 1, 3, 1, 4, 0.25, 0)

	idx01 := densePairIndex(0, 1, nGenomes)
	idx23 := densePairIndex(2, 3, nGenomes)
	if stats.matches[idx01] != 1 || stats.sumPrefixes[idx01] != 25 {
		t.Fatalf("pair (0,1): got matches=%d sum=%d", stats.matches[idx01], stats.sumPrefixes[idx01])
	}
	if stats.matches[idx23] != 1 || stats.sumPrefixes[idx23] != 22 {
		t.Fatalf("pair (2,3): got matches=%d sum=%d", stats.matches[idx23], stats.sumPrefixes[idx23])
	}
	if counts.prefixes[idx01] != 0 || counts.prefixes[idx23] != 0 || counts.touched[0] != 0 {
		t.Fatalf("mask counts were not reset: prefixes=%v touched=%v", counts.prefixes, counts.touched)
	}

	setDensePairPrefix(counts, pairStarts, 1, 0, 23)
	mergeDenseMaskCounts(stats, counts, 2, 2, 1, 4, 0.25, 0)
	if stats.matches[idx01] != 2 || stats.sumPrefixes[idx01] != 48 {
		t.Fatalf("pair (0,1) after second mask: got matches=%d sum=%d", stats.matches[idx01], stats.sumPrefixes[idx01])
	}
}

func TestBuildGenomeOrdinals(t *testing.T) {
	code0 := uint64(3)
	code1 := uint64(1)<<BITS_GENOME_IDX | 2
	code2 := uint64(1)<<BITS_GENOME_IDX | 7
	id2name := map[uint64][]byte{code2: []byte("c"), code0: []byte("a"), code1: []byte("b")}

	codes, ordinals, err := buildGenomeOrdinals(id2name, 2)
	if err != nil {
		t.Fatal(err)
	}
	wantCodes := []uint32{uint32(code0), uint32(code1), uint32(code2)}
	if !reflect.DeepEqual(codes, wantCodes) {
		t.Fatalf("codes: got %v, want %v", codes, wantCodes)
	}
	for ordinal, code := range codes {
		if got := genomeOrdinal(code, ordinals); got != uint32(ordinal) {
			t.Fatalf("ordinal for code %d: got %d, want %d", code, got, ordinal)
		}
	}
}
