package cmd

import (
	"bufio"
	"errors"
	"io"
	"math/rand/v2"
	"os"
	"path/filepath"
	"reflect"
	"slices"
	"strconv"
	"strings"
	"sync"
	"testing"
)

// Include every identifier byte, sparse byte ranges, repeated IDs, and partial
// final runs. Coordinates/flags distinguish ties so stable ordering is checked.
func TestSeedSpillRadixStableAndBudgeted(t *testing.T) {
	for _, n := range []int{1023, 1024, 11000} {
		for _, bits := range []int{0, 8, 34, 64} {
			for _, share := range []int{1023, 1024, 1500, 2047, 2048, 8192, 32768} {
				b := newSeedMemoryBudget(int64(share)*seedAnchorBytes, 1)
				b.tempDir = t.TempDir()
				s := newSeedSpillStore(b)
				// Exercise reserved producer space as well as the radix arrays.
				s.transferReserve = share / 8
				rng := rand.New(rand.NewPCG(17, 29))
				want := make([]seedAnchor, n)
				for i := range want {
					id := rng.Uint64()
					if bits < 64 {
						id &= (uint64(1) << bits) - 1
					}
					if bits > 0 && i%3 == 0 {
						id = 255 // many equal IDs, including the last byte bucket
					}
					want[i] = seedAnchor{batchGenomeIndex: id, qBegin: int32(i), length: 31, qrc: i%2 == 0}
				}
				for start := 0; start < n; start += 257 {
					if err := s.addBatch(want[start:min(n, start+257)]); err != nil {
						t.Fatal(err)
					}
				}
				slices.SortStableFunc(want, compareSeedAnchorGenome)
				var got []seedAnchor
				if err := s.consume(func(a seedAnchor) error { got = append(got, a); return nil }); err != nil {
					t.Fatal(err)
				}
				if !slices.Equal(got, want) || s.peakBufferBytes > int64(share)*seedAnchorBytes {
					t.Fatalf("n=%d bits=%d share=%d: unstable grouping or budget exceeded (%d)", n, bits, share, s.peakBufferBytes)
				}
				if err := s.close(); err != nil {
					t.Fatal(err)
				}
			}
		}
	}
}

func TestSeedSpillParallelReplayBoundaries(t *testing.T) {
	for _, limit := range []int64{128 << 10, 4 << 20} {
		for _, oversized := range []bool{false, true} {
			idx := seedChainingTestIndex(4)
			b := newSeedMemoryBudget(limit, 1)
			b.tempDir = t.TempDir()
			s := newSeedSpillStore(b)
			defer s.close()
			original := seedChainingTestResults(9000) // crosses descriptor and buffer boundaries
			for _, targets := range original {
				for _, r := range targets {
					// Repeated anchors make the small replay buffer fill between
					// genome boundaries, while preserving the reference chain.
					r.Subs = slices.Repeat(r.Subs, 3)
				}
			}
			if oversized {
				// More than the small replay buffer, but duplicates compact to one anchor
				// before chaining. Normal groups occur both before and after it.
				r := original[0][1000]
				r.Subs = make([]SubstrPair, 20000)
				for i := range r.Subs {
					r.Subs[i] = SubstrPair{Len: 31}
				}
			}
			for _, targets := range original {
				for _, r := range targets {
					for _, p := range r.Subs {
						if err := s.add(seedAnchor{batchGenomeIndex: r.BatchGenomeIndex, qBegin: p.QBegin, tBegin: p.TBegin, length: p.Len, qrc: p.QRC, trc: p.TRC}); err != nil {
							t.Fatal(err)
						}
					}
				}
			}
			var arena seedSearchResultArena
			got, err := idx.chainSeedSpillStore(s, &arena)
			if err != nil {
				t.Fatal(err)
			}
			want := idx.chainSeedResults(original, 9000)
			slices.SortFunc(*want, func(a, b *seedSearchResult) int {
				return compareSeedAnchorGenome(seedAnchor{batchGenomeIndex: a.BatchGenomeIndex}, seedAnchor{batchGenomeIndex: b.BatchGenomeIndex})
			})
			if len(*got) != len(*want) {
				t.Fatalf("got %d genomes, want %d", len(*got), len(*want))
			}
			for i, r := range *got {
				w := (*want)[i]
				if r.BatchGenomeIndex != w.BatchGenomeIndex || r.Score != w.Score || !slices.Equal(r.chainRegions, w.chainRegions) || len(r.Subs) != 0 {
					t.Fatalf("genome %d: parallel replay changed order, score, bounds or retained anchors", r.BatchGenomeIndex)
				}
			}
			if s.peakBufferBytes > limit || s.replayBufferBytes != 0 {
				t.Fatalf("replay exceeded charged budget or retained charge: peak=%d", s.peakBufferBytes)
			}
			recycleSeedSearchResultSlice(got)
			recycleSeedSearchResultSlice(want)
			arena.recycle()
		}
	}
}

func TestSeedSpillParallelReplayReadError(t *testing.T) {
	b := newSeedMemoryBudget(128<<10, 1)
	b.tempDir = t.TempDir()
	s := newSeedSpillStore(b)
	defer s.close()
	for i := range 12000 {
		if err := s.add(seedAnchor{batchGenomeIndex: uint64(i), length: 31}); err != nil {
			t.Fatal(err)
		}
	}
	// Break the final run after complete batches have already been delivered.
	if err := s.flush(); err != nil {
		t.Fatal(err)
	}
	last := s.runs[len(s.runs)-1]
	info, err := os.Stat(last.path)
	if err != nil {
		t.Fatal(err)
	}
	if err := os.Truncate(last.path, info.Size()-1); err != nil {
		t.Fatal(err)
	}
	var arena seedSearchResultArena
	defer arena.recycle()
	got, err := seedChainingTestIndex(4).chainSeedSpillStore(s, &arena)
	if got != nil || !errors.Is(err, io.ErrUnexpectedEOF) || s.replayBufferBytes != 0 {
		t.Fatalf("failed replay retained results/buffer or lost read error: %v", err)
	}
	if err := s.close(); err != nil {
		t.Fatal(err)
	}
	entries, err := os.ReadDir(b.tempDir)
	if err != nil || len(entries) != 0 || len(b.slots) != 0 {
		t.Fatal("failed parallel replay retained temporary files or a query slot")
	}
}

func TestSeedSpillSmallReplayDoesNotReserveWholeBatch(t *testing.T) {
	s := newSeedSpillStore(newSeedMemoryBudget(64<<20, 1))
	defer s.close()
	for i := range 3 {
		if err := s.add(seedAnchor{batchGenomeIndex: uint64(i), length: 31}); err != nil {
			t.Fatal(err)
		}
	}
	var arena seedSearchResultArena
	defer arena.recycle()
	got, err := seedChainingTestIndex(4).chainSeedSpillStore(s, &arena)
	if err != nil {
		t.Fatal(err)
	}
	defer recycleSeedSearchResultSlice(got)
	if len(*got) != 3 || s.peakBufferBytes > 16<<10 {
		t.Fatalf("three-anchor query reserved a large replay batch: genomes=%d, bytes=%d", len(*got), s.peakBufferBytes)
	}
}

func BenchmarkSeedSpillGrouping(b *testing.B) {
	rng := rand.New(rand.NewPCG(17, 29))
	input := make([]seedAnchor, 1<<20)
	for i := range input {
		input[i] = seedAnchor{batchGenomeIndex: rng.Uint64() & ((1 << 34) - 1), qBegin: int32(i), length: 31}
	}
	for _, radix := range []bool{false, true} {
		name := "comparison"
		if radix {
			name = "radix"
		}
		b.Run(name, func(b *testing.B) {
			s := &seedSpillStore{anchors: slices.Clone(input), sortScratch: make([]seedAnchor, len(input))}
			b.ReportAllocs()
			b.ResetTimer()
			for range b.N {
				b.StopTimer()
				copy(s.anchors, input)
				b.StartTimer()
				if radix {
					s.sortAnchors()
				} else {
					slices.SortStableFunc(s.anchors, compareSeedAnchorGenome)
				}
			}
		})
	}
}

func TestSeedSpillStoreRoundTripAndCleanup(t *testing.T) {
	for _, limit := range []int64{24, 96, 24 * 512, 1 << 20} {
		t.Run(strconv.FormatInt(limit, 10), func(t *testing.T) {
			b := newSeedMemoryBudget(limit, 1)
			b.tempDir = t.TempDir()
			s := newSeedSpillStore(b)
			defer s.close()
			var spillStarts int
			s.onFirstSpill = func() {
				spillStarts++
				if s.spilled || len(s.anchors) == 0 {
					t.Fatal("spill start reported after writing or without buffered anchors")
				}
			}
			var want []seedAnchor
			for i := 0; i < 513; i++ {
				a := seedAnchor{batchGenomeIndex: uint64((i * 19) % 7), qBegin: int32(i), tBegin: int32(1000 - i), length: uint8(15 + i%17), qrc: i&1 != 0, trc: i&2 != 0}
				want = append(want, a)
				if err := s.add(a); err != nil {
					t.Fatal(err)
				}
				if s.peakBufferBytes > limit {
					t.Fatalf("collection buffer %d exceeds %d", s.peakBufferBytes, limit)
				}
			}
			slices.SortStableFunc(want, compareSeedAnchorGenome)
			var got []seedAnchor
			if err := s.consume(func(a seedAnchor) error { got = append(got, a); return nil }); err != nil {
				t.Fatal(err)
			}
			if !reflect.DeepEqual(got, want) {
				t.Fatal("run merge lost, reordered or changed anchors")
			}
			if (limit < 1<<20) != s.spilled {
				t.Fatal("unexpected spill mode")
			}
			if s.spilled && spillStarts != 1 || !s.spilled && spillStarts != 0 {
				t.Fatalf("first-spill notification count: %d (spilled=%t)", spillStarts, s.spilled)
			}
			if !s.spilled && s.dir != "" {
				t.Fatal("below threshold created temporary files")
			}
			dir := s.dir
			if err := s.close(); err != nil {
				t.Fatal(err)
			}
			if dir != "" {
				if _, err := os.Stat(dir); !os.IsNotExist(err) {
					t.Fatal("temporary directory remains")
				}
			}
			if len(b.slots) != 0 {
				t.Fatal("query slot remains held")
			}
			if s.onFirstSpill != nil {
				t.Fatal("closed store retains its spill notification callback")
			}
		})
	}
}

func TestSeedSpillStoreErrorsReleaseSlot(t *testing.T) {
	for _, kind := range []string{"create", "read", "missing-record", "callback"} {
		t.Run(kind, func(t *testing.T) {
			b := newSeedMemoryBudget(24, 1)
			b.tempDir = t.TempDir()
			s := newSeedSpillStore(b)
			if kind == "create" {
				b.tempDir = filepath.Join(b.tempDir, "missing")
			}
			if err := s.add(seedAnchor{length: 31}); err != nil {
				t.Fatal(err)
			}
			err := s.add(seedAnchor{qBegin: 1, length: 31})
			sentinel := errors.New("consumer failure")
			if kind == "read" || kind == "missing-record" {
				if err != nil {
					t.Fatal(err)
				}
				if err = os.Truncate(s.runs[0].path, func() int64 {
					if kind == "missing-record" {
						return 0
					}
					return 1
				}()); err != nil {
					t.Fatal(err)
				}
			}
			if kind != "create" {
				err = s.consume(func(seedAnchor) error {
					if kind == "callback" {
						return sentinel
					}
					return nil
				})
			}
			if err == nil {
				t.Fatal("expected I/O/consumer failure")
			}
			if (kind == "read" || kind == "missing-record") && !errors.Is(err, io.ErrUnexpectedEOF) {
				t.Fatalf("wrong read error: %v", err)
			}
			if kind == "callback" && !errors.Is(err, sentinel) {
				t.Fatalf("wrong callback error: %v", err)
			}
			if err = s.close(); err != nil {
				t.Fatal(err)
			}
			entries, _ := os.ReadDir(b.tempDir)
			if len(entries) != 0 {
				t.Fatal("failed query left spill files")
			}
			if len(b.slots) != 0 {
				t.Fatal("failed query retained a slot")
			}
		})
	}
}

func TestSeedMemoryBudgetConcurrentShares(t *testing.T) {
	b := newSeedMemoryBudget(24*400, 4)
	b.tempDir = t.TempDir()
	var wg sync.WaitGroup
	for range 8 {
		wg.Add(1)
		go func() {
			defer wg.Done()
			s := newSeedSpillStore(b)
			defer s.close()
			for range 120 {
				if err := s.add(seedAnchor{length: 31}); err != nil {
					t.Error(err)
					return
				}
			}
			if s.peakBufferBytes > 24*100 {
				t.Error("query exceeded its share")
			}
			if err := s.consume(func(seedAnchor) error { return nil }); err != nil {
				t.Error(err)
			}
		}()
	}
	wg.Wait()
	if len(b.slots) != 0 {
		t.Fatal("query slot leak")
	}
}

func TestSeedSpillChainingAndCutoffTies(t *testing.T) {
	for _, limit := range []int64{24, 24 * 17, 1 << 20} {
		idx := seedChainingTestIndex(4)
		original := seedChainingTestResults(100)
		b := newSeedMemoryBudget(limit, 1)
		b.tempDir = t.TempDir()
		s := newSeedSpillStore(b)
		// Interleave genomes so every candidate's anchors can cross many runs.
		for anchor := 0; anchor < 4; anchor++ {
			for _, batch := range original {
				for _, r := range batch {
					if anchor >= len(r.Subs) {
						continue
					}
					p := r.Subs[anchor]
					if err := s.add(seedAnchor{batchGenomeIndex: r.BatchGenomeIndex, qBegin: p.QBegin, tBegin: p.TBegin, length: p.Len, qrc: p.QRC, trc: p.TRC}); err != nil {
						t.Fatal(err)
					}
				}
			}
		}
		var arena seedSearchResultArena
		got, err := idx.chainSeedSpillStore(s, &arena)
		if err != nil {
			t.Fatal(err)
		}
		want := idx.chainSeedResults(original, 100)
		for _, topN := range []int{0, 1, 10} {
			w, g := slices.Clone(*want), slices.Clone(*got)
			trimSeedSearchResults(&w, topN)
			trimSeedSearchResults(&g, topN)
			snapshots := func(rs []*seedSearchResult) map[uint64]any {
				m := make(map[uint64]any)
				for _, r := range rs {
					m[r.BatchGenomeIndex] = struct {
						Score   float32
						Regions []seedChainRegion
					}{r.Score, slices.Clone(r.chainRegions)}
				}
				return m
			}
			if !reflect.DeepEqual(snapshots(w), snapshots(g)) {
				t.Fatalf("limit %d topN %d changed chaining bounds/scores/cutoff ties", limit, topN)
			}
		}
		recycleSeedSearchResultSlice(got)
		recycleSeedSearchResultSlice(want)
		arena.recycle()
		if err := s.close(); err != nil {
			t.Fatal(err)
		}
	}
}

func TestSeedSpillEncodingUsesReusableIOBuffers(t *testing.T) {
	w := bufio.NewWriterSize(io.Discard, 64<<10)
	a := seedAnchor{length: 31, qrc: true}
	if n := testing.AllocsPerRun(100, func() {
		for range 100 {
			if err := writeSeedSpillAnchor(w, a); err != nil {
				t.Fatal(err)
			}
		}
		if err := w.Flush(); err != nil {
			t.Fatal(err)
		}
	}); n != 0 {
		t.Fatalf("spill encoding allocates %.1f objects per batch", n)
	}
}

func TestSeedMemoryRejectsInvalidLibraryLimits(t *testing.T) {
	for _, limit := range []int64{-1, 1, seedAnchorBytes - 1} {
		index, err := NewIndexSearcher("", &IndexSearchingOptions{MaxSeedMemory: limit})
		if index != nil || err == nil || !strings.Contains(err.Error(), "MaxSeedMemory") {
			t.Fatalf("invalid seed memory limit %d was not rejected: %v", limit, err)
		}
	}
}
