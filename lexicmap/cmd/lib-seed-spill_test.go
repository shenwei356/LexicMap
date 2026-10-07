package cmd

import (
	"bufio"
	"errors"
	"io"
	"os"
	"path/filepath"
	"reflect"
	"slices"
	"strconv"
	"strings"
	"sync"
	"testing"
)

func TestSeedSpillStoreRoundTripAndCleanup(t *testing.T) {
	for _, limit := range []int64{24, 96, 24 * 512, 1 << 20} {
		t.Run(strconv.FormatInt(limit, 10), func(t *testing.T) {
			b := newSeedMemoryBudget(limit, 1)
			b.tempDir = t.TempDir()
			s := newSeedSpillStore(b)
			defer s.close()
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
