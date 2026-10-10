package kv

import (
	"bytes"
	"errors"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"reflect"
	"slices"
	"testing"
)

// readComparableIndex15File removes the independently generated header pairing token.
func readComparableIndex15File(t *testing.T, file, ext string) []byte {
	t.Helper()
	data, err := os.ReadFile(file + ext)
	if err != nil {
		t.Fatal(err)
	}
	if ext == KVIndexFileExt {
		clear(data[14:16])
	} else if ext == KVIndex15FileExt {
		clear(data[46:48])
	}
	return data
}

// TestIndex15WriterCloseError propagates secondary flush failures and closes every output.
func TestIndex15WriterCloseError(t *testing.T) {
	file := filepath.Join(t.TempDir(), "seeds")
	writer, err := NewWriterWithIndex15(16, 0, 1, file, 7, 6, false, MinIndex15Threshold)
	if err != nil {
		t.Fatal(err)
	}
	if err := writer.WriteDataOfAMask(map[uint64]*[]uint64{0: &[]uint64{1}}); err != nil {
		t.Fatal(err)
	}
	if err := writer.fh15.Close(); err != nil {
		t.Fatal(err)
	}
	if err := writer.Close(); !errors.Is(err, os.ErrClosed) {
		t.Fatalf("secondary flush failure was not returned: %v", err)
	}
	for _, file := range []*os.File{writer.fh, writer.fhi, writer.fh15} {
		if _, err := file.Stat(); !errors.Is(err, os.ErrClosed) {
			t.Fatalf("output descriptor remains open after failure: %v", err)
		}
	}
}

// TestLegacyWriterRevisitedAnchor preserves the original prefix-transition rule.
func TestLegacyWriterRevisitedAnchor(t *testing.T) {
	file := filepath.Join(t.TempDir(), "seeds")
	writer, err := NewWriter(16, 0, 1, file, 1, 1, false)
	if err != nil {
		t.Fatal(err)
	}
	// The leading mask prefix changes, revisiting anchor zero in the next pair.
	const revisited = uint64(1) << 30
	data := map[uint64]*[]uint64{0: &[]uint64{1}, 1 << 28: &[]uint64{2},
		revisited: &[]uint64{3}, revisited | 2<<28: &[]uint64{4}}
	if err := writer.WriteDataOfAMask(data); err != nil {
		t.Fatal(err)
	}
	if err := writer.Close(); err != nil {
		t.Fatal(err)
	}
	_, _, indexes, _, _, _, err := ReadKVIndex(file + KVIndexFileExt)
	if err != nil {
		t.Fatal(err)
	}
	if indexes[0][2] != revisited || indexes[0][3]&1 != 0 {
		t.Fatalf("revisited anchor lost its last prefix transition: %v", indexes[0])
	}
}

// TestIndex15WriterMatchesReindex covers empty masks, odd pairs, prefix crossings,
// large postings, both position widths, and blocks above and below the threshold.
func TestIndex15WriterMatchesReindex(t *testing.T) {
	const k, maskPrefix, anchorPrefix = uint8(16), uint8(7), uint8(6)
	prefix := uint64(0b0110110) << ((k - maskPrefix) << 1)
	masks := make([]map[uint64]*[]uint64, 6)
	for i := range masks {
		masks[i] = make(map[uint64]*[]uint64)
	}
	masks[1][prefix] = &[]uint64{}
	for _, mask := range []int{2, 3} {
		for _, suffix := range []uint64{0, 1, 4, 20, 63, 64, 68, 128, 132, 191} {
			values := make([]uint64, 700)
			for i := range values {
				values[i] = uint64(i % 257) // preserve duplicate and unsorted postings
			}
			masks[mask][prefix|suffix] = &values
		}
	}
	masks[3][prefix|192] = &[]uint64{17}
	large := make([]uint64, 70000)
	masks[4][prefix] = &large
	queries := make([][]uint64, len(masks))
	for i, data := range masks {
		for key := range data {
			queries[i] = append(queries[i], key)
		}
		queries[i] = append(queries[i], prefix|8, prefix|72)
		slices.Sort(queries[i])
	}
	for _, narrow := range []bool{false, true} {
		for _, threshold := range []uint64{MinIndex15Threshold, math.MaxUint64} {
			t.Run(fmt.Sprintf("narrow=%t/threshold=%d", narrow, threshold), func(t *testing.T) {
				dir := t.TempDir()
				plain, direct, merged := filepath.Join(dir, "plain"), filepath.Join(dir, "direct"), filepath.Join(dir, "merged")
				for _, file := range []string{plain, direct} {
					limit := uint64(0)
					if file == direct {
						limit = threshold
					}
					writer, err := NewWriterWithIndex15(k, 7, len(masks), file, maskPrefix, anchorPrefix, narrow, limit)
					if err != nil {
						t.Fatal(err)
					}
					for _, data := range masks {
						if err := writer.WriteDataOfAMask(data); err != nil {
							t.Fatal(err)
						}
					}
					if err := writer.Close(); err != nil {
						t.Fatal(err)
					}
				}
				search := func(file string) [][]index15ResultSnapshot {
					searcher, err := NewSearcher(file, 1)
					if err != nil {
						t.Fatal(err)
					}
					defer searcher.Close()
					var snapshots [][]index15ResultSnapshot
					for _, length := range []uint8{13, 15, 16} {
						results, err := searcher.Search2(queries, length, false, false)
						if err != nil {
							t.Fatal(err)
						}
						snapshots = append(snapshots, snapshotIndex15Results(results, true))
					}
					return snapshots
				}
				wantSearch := search(plain)
				reader, err := NewReader(plain)
				if err != nil {
					t.Fatal(err)
				}
				writer, err := NewWriterWithIndex15(k, 7, len(masks), merged, maskPrefix, anchorPrefix, narrow, threshold)
				if err != nil {
					t.Fatal(err)
				}
				merger, err := NewMerger([]*Reader{reader}, writer)
				if err != nil {
					t.Fatal(err)
				}
				for range masks {
					if err := merger.WriteMask(); err != nil {
						t.Fatal(err)
					}
				}
				if err := reader.Close(); err != nil {
					t.Fatal(err)
				}
				if err := writer.Close(); err != nil {
					t.Fatal(err)
				}
				stats, err := CreateKVIndex15(plain, 4096, threshold)
				if err != nil {
					t.Fatal(err)
				}
				if threshold == MinIndex15Threshold && stats.IndexedBlocks == 0 {
					t.Fatal("fixture produced no indexed blocks")
				}
				for _, file := range []string{direct, merged} {
					for _, ext := range []string{"", KVIndexFileExt, KVIndex15FileExt} {
						if !bytes.Equal(readComparableIndex15File(t, plain, ext), readComparableIndex15File(t, file, ext)) {
							t.Fatalf("%s%s differs from standalone reindexing", file, ext)
						}
					}
					if got := search(file); !reflect.DeepEqual(got, wantSearch) {
						t.Fatalf("%s search differs from the ordinary index", file)
					}
				}
			})
		}
	}
}

// TestIndex15WriterDoesNotReopenSeeds removes the output path while it is open.
// Both indexes must finish using only the metadata supplied during seed writes.
func TestIndex15WriterDoesNotReopenSeeds(t *testing.T) {
	file := filepath.Join(t.TempDir(), "seeds")
	writer, err := NewWriterWithIndex15(16, 0, 1, file, 7, 6, false, MinIndex15Threshold)
	if err != nil {
		t.Fatal(err)
	}
	moved := file + ".moved"
	if err := os.Rename(file, moved); err != nil {
		t.Fatal(err)
	}
	values := make([]uint64, 700)
	if err := writer.WriteDataOfAMask(map[uint64]*[]uint64{0: &values}); err != nil {
		t.Fatal(err)
	}
	if err := writer.Close(); err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(file); !os.IsNotExist(err) {
		t.Fatalf("recreated seeds path: %v", err)
	}
	data, err := os.ReadFile(moved)
	if err != nil {
		t.Fatal(err)
	}
	index, err := os.ReadFile(file + KVIndex15FileExt)
	if err != nil || len(index) <= index15HeaderSize || be.Uint64(index[32:40]) != uint64(len(data)) {
		t.Fatalf("secondary index was not completed without reopening seeds: %v", err)
	}
}
