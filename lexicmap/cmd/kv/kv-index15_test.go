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

package kv

import (
	"bytes"
	"math"
	"os"
	"path/filepath"
	"reflect"
	"testing"
)

func TestIndex15Encoding(t *testing.T) {
	if _, err := anchorPrefixForPartitions(1); err == nil {
		t.Fatal("expected partitions below 4 to fail")
	}
	if prefix, err := anchorPrefixForPartitions(4096); err != nil || prefix != 6 {
		t.Fatalf("unexpected anchor prefix for 4096 partitions: %d, %v", prefix, err)
	}
	if primary, secondary, err := index15PrefixLengths(16, 7, 5); err != nil || primary != 12 || secondary != 14 {
		t.Fatalf("unexpected dynamic prefix lengths: %d/%d, %v", primary, secondary, err)
	}
	if _, _, err := index15PrefixLengths(31, 7, 1); err == nil {
		t.Fatal("expected a 5-byte checkpoint-suffix overflow")
	}
	if _, _, err := index15PrefixLengths(13, 7, 6); err == nil {
		t.Fatal("expected a missing secondary-prefix error")
	}
	for _, offset := range []uint64{0, 1, 255, 256, 1<<32 + 7, index15RawOffsetMask >> 1} {
		seedsOffset, err := MakeSeedsOffset(offset, false)
		if err != nil || seedsOffset>>1 != offset || IsIndex15Offset(seedsOffset) {
			t.Fatalf("invalid seeds offset round trip for %d: %d, %v", offset, seedsOffset, err)
		}
		secondOffset, err := MakeSeedsOffset(offset, true)
		if err != nil || secondOffset>>1 != offset || secondOffset&1 == 0 {
			t.Fatalf("invalid second-k-mer offset round trip for %d: %d, %v", offset, secondOffset, err)
		}
		idx15Offset, err := MakeIndex15Offset(offset)
		if err != nil || !IsIndex15Offset(idx15Offset) || RawIndexOffset(idx15Offset) != offset {
			t.Fatalf("invalid idx15 offset round trip for %d: %d, %v", offset, idx15Offset, err)
		}
	}
	if _, err := MakeSeedsOffset(index15RawOffsetMask>>1+1, false); err == nil {
		t.Fatal("expected an error for an overflowing seeds offset")
	}
	if _, err := MakeIndex15Offset(index15OffsetTag); err == nil {
		t.Fatal("expected an error for an overflowing idx15 offset")
	}

	for width := uint8(1); width <= 8; width++ {
		value := uint64(1)<<(width*8-1) | uint64(width)
		buf := make([]byte, width)
		putUintWidth(buf, value, width)
		if got := uintWidth(buf, width); got != value {
			t.Fatalf("%d-byte integer round trip: got %d, expected %d", width, got, value)
		}
	}

	for _, value := range []uint64{0, 1, 1<<38 - 1} {
		buf := make([]byte, index15SuffixBytes)
		putUint40(buf, value)
		if got := uint40(buf); got != value {
			t.Fatalf("40-bit integer round trip: got %d, expected %d", got, value)
		}
	}
	for _, test := range []struct {
		value uint64
		width uint8
	}{{0, 1}, {255, 1}, {256, 2}, {1<<24 - 1, 3}, {1 << 24, 4}, {math.MaxUint64, 8}} {
		if got := index15OffsetWidth(test.value); got != test.width {
			t.Fatalf("offset width for %d: got %d, expected %d", test.value, got, test.width)
		}
	}

	var encoded bytes.Buffer
	block := index15BuildBlock{pairFirstCheckpoint: 42, maxRelOffset: 256}
	width, size, err := writeIndex15Block(&encoded, &block, make([]byte, index15BlockFixedBytes+15*8))
	if err != nil {
		t.Fatal(err)
	}
	if width != 2 || size != 124 || encoded.Len() != 124 {
		t.Fatalf("unexpected block encoding size: width=%d size=%d bytes=%d", width, size, encoded.Len())
	}
	if checkpoint := be.Uint64(encoded.Bytes()[11:19]); checkpoint != block.pairFirstCheckpoint {
		t.Fatalf("unexpected pair-first checkpoint: %d", checkpoint)
	}
}

func TestIndex15LookupValidation(t *testing.T) {
	data := make([]byte, index15BlockFixedBytes+15)
	idx := &index15{data: data, k: 16, primaryPrefix: 13, secondaryPrefix: 15, seedsSize: 1000}
	tagged, err := MakeIndex15Offset(0)
	if err != nil {
		t.Fatal(err)
	}
	if _, _, _, err = idx.lookup(tagged, 0, 15); err == nil {
		t.Fatal("expected an invalid offset-width error")
	}
	data[0] = 1
	be.PutUint64(data[3:11], 100)
	be.PutUint64(data[11:19], 42)
	be.PutUint16(data[1:3], 1)
	checkpoint, offset, empty, err := idx.lookup(tagged, 0, 15)
	if err != nil || empty || checkpoint != 42 || offset != 100 {
		t.Fatalf("unexpected base lookup: checkpoint=%d offset=%d empty=%v err=%v", checkpoint, offset, empty, err)
	}
	be.PutUint16(data[1:3], 1<<1)
	data[index15BlockFixedBytes] = 255
	be.PutUint64(data[3:11], 900)
	if _, _, _, err = idx.lookup(tagged, 4, 15); err == nil {
		t.Fatal("expected an out-of-range relative-offset error")
	}
	truncated, err := MakeIndex15Offset(1)
	if err != nil {
		t.Fatal(err)
	}
	if _, _, _, err = idx.lookup(truncated, 0, 15); err == nil {
		t.Fatal("expected a truncated-block error")
	}
}

type index15ResultSnapshot struct {
	IQuery   int
	IQuery2  int
	Len      uint8
	IsSuffix bool
	Values   []uint64
}

func snapshotIndex15Results(results *[]SearchResult, includeQuery2 bool) []index15ResultSnapshot {
	snapshots := make([]index15ResultSnapshot, len(*results))
	for i, result := range *results {
		snapshots[i] = index15ResultSnapshot{
			IQuery:   result.IQuery,
			Len:      result.Len,
			IsSuffix: result.IsSuffix,
			Values:   append([]uint64(nil), result.Values...),
		}
		if includeQuery2 {
			snapshots[i].IQuery2 = result.IQuery2
		}
	}
	RecycleSearchResults(results)
	return snapshots
}

func collectIndex15Searches(t *testing.T, searcher *Searcher, queries []uint64) [][]index15ResultSnapshot {
	t.Helper()
	all := make([][]index15ResultSnapshot, 0, len(queries)*3+1)
	for _, prefixLength := range []uint8{13, 15, 16} {
		for _, query := range queries {
			results, err := searcher.Search([]uint64{query}, prefixLength, false, false)
			if err != nil {
				t.Fatal(err)
			}
			all = append(all, snapshotIndex15Results(results, false))
		}
	}
	queryList := append([]uint64(nil), queries...)
	queryLists := [][]uint64{queryList}
	results, err := searcher.Search2(queryLists, 15, false, false)
	if err != nil {
		t.Fatal(err)
	}
	all = append(all, snapshotIndex15Results(results, true))
	return all
}

func TestCreateKVIndex15PreservesSearchResults(t *testing.T) {
	dir := t.TempDir()
	file := filepath.Join(dir, "seeds")
	const k = uint8(16)
	const maskPrefix = uint8(7)
	const anchorPrefix = uint8(6)

	prefix := uint64(0b0110110) << ((k - maskPrefix) << 1)
	suffixes := []uint64{0, 1, 4, 20, 63, 64, 68, 128, 132, 191}
	dataOfMask := make(map[uint64]*[]uint64, len(suffixes))
	queries := make([]uint64, 0, len(suffixes)+2)
	for i, suffix := range suffixes {
		kmer := prefix | suffix
		values := make([]uint64, 700)
		for j := range values {
			values[j] = uint64(i + 1)
		}
		dataOfMask[kmer] = &values
		queries = append(queries, kmer)
	}
	queries = append(queries, prefix|8, prefix|72) // empty 15-bp subprefixes
	data := []*map[uint64]*[]uint64{&dataOfMask}
	if _, err := WriteKVData(k, 0, data, file, maskPrefix, anchorPrefix, 512, false); err != nil {
		t.Fatal(err)
	}
	seedsBefore, err := os.ReadFile(file)
	if err != nil {
		t.Fatal(err)
	}

	legacy, err := NewSearcher(file, 1)
	if err != nil {
		t.Fatal(err)
	}
	want := collectIndex15Searches(t, legacy, queries)
	if err = legacy.Close(); err != nil {
		t.Fatal(err)
	}
	if _, err = CreateKVIndex15(file, 4096, MinIndex15Threshold-1); err == nil {
		t.Fatal("expected a threshold below 4 KiB to fail")
	}

	var progress [][2]uint64
	stats, err := CreateKVIndex15WithProgress(file, 4096, math.MaxUint64, func(processedMasks, totalMasks uint64) {
		progress = append(progress, [2]uint64{processedMasks, totalMasks})
	})
	if err != nil {
		t.Fatal(err)
	}
	if expected := [][2]uint64{{0, 1}, {1, 1}}; !reflect.DeepEqual(progress, expected) {
		t.Fatalf("unexpected mask progress: got %v, expected %v", progress, expected)
	}
	if stats.IndexedBlocks != 0 || stats.Index15Bytes != index15HeaderSize {
		t.Fatalf("unexpected high-threshold stats: %+v", stats)
	}
	unindexed, err := NewSearcher(file, 1)
	if err != nil {
		t.Fatal(err)
	}
	if got := collectIndex15Searches(t, unindexed, queries); !reflect.DeepEqual(got, want) {
		t.Fatalf("unindexed tagged search differs from legacy:\ngot:  %#v\nwant: %#v", got, want)
	}
	if err = unindexed.Close(); err != nil {
		t.Fatal(err)
	}

	stats, err = CreateKVIndex15(file, 4096, MinIndex15Threshold)
	if err != nil {
		t.Fatal(err)
	}
	if stats.IndexedBlocks == 0 || stats.Index15Bytes <= index15HeaderSize {
		t.Fatalf("unexpected minimum-threshold stats: %+v", stats)
	}
	indexed, err := NewSearcher(file, 1)
	if err != nil {
		t.Fatal(err)
	}
	if got := collectIndex15Searches(t, indexed, queries); !reflect.DeepEqual(got, want) {
		t.Fatalf("idx15 search differs from legacy:\ngot:  %#v\nwant: %#v", got, want)
	}
	// Streaming must also follow tagged offsets and empty secondary buckets.
	for _, prefixLength := range []uint8{13, 15, 16} {
		queries2 := [][]uint64{queries}
		collected, err := indexed.Search2(queries2, prefixLength, false, false)
		if err != nil {
			t.Fatal(err)
		}
		expected := seedDataSnapshot(*collected)
		RecycleSearchResults(collected)
		var streamed []SearchResult
		err = indexed.Search2Stream(queries2, prefixLength, false, false, func(r SearchResult) error {
			if len(r.Values) > seedPosBatchSize+1 {
				t.Fatal("unbounded idx15 seed data batch")
			}
			r.Values = append([]uint64(nil), r.Values...)
			streamed = append(streamed, r)
			return nil
		})
		if err != nil || !reflect.DeepEqual(expected, seedDataSnapshot(streamed)) {
			t.Fatalf("idx15 stream differs at prefix %d: %v", prefixLength, err)
		}
	}
	if err = indexed.Close(); err != nil {
		t.Fatal(err)
	}

	stats, err = CreateKVIndex15(file, 1024, MinIndex15Threshold)
	if err != nil {
		t.Fatal(err)
	}
	dynamic, err := NewSearcher(file, 1)
	if err != nil {
		t.Fatal(err)
	}
	if dynamic.index15.primaryPrefix != 12 || dynamic.index15.secondaryPrefix != 14 {
		t.Fatalf("unexpected prefixes for 1024 partitions: %d/%d", dynamic.index15.primaryPrefix, dynamic.index15.secondaryPrefix)
	}
	if got := collectIndex15Searches(t, dynamic, queries); !reflect.DeepEqual(got, want) {
		t.Fatalf("dynamic-prefix idx15 search differs from legacy:\ngot:  %#v\nwant: %#v", got, want)
	}
	if err = dynamic.Close(); err != nil {
		t.Fatal(err)
	}

	inMemory, err := NewInMemomrySearcher(file)
	if err != nil {
		t.Fatal(err)
	}
	queryList := append([]uint64(nil), queries...)
	queryLists := [][]uint64{queryList}
	inMemoryResults, err := inMemory.Search2(queryLists, 15, false, false)
	if err != nil {
		t.Fatal(err)
	}
	if got := snapshotIndex15Results(inMemoryResults, true); !reflect.DeepEqual(got, want[len(want)-1]) {
		t.Fatalf("in-memory search with the tagged primary index differs from legacy:\ngot:  %#v\nwant: %#v", got, want[len(want)-1])
	}
	if err = inMemory.Close(); err != nil {
		t.Fatal(err)
	}

	seedsAfter, err := os.ReadFile(file)
	if err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(seedsAfter, seedsBefore) {
		t.Fatal("reindexing changed the seeds file")
	}
	idx15Info, err := os.Stat(file + KVIndex15FileExt)
	if err != nil {
		t.Fatal(err)
	}
	if uint64(idx15Info.Size()) != stats.Index15Bytes {
		t.Fatalf("idx15 size mismatch: got %d, expected %d", idx15Info.Size(), stats.Index15Bytes)
	}
	idx15File := file + KVIndex15FileExt
	missingFile := idx15File + ".missing"
	if err = os.Rename(idx15File, missingFile); err != nil {
		t.Fatal(err)
	}
	if _, err = NewSearcher(file, 1); err == nil {
		t.Fatal("expected opening a tagged index without idx15 to fail")
	}
	if err = os.Rename(missingFile, idx15File); err != nil {
		t.Fatal(err)
	}
	starts, err := ReadKVIndexStarts(file + KVIndexFileExt)
	if err != nil {
		t.Fatal(err)
	}
	if len(starts) != 1 || starts[0][1] == 0 || starts[0][1]&1 != 0 {
		t.Fatalf("unexpected normalized mask start: %v", starts)
	}

	progress = progress[:0]
	if err = CreateKVIndexWithProgress(file, 4096, func(processedMasks, totalMasks uint64) {
		progress = append(progress, [2]uint64{processedMasks, totalMasks})
	}); err != nil {
		t.Fatal(err)
	}
	if expected := [][2]uint64{{0, 1}, {1, 1}}; !reflect.DeepEqual(progress, expected) {
		t.Fatalf("unexpected primary-index mask progress: got %v, expected %v", progress, expected)
	}
	if _, err := os.Stat(idx15File); !os.IsNotExist(err) {
		t.Fatalf("obsolete secondary index was not removed: %v", err)
	}
	// Rebuilding a single-level index again must also work without a sidecar.
	if err := CreateKVIndex(file, 4096); err != nil {
		t.Fatal(err)
	}
	legacyReindexed, err := NewSearcher(file, 1)
	if err != nil {
		t.Fatal(err)
	}
	if got := collectIndex15Searches(t, legacyReindexed, queries); !reflect.DeepEqual(got, want) {
		t.Fatalf("reindexed primary search differs from legacy:\ngot:  %#v\nwant: %#v", got, want)
	}
	if err = legacyReindexed.Close(); err != nil {
		t.Fatal(err)
	}
}
