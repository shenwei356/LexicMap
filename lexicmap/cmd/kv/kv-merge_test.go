package kv

import (
	"bytes"
	"math"
	"math/rand/v2"
	"os"
	"path/filepath"
	"testing"
)

func mergeTestInputs() [][]map[uint64]*[]uint64 {
	random := rand.New(rand.NewPCG(17, 29))
	inputs := make([][]map[uint64]*[]uint64, 5)
	for file := range inputs {
		inputs[file] = make([]map[uint64]*[]uint64, 40)
		for mask := range inputs[file] {
			data := make(map[uint64]*[]uint64)
			inputs[file][mask] = data
			switch mask {
			case 0: // all empty
				continue
			case 1: // one key, repeated and deliberately unsorted positions
				data[0] = &[]uint64{uint64(file + 9), 1, 1, 0}
				continue
			case 2: // even key count, anchors in both halves of a pair
				data[1<<58] = &[]uint64{}
				data[2<<58] = &[]uint64{uint64(file)}
				continue
			case 3: // odd key count and large deltas
				data[0] = &[]uint64{0}
				data[math.MaxUint64-1] = &[]uint64{math.MaxUint64}
				data[math.MaxUint64] = &[]uint64{math.MaxUint64 - 1}
				continue
			case 4: // keys only present in one input, empty inputs mixed in
				if file%2 == 0 {
					data[uint64(file)] = &[]uint64{uint64(file)}
				}
				continue
			case 5: // postings exceed both reader and writer buffers
				values := make([]uint64, 6000)
				for i := range values {
					values[i] = uint64((6000 - i) % 257)
				}
				data[99] = &values
				continue
			}
			for i := 0; i < 50; i++ {
				// Shared keys exercise stable merging; others vary by input.
				key := random.Uint64() >> 2
				if i%3 == 0 {
					key = uint64(i)<<58 | uint64(i)
				}
				if _, exists := data[key]; exists {
					continue
				}
				values := make([]uint64, random.IntN(12))
				for j := range values {
					values[j] = random.Uint64()
				}
				data[key] = &values
			}
		}
	}
	return inputs
}

func writeMergeTestFile(t *testing.T, file string, masks []map[uint64]*[]uint64, threeBytes bool) {
	t.Helper()
	writer, err := NewWriter(32, 7, len(masks), file, 1, 2, threeBytes)
	if err != nil {
		t.Fatal(err)
	}
	for _, data := range masks {
		if err = writer.WriteDataOfAMask(data); err != nil {
			t.Fatal(err)
		}
	}
	if err = writer.Close(); err != nil {
		t.Fatal(err)
	}
}

func mergeTestFiles(t *testing.T, files []string, output string, masks int, threeBytes, streaming bool) {
	t.Helper()
	readers := make([]*Reader, len(files))
	for i, file := range files {
		var err error
		readers[i], err = NewReader(file)
		if err != nil {
			t.Fatal(err)
		}
		defer readers[i].Close()
	}
	writer, err := NewWriter(32, 7, masks, output, 1, 2, threeBytes)
	if err != nil {
		t.Fatal(err)
	}
	if streaming {
		merger, err := NewMerger(readers, writer)
		if err != nil {
			t.Fatal(err)
		}
		for i := 0; i < masks; i++ {
			if err = merger.WriteMask(); err != nil {
				t.Fatalf("mask %d: %v", i, err)
			}
		}
	} else {
		for i := 0; i < masks; i++ {
			data := make(map[uint64]*[]uint64)
			for _, reader := range readers {
				if err = reader.ReadDataOfAMaskAndAppendToMap(&data); err != nil {
					t.Fatal(err)
				}
			}
			if err = writer.WriteDataOfAMask(data); err != nil {
				t.Fatal(err)
			}
		}
	}
	if err = writer.Close(); err != nil {
		t.Fatal(err)
	}
}

func equalMergeTestFiles(t *testing.T, want, got string) {
	t.Helper()
	for _, ext := range []string{"", KVIndexFileExt} {
		expected, err := os.ReadFile(want + ext)
		if err != nil {
			t.Fatal(err)
		}
		actual, err := os.ReadFile(got + ext)
		if err != nil {
			t.Fatal(err)
		}
		if !bytes.Equal(expected, actual) {
			t.Fatalf("merged file differs: %s (%d vs %d bytes)", ext, len(expected), len(actual))
		}
	}
}

func TestMergerMatchesMap(t *testing.T) {
	inputs := mergeTestInputs()
	for _, width := range []struct {
		name   string
		output bool
		mixed  bool
	}{
		{"eight", false, false}, {"seven", true, false},
		{"mixed-to-eight", false, true}, {"mixed-to-seven", true, true},
	} {
		t.Run(width.name, func(t *testing.T) {
			dir := t.TempDir()
			files := make([]string, len(inputs))
			for i, masks := range inputs {
				files[i] = filepath.Join(dir, string(rune('a'+i)))
				narrow := width.output
				if width.mixed && i%2 == 0 {
					narrow = !narrow
				}
				writeMergeTestFile(t, files[i], masks, narrow)
			}
			want, got := filepath.Join(dir, "map"), filepath.Join(dir, "stream")
			mergeTestFiles(t, files, want, len(inputs[0]), width.output, false)
			mergeTestFiles(t, files, got, len(inputs[0]), width.output, true)
			equalMergeTestFiles(t, want, got)
			// Read every mask through the existing decoder after header patches.
			reader, err := NewReader(got)
			if err != nil {
				t.Fatal(err)
			}
			defer reader.Close()
			for range inputs[0] {
				data, err := reader.ReadDataOfAMaskAsMap()
				if err != nil {
					t.Fatal(err)
				}
				RecycleKmerData(data)
			}
			// Recursive merge must concatenate positions in the same order.
			if !width.mixed {
				group1, group2 := filepath.Join(dir, "group1"), filepath.Join(dir, "group2")
				mergeTestFiles(t, files[:2], group1, len(inputs[0]), width.output, true)
				mergeTestFiles(t, files[2:4], group2, len(inputs[0]), width.output, true)
				final := filepath.Join(dir, "recursive")
				mergeTestFiles(t, []string{group1, group2, files[4]}, final, len(inputs[0]), width.output, true)
				equalMergeTestFiles(t, want, final)
			}
		})
	}
}

func TestMergerRejectsBrokenMask(t *testing.T) {
	for _, kind := range []string{"truncated-header", "truncated-postings", "wrong-count", "missing-last", "huge-postings"} {
		t.Run(kind, func(t *testing.T) {
			dir := t.TempDir()
			input := filepath.Join(dir, "input")
			writeMergeTestFile(t, input, []map[uint64]*[]uint64{{0: &[]uint64{1}}}, false)
			data, err := os.ReadFile(input)
			if err != nil {
				t.Fatal(err)
			}
			switch kind {
			case "truncated-header":
				data = data[:35]
			case "truncated-postings":
				data = data[:len(data)-1]
			case "wrong-count":
				be.PutUint64(data[32:40], 2)
			case "missing-last":
				data[40] &^= 128
			case "huge-postings":
				// Pair header for key 0, followed by an oversized value count.
				data = append(data[:40], 192, 0, 0, 7)
				data = be.AppendUint64(data, math.MaxUint64)
				data = append(data, 0)
			}
			if err = os.WriteFile(input, data, 0600); err != nil {
				t.Fatal(err)
			}
			reader, err := NewReader(input)
			if err != nil {
				t.Fatal(err)
			}
			defer reader.Close()
			writer, err := NewWriter(32, 7, 1, filepath.Join(dir, "output"), 1, 2, false)
			if err != nil {
				t.Fatal(err)
			}
			defer writer.Close()
			merger, err := NewMerger([]*Reader{reader}, writer)
			if err != nil {
				t.Fatal(err)
			}
			if err = merger.WriteMask(); err == nil {
				t.Fatal("accepted broken mask")
			}
		})
	}
}
