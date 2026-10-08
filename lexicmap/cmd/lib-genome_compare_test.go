package cmd

import (
	"math/rand"
	"slices"
	"sync"
	"testing"

	rtree "github.com/shenwei356/LexicMap/lexicmap/cmd/tree"
)

// legacyFragmentPairs is an independent oracle for the original uint16
// cross-product counts and count-only, unstable Top-N sorting. Row order is
// normalized by callers, just as both alignment paths sort their returned pairs.
func legacyFragmentPairs(a, b []rtree.BatchEntry, options *FragmentComparatorOptions) []uint64 {
	counts := make(map[uint64]uint16)
	for ia, ib := 0, 0; ia < len(a) && ib < len(b); {
		if a[ia].Key < b[ib].Key {
			ia++
			continue
		}
		if b[ib].Key < a[ia].Key {
			ib++
			continue
		}
		ja, jb := ia+1, ib+1
		for ja < len(a) && a[ja].Key == a[ia].Key {
			ja++
		}
		for jb < len(b) && b[jb].Key == b[ib].Key {
			jb++
		}
		for _, ea := range a[ia:ja] {
			for _, eb := range b[ib:jb] {
				counts[uint64(ea.Val>>1)<<32|uint64(eb.Val>>1)]++
			}
		}
		ia, ib = ja, jb
	}
	pairs := make([]uint64, 0)
	for key, count := range counts {
		if count >= options.MinSharedKmers {
			pairs = append(pairs, key)
		}
	}
	if options.TopNFragments > 0 {
		slices.Sort(pairs)
		rows := make(map[uint64][]uint64)
		for _, key := range pairs {
			rows[key>>32] = append(rows[key>>32], key&0xffffffff)
		}
		pairs = pairs[:0]
		for row, subjects := range rows {
			if len(subjects) > options.TopNFragments {
				slices.SortFunc(subjects, func(a, b uint64) int {
					return int(counts[row<<32|b]) - int(counts[row<<32|a])
				})
				subjects = subjects[:options.TopNFragments]
			}
			for _, subject := range subjects {
				pairs = append(pairs, row<<32|subject)
			}
		}
	}
	slices.Sort(pairs)
	return pairs
}

func TestFragmentPairsMatchLegacy(t *testing.T) {
	var a, b []rtree.BatchEntry
	rng := rand.New(rand.NewSource(61))
	// More than 12 subjects exercises quicksort as well as small insertion sorts.
	// Rows contain ties, unequal counts and repeated positions for the same key.
	for row := uint32(0); row < 7; row++ {
		for subject := uint32(0); subject < 97; subject++ {
			key := uint64(row)*97 + uint64(subject)
			a = append(a, rtree.BatchEntry{Key: key, Val: row << 1})
			count := 1 + rng.Intn(5)
			if row == 0 {
				count = 3 // every subject ties in this row
			}
			for range count {
				b = append(b, rtree.BatchEntry{Key: key, Val: subject<<1 | 1})
			}
		}
	}
	// Repeated occurrences on both sides must contribute the full cross product.
	for range 3 {
		a = append(a, rtree.BatchEntry{Key: 9000, Val: 20 << 1})
	}
	for range 4 {
		b = append(b, rtree.BatchEntry{Key: 9000, Val: 30<<1 | 1})
	}
	// Preserve uint16 wraparound, including a touched pair whose count becomes zero.
	for range 256 {
		a = append(a, rtree.BatchEntry{Key: 9001, Val: 21 << 1})
		b = append(b, rtree.BatchEntry{Key: 9001, Val: 31<<1 | 1})
	}
	for _, threshold := range []uint16{0, 1, 3, 6, 12, 13} {
		for _, topN := range []int{0, 1, 5, 17, 96, 97, 100} {
			options := &FragmentComparatorOptions{K: 11, MinSharedKmers: threshold, TopNFragments: topN}
			cpr := &FragmentComparator{options: options}
			got := cpr.scanPairsMerged(a, b)
			slices.Sort(*got)
			want := legacyFragmentPairs(a, b, options)
			if !slices.Equal(*got, want) {
				t.Fatalf("threshold=%d topN=%d: pairs differ from original selection", threshold, topN)
			}
			RecycleFragmentCompareResult(got)
		}
	}
}

func TestFragmentPairsEmptyAndDisjoint(t *testing.T) {
	cpr := NewFragmentComparator(&FragmentComparatorOptions{K: 11, TopNFragments: 5}, nil)
	a := []rtree.BatchEntry{{Key: 1, Val: 0}}
	b := []rtree.BatchEntry{{Key: 2, Val: 1}}
	for _, inputs := range [][2][]rtree.BatchEntry{{nil, nil}, {a, nil}, {nil, b}, {a, b}} {
		pairs := cpr.scanPairsMerged(inputs[0], inputs[1])
		if len(*pairs) != 0 {
			t.Fatal("empty or disjoint inputs produced fragment pairs")
		}
		RecycleFragmentCompareResult(pairs)
	}
}

func TestFragmentPairsBlocksReuseAndFallback(t *testing.T) {
	cpr := &FragmentComparator{options: &FragmentComparatorOptions{K: 11, MinSharedKmers: 1, TopNFragments: 1}}
	// These dimensions exceed one dense block while containing just four pairs.
	// The last query row reuses the same counter cells as the first block.
	const lastRow = maxFragmentPairCounterCells / 1024
	a := []rtree.BatchEntry{{Key: 1, Val: 0}, {Key: 1, Val: uint32(lastRow) << 1}}
	b := []rtree.BatchEntry{{Key: 1, Val: 1}, {Key: 1, Val: 1023<<1 | 1}}
	_, _, blockRows, ok := fragmentPairCounterShape(a, b)
	if !ok || blockRows != lastRow {
		t.Fatal("fixture must exercise multiple dense row blocks")
	}
	// Very large sparse IDs and more than 16 blocks must use hash counting.
	sparseA := []rtree.BatchEntry{{Key: 1, Val: 0xfffffffe}}
	wideB := []rtree.BatchEntry{{Key: 1, Val: 0xffffffff}}
	manyRows := []rtree.BatchEntry{{Key: 1, Val: uint32(maxFragmentPairCounterCells*maxFragmentPairCounterBlocks) << 1}}
	for _, input := range [][2][]rtree.BatchEntry{{sparseA, b}, {a, wideB}, {manyRows, {{Key: 1, Val: 1}}}} {
		if _, _, _, ok := fragmentPairCounterShape(input[0], input[1]); ok {
			t.Fatal("large or sparse ID dimensions selected dense rescanning")
		}
	}
	for _, input := range [][2][]rtree.BatchEntry{{a, b}, {sparseA, b}, {a, wideB}, {a[:1], b[:1]}, {a, b}, {nil, b}, {manyRows, {{Key: 1, Val: 1}}}} {
		originalA, originalB := slices.Clone(input[0]), slices.Clone(input[1])
		pairs := cpr.scanPairsMerged(input[0], input[1])
		slices.Sort(*pairs)
		if !slices.Equal(*pairs, legacyFragmentPairs(input[0], input[1], cpr.options)) {
			t.Fatal("block boundaries, scratch reuse or hash fallback changed pair selection")
		}
		RecycleFragmentCompareResult(pairs)
		if !slices.Equal(input[0], originalA) || !slices.Equal(input[1], originalB) {
			t.Fatal("pair scanning modified shared sorted entries")
		}
		bytes := cap(cpr.pairCounter)*2 + cap(cpr.pairVisited)*8
		if bytes > maxFragmentPairCounterBytes {
			t.Fatalf("retained counter scratch %d exceeds budget", bytes)
		}
		for _, mask := range cpr.pairVisited {
			if mask != 0 {
				t.Fatal("visited bitmap word was not reset")
			}
		}
	}
	// Check the entire retained backing array after reuse, including tails beyond
	// the current length after a smaller comparison, without rescanning it per case.
	for _, count := range cpr.pairCounter[:cap(cpr.pairCounter)] {
		if count != 0 {
			t.Fatal("visited counter cell was not reset")
		}
	}
}

// Default-size fragments from two 20 Mb genomes must stay on the dense path
// within both the counter memory budget and the maximum number of rescans.
func TestFragmentPairCounterShape20Mb(t *testing.T) {
	const fragments = (20_000_000 + 1020 - 1) / 1020
	entries := []rtree.BatchEntry{{Key: 1, Val: uint32(fragments-1) << 1}}
	rows, cols, blockRows, ok := fragmentPairCounterShape(entries, entries)
	if !ok || rows != fragments || cols != fragments {
		t.Fatal("20 Mb genome pair did not select dense counting")
	}
	if blocks := (rows + blockRows - 1) / blockRows; blocks > maxFragmentPairCounterBlocks {
		t.Fatal("20 Mb genome pair exceeds the rescan limit")
	}
	cells := blockRows * cols
	if cells*2+(cells+63)/64*8 > maxFragmentPairCounterBytes {
		t.Fatal("20 Mb genome pair exceeds the counter memory budget")
	}
}

func TestFragmentPairsConcurrentImmutableEntries(t *testing.T) {
	a := []rtree.BatchEntry{{Key: 1, Val: 0}, {Key: 1, Val: 2}, {Key: 3, Val: 4}}
	b := []rtree.BatchEntry{{Key: 1, Val: 1}, {Key: 1, Val: 1}, {Key: 1, Val: 3}, {Key: 3, Val: 5}}
	originalA, originalB := slices.Clone(a), slices.Clone(b)
	options := &FragmentComparatorOptions{K: 11, MinSharedKmers: 1, TopNFragments: 1}
	want := legacyFragmentPairs(a, b, options)
	var wg sync.WaitGroup
	for range 8 {
		wg.Add(1)
		go func() {
			defer wg.Done()
			cpr := &FragmentComparator{options: options}
			for range 30 {
				pairs := cpr.scanPairsMerged(a, b)
				slices.Sort(*pairs)
				if !slices.Equal(*pairs, want) {
					t.Error("concurrent pair counting changed the original selection")
				}
				RecycleFragmentCompareResult(pairs)
			}
		}()
	}
	wg.Wait()
	if !slices.Equal(a, originalA) || !slices.Equal(b, originalB) {
		t.Fatal("concurrent pair scanning modified shared entries")
	}
}

// BenchmarkFragmentPairs isolates counting and Top-N selection with real
// contig-aware, sampled OrthoANI entries; FASTA reading and indexing are excluded.
func BenchmarkFragmentPairs(b *testing.B) {
	options := &FragmentComparatorOptions{K: 11, Scaled: 4, TopNFragments: 5,
		MinSharedKmers: max(3, MinSharedKmersThresholdExact(1020, 11, 4, 0.80, 0.99))}
	cpr := NewFragmentComparator(options, nil)
	load := func(file string) []rtree.BatchEntry {
		q, err := ReadGenomeFromFile("../t.demo/"+file+".fa", nil, false)
		if err != nil {
			b.Fatal(err)
		}
		defer RecycleGQuery(q)
		fragments, _ := seqs2fragments(&q.seqs, 1020, 1020)
		defer recycleFragments(fragments)
		entries, err := cpr.IndexA(fragments)
		if err != nil {
			b.Fatal(err)
		}
		b.Cleanup(func() { RecycleResultOfIndexA(entries) })
		return entries
	}
	a, other := load("GCF_002949675.1"), load("GCF_002950215.1")
	for _, input := range []struct {
		name string
		b    []rtree.BatchEntry
	}{{"self", a}, {"related", other}} {
		b.Run(input.name, func(b *testing.B) {
			pairs := cpr.scanPairsMerged(a, input.b)
			b.ReportMetric(float64(len(*pairs)), "pairs/op")
			RecycleFragmentCompareResult(pairs)
			b.ReportAllocs()
			for b.Loop() {
				pairs := cpr.scanPairsMerged(a, input.b)
				RecycleFragmentCompareResult(pairs)
			}
		})
	}
}
