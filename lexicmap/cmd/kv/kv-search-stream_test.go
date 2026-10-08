package kv

import (
	"errors"
	"fmt"
	"path/filepath"
	"reflect"
	"testing"
)

// seedDataKey identifies the query/match metadata used to compare collected
// and streamed output. Batch boundaries are intentionally absent: one collected
// result may be split across many callbacks without changing its seed data.
type seedDataKey struct {
	// mask is the global mask index (IQuery); query is the within-mask query
	// index (IQuery2) used by Search2.
	mask, query int
	// length is the common prefix/suffix length carried by every batch.
	length uint8
	// suffix distinguishes prefix matches from reversed-key suffix matches.
	suffix bool
}

// seedDataSnapshot joins batches with identical metadata, keeping value order
// within each group. Appending copies the borrowed Values, so the snapshot owns
// its data even after the streaming decoder reuses its buffer.
func seedDataSnapshot(rs []SearchResult) map[seedDataKey][]uint64 {
	m := make(map[seedDataKey][]uint64)
	for _, r := range rs {
		if len(r.Values) == 0 {
			continue
		}
		key := seedDataKey{r.IQuery, r.IQuery2, r.Len, r.IsSuffix}
		m[key] = append(m[key], r.Values...)
	}
	return m
}

func TestSeedDataStreamsMatchCollectedResults(t *testing.T) {
	for _, threeBytes := range []bool{false, true} {
		t.Run(fmt.Sprintf("threeBytes=%t", threeBytes), func(t *testing.T) {
			const k uint8 = 31
			const base uint64 = (0b0111 << (62 - 4)) | (0b0011 << 54) | 1024
			data := map[uint64]*[]uint64{}
			for i := uint64(0); i < 3; i++ {
				posts := make([]uint64, 1027+int(i))
				for j := range posts {
					posts[j] = uint64(j)<<2 | i&1
				}
				data[base+i] = &posts
			}
			file := filepath.Join(t.TempDir(), "seeds")
			if _, err := WriteKVData(k, 3, []*map[uint64]*[]uint64{&data}, file, 2, 2, 512, threeBytes); err != nil {
				t.Fatal(err)
			}
			disk, err := NewSearcher(file, 1)
			if err != nil {
				t.Fatal(err)
			}
			defer disk.Close()
			memory, err := NewInMemomrySearcher(file)
			if err != nil {
				t.Fatal(err)
			}
			defer memory.Close()
			for _, searcher := range []interface {
				Search([]uint64, uint8, bool, bool) (*[]SearchResult, error)
				Search2([][]uint64, uint8, bool, bool) (*[]SearchResult, error)
				SearchStream([]uint64, uint8, bool, bool, func(SearchResult) error) error
				Search2Stream([][]uint64, uint8, bool, bool, func(SearchResult) error) error
			}{disk, memory} {
				for _, check := range []bool{false, true} {
					for _, reverse := range []bool{false, true} {
						for _, multi := range []bool{false, true} {
							var want *[]SearchResult
							var err error
							if multi {
								want, err = searcher.Search2([][]uint64{{base + 1, base + 2}}, 28, check, reverse)
							} else {
								want, err = searcher.Search([]uint64{base}, 28, check, reverse)
							}
							if err != nil {
								t.Fatal(err)
							}
							var got []SearchResult
							consume := func(r SearchResult) error {
								if len(r.Values) > seedPosBatchSize+1 {
									t.Fatal("seed data batch grew without a bound")
								}
								r.Values = append([]uint64(nil), r.Values...)
								got = append(got, r)
								return nil
							}
							if multi {
								err = searcher.Search2Stream([][]uint64{{base + 1, base + 2}}, 28, check, reverse, consume)
							} else {
								err = searcher.SearchStream([]uint64{base}, 28, check, reverse, consume)
							}
							if err != nil {
								t.Fatal(err)
							}
							if !reflect.DeepEqual(seedDataSnapshot(*want), seedDataSnapshot(got)) {
								t.Fatalf("changed decoded seed data: %T check=%t reverse=%t multi=%t", searcher, check, reverse, multi)
							}
							RecycleSearchResults(want)
						}
					}
				}
				sentinel := errors.New("stop consumer")
				if err := searcher.SearchStream([]uint64{base}, 28, false, false, func(SearchResult) error { return sentinel }); !errors.Is(err, sentinel) {
					t.Fatalf("did not propagate consumer error: %v", err)
				}
				// The failed streaming call must return its decoding kit.
				rs, err := searcher.Search([]uint64{base}, 28, false, false)
				if err != nil {
					t.Fatal(err)
				}
				RecycleSearchResults(rs)
			}
		})
	}
}
