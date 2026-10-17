package cmd

import (
	"errors"
	"reflect"
	"slices"
	"testing"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
)

// borrowedSeedStreamer splits keys across callbacks and overwrites each borrowed
// buffer immediately afterward, as a real decoder is allowed to do.
type borrowedSeedStreamer struct {
	prefix, suffix []kv.SearchResult
	prefixErr      error
	suffixCalled   bool
}

func (s *borrowedSeedStreamer) send(results []kv.SearchResult, consume func(kv.SearchResult) error) error {
	var buffer [7]uint64
	for _, r := range results {
		for values := r.Values; len(values) > 0; {
			n := min(len(buffer), len(values))
			copy(buffer[:], values[:n])
			r.Values = buffer[:n]
			if err := consume(r); err != nil {
				return err
			}
			clear(buffer[:])
			values = values[n:]
		}
	}
	return nil
}

func (s *borrowedSeedStreamer) SearchStream(_ []uint64, _ uint8, _, _ bool, consume func(kv.SearchResult) error) error {
	if err := s.send(s.prefix, consume); err != nil {
		return err
	}
	return s.prefixErr
}

func (s *borrowedSeedStreamer) Search2Stream(_ [][]uint64, _ uint8, _, _ bool, consume func(kv.SearchResult) error) error {
	s.suffixCalled = true
	return s.send(s.suffix, consume)
}

func TestStreamSeedAnchorsMatchesCollected(t *testing.T) {
	for _, filtered := range []bool{false, true} {
		idx := &Index{k: 31, info: &IndexInfo{GenomeBatches: 3}}
		locses := [][]int{{10 << 1, 80<<1 | 1}, {50<<1 | 1, 90 << 1}}
		reverseLocses := [][]int{{1, 0}, {0, 1}}
		var prefix, suffix []kv.SearchResult
		for mask := range 2 {
			for pass := range 2 {
				r := kv.SearchResult{IQuery: mask, IQuery2: 1, Len: 21, IsSuffix: pass == 1}
				for i := range 600 { // several decoder and collector batches
					genome := uint64(i%3)<<BITS_GENOME_IDX | uint64(i%5)
					flags := uint64(i%2)<<BITS_REVERSE | uint64(pass)
					r.Values = append(r.Values, genome<<BITS_NONE_IDX|uint64(i*3)<<BITS_FLAGS|flags)
				}
				if pass == 0 {
					prefix = append(prefix, r)
				} else {
					suffix = append(suffix, r)
				}
			}
		}
		var filter *map[uint64]*[]uint64
		if filtered {
			ids := map[uint64]*[]uint64{0: nil, uint64(1)<<BITS_GENOME_IDX | 1: nil}
			filter = &ids
		}
		collected := append(slices.Clone(prefix), suffix...)
		// The serial collector recycles Values, so give it independent arrays.
		for i := range collected {
			collected[i].Values = slices.Clone(collected[i].Values)
		}
		ch := make(chan *[]kv.SearchResult, 1)
		ch <- &collected
		close(ch)
		var arena seedSearchResultArena
		wantResults := collectSeedAnchorsSerial(idx, ch, &locses, &reverseLocses, filter, &arena)
		counts := func(groups [][]*seedSearchResult) map[uint64]map[SubstrPair]int {
			out := make(map[uint64]map[SubstrPair]int)
			for _, group := range groups {
				for _, r := range group {
					m := make(map[SubstrPair]int)
					for _, sub := range r.Subs {
						m[sub]++
					}
					out[r.BatchGenomeIndex] = m
				}
			}
			return out
		}
		want := counts([][]*seedSearchResult{wantResults})
		arena.recycle()
		streamer := &borrowedSeedStreamer{prefix: prefix, suffix: suffix}
		collector := newSeedAnchorCollector(idx, 2)
		batches, values, err := streamSeedAnchors(streamer, collector, nil, nil, locses, reverseLocses, 15, filter)
		collector.finish()
		if err != nil || values != 2400 || batches != 344 || !streamer.suffixCalled {
			t.Fatalf("stream stats: batches=%d values=%d err=%v", batches, values, err)
		}
		if !reflect.DeepEqual(counts(collector.results), want) {
			t.Fatalf("filtered=%t: streamed anchors differ from collected anchors", filtered)
		}
		collector.clearResults()
		for i := range collector.arenas {
			collector.arenas[i].recycle()
		}
	}
}

func TestStreamSeedAnchorsDrainsOnReadError(t *testing.T) {
	sentinel := errors.New("truncated seed data")
	streamer := &borrowedSeedStreamer{prefixErr: sentinel, prefix: []kv.SearchResult{{
		IQuery: 0, Len: 21, Values: []uint64{uint64(7)<<BITS_NONE_IDX | uint64(10)<<BITS_FLAGS},
	}}}
	idx := &Index{k: 31, info: &IndexInfo{GenomeBatches: 2}}
	collector := newSeedAnchorCollector(idx, 2)
	_, _, err := streamSeedAnchors(streamer, collector, nil, nil, [][]int{{0}}, nil, 15, nil)
	collector.finish()
	if !errors.Is(err, sentinel) || streamer.suffixCalled || collector.anchorCount() != 1 {
		t.Fatalf("partial collection: err=%v suffix=%t anchors=%d", err, streamer.suffixCalled, collector.anchorCount())
	}
	collector.clearResults()
	for i := range collector.arenas {
		collector.arenas[i].recycle()
	}
}
