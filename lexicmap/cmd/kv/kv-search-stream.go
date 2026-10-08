package kv

// searchResultSink lets the collected and streaming APIs share exactly the same
// matching/decoding code. With consume == nil, every matched key is appended to
// a pooled result slice and its Values can grow to hold all matching locations.
// With a callback, a single fixed seedDataBuffer is reused and flush delivers
// small borrowed batches instead of retaining a query's complete KV results.
//
// One search/search2 call owns a sink. It is not shared across decoding goroutines.
// The decoder starts each matched key with start, fills metadata and Values,
// flushes full batches in streaming mode, and flushes the final tail on success.
type searchResultSink struct {
	// results is non-nil only in collected mode. The decoder returns this
	// pooled container to its caller on success and recycles it on error.
	results *[]SearchResult
	// consume is non-nil only in streaming mode. It runs synchronously and
	// must finish reading/copying Values before returning; a copy of SearchResult
	// alone does not own the Values backing array. Its error aborts decoding.
	consume func(SearchResult) error
	// stream is non-nil only in streaming mode and owns the reused result and
	// location array. It is separate from searchResultSink so collected calls
	// do not allocate the fixed streaming buffer.
	stream *seedDataBuffer
}

// seedDataBuffer holds encoded reference seed locations and query/match metadata.
// The matched reference k-mer is not retained. Values borrows this fixed array
// and must be consumed before the callback returns and the buffer is reused.
type seedDataBuffer struct {
	// current carries IQuery/IQuery2, matched length, suffix mode, and a Values
	// view into values. A flush resets only Values length so the next batch for
	// the same matched key keeps its metadata. start resets metadata for a new key.
	current SearchResult
	// values contains encoded reference genome identifiers, coordinates, and
	// strand/reverse flags. Disk decoding first reads one location to check its
	// reverse flag and, if accepted, appends it to Values. It then appends a full
	// 256-item block before flushing, so the first batch can contain 1+256=257
	// locations. With only 256 slots, that last append would grow the slice and
	// allocate a replacement array. The extra slot costs just 8 bytes and keeps
	// this buffer fixed while preserving the existing 256-item decoding blocks.
	// In-memory decoding has no separate pre-read and flushes at 256 locations.
	// Neither path retains the matched reference k-mer here.
	values [seedPosBatchSize + 1]uint64
}

// newSearchResultSink chooses the mode once per search call. Collected mode
// borrows an empty pooled result slice; streaming mode allocates one fixed
// buffer and initializes current.Values to an empty view of its array.
func newSearchResultSink(consume func(SearchResult) error) searchResultSink {
	s := searchResultSink{consume: consume}
	if consume == nil {
		s.results = poolSearchResults.Get().(*[]SearchResult)
		*s.results = (*s.results)[:0]
	} else {
		s.stream = &seedDataBuffer{}
		s.stream.current.Values = s.stream.values[:0]
	}
	return s
}

// start begins a new matched-key/query result and returns writable storage.
// In streaming mode it first delivers the preceding key's remaining Values,
// then clears all metadata while retaining the fixed array. The decoder must
// set IQuery/IQuery2, Len, and IsSuffix before adding locations. The returned
// address must not be retained across another start or the end of the search.
// A flush error prevents the new result from being started.
func (s *searchResultSink) start() (*SearchResult, error) {
	if s.consume == nil {
		return appendSearchResult(s.results), nil
	}
	if err := s.flush(); err != nil {
		return nil, err
	}
	s.stream.current = SearchResult{Values: s.stream.current.Values[:0]}
	return &s.stream.current, nil
}

// flush delivers a nonempty streaming batch and resets Values to length zero
// after the callback succeeds. Metadata and array capacity are retained so
// decoding the same key can continue without allocation. Collected mode and
// empty batches are no-ops. On callback error the decoder must abort immediately;
// the batch is left intact and the callback error is returned unchanged.
func (s *searchResultSink) flush() error {
	if s.consume == nil || len(s.stream.current.Values) == 0 {
		return nil
	}
	if err := s.consume(s.stream.current); err != nil {
		return err
	}
	s.stream.current.Values = s.stream.current.Values[:0]
	return nil
}
