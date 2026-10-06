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

package cmd

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"os"
	"sync"
)

// alignmentPayloadBudget limits retained final-alignment output buffers across
// all concurrent queries. It does not represent total process memory.
type alignmentPayloadBudget struct {
	mu      sync.Mutex
	limit   int64
	used    int64
	tempDir string
}

func newAlignmentPayloadBudget(limit int64) *alignmentPayloadBudget {
	if limit <= 0 {
		return nil
	}
	return &alignmentPayloadBudget{limit: limit}
}

func (b *alignmentPayloadBudget) tryReserve(n int64) bool {
	if n <= 0 {
		return true
	}
	b.mu.Lock()
	defer b.mu.Unlock()
	if n > b.limit-b.used {
		return false
	}
	b.used += n
	return true
}

func (b *alignmentPayloadBudget) release(n int64) {
	if n <= 0 {
		return
	}
	b.mu.Lock()
	b.used -= n
	if b.used < 0 {
		panic("negative alignment payload memory usage")
	}
	b.mu.Unlock()
}

func (b *alignmentPayloadBudget) usedBytes() int64 {
	b.mu.Lock()
	defer b.mu.Unlock()
	return b.used
}

type alignmentPayloadRef struct {
	offset  int64
	lengths [4]uint32
}

// alignmentPayloadStore hides the in-memory/spilled payload choice from search
// ordering, chunk merging, coverage filtering, and output formatting.
type alignmentPayloadStore struct {
	budget   *alignmentPayloadBudget
	reserved int64

	file    *os.File
	writer  *bufio.Writer
	written int64
	refs    map[*AlignmentResult]alignmentPayloadRef

	spilled bool
	closed  bool
}

func newAlignmentPayloadStore(budget *alignmentPayloadBudget) *alignmentPayloadStore {
	if budget == nil {
		return nil
	}
	return &alignmentPayloadStore{budget: budget}
}

func alignmentPayloadBytes(a *AlignmentResult) int64 {
	return int64(cap(a.CIGAR)) + int64(cap(a.QSeq)) + int64(cap(a.TSeq)) + int64(cap(a.Alignment))
}

func searchResultPayloadBytes(r *SearchResult) (n int64) {
	for _, sd := range *r.SimilarityDetails {
		for i := range sd.Alignments {
			n += alignmentPayloadBytes(&sd.Alignments[i])
		}
	}
	return n
}

func searchResultsPayloadBytes(rs *[]*SearchResult) (n int64) {
	for _, r := range *rs {
		if r != nil {
			n += searchResultPayloadBytes(r)
		}
	}
	return n
}

// retain keeps a completed genome result in memory while the shared budget has
// room. On the first failed reservation, this query spills all earlier payloads
// and remains in spill mode for the rest of its lifetime.
func (s *alignmentPayloadStore) retain(r *SearchResult, previous *[]*SearchResult) error {
	if s.spilled {
		return s.spillSearchResult(r)
	}

	n := searchResultPayloadBytes(r)
	if s.budget.tryReserve(n) {
		s.reserved += n
		return nil
	}

	if err := s.openSpillFile(); err != nil {
		return err
	}
	for _, previousResult := range *previous {
		if err := s.spillSearchResult(previousResult); err != nil {
			return err
		}
	}
	if err := s.spillSearchResult(r); err != nil {
		return err
	}

	s.budget.release(s.reserved)
	s.reserved = 0
	s.spilled = true
	return nil
}

func (s *alignmentPayloadStore) openSpillFile() error {
	if s.file != nil {
		return nil
	}
	f, err := os.CreateTemp(s.budget.tempDir, "lexicmap-alignment-*.bin")
	if err != nil {
		return fmt.Errorf("create alignment-result spill file: %w", err)
	}
	s.file = f
	s.writer = bufio.NewWriterSize(f, 64<<10)
	s.refs = make(map[*AlignmentResult]alignmentPayloadRef)
	return nil
}

func (s *alignmentPayloadStore) spillSearchResult(r *SearchResult) error {
	for _, sd := range *r.SimilarityDetails {
		for i := range sd.Alignments {
			if err := s.spillAlignment(&sd.Alignments[i]); err != nil {
				return err
			}
		}
	}
	return nil
}

func (s *alignmentPayloadStore) spillAlignment(a *AlignmentResult) error {
	if _, ok := s.refs[a]; ok {
		return nil
	}
	if len(s.refs) == math.MaxUint32 {
		return fmt.Errorf("too many spilled alignment payloads")
	}

	payloads := [4][]byte{a.CIGAR, a.QSeq, a.TSeq, a.Alignment}
	ref := alignmentPayloadRef{offset: s.written}
	for i, payload := range payloads {
		if uint64(len(payload)) > math.MaxUint32 {
			return fmt.Errorf("alignment payload is too large: %d bytes", len(payload))
		}
		ref.lengths[i] = uint32(len(payload))
		n, err := s.writer.Write(payload)
		s.written += int64(n)
		if err != nil {
			return fmt.Errorf("write alignment-result spill file: %w", err)
		}
		if n != len(payload) {
			return fmt.Errorf("write alignment-result spill file: %w", io.ErrShortWrite)
		}
	}

	s.refs[a] = ref
	a.CIGAR = nil
	a.QSeq = nil
	a.TSeq = nil
	a.Alignment = nil
	return nil
}

// reconcile releases budget for results removed during chunk merging/filtering.
func (s *alignmentPayloadStore) reconcile(results *[]*SearchResult) {
	if s == nil {
		return
	}
	if s.spilled {
		live := 0
		for _, r := range *results {
			if r == nil {
				continue
			}
			for _, sd := range *r.SimilarityDetails {
				for i := range sd.Alignments {
					if _, ok := s.refs[&sd.Alignments[i]]; ok {
						live++
					}
				}
			}
		}
		if live == len(s.refs) {
			return
		}
		refs := make(map[*AlignmentResult]alignmentPayloadRef, live)
		for _, r := range *results {
			if r == nil {
				continue
			}
			for _, sd := range *r.SimilarityDetails {
				for i := range sd.Alignments {
					a := &sd.Alignments[i]
					if ref, ok := s.refs[a]; ok {
						refs[a] = ref
					}
				}
			}
		}
		clear(s.refs)
		s.refs = refs
		return
	}
	live := searchResultsPayloadBytes(results)
	if live < s.reserved {
		s.budget.release(s.reserved - live)
		s.reserved = live
	}
}

// finalize makes spilled payloads readable after alignment is complete.
func (s *alignmentPayloadStore) finalize() error {
	if s == nil || s.writer == nil {
		return nil
	}
	if err := s.writer.Flush(); err != nil {
		return fmt.Errorf("flush alignment-result spill file: %w", err)
	}
	s.writer = nil
	return nil
}

// payload returns the four output fields. scratch is reused by the caller and
// only grows to the largest single spilled HSP being printed.
func (s *alignmentPayloadStore) payload(a *AlignmentResult, scratch []byte) (
	cigar, qseq, tseq, alignment, nextScratch []byte, err error,
) {
	if s == nil || !s.spilled {
		return a.CIGAR, a.QSeq, a.TSeq, a.Alignment, scratch, nil
	}
	ref, ok := s.refs[a]
	if s.file == nil || !ok {
		return nil, nil, nil, nil, scratch, fmt.Errorf("alignment payload reference not found")
	}
	total := 0
	for _, n := range ref.lengths {
		total += int(n)
	}
	if cap(scratch) < total {
		scratch = make([]byte, total)
	} else {
		scratch = scratch[:total]
	}
	if total > 0 {
		n, readErr := s.file.ReadAt(scratch, ref.offset)
		if readErr != nil {
			return nil, nil, nil, nil, scratch, fmt.Errorf("read alignment-result spill file: %w", readErr)
		}
		if n != total {
			return nil, nil, nil, nil, scratch, fmt.Errorf("read alignment-result spill file: %w", io.ErrUnexpectedEOF)
		}
	}

	i := int(ref.lengths[0])
	j := i + int(ref.lengths[1])
	k := j + int(ref.lengths[2])
	return scratch[:i], scratch[i:j], scratch[j:k], scratch[k:], scratch, nil
}

// close releases the shared budget and removes any spill file. It is idempotent.
func (s *alignmentPayloadStore) close() error {
	if s == nil || s.closed {
		return nil
	}
	s.closed = true

	s.budget.release(s.reserved)
	s.reserved = 0

	var firstErr error
	if s.writer != nil {
		if err := s.writer.Flush(); err != nil {
			firstErr = err
		}
		s.writer = nil
	}
	var name string
	if s.file != nil {
		name = s.file.Name()
		if err := s.file.Close(); err != nil && firstErr == nil {
			firstErr = err
		}
		s.file = nil
	}
	if name != "" {
		if err := os.Remove(name); err != nil && !os.IsNotExist(err) && firstErr == nil {
			firstErr = err
		}
	}
	clear(s.refs)
	s.refs = nil
	return firstErr
}
