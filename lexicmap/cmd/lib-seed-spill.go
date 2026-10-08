package cmd

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"slices"
	"unsafe"
)

/*
#Seed spilling: background, storage, and readback

Queries with very many seed matches in large indexes can run out of memory
during seed collection. A single decoded reference seed location can produce
multiple anchors when the matching seed occurs at multiple query positions.
Keeping all decoded KV data and expanded anchors in memory can therefore be
expensive, even when only a small number of genomes will pass global Top-N.
Top-N cannot be decided yet: each genome's chaining score needs its complete
set of anchors from all seed partitions.

--max-seed-memory enables streaming seed data and spilling collected anchors
to temporary files, using the existing index without rebuilding it. It is
disabled by default (0). Reducing -J/--max-query-conc also reduces the memory
used by concurrent queries. Spilling trades sorting and disk I/O for smaller
anchor collection buffers; a very small budget causes frequent writes/merges.

## Collection and the first spill

The budget is divided into fixed shares for concurrent query slots. Each query
owns one seedSpillStore, and idle slots do not lend their shares to active ones.
Decoded KV data arrive in small borrowed batches. Searchers expand those values
into query/target anchors, apply genome/TaxId filters, and submit accepted anchors
to the store through bounded producer buffers.

	Existing seed index -> batched KV decoding -> expanded anchors
	                                                  |
	                                                  v
	                                          In-memory collection
	                                                  |
	                       buffer growth fits? -------+-------+
	                             yes: continue                | no
	                                                          v
	                                               Stable sort by genome ID
	                                                          |
	                                                          v
	                                                 Write a sorted run
	                                                          |
	                                                          v
	                                              Reset length; reuse buffer

Before growing the collection array, the store accounts for producer buffers,
old/new collection arrays during copying, and radix-sort workspace. If growth
would exceed the query share, flush writes the ENTIRE current anchor buffer,
then resets its length to zero while retaining capacity. It does not keep the
first budget-sized portion permanently in RAM and write only later anchors.
Sorting/growth reserves mean the spill can start before anchor payload alone
reaches the configured share. Stable sorting groups genomes while preserving
the accepted arrival order of anchors within each genome; it does not deduplicate
anchors, calculate chaining scores, or apply Top-N.

## Temporary files and the 20-byte record format

Each query gets a separate directory, created lazily on its first spill under
the system temporary directory (TMPDIR). Choose a disk-backed directory with
enough free space: tmpfs stores the files in RAM. A run contains many genomes,
and one genome can occur in many runs. Run files have no header or per-genome
delimiter; they contain consecutive, explicitly encoded 20-byte records:

	Byte range  [0:8]      [8:12]       [12:16]       [16]    [17]   [18:20]
	Field       genome ID  query begin  target begin  length  flags  reserved
	Size        8 B        4 B          4 B           1 B     1 B    2 B

	8 + 4 + 4 + 1 + 1 + 2 = 20 bytes per anchor

Multi-byte fields are little-endian. The genome ID packs the genome-batch and
genome identifiers. Query/target coordinates retain their int32 bit patterns.
Flags pack query/target reverse-complement booleans into bits 0/1; the reserved
bytes are zero. These are 18 useful bytes plus two explicit reserved bytes,
independent of Go struct padding (seedAnchor occupies 24 bytes on amd64).
Records store expanded anchors, not raw uint64 KV data, k-mers, sequences, or
alignment results. Query identity comes from the owning store, so it is not
repeated in each record. writeSeedSpillAnchor/readSeedSpillAnchor encode/decode
this layout.

For example, letters below identify genomes, and each letter is one anchor:

	First buffer:  B A C A B -> sort -> run-000000000000: A A B B C
	Second buffer: D B A D B -> sort -> run-000000000001: A B B D D
	Final tail:    A C       -> sort -> run-000000000002: A C

During collection, each group of 16 chronological runs at the same merge level
is merged into one run at the next level. Inputs are deleted only after the
new output succeeds. This keeps at most 15 runs per level and avoids retaining
metadata for every flush. Each merge opens at most 16 readers and one writer.

## Readback, chaining, and global Top-N

consume runs only after every seed searcher has finished, since any partition
can still contribute anchors to a genome. If no spill occurred, it sorts and
iterates the in-memory anchors without creating files. If any spill occurred,
it also writes the final in-memory tail, releases the collection/sort arrays,
and compacts to at most 16 remaining runs before streaming their final merge.

The temporary files are not another searchable seed index. Matching has already
finished against the original index. Readback sequentially merges sorted runs,
holding one head record per input in a value-slice min-heap, plus I/O buffers.
Equal genome IDs become contiguous, even when their anchors span many files:

	run 0: A A B B C
	run 1: A B B D D  -- merge --> A A A A | B B B B | C C | D D
	run 2: A C                      genome A  genome B   C     D

	Complete genome anchors -> deduplicate and chain
	                        -> retain score and seedChainRegion boundaries
	All genome scores       -> global Top-N with cutoff ties
	                        -> read target sequences from the original index
	                        -> existing alignment, chunk merging, and output

The chaining caller reconstructs complete genomes into bounded batches and
uses the ordinary parallel chaining implementation. Tiny shares and single-
thread searches process one complete genome at a time. A genome is never split
into independently scored partial groups, and no per-run Top-N is applied.
This preserves global selection by chaining score, including cutoff ties.
After chaining, the store's deferred close removes its directory on normal
returns and returned errors, and releases its query slot.

## Budget scope

The budget covers producer/collection buffers, radix workspace, and normal
chaining batches. It excludes a single oversized genome or the serial fallback's
complete anchor array, worker chaining scratch, decoding/I/O buffers, candidate
metadata, loaded index data, and alignment/output memory. Unreachable arrays
may also await GC. This limits accounted anchor buffers, not total heap/RSS;
a genome with unusually many anchors can still exceed available RAM.
*/

// seedSpillRecordBytes is the fixed on-disk size of one expanded seedAnchor:
// 8 (genome ID) + 4 (query begin) + 4 (target begin) + 1 (match length)
// + 1 (two strand flags packed into one byte) + 2 (reserved zeros) = 20 bytes.
// The 18 useful bytes and two reserved bytes are an explicit file format, not
// Go struct padding. The in-memory seedAnchor is 24 bytes on amd64; its fields
// are encoded individually, so platform-dependent padding is never written.
const seedSpillRecordBytes = 20

// seedSpillMergeFanIn limits one merge to 16 input files and one output file.
// Equal-level compaction avoids accumulating unbounded file metadata when a
// small collection budget causes frequent spills.
const seedSpillMergeFanIn = 16

// Small groups avoid a second array and use the existing comparison sort.
const seedSpillRadixMin = 1024

// seedAnchorBytes includes struct padding (24 bytes on amd64). Collection budget
// calculations use this size, not the smaller serialized record size.
const seedAnchorBytes = int64(unsafe.Sizeof(seedAnchor{}))

// seedMemoryBudget divides collection space between query slots. Acquiring a
// slot before using searchers prevents queries waiting for space from holding
// decoding resources. Each query receives the same fixed share; idle slots do
// not lend their space. A store holds its slot through chaining until close.
//
// This accounts for producer batches, collection arrays (including old/new arrays
// during growth), radix workspace, and normal replay batches. It excludes a
// single oversized genome or the serial fallback's complete anchor array, worker
// chaining scratch, decoding/I/O buffers, and candidate metadata. Unreachable
// arrays may await GC, so this is not a total heap/RSS limit. The channel is shared
// across queries; the other fields are fixed before use.
type seedMemoryBudget struct {
	// slots is a counting semaphore: send to acquire a query share, receive to
	// release it. Capacity is the effective number of simultaneously active stores.
	slots chan struct{}
	// anchorsPerQuery limits collection, transfer, and radix arrays in seedAnchor
	// elements. Normal replay batches also use the unused bytes of this share.
	anchorsPerQuery int
	// tempDir is the parent of per-query directories. Empty uses os.MkdirTemp's
	// system default, including TMPDIR. Tests override it before creating stores.
	tempDir string
}

// newSeedMemoryBudget returns nil for limit <= 0 (spilling disabled). Callers
// must validate positive limits to be >= seedAnchorBytes before calling it.
// queries <= 0 means one slot; a small limit reduces the slot count so every
// share fits at least one anchor. Integer-division remainders are left unused,
// and the per-query element count is capped to fit an int.
func newSeedMemoryBudget(limit int64, queries int) *seedMemoryBudget {
	if limit <= 0 {
		return nil
	}
	queries = int(min(int64(max(queries, 1)), max(int64(1), limit/seedAnchorBytes)))
	return &seedMemoryBudget{slots: make(chan struct{}, queries), anchorsPerQuery: int(min(limit/int64(queries)/seedAnchorBytes, int64(^uint(0)>>1)))}
}

// seedSpillRun describes a temporary file sorted by packed genome identifier.
// Anchors within each genome retain their accepted arrival order. A merged run
// represents a contiguous range of earlier runs and inherits the oldest order
// in that range, preserving chronology across subsequent merge passes.
type seedSpillRun struct {
	// path belongs to the owning store's directory. combine deletes replaced
	// inputs; close deletes all files still present, including partial files.
	path string
	// level is the merge generation: zero for a collection flush, otherwise
	// one above the highest input level. Collection compacts equal-level groups:
	// 16 level-zero runs become one level-one run; 16 of those become level two.
	level int
	// order breaks ties between equal genome identifiers. It identifies the
	// oldest input run, not the file number of a newly created merged output.
	order uint64
}

// seedSpillStore owns one query's anchors and temporary files. Its lifecycle is
// newSeedSpillStore -> add -> consume -> close. consume is called once after
// all seed searchers finish. Defer close immediately after construction so
// returned errors also remove files and release the budget slot.
//
// Collection stays in anchors until growth would exceed the query share. After
// the first spill, sorted runs hold earlier anchors and anchors holds the newest
// tail. consume sorts the in-memory array or merges files, emitting one genome
// group at a time. It does not deduplicate, chain, rank, or discard anchors.
//
// Methods are not safe for concurrent use. collectAndChainSpilledSeeds protects
// add with a per-query mutex and calls consume only after joining all producers.
type seedSpillStore struct {
	// budget provides the fixed collection share and semaphore. This store
	// acquires exactly one slot at construction and releases it in close.
	budget *seedMemoryBudget
	// anchors is the flat, pointer-free array charged to the budget. flush
	// resets its length and reuses capacity; spilled consume drops the array.
	anchors []seedAnchor
	// sortScratch is the second, pointer-free array for stable radix grouping.
	// Both capacities count against the query share. It is reused between flushes
	// and released before replay, so replay can use the resulting budget headroom.
	sortScratch []seedAnchor
	// transferReserve charges the fixed producer batches before producers start.
	// It is zero for direct store users and released after all searchers join.
	transferReserve int
	// replayBufferBytes charges the one flat buffer for normal chaining batches.
	// A single genome larger than this buffer still follows the documented
	// unbudgeted replay fallback; the normal buffer remains charged during it.
	replayBufferBytes int64
	// onFirstSpill optionally reports the transition before the first nonempty
	// flush starts sorting/writing. flush clears it before calling it, so even a
	// failed write cannot report the start twice. The callback runs synchronously
	// under the collector's store lock; direct store users need no logger.
	onFirstSpill func()
	// runs lists live sorted files in chronological order. After collection
	// compaction, at most 15 runs remain at any one merge level.
	runs []seedSpillRun
	// dir is created lazily at the first spill and remains empty for queries
	// that fit in their collection array, avoiding unnecessary temporary files.
	dir string
	// nextRun assigns unique file names, including merged outputs. A newly
	// flushed run also uses this number as its chronological order.
	nextRun uint64
	// count counts every successfully accepted anchor, including duplicates.
	// consume compares its emitted count to detect whole-record truncation.
	count uint64
	// peakBufferBytes is debug telemetry for accounted anchor-buffer allocations,
	// including simultaneous old/replacement arrays, rather than observed RSS.
	peakBufferBytes int64
	// spilled stays true after the first successful flush and selects the file
	// path in consume even if the current in-memory tail is empty.
	spilled bool
	// closed makes cleanup idempotent, even if directory removal reports an error.
	closed bool
	// writer reuses one 64 KiB output buffer across all run files. writeRun
	// detaches its file reference on return. This buffer is outside the budget.
	writer *bufio.Writer
	// readers reuses up to 16 input buffers of 32 KiB each. mergeRuns detaches
	// file references on return. Buffers are allocated lazily, outside the budget.
	readers []*bufio.Reader
}

// newSeedSpillStore blocks until a slot is available in the non-nil budget.
// It allocates no anchor arrays or files yet. Call it before taking searcher
// tokens, and defer close immediately so failure cannot leak a query slot.
func newSeedSpillStore(b *seedMemoryBudget) *seedSpillStore {
	b.slots <- struct{}{}
	return &seedSpillStore{budget: b}
}

// add accepts one expanded, already-filtered anchor. Before growing the array,
// it accounts for both old/new capacities; if their sum would exceed the share,
// it flushes and reuses the existing array. Sorting/writing/compaction happen
// synchronously, making the producer wait instead of queuing extra anchors.
// For example, replacing a 512-element array with a 1024-element array requires
// room for 1536 elements during copying, not just the final 1024 elements.
// On an I/O error the new anchor is not accepted; the caller must abort and close.
func (s *seedSpillStore) add(anchor seedAnchor) error {
	if len(s.anchors) == cap(s.anchors) {
		if err := s.makeRoom(); err != nil {
			return err
		}
	}
	s.anchors = append(s.anchors, anchor)
	s.count++
	return nil
}

// addBatch copies a borrowed producer batch, flushing only at array boundaries.
// No producer-owned slice is retained. Copying spans amortizes the per-anchor
// append/capacity checks, while each successfully copied span updates count.
func (s *seedSpillStore) addBatch(anchors []seedAnchor) error {
	for len(anchors) > 0 {
		if len(s.anchors) == cap(s.anchors) {
			if err := s.makeRoom(); err != nil {
				return err
			}
		}
		n := min(len(anchors), cap(s.anchors)-len(s.anchors))
		s.anchors = append(s.anchors, anchors[:n]...)
		s.count += uint64(n)
		anchors = anchors[n:]
	}
	return nil
}

// makeRoom reserves radix workspace before growing a large collection array.
// Tiny shares keep the comparison-sort path, which needs no second array.
// Growth accounts for the old array too, even though it becomes unreachable.
func (s *seedSpillStore) makeRoom() error {
	available := s.budget.anchorsPerQuery - s.transferReserve
	maximum := available
	if available >= seedSpillRadixMin {
		// Shares smaller than two radix arrays stop below the radix threshold.
		maximum = max(seedSpillRadixMin-1, available/2)
	}
	capacity := cap(s.anchors)
	next := min(maximum, max(512, capacity*2))
	if capacity > 0 && (next <= capacity || capacity+next+cap(s.sortScratch) > available) {
		return s.flush()
	}
	buf := make([]seedAnchor, len(s.anchors), next)
	copy(buf, s.anchors)
	s.anchors = buf
	s.recordBufferPeak(int64(capacity) * seedAnchorBytes)
	return nil
}

// recordBufferPeak includes all live, explicitly budgeted buffers. extra is the
// old collection array during growth. Unreachable arrays awaiting GC, fixed I/O
// buffers, candidate pages, and chaining scratch retain their existing exclusions.
func (s *seedSpillStore) recordBufferPeak(extra int64) {
	n := int64(s.transferReserve+cap(s.anchors)+cap(s.sortScratch))*seedAnchorBytes + s.replayBufferBytes + extra
	s.peakBufferBytes = max(s.peakBufferBytes, n)
}

// sortAnchors stably groups by the entire uint64 genome identifier. Byte-wise
// counting passes avoid comparison sorting billions of records; bytes shared by
// every identifier are skipped. Equal identifiers preserve accepted arrival order,
// exactly as SortStableFunc, including across run compaction. The two array
// capacities have already been reserved by makeRoom and are charged together.
func (s *seedSpillStore) sortAnchors() {
	if len(s.anchors) < seedSpillRadixMin {
		slices.SortStableFunc(s.anchors, compareSeedAnchorGenome)
		return
	}
	first := s.anchors[0].batchGenomeIndex
	var varying uint64
	for _, a := range s.anchors {
		varying |= a.batchGenomeIndex ^ first
	}
	if varying == 0 {
		return
	}
	if cap(s.sortScratch) < cap(s.anchors) {
		// Equal capacities keep odd-pass swaps from shrinking the reusable
		// collection array when the final run is only partially filled.
		s.sortScratch = make([]seedAnchor, len(s.anchors), cap(s.anchors))
	}
	s.recordBufferPeak(0)
	src, dst := s.anchors, s.sortScratch[:len(s.anchors)]
	for shift := uint(0); shift < 64; shift += 8 {
		if byte(varying>>shift) == 0 {
			continue
		}
		var counts [256]int
		for _, a := range src {
			counts[byte(a.batchGenomeIndex>>shift)]++
		}
		// 257 prefix boundaries cover all 256 byte values, including the end
		// of bucket 255. Entry 256 is the total count, rather than a bucket.
		var offsets [257]int
		for i, n := range counts {
			offsets[i+1] = offsets[i] + n
		}
		for _, a := range src {
			bucket := byte(a.batchGenomeIndex >> shift)
			dst[offsets[bucket]] = a
			offsets[bucket]++
		}
		src, dst = dst, src
	}
	// Keep the sorted array directly; an odd pass count needs no copy-back.
	s.anchors, s.sortScratch = src, dst[:0]
}

// compareSeedAnchorGenome compares only packed genome-batch/genome identifiers.
// Equal identifiers deliberately compare equal: stable sorting preserves arrival
// order and leaves coordinate sorting/deduplication to ClearSubstrPairs.
func compareSeedAnchorGenome(a, b seedAnchor) int {
	if a.batchGenomeIndex < b.batchGenomeIndex {
		return -1
	}
	if a.batchGenomeIndex > b.batchGenomeIndex {
		return 1
	}
	return 0
}

// writeSeedSpillAnchor appends one expanded anchor as a fixed 20-byte record.
// Multi-byte integers use little-endian encoding. Byte offsets are:
//
//	[0:8]   packed genome-batch/genome identifier (uint64)
//	[8:12]  query begin (int32 bit pattern)
//	[12:16] target begin (int32 bit pattern)
//	[16]    match length (uint8)
//	[17]    flags: bit 0 = query reverse complement, bit 1 = target reverse complement
//	[18:20] reserved zeros, included in the 20-byte record size
//
// This is not the raw uint64 KV data: matching has already expanded it into
// query/target coordinates. No k-mer or sequence is stored. No query ID is
// needed because each seedSpillStore and its run files belong to one query.
// Converting coordinates to uint32 preserves their signed int32 bit patterns.
// AvailableBuffer avoids allocating a separate record slice. The caller owns
// the writer and must flush it before closing the file.
func writeSeedSpillAnchor(w *bufio.Writer, a seedAnchor) error {
	if w.Available() < seedSpillRecordBytes {
		if err := w.Flush(); err != nil {
			return err
		}
	}
	record := w.AvailableBuffer()
	// Bytes 0..7: identify the genome to group anchors across all run files.
	record = binary.LittleEndian.AppendUint64(record, a.batchGenomeIndex)
	// Bytes 8..11 and 12..15: the two coordinates needed by chaining.
	record = binary.LittleEndian.AppendUint32(record, uint32(a.qBegin))
	record = binary.LittleEndian.AppendUint32(record, uint32(a.tBegin))
	// Pack the two bool fields into bits 0 and 1 of byte 17.
	var flags byte
	if a.qrc {
		flags |= 1
	}
	if a.trc {
		flags |= 2
	}
	// Bytes 16..19: length, packed flags, and two explicit reserved bytes.
	record = append(record, a.length, flags, 0, 0)
	_, err := w.Write(record)
	return err
}

// readSeedSpillAnchor consumes one record in writeSeedSpillAnchor's format.
// io.EOF means a record boundary; a partial trailing record returns
// io.ErrUnexpectedEOF. The returned value keeps no reference to the borrowed
// Peek buffer, which the reader can reuse as soon as the record is discarded.
// Decode the offsets documented above; bytes 18..19 are reserved and ignored.
// uint32-to-int32 conversion restores the original signed coordinate bit pattern.
func readSeedSpillAnchor(r *bufio.Reader) (seedAnchor, error) {
	record, err := r.Peek(seedSpillRecordBytes)
	if err != nil {
		if err == io.EOF && len(record) > 0 {
			err = io.ErrUnexpectedEOF
		}
		return seedAnchor{}, err
	}
	a := seedAnchor{
		batchGenomeIndex: binary.LittleEndian.Uint64(record[:8]),
		qBegin:           int32(binary.LittleEndian.Uint32(record[8:12])),
		tBegin:           int32(binary.LittleEndian.Uint32(record[12:16])),
		length:           record[16],
		qrc:              record[17]&1 != 0,
		trc:              record[17]&2 != 0,
	}
	_, err = r.Discard(seedSpillRecordBytes)
	return a, err
}

// writeRun creates a file with the supplied level/order. write must produce
// genome-sorted records and must not retain the borrowed writer. A run is
// returned only after writing, flushing, and closing succeed. On error any
// partial file stays inside s.dir for the deferred close to remove.
// The writer buffer is reused across calls and detached from the file on return.
// A file is just concatenated 20-byte records, with no header or per-genome
// delimiter. It can contain many genomes, and one genome can span many runs.
// The genome ID in each record supplies the grouping boundary during merging.
func (s *seedSpillStore) writeRun(level int, order uint64, write func(*bufio.Writer) error) (seedSpillRun, error) {
	if s.dir == "" {
		var err error
		s.dir, err = os.MkdirTemp(s.budget.tempDir, "lexicmap-seeds-*")
		if err != nil {
			return seedSpillRun{}, fmt.Errorf("create seed temporary directory: %w", err)
		}
	}
	path := filepath.Join(s.dir, fmt.Sprintf("run-%012d", s.nextRun))
	s.nextRun++
	f, err := os.Create(path)
	if err != nil {
		return seedSpillRun{}, err
	}
	if s.writer == nil {
		s.writer = bufio.NewWriterSize(nil, 64<<10)
	}
	w := s.writer
	w.Reset(f)
	defer w.Reset(nil)
	err = write(w)
	if err == nil {
		err = w.Flush()
	}
	closeErr := f.Close()
	if err == nil {
		err = closeErr
	}
	if err != nil {
		return seedSpillRun{}, fmt.Errorf("write seed spill: %w", err)
	}
	return seedSpillRun{path: path, level: level, order: order}, nil
}

// flush stably groups the current anchors by genome, writes a level-zero run,
// and resets the collection slice for reuse. An empty slice is a no-op.
// It then repeatedly replaces the newest 16 equal-level runs with one higher-
// level run. This bounds metadata to O(log(number of flushed runs)) and limits
// each merge's open input files without reading complete runs into memory.
func (s *seedSpillStore) flush() error {
	if len(s.anchors) == 0 {
		return nil
	}
	if notify := s.onFirstSpill; notify != nil {
		s.onFirstSpill = nil
		notify()
	}
	// Stable grouping and chronological run tie breaks preserve the incoming
	// anchor order within each genome before the existing deduplication step.
	s.sortAnchors()
	run, err := s.writeRun(0, s.nextRun, func(w *bufio.Writer) error {
		for _, a := range s.anchors {
			if err := writeSeedSpillAnchor(w, a); err != nil {
				return err
			}
		}
		return nil
	})
	if err != nil {
		return err
	}
	s.spilled = true
	s.runs = append(s.runs, run)
	s.anchors = s.anchors[:0]
	// Compact equal-level chronological runs, bounding open files and run
	// metadata even when a very small budget creates millions of runs.
	// A merge can complete another group at the next level, so repeat until
	// the newest group has fewer than 16 runs or contains mixed levels.
	for len(s.runs) >= seedSpillMergeFanIn {
		start := len(s.runs) - seedSpillMergeFanIn
		group := s.runs[start:]
		same := true
		for _, r := range group {
			if r.level != run.level {
				same = false
				break
			}
		}
		if !same {
			break
		}
		run, err = s.combine(group)
		if err != nil {
			return err
		}
		clear(s.runs[start:])
		s.runs = append(s.runs[:start], run)
	}
	return nil
}

// combine merges a nonempty chronological group of at most 16 runs into a
// new file. The output inherits the oldest input's order and advances the level.
// Inputs are removed only after the output is successfully closed. Callers
// replace their metadata only on success; close removes all remaining files
// after either a merge failure or a failure to remove an input.
func (s *seedSpillStore) combine(runs []seedSpillRun) (seedSpillRun, error) {
	level := 0
	for _, r := range runs {
		level = max(level, r.level+1)
	}
	run, err := s.writeRun(level, runs[0].order, func(w *bufio.Writer) error {
		return s.mergeRuns(runs, func(a seedAnchor) error { return writeSeedSpillAnchor(w, a) })
	})
	if err != nil {
		return seedSpillRun{}, err
	}
	for _, r := range runs {
		if err := os.Remove(r.path); err != nil {
			return seedSpillRun{}, err
		}
	}
	return run, nil
}

// consume delivers all accepted anchors grouped by genome, preserving arrival
// order within each genome. Call it once after the last add. emit receives an
// independent value synchronously; its first error stops delivery unchanged.
//
// Without a spill, it sorts/iterates the collection array. With a spill, it writes
// the tail, drops that array, compacts to at most 16 runs, and streams their merge.
// The chaining caller buffers complete genomes in bounded batches (or uses the
// one-genome fallback); this method never reconstructs a whole query in memory.
// Successful delivery must match count, detecting missing whole records that a
// normal EOF cannot reveal. Files and the slot stay owned by the store until close.
func (s *seedSpillStore) consume(emit func(seedAnchor) error) (err error) {
	var consumed uint64
	originalEmit := emit
	emit = func(a seedAnchor) error { consumed++; return originalEmit(a) }
	// Truncation at a multiple of 20 bytes produces a clean record-boundary EOF.
	// Validate totals only after success, preserving genuine I/O/callback errors.
	defer func() {
		if err == nil && consumed != s.count {
			err = fmt.Errorf("seed spill record count %d != %d: %w", consumed, s.count, io.ErrUnexpectedEOF)
		}
	}()
	if !s.spilled {
		s.sortAnchors()
		s.sortScratch = nil
		for _, a := range s.anchors {
			if err := emit(a); err != nil {
				return err
			}
		}
		return nil
	}
	if err := s.flush(); err != nil {
		return err
	}
	// The anchor array is no longer needed while reading/merging spilled runs.
	s.anchors = nil
	s.sortScratch = nil
	// Collection can leave mixed levels. Combine chronological prefixes until
	// the final merge can open every remaining input within the 16-reader limit.
	for len(s.runs) > seedSpillMergeFanIn {
		run, err := s.combine(s.runs[:seedSpillMergeFanIn])
		if err != nil {
			return err
		}
		tail := slices.Clone(s.runs[seedSpillMergeFanIn:])
		clear(s.runs)
		s.runs = append(s.runs[:0], run)
		s.runs = append(s.runs, tail...)
	}
	return s.mergeRuns(s.runs, emit)
}

// close drops buffer/metadata references, recursively removes the temporary
// directory, and releases the slot. It does not flush or consume pending anchors:
// this is cleanup for both success and failure. The slot is released even if
// removal fails, and a second close is a no-op.
func (s *seedSpillStore) close() error {
	if s.closed {
		return nil
	}
	s.closed = true
	s.anchors = nil
	s.sortScratch = nil
	s.writer = nil
	s.onFirstSpill = nil
	clear(s.readers)
	s.readers = nil
	clear(s.runs)
	s.runs = nil
	var err error
	if s.dir != "" {
		err = os.RemoveAll(s.dir)
	}
	<-s.budget.slots
	return err
}

// seedRunReader is one min-heap entry in mergeRuns: a run's head record and the
// input needed to advance it. There is one entry per input run (at most 16),
// rather than a pointer/object allocation for every anchor.
type seedRunReader struct {
	// reader is borrowed from the store's reusable buffers. mergeRuns owns the
	// associated file and detaches/closes it before returning.
	reader *bufio.Reader
	// anchor is the next record to emit, replaced by this run's next record
	// after emission. At EOF the entire heap entry is removed.
	anchor seedAnchor
	// order comes from seedSpillRun and preserves chronology when several
	// heap heads belong to the same genome.
	order uint64
}

// seedRunLess orders first by genome and then by run chronology. Runs are already
// stable internally, so this preserves same-genome order across file boundaries
// and across runs produced by earlier compaction passes.
func seedRunLess(a, b seedRunReader) bool {
	if a.anchor.batchGenomeIndex != b.anchor.batchGenomeIndex {
		return a.anchor.batchGenomeIndex < b.anchor.batchGenomeIndex
	}
	return a.order < b.order
}

// mergeRuns performs a stable k-way merge of at most 16 genome-sorted runs.
// Only one head record per input is held in a value-slice min-heap; reader
// buffers are reused. emit synchronously writes a new run or passes anchors to
// chaining, and its errors propagate unchanged. All opened files are closed
// and readers detached on every exit. combine or close owns input-file deletion.
func (s *seedSpillStore) mergeRuns(runs []seedSpillRun, emit func(seedAnchor) error) error {
	readers := make([]seedRunReader, 0, len(runs))
	files := make([]*os.File, 0, len(runs))
	defer func() {
		for _, r := range s.readers {
			r.Reset(nil)
		}
		for _, f := range files {
			_ = f.Close()
		}
	}()
	// Read one head per input. Empty files have no heap entry; consume's final
	// total-count check detects missing whole records on the successful path.
	for i, run := range runs {
		f, err := os.Open(run.path)
		if err != nil {
			return err
		}
		files = append(files, f)
		if i == len(s.readers) {
			s.readers = append(s.readers, bufio.NewReaderSize(nil, 32<<10))
		}
		r := s.readers[i]
		r.Reset(f)
		a, err := readSeedSpillAnchor(r)
		if err == io.EOF {
			continue
		}
		if err != nil {
			return fmt.Errorf("read seed spill: %w", err)
		}
		readers = append(readers, seedRunReader{reader: r, anchor: a, order: run.order})
	}
	// A sorted array satisfies the min-heap property. Each iteration below only
	// changes the root, so maintaining the heap needs just one sift-down.
	slices.SortFunc(readers, func(a, b seedRunReader) int {
		if seedRunLess(a, b) {
			return -1
		}
		if seedRunLess(b, a) {
			return 1
		}
		return 0
	})
	// Emit the smallest head, advance that input, and restore the heap. At EOF
	// replace the root with the last entry instead of retaining an exhausted run.
	for len(readers) > 0 {
		if err := emit(readers[0].anchor); err != nil {
			return err
		}
		a, err := readSeedSpillAnchor(readers[0].reader)
		if err == io.EOF {
			readers[0] = readers[len(readers)-1]
			readers = readers[:len(readers)-1]
		} else if err != nil {
			return fmt.Errorf("read seed spill: %w", err)
		} else {
			readers[0].anchor = a
		}
		// Only the root changed, so restore the min-heap without per-item objects.
		for i := 0; i < len(readers); {
			child := 2*i + 1
			if child >= len(readers) {
				break
			}
			if child+1 < len(readers) && seedRunLess(readers[child+1], readers[child]) {
				child++
			}
			if !seedRunLess(readers[child], readers[i]) {
				break
			}
			readers[i], readers[child] = readers[child], readers[i]
			i = child
		}
	}
	return nil
}
