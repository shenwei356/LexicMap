package kv

import "io"

// index15Entry is one sparse primary checkpoint, including the mask's start entry.
type index15Entry struct {
	checkpoint uint64 // first target k-mer identifying the primary anchor
	offset     uint64 // tagged secondary offset or encoded seeds offset
}

// index15Builder consumes encoded pair boundaries from either a writer or a
// reindexing reader. It adds longer-prefix checkpoints to large blocks in large
// indexes to reduce scanning during prefix matching. Batch search gains are most
// noticeable when -n/--top-n-genomes limits downstream alignments; without this
// limit, alignment often dominates runtime. It retains one block and one mask's
// sparse primary entries; postings and completed secondary blocks are never
// retained.
type index15Builder struct {
	wi, w15                           io.Writer                           // primary and secondary outputs
	k, primaryPrefix, secondaryPrefix uint8                               // validated prefix lengths
	threshold                         uint64                              // minimum encoded block size
	getPrimary                        func(uint64) uint64                 // extracts the anchor within a mask
	block                             index15BuildBlock                   // current primary-prefix block
	hasBlock                          bool                                // whether the current mask has a block
	entries                           []index15Entry                      // reused across masks
	stats                             Index15Stats                        // per-chunk block and size totals
	scratch                           [index15BlockFixedBytes + 15*8]byte // bounded encoding scratch
}

// newIndex15Builder prepares a chunk whose headers have already been written.
// Its prefix lengths and threshold must have been validated by the caller.
func newIndex15Builder(wi, w15 io.Writer, k, maskPrefix, anchorPrefix uint8, chunkSize, threshold uint64) *index15Builder {
	return &index15Builder{
		wi: wi, w15: w15, k: k, primaryPrefix: maskPrefix + anchorPrefix,
		secondaryPrefix: maskPrefix + anchorPrefix + 2, threshold: threshold,
		getPrimary: AnchorExtracter(k, maskPrefix, anchorPrefix),
		stats:      Index15Stats{TotalBlocks: chunkSize * (uint64(1) << (anchorPrefix << 1)), Index15Bytes: index15HeaderSize},
	}
}

// beginMask reuses the primary directory without retaining the previous mask.
func (b *index15Builder) beginMask() {
	b.entries = b.entries[:0]
	b.hasBlock = false
}

// addPair records the encoded extent and pair-first checkpoint for one pair.
// A prefix starting at the second key shares this pair's seeds offset.
func (b *index15Builder) addPair(key1, key2 uint64, hasSecond bool, pairOffset, pairEnd uint64) error {
	if len(b.entries) == 0 && !b.hasBlock {
		offset, err := MakeSeedsOffset(pairOffset, false)
		if err != nil {
			return err
		}
		b.entries = append(b.entries, index15Entry{key1, offset})
	}
	if err := b.addKmer(key1, key1, pairOffset, pairEnd, pairOffset, false); err != nil {
		return err
	}
	if hasSecond {
		return b.addKmer(key2, key1, pairOffset, pairEnd, pairEnd, true)
	}
	return nil
}

// addKmer closes a preceding primary block at the same boundary as reindexing.
func (b *index15Builder) addKmer(key, pairKey, pairOffset, pairEnd, previousEnd uint64, isSecond bool) error {
	prefix := b.getPrimary(key)
	if !b.hasBlock || prefix != b.block.prefix {
		if err := b.finishBlock(previousEnd); err != nil {
			return err
		}
		b.block = index15BuildBlock{prefix: prefix, indexKmer: key,
			pairFirstCheckpoint: pairKey, baseOffset: pairOffset,
			endOffset: pairEnd, baseIsSecond: isSecond}
		b.hasBlock = true
	}
	b.block.add(key, pairKey, pairOffset, pairEnd, b.k, b.primaryPrefix, b.secondaryPrefix)
	return nil
}

// finishBlock writes a secondary block only when its encoded extent meets the threshold.
func (b *index15Builder) finishBlock(endOffset uint64) error {
	if !b.hasBlock {
		return nil
	}
	b.block.endOffset = endOffset
	b.stats.NonEmptyBlocks++
	offset, err := MakeSeedsOffset(b.block.baseOffset, b.block.baseIsSecond)
	if err != nil {
		return err
	}
	if endOffset-b.block.baseOffset >= b.threshold {
		offset, err = MakeIndex15Offset(b.stats.Index15Bytes)
		if err != nil {
			return err
		}
		width, size, err := writeIndex15Block(b.w15, &b.block, b.scratch[:])
		if err != nil {
			return err
		}
		b.stats.IndexedBlocks++
		b.stats.WidthCounts[width-1]++
		b.stats.Index15Bytes += uint64(size)
	}
	b.entries = append(b.entries, index15Entry{b.block.indexKmer, offset})
	return nil
}

// endMask flushes the final block and its sparse primary directory, including empty masks.
func (b *index15Builder) endMask(endOffset uint64) error {
	if err := b.finishBlock(endOffset); err != nil {
		return err
	}
	be.PutUint64(b.scratch[:8], uint64(len(b.entries)))
	if _, err := b.wi.Write(b.scratch[:8]); err != nil {
		return err
	}
	for _, entry := range b.entries {
		be.PutUint64(b.scratch[:8], entry.checkpoint)
		be.PutUint64(b.scratch[8:16], entry.offset)
		if _, err := b.wi.Write(b.scratch[:16]); err != nil {
			return err
		}
	}
	return nil
}
