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
	"fmt"
	"path/filepath"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
)

// validateSeedIndexInfo bounds seed work queues before allocation or publication.
// Genome-only commands do not require seed files or their summary fields.
func validateSeedIndexInfo(info *IndexInfo) error {
	if info.MainVersion != MainVersion {
		return fmt.Errorf("index main versions do not match: %d (index) != %d (tool)", info.MainVersion, MainVersion)
	}
	// Index construction supports at most 128 nonempty seed chunks.
	if info.K == 0 || info.K > 32 || info.Masks < 1 || info.Chunks < 1 || info.Chunks > 128 || info.Chunks > info.Masks {
		return fmt.Errorf("invalid index summary: k=%d, masks=%d, chunks=%d (expected k in [1, 32], positive masks, and chunks in [1, min(128, masks)])", info.K, info.Masks, info.Chunks)
	}
	return nil
}

// validateGenomeIndexInfo bounds genome pools and rejects invalid genome totals.
// Empty intermediate batches and old indexes without InputBases remain valid.
func validateGenomeIndexInfo(info *IndexInfo) error {
	if info.MainVersion != MainVersion {
		return fmt.Errorf("index main versions do not match: %d (index) != %d (tool)", info.MainVersion, MainVersion)
	}
	if info.GenomeBatches < 1 || info.GenomeBatches > 1<<BITS_BATCH_IDX {
		return fmt.Errorf("invalid index summary: genome batches must be in [1, %d], got %d", 1<<BITS_BATCH_IDX, info.GenomeBatches)
	}
	if info.Genomes < 0 || uint64(info.Genomes) > uint64(info.GenomeBatches)<<BITS_GENOME_IDX ||
		info.InputGenomes < 0 || info.InputBases < 0 || info.ContigInterval < 0 {
		return fmt.Errorf("invalid index summary: genomes=%d, input genomes=%d, input bases=%d, contig interval=%d", info.Genomes, info.InputGenomes, info.InputBases, info.ContigInterval)
	}
	return nil
}

// seedChunkStarts keeps only the mask starts needed by sequential seed scans.
// Sequential scans do not need a full prefix lookup table for every mask.
type seedChunkStarts struct {
	starts [][2]uint64 // first k-mer and encoded file offset; zero offset means empty
	config uint8       // seed position width from the data header
}

// readSeedChunkStarts validates the selected data/index layouts before output opens.
// Retaining starts costs O(masks), independent of the number of prefix partitions.
func readSeedChunkStarts(dir string, info *IndexInfo, k, masks, selected int) (int, map[int]seedChunkStarts, error) {
	if err := validateSeedIndexInfo(info); err != nil {
		return 0, nil, err
	}
	if info.Masks != masks || int(info.K) != k {
		return 0, nil, fmt.Errorf("invalid index summary: chunks=%d, masks=%d, k=%d; mask file has %d masks with k=%d",
			info.Chunks, info.Masks, info.K, masks, k)
	}
	chunkSize := 1 + (masks-1)/info.Chunks
	first, last := 0, (masks-1)/chunkSize
	if selected > 0 {
		first, last = (selected-1)/chunkSize, (selected-1)/chunkSize
	}
	chunks := make(map[int]seedChunkStarts, last-first+1)
	for chunk := first; chunk <= last; chunk++ {
		file := filepath.Join(dir, DirSeeds, chunkFile(chunk))
		begin := chunk * chunkSize
		count := min(chunkSize, masks-begin)
		config, err := checkSeedChunkLayout(file, k, begin, count)
		if err != nil {
			return 0, nil, err
		}
		indexFile := file + kv.KVIndexFileExt
		starts, err := kv.ReadKVIndexStarts(indexFile)
		if err != nil {
			return 0, nil, fmt.Errorf("read seed index starts %s: %w", indexFile, err)
		}
		chunks[chunk] = seedChunkStarts{starts: starts, config: config}
	}
	return chunkSize, chunks, nil
}

// checkSeedChunkLayout compares summary-derived mask ranges with both file headers.
// Prefix partition counts may differ after an interrupted rebuild; they are not compared.
func checkSeedChunkLayout(file string, k, begin, count int) (uint8, error) {
	dataK, dataBegin, dataCount, config, err := kv.ReadKVDataHeader(file)
	if err != nil {
		return 0, fmt.Errorf("read seed header %s: %w", file, err)
	}
	indexFile := file + kv.KVIndexFileExt
	indexK, indexBegin, indexCount, _, _, err := kv.ReadKVIndexInfo(indexFile)
	if err != nil {
		return 0, fmt.Errorf("read seed index header %s: %w", indexFile, err)
	}
	if int(dataK) != k || int(indexK) != k || dataBegin != begin || indexBegin != begin || dataCount != count || indexCount != count {
		return 0, fmt.Errorf("invalid seed chunk layout in %s: expected k=%d, first mask=%d, masks=%d; data=(%d,%d,%d), index=(%d,%d,%d)",
			file, k, begin, count, dataK, dataBegin, dataCount, indexK, indexBegin, indexCount)
	}
	return config, nil
}

// checkSeedIndexLayout validates all chunk headers before an in-place rebuild.
// It reads only metadata; seed contents are checked while rebuilding each chunk.
func checkSeedIndexLayout(dir string, info *IndexInfo) error {
	if err := validateSeedIndexInfo(info); err != nil {
		return err
	}
	width := 1 + (info.Masks-1)/info.Chunks
	for chunk := 0; chunk < info.Chunks; chunk++ {
		begin := chunk * width
		count := max(0, min(width, info.Masks-begin)) // old indexes can have empty tails
		if _, err := checkSeedChunkLayout(filepath.Join(dir, DirSeeds, chunkFile(chunk)), int(info.K), begin, count); err != nil {
			return err
		}
	}
	return nil
}
