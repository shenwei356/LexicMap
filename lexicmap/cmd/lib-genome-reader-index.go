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
	"path/filepath"
	"sync"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
)

// sharedGenomeReaderIndex publishes one immutable directory per batch. Failed
// loads are retryable; only successful metadata is retained by the Index.
type sharedGenomeReaderIndex struct {
	mu    sync.Mutex
	index *genome.ReaderIndex
}

// newGenomeReader shares metadata in both pooled and file-budget fallback paths.
// File descriptors and mutable decoding buffers still belong to each reader.
func (idx *Index) newGenomeReader(batch int) (*genome.Reader, error) {
	value, ok := idx.genomeReaderIndexes.Load(batch)
	if !ok {
		value, _ = idx.genomeReaderIndexes.LoadOrStore(batch, &sharedGenomeReaderIndex{})
	}
	slot := value.(*sharedGenomeReaderIndex)
	slot.mu.Lock()
	index := slot.index
	if index == nil {
		file := filepath.Join(idx.path, DirGenomes, batchDir(batch), FileGenomes)
		var err error
		index, err = genome.ReadIndex(file)
		if err != nil {
			slot.mu.Unlock()
			return nil, err
		}
		slot.index = index
	}
	slot.mu.Unlock()
	return index.NewReader()
}
