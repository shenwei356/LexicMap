// Copyright © 2023-2024 Wei Shen <shenwei356@gmail.com>
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
	"io/fs"
	"os"
	"path/filepath"
	"strings"
	"sync"
	"time"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
	"github.com/shenwei356/lexichash"
)

// mergeIndexes merge multiple indexes to a big one
func mergeIndexes(lh *lexichash.LexicHash, opt *IndexBuildingOptions, kvChunks int,
	outdir string, paths []string, tmpDir string, round int) error {
	timeStart := time.Now()
	if opt.Verbose || opt.Log2File {
		log.Infof("  [round %d]", round)
	}

	nIndexes := len(paths)
	chunkSize := opt.MaxOpenFiles // chunk of indexes
	if chunkSize > nIndexes {
		chunkSize = nIndexes
	}
	batches := (len(paths) + chunkSize - 1) / chunkSize

	var j, begin, end int
	tmpIndexes := make([]string, 0, 8)

	var mergeThreads int
	var wg sync.WaitGroup

	var pathB []string

	for j = 0; j < batches; j++ { // each chunk for storing kmer-value data
		begin = j * chunkSize
		end = begin + chunkSize
		if end > nIndexes {
			end = nIndexes
		}
		pathB = paths[begin:end]

		mergeThreads = opt.MergeThreads
		for mergeThreads*(len(pathB)+2) > opt.MaxOpenFiles { // 2 is for output file and index file
			mergeThreads--
		}
		if mergeThreads < 1 {
			mergeThreads = 1
		}
		tokens := make(chan int, mergeThreads)

		outdir1 := filepath.Join(tmpDir, fmt.Sprintf("r%d_b%d", round, j+1))

		if opt.Verbose || opt.Log2File {
			log.Infof("    batch %d/%d, merging %d indexes to %s with %d threads...", j+1, batches, end-begin, outdir1, mergeThreads)
		}

		tmpIndexes = append(tmpIndexes, outdir1)

		// --------------------------------------------------------------------
		// data structure

		err := os.MkdirAll(outdir1, 0755)
		if err != nil {
			checkError(fmt.Errorf("failed to create dir: %s", err))
		}

		// seeds
		dirSeeds := filepath.Join(outdir1, DirSeeds)
		err = os.MkdirAll(dirSeeds, 0755)
		if err != nil {
			checkError(fmt.Errorf("failed to create dir: %s", err))
		}

		// genomes
		dirGenomes := filepath.Join(outdir1, DirGenomes)
		err = os.MkdirAll(dirGenomes, 0755)
		if err != nil {
			checkError(fmt.Errorf("failed to create dir: %s", err))
		}

		// --------------------------------------------------------------------
		// kmer-value data

		for chunk := 0; chunk < kvChunks; chunk++ {
			tokens <- 1
			wg.Add(1)

			go func(chunk int) {
				defer func() {
					wg.Done()
					<-tokens
				}()

				var rdr *kv.Reader
				var i int
				var kmer uint64
				var values, values1 *[]uint64
				var ok bool

				// read information from an existing index file
				fileIdx := filepath.Join(pathB[0], DirSeeds, chunkFile(chunk)+kv.KVIndexFileExt)
				rdrIdx, err := kv.NewIndexReader(fileIdx)
				if err != nil {
					checkError(fmt.Errorf("failed to read info from an index file: %s", err))
				}

				// outfile
				file := filepath.Join(dirSeeds, chunkFile(chunk))
				wtr, err := kv.NewWriter(rdrIdx.K, rdrIdx.ChunkIndex, rdrIdx.ChunkSize, file)
				if err != nil {
					checkError(fmt.Errorf("failed to write a k-mer data file: %s", err))
				}

				rdrs := make([]*kv.Reader, len(pathB))
				for i, db := range pathB {
					rdrs[i], err = kv.NewReader(filepath.Join(db, DirSeeds, chunkFile(chunk)))
					if err != nil {
						checkError(fmt.Errorf("failed to read kv-data file: %s", err))
					}
				}

				m := kv.PoolKmerData.Get().(*map[uint64]*[]uint64)
				for c := 0; c < rdrIdx.ChunkSize; c++ { // for all mask
					clear(*m)

					for i, rdr = range rdrs {
						m1, err := rdr.ReadDataOfAMaskAsMap()
						if err != nil {
							checkError(fmt.Errorf("failed to read data of mask %d from file %s: %s",
								c+rdr.ChunkIndex, pathB[i], err))
						}

						for kmer, values1 = range *m1 {
							if values, ok = (*m)[kmer]; !ok {
								tmp := make([]uint64, 0, len(*values1))
								values = &tmp
								(*m)[kmer] = values
							}
							*values = append(*values, (*values1)...)
						}
						kv.RecycleKmerData(m1)
					}

					err = wtr.WriteDataOfAMask(*m, opt.Partitions)
					if err != nil {
						checkError(fmt.Errorf("failed to write to k-mer data file: %s", err))
					}
				}
				kv.RecycleKmerData(m)

				for _, rdr = range rdrs {
					err = rdr.Close()
					if err != nil {
						checkError(fmt.Errorf("failed to close kv-data file: %s", err))
					}
				}

				rdrIdx.Close()

				err = wtr.Close()
				if err != nil {
					checkError(fmt.Errorf("failed to close kv-data file: %s", err))
				}

			}(chunk)
		}
		wg.Wait()

		// -------------------------------------------------------------------
		// genomes/, just move
		var dirGenomesIn, dirG string
		var files []fs.DirEntry
		var file fs.DirEntry
		for _, db := range pathB {
			dirGenomesIn = filepath.Join(db, DirGenomes)
			files, err = os.ReadDir(dirGenomesIn)
			if err != nil {
				checkError(fmt.Errorf("failed to read genome dir: %s", err))
			}
			for _, file = range files {
				dirG = file.Name()
				if file.IsDir() && strings.HasPrefix(dirG, "batch_") {
					err = os.Rename(filepath.Join(dirGenomesIn, dirG), filepath.Join(dirGenomes, dirG))
					if err != nil {
						checkError(fmt.Errorf("failed to move genome data"))
					}
				}
			}
		}

		// -------------------------------------------------------------------
		// genomes.map.bin, just concatenate them
		fh, err := os.Create(filepath.Join(outdir1, FileGenomeIndex))
		if err != nil {
			checkError(fmt.Errorf("failed to write genome index mapping file: %s", err))
		}
		bw := bufio.NewWriter(fh)
		for _, db := range pathB {
			fh1, err := os.Open(filepath.Join(db, FileGenomeIndex))
			if err != nil {
				checkError(fmt.Errorf("failed to open genome index mapping file: %s", err))
			}
			br := bufio.NewReader(fh1)
			_, err = io.Copy(bw, br)
			if err != nil {
				checkError(fmt.Errorf("failed to copy genome index mapping data: %s", err))
			}
			err = fh1.Close()
			if err != nil {
				checkError(fmt.Errorf("failed to close genome index mapping file: %s", err))
			}
		}
		bw.Flush()
		err = fh.Close()
		if err != nil {
			checkError(fmt.Errorf("failed to close genome index mapping file: %s", err))
		}

		// -------------------------------------------------------------------
		// info.toml, copy one and update the genome number
		info, err := readIndexInfo(filepath.Join(pathB[0], FileInfo))
		if err != nil {
			checkError(fmt.Errorf("failed to open info file: %s", err))
		}

		for _, db := range pathB[1:] {
			info2, err := readIndexInfo(filepath.Join(db, FileInfo))
			if err != nil {
				checkError(fmt.Errorf("failed to open info file: %s", err))
			}

			info.Genomes += info2.Genomes
			info.GenomeBatches += info2.GenomeBatches
		}

		err = writeIndexInfo(filepath.Join(outdir1, FileInfo), info)
		if err != nil {
			checkError(fmt.Errorf("failed to write info file: %s", err))
		}

		// -------------------------------------------------------------------
		// masks.bin, just copy one
		err = os.Rename(filepath.Join(pathB[0], FileMasks), filepath.Join(outdir1, FileMasks))
		if err != nil {
			checkError(fmt.Errorf("failed to move genome data"))
		}

	}

	if opt.Verbose || opt.Log2File {
		log.Infof("  [round %d] finished in %s", round, time.Since(timeStart))
	}

	if len(tmpIndexes) == 1 {
		// move it to outdir
		if opt.Verbose || opt.Log2File {
			log.Infof("rename %s to %s", tmpIndexes[0], outdir)
		}

		// delete old one, actually it's empty
		err := os.RemoveAll(outdir)
		if err != nil {
			checkError(fmt.Errorf("failed to remove empty directory: %s", err))
		}

		err = os.Rename(tmpIndexes[0], outdir)
		if err != nil {
			checkError(fmt.Errorf("failed to move index directory: %s", err))
		}

		return nil
	}

	mergeIndexes(lh, opt, kvChunks, outdir, tmpIndexes, tmpDir, round+1)
	return nil
}
