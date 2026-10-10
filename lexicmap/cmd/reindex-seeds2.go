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
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/dustin/go-humanize"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
	"github.com/shenwei356/bio/seq"
	"github.com/spf13/cobra"
	"github.com/vbauerster/mpb/v8"
	"github.com/vbauerster/mpb/v8/decor"
)

var reindexSeeds2Cmd = &cobra.Command{
	Use:   "reindex-seeds2",
	Short: "Recreate adaptive two-level indexes of seeds data",
	Long: `Recreate adaptive two-level indexes of seeds data.

These indexes speed up prefix matching in large indexes by using longer prefixes
to locate seeds within large data blocks.

The seeds files are not changed. The primary prefix length is the mask prefix
plus log4(partitions), and the secondary prefix adds 2 bases. With the standard
7-bp mask prefix and 4096 partitions, these lengths are 13 and 15 bp. Large
primary-prefix blocks receive a secondary index, while smaller blocks retain
the original linear-scan path. The command rewrites each primary .idx file and
creates its matching .idx15 file. Different seeds files are processed in
parallel using -j/--threads.

Overall gains are most noticeable for batch queries with "lexicmap search" when
-n/--top-n-genomes limits downstream alignments, making seed matching a larger
share of runtime. Without this limit, sequence alignment often dominates runtime.
`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		seq.ValidateSeq = false

		defer setupCommandLog(opt, "")()

		dbDir := getFlagString(cmd, "index")
		if dbDir == "" {
			checkError(fmt.Errorf("flag -d/--index needed"))
		}
		partitions := getFlagPositiveInt(cmd, "partitions")
		threshold, err := parseIndex15Threshold(getFlagString(cmd, "threshold"))
		checkError(err)
		checkError(reindexSeeds2(dbDir, partitions, threshold, opt))
	},
}

// parseIndex15Threshold parses the minimum block size shared by index and reindex-seeds2.
func parseIndex15Threshold(value string) (uint64, error) {
	if strings.HasPrefix(strings.TrimSpace(value), "-") {
		return 0, fmt.Errorf("threshold should not be negative")
	}
	threshold, err := ParseByteSize(value)
	if err != nil {
		return 0, err
	}
	if threshold < int64(kv.MinIndex15Threshold) {
		return 0, fmt.Errorf("threshold should be at least 4K (4096 bytes)")
	}
	return uint64(threshold), nil
}

// reindexSeeds2 creates adaptive seed indexes for the final chunks of a database.
// Longer secondary prefixes reduce scanning within large seed blocks in large
// indexes. Batch search gains are most noticeable when -n/--top-n-genomes limits
// downstream alignments; without this limit, alignment often dominates runtime.
func reindexSeeds2(dbDir string, partitions int, threshold uint64, opt *Options) error {
	outputLog := opt.Verbose || opt.Log2File

	fileInfo := filepath.Join(dbDir, FileInfo)
	info, err := readIndexInfo(fileInfo)
	if err != nil {
		return fmt.Errorf("failed to read info file: %s", err)
	}

	if outputLog {
		log.Infof("recreating adaptive seed indexes with %d partitions and a %s block threshold for: %s",
			partitions, humanize.IBytes(threshold), dbDir)
	}
	showProgressBar := opt.Verbose
	var pbs *mpb.Progress
	var chunksBar *mpb.Bar
	if showProgressBar {
		pbs = mpb.New(mpb.WithWidth(40), mpb.WithOutput(os.Stderr))
		chunksBar = pbs.AddBar(int64(info.Chunks),
			mpb.BarPriority(0),
			mpb.PrependDecorators(
				decor.Name("processed files: ", decor.WC{W: len("processed files: "), C: decor.DindentRight}),
				decor.Name("", decor.WCSyncSpaceR),
				decor.CountersNoUnit("%d / %d", decor.WCSyncWidth),
			),
			mpb.AppendDecorators(
				decor.Name("ETA: ", decor.WC{W: len("ETA: ")}),
				decor.EwmaETA(decor.ET_STYLE_GO, 3),
				decor.OnComplete(decor.Name(""), ". done"),
			),
		)
	}

	// reindexResult carries one chunk's summary and completion status.
	type reindexResult struct {
		stats    kv.Index15Stats // block counts and secondary index size
		duration time.Duration   // elapsed time for progress estimates
		err      error           // chunk-specific indexing failure
	}
	workers := opt.NumCPUs
	if workers > info.Chunks {
		workers = info.Chunks
	}
	if workers < 1 {
		workers = 1
	}
	jobs := make(chan int, info.Chunks)
	results := make(chan reindexResult, info.Chunks)
	// Reindex independent seeds files in parallel. CreateKVIndex15 performs
	// one sequential scan within each file.
	for i := 0; i < workers; i++ {
		go func() {
			for chunk := range jobs {
				start := time.Now()
				file := filepath.Join(dbDir, DirSeeds, chunkFile(chunk))
				var masksBar *mpb.Bar
				var progress kv.Index15ProgressFunc
				if showProgressBar {
					progress = func(processedMasks, totalMasks uint64) {
						if masksBar == nil {
							name := fmt.Sprintf("  chunk %03d masks: ", chunk)
							masksBar = pbs.AddBar(int64(totalMasks),
								mpb.BarPriority(chunk+1),
								mpb.BarRemoveOnComplete(),
								mpb.PrependDecorators(
									decor.Name(name, decor.WC{W: len(name), C: decor.DindentRight}),
									decor.CountersNoUnit("%d / %d"),
								),
								mpb.AppendDecorators(decor.Percentage()),
							)
						}
						masksBar.SetCurrent(int64(processedMasks))
					}
				}
				stats, err := kv.CreateKVIndex15WithProgress(file, partitions, threshold, progress)
				if masksBar != nil {
					if err != nil {
						masksBar.Abort(true)
					} else {
						masksBar.SetTotal(-1, true)
					}
				}
				if err != nil {
					err = fmt.Errorf("%s: %w", file, err)
				}
				results <- reindexResult{stats: stats, duration: time.Since(start), err: err}
			}
		}()
	}
	for chunk := 0; chunk < info.Chunks; chunk++ {
		jobs <- chunk
	}
	close(jobs)

	var total kv.Index15Stats
	var firstErr error
	threadsFloat := float64(workers)
	for i := 0; i < info.Chunks; i++ {
		result := <-results
		if result.err != nil && firstErr == nil {
			firstErr = result.err
		}
		total.Add(result.stats)
		if showProgressBar {
			chunksBar.EwmaIncrBy(1, time.Duration(float64(result.duration)/threadsFloat))
		}
	}
	if showProgressBar {
		pbs.Wait()
	}
	if firstErr != nil {
		return firstErr
	}

	info.Partitions = partitions
	if err = writeIndexInfo(fileInfo, info); err != nil {
		return fmt.Errorf("failed to update index information file: %s", err)
	}

	if outputLog {
		percentage := float64(0)
		if total.NonEmptyBlocks > 0 {
			percentage = float64(total.IndexedBlocks) / float64(total.NonEmptyBlocks) * 100
		}
		log.Infof("primary-prefix blocks: %s total, %s non-empty, %s indexed (%.2f%% of non-empty)",
			humanize.Comma(int64(total.TotalBlocks)), humanize.Comma(int64(total.NonEmptyBlocks)),
			humanize.Comma(int64(total.IndexedBlocks)), percentage)
		log.Infof("idx15 size: %s", humanize.IBytes(total.Index15Bytes))
		for width, count := range total.WidthCounts {
			if count > 0 {
				log.Infof("  offset width %d byte(s): %d blocks", width+1, count)
			}
		}
	}
	return nil
}

func init() {
	utilsCmd.AddCommand(reindexSeeds2Cmd)

	reindexSeeds2Cmd.Flags().StringP("index", "d", "",
		formatFlagUsage(`Index directory created by "lexicmap index".`))
	reindexSeeds2Cmd.Flags().Int("partitions", 4096,
		formatFlagUsage(`Number of anchor partitions. The value needs to be a power of 4 and at least 4. The primary prefix is mask-prefix + log4(partitions), and the secondary prefix adds 2 bases; the 5-byte checkpoint suffix also requires K - primary-prefix <= 20. The default 4096 produces 13/15-bp prefixes with the standard 7-bp mask prefix.`))
	reindexSeeds2Cmd.Flags().String("threshold", "8K",
		formatFlagUsage(`Minimum size of a primary-prefix seeds block for creating its secondary index. Higher values reduce .idx15 size but limit the speedup. The minimum is 4K because seed readers use 4 KiB buffers. The default 8K is recommended.`))

	reindexSeeds2Cmd.SetUsageTemplate(usageTemplate(""))
}
