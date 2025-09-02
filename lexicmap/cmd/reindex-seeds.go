// Copyright © 2023-2025 Wei Shen <shenwei356@gmail.com>
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
	"sync"
	"time"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
	"github.com/shenwei356/bio/seq"
	"github.com/spf13/cobra"
	"github.com/vbauerster/mpb/v8"
	"github.com/vbauerster/mpb/v8/decor"
)

var reindexSeedsCmd = &cobra.Command{
	Use:   "reindex-seeds",
	Short: "Recreate indexes of k-mer-value (seeds) data",
	Long: `Recreate indexes of k-mer-value (seeds) data

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		seq.ValidateSeq = false

		// ------------------------------

		dbDir := getFlagString(cmd, "index")
		if dbDir == "" {
			checkError(fmt.Errorf("flag -d/--index needed"))
		}

		partitions := getFlagPositiveInt(cmd, "partitions")

		// ---------------------------------------------------------------

		if opt.Verbose {
			log.Infof("recreating seed indexes with %d partitions for: %s", partitions, dbDir)
		}

		// info file for the number of genome batches
		fileInfo := filepath.Join(dbDir, FileInfo)
		info, err := readIndexInfo(fileInfo)
		if err != nil {
			checkError(fmt.Errorf("failed to read info file: %s", err))
		}

		// ---------------------------------------------------------------

		timeStart := time.Now()
		defer func() {
			if opt.Verbose {
				log.Info()
				log.Infof("elapsed time: %s", time.Since(timeStart))
				log.Info()
			}
		}()

		showProgressBar := opt.Verbose

		// process bar
		var pbs *mpb.Progress
		var bar *mpb.Bar
		var chDuration chan time.Duration
		var doneDuration chan int
		if showProgressBar {
			pbs = mpb.New(mpb.WithWidth(40), mpb.WithOutput(os.Stderr))
			bar = pbs.AddBar(int64(info.Chunks),
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

			chDuration = make(chan time.Duration, opt.NumCPUs)
			doneDuration = make(chan int)
			go func() {
				for t := range chDuration {
					bar.EwmaIncrBy(1, t)
				}
				doneDuration <- 1
			}()
		}

		var wg sync.WaitGroup
		tokens := make(chan int, opt.NumCPUs)
		threadsFloat := float64(opt.NumCPUs)
		for chunk := 0; chunk < info.Chunks; chunk++ {
			file := filepath.Join(dbDir, DirSeeds, chunkFile(chunk))
			wg.Add(1)
			tokens <- 1

			go func(file string) {
				timeStart := time.Now()
				err := kv.CreateKVIndex(file, partitions)
				checkError(err)
				if showProgressBar {
					chDuration <- time.Duration(float64(time.Since(timeStart)) / threadsFloat)
				}
				<-tokens
				wg.Done()
			}(file)
		}
		wg.Wait()

		if showProgressBar {
			close(chDuration)
			<-doneDuration
			pbs.Wait()
		}

		if opt.Verbose {
			log.Infof("update index information file: %s", fileInfo)
		}
		info.Partitions = partitions
		err = writeIndexInfo(fileInfo, info)
		if err != nil {
			checkError(fmt.Errorf("failed to write info file: %s", err))
		}
		if opt.Verbose {
			log.Infof("  finished updating the index information file: %s", fileInfo)
		}
	},
}

func init() {
	utilsCmd.AddCommand(reindexSeedsCmd)

	reindexSeedsCmd.Flags().StringP("index", "d", "",
		formatFlagUsage(`Index directory created by "lexicmap index".`))
	reindexSeedsCmd.Flags().IntP("partitions", "", 4096,
		formatFlagUsage(`Number of partitions for re-indexing seeds (k-mer-value data) files. The value needs to be the power of 4.`))

	reindexSeedsCmd.SetUsageTemplate(usageTemplate(""))
}
