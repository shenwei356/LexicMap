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
	"errors"
	"fmt"
	"path/filepath"
	"sync"
	"time"

	"github.com/dustin/go-humanize"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
	"github.com/shenwei356/bio/seq"
	"github.com/spf13/cobra"
)

var countbasesCmd = &cobra.Command{
	Use:   "recount-bases",
	Short: "Recount bases for index version <=3.2",
	Long: `Recount bases for index version <=3.2

This command is only needed for indexes created by LexicMap v0.6.0 (3c257ca) or before versions.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		seq.ValidateSeq = false

		defer setupCommandLog(opt, "")()
		outputLog := opt.Verbose || opt.Log2File

		// ------------------------------

		dbDir := getFlagString(cmd, "index")
		if dbDir == "" {
			checkError(fmt.Errorf("flag -d/--index needed"))
		}

		// info file
		fileInfo := filepath.Join(dbDir, FileInfo)
		info, err := readIndexInfo(fileInfo)
		if err != nil {
			checkError(fmt.Errorf("failed to read info file: %s", err))
		}
		if info.MainVersion != MainVersion {
			checkError(fmt.Errorf("index main versions do not match: %d (index) != %d (tool). please re-create the index", info.MainVersion, MainVersion))
		}

		startTime := time.Now()

		old := info.InputBases
		totalBases, err := updateInputBases(info, dbDir, opt.NumCPUs)
		checkError(err)

		if outputLog {
			log.Infof("update input bases from %d to %s in %s", old, humanize.Comma(totalBases), time.Since(startTime))
		}
	},
}

func init() {
	// utilsCmd.AddCommand(countbasesCmd)

	countbasesCmd.Flags().StringP("index", "d", "",
		formatFlagUsage(`Index directory created by "lexicmap index".`))

	countbasesCmd.SetUsageTemplate(usageTemplate(""))
}

func updateInputBases(info *IndexInfo, dbDir string, threads int) (int64, error) {
	// Publish totals only after every batch has been read and closed successfully.
	var totalBases int64
	var readErr error
	var mu sync.Mutex // guards the total and first error from concurrent batches

	// extract genome sizes
	var wg sync.WaitGroup
	tokens := make(chan int, threads)
	for i := 0; i < info.GenomeBatches; i++ {
		wg.Add(1)
		tokens <- 1
		go func(i int) {
			defer func() {
				<-tokens
				wg.Done()
			}()
			fileGenomes := filepath.Join(dbDir, DirGenomes, batchDir(i), FileGenomes)
			var bases int64
			rdr, err := genome.NewReader(fileGenomes)
			if err == nil {
				bases, err = rdr.TotalBases()
				// Close before releasing the slot, including failed reads. Keep both errors.
				err = errors.Join(err, rdr.Close())
			}

			mu.Lock()
			defer mu.Unlock()
			if err != nil {
				if readErr == nil {
					readErr = fmt.Errorf("failed to count bases in %s: %w", fileGenomes, err)
				}
				return
			}
			totalBases += bases
		}(i)
	}
	wg.Wait()
	if readErr != nil {
		return 0, readErr
	}

	// update info file
	info.InputBases = totalBases

	err := writeIndexInfo(filepath.Join(dbDir, FileInfo), info)
	if err != nil {
		return 0, (fmt.Errorf("failed to write info file: %s", err))
	}

	return totalBases, err
}
