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
	"regexp"
	"sort"
	"time"

	"github.com/pkg/errors"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/lexichash"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
)

var remergeCmd = &cobra.Command{
	Use:   "remerge",
	Short: "Rerun the merging step for an unfinished index",
	Long: `Rerun the merging step for an unfinished index

When to use this command?

- Only one thread is used for merging indexes, which happens when there are
  a lot (>200 batches) of batches ($inpu_files / --batch-size) and the value
  of --max-open-files is not big enough. E.g.,

  22:54:24.420 [INFO] merging 297 indexes...
  22:54:24.455 [INFO]   [round 1]
  22:54:24.455 [INFO]     batch 1/1, merging 297 indexes to xxx.lmi.tmp/r1_b1 with 1 threads...

  ► Then you can run this command with a bigger --max-open-files (e.g., 4096) and 
  -J/--seed-data-threads (e.g., 12. 12 needs be <= 4096/(297+2)=13.7).
  And you need to set a bigger 'ulimit -n' if the value of --max-open-files is bigger than 1024.

- The Slurm/PBS job time limit is almost reached and the merging step won't be finished before that.

- Disk quota is reached in the merging step.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		seq.ValidateSeq = false

		var fhLog *os.File
		if opt.Log2File {
			fhLog = addLog(opt.LogFile, opt.Verbose)
		}

		outputLog := opt.Verbose || opt.Log2File

		timeStart := time.Now()
		defer func() {
			if outputLog {
				log.Info()
				log.Infof("elapsed time: %s", time.Since(timeStart))
				log.Info()
			}
			if opt.Log2File {
				fhLog.Close()
			}
		}()

		// ---------------------------------------------------------------

		dbDir := getFlagString(cmd, "index")
		if dbDir == "" {
			checkError(fmt.Errorf("index directory is need"))
		}

		tmpDir := filepath.Clean(dbDir) + ExtTmpDir
		ok, err := pathutil.DirExists((tmpDir))
		if err != nil {
			checkError(fmt.Errorf("index directory is need"))
		}
		if !ok {
			checkError(fmt.Errorf("tmp directory is not found: %s", tmpDir))
		}

		mergeThreads := getFlagPositiveInt(cmd, "seed-data-threads")

		maxOpenFiles := getFlagPositiveInt(cmd, "max-open-files")

		// ---------------------------------------------------------------
		// check indexes of all batches

		if opt.Verbose || opt.Log2File {
			log.Infof("checking indexes ...")
		}

		// batch dirs
		batchDirs := make([]string, 0, 512)
		pattern := regexp.MustCompile(`^batch_\d+$`)
		files, err := os.ReadDir(tmpDir)
		if err != nil {
			checkError(errors.Errorf("failed to read dir: %s", err))
		}
		for _, file := range files {
			if file.Name() == "." || file.Name() == ".." {
				continue
			}
			if file.IsDir() && pattern.MatchString(file.Name()) {
				batchDirs = append(batchDirs, filepath.Join(tmpDir, file.Name()))
			}
		}

		if len(batchDirs) == 0 {
			checkError(fmt.Errorf("no indexes found in %s", tmpDir))
		} else if opt.Verbose || opt.Log2File {
			log.Infof("  %d index directries found in %s", len(batchDirs), tmpDir)
		}

		// ---------------------------------------------------------------
		// prepare arguments for mergeIndexes

		sort.Strings(batchDirs)
		OneIndex := batchDirs[0]

		// lh *lexichash.LexicHash,                 read from one batch
		fileMask := filepath.Join(OneIndex, FileMasks)
		ok, err = pathutil.Exists(fileMask)
		if err != nil || !ok {
			checkError(fmt.Errorf("mask file not found: %s. Was the index merged?", fileMask))
		}
		var lh *lexichash.LexicHash
		lh, err = lexichash.NewFromFile(fileMask)
		if err != nil {
			checkError(fmt.Errorf("checking mask file: %s", err))
		}
		// fmt.Println(len(lh.Masks))

		// maskPrefix uint8, anchorPrefix uint8,    read from one batch with ReadKVIndexInfo
		var maskPrefix, anchorPrefix uint8
		fileSeedChunk := filepath.Join(OneIndex, DirSeeds, chunkFile(0))
		_, _, _, maskPrefix, anchorPrefix, err = kv.ReadKVIndexInfo(filepath.Clean(fileSeedChunk) + kv.KVIndexFileExt)
		if err != nil {
			checkError(fmt.Errorf("checking seed information: %s", err))
		}
		// fmt.Println(maskPrefix, anchorPrefix)

		// kvChunks int,                            read from one batch, info file
		var info *IndexInfo
		info, err = readIndexInfo(filepath.Join(OneIndex, FileInfo))
		if err != nil {
			checkError(fmt.Errorf("failed to open info file: %s", err))
		}
		kvChunks := info.Chunks
		if mergeThreads > kvChunks {
			mergeThreads = kvChunks
		}

		// opt *IndexBuildingOptions,               create one, used: opt.Verbose, opt.Log2File, opt.MaxOpenFiles, opt.MergeThreads
		bopt := &IndexBuildingOptions{
			// general
			NumCPUs:      opt.NumCPUs,
			Verbose:      opt.Verbose,
			Log2File:     opt.Log2File,
			MaxOpenFiles: maxOpenFiles,
			MergeThreads: mergeThreads,
		}

		// outdir string,                           dbDir
		// paths []string,                          batchDirs
		// tmpDir string,                           tmpDir
		// round int                                1

		err = mergeIndexes(lh, maskPrefix, anchorPrefix, bopt, kvChunks, dbDir, batchDirs, tmpDir, 1)
		if err != nil {
			checkError(fmt.Errorf("failed to merge indexes: %s", err))
		}

		// clean tmp dir
		err = os.RemoveAll(tmpDir)
		if err != nil {
			checkError(fmt.Errorf("failed to remove tmp directory: %s", err))
		}
	},
}

func init() {
	utilsCmd.AddCommand(remergeCmd)

	remergeCmd.Flags().StringP("index", "d", "",
		formatFlagUsage(`Index directory created by "lexicmap index".`))

	remergeCmd.Flags().IntP("seed-data-threads", "J", 8,
		formatFlagUsage(`Number of threads for writing seed data and merging seed chunks from all batches, the value should be in range of [1, -c/--chunks]. If there are >100 batches, please also increase the value of --max-open-files and set a bigger "ulimit -n" in shell.`))

	remergeCmd.Flags().IntP("max-open-files", "", 1024,
		formatFlagUsage(`Maximum opened files, used in merging indexes. If there are >100 batches, please increase this value and set a bigger "ulimit -n" in shell.`))

	remergeCmd.SetUsageTemplate(usageTemplate("[flags] -d <index path>"))
}
