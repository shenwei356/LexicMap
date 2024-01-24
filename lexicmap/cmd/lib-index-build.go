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
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OFTestSerializationTestSerialization ANY KIND, EXPRESS OR
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
	"time"

	"github.com/shenwei356/lexichash"
	"github.com/vbauerster/mpb/v8"
	"github.com/vbauerster/mpb/v8/decor"
)

// TmpDirExt is the path extension for temporary files
const TmpDirExt = ".tmp"

type IndexBuildingOptions struct {
	// general
	NumCPUs      int
	Verbose      bool // show log
	Force        bool // force overwrite existed index
	MaxOpenFiles int  // maximum opened files, used in merging indexes

	// LexicHash

	K                int   // k-mer size
	Masks            int   // number of masks
	RandSeed         int64 // random seed
	PrefixForCheckLC int   // length of prefix for checking low-complexity

	// k-mer index

	Chunks int // the number of chunks for storing k-mer index

	// genome batches

	GenomeBatchSize int // the maximum number of genomes of a batch

	// genome

	ReRefName    *regexp.Regexp
	ReSeqExclude []*regexp.Regexp
}

// CheckIndexBuildingOptions check the options
func CheckIndexBuildingOptions(opt *IndexBuildingOptions) error {
	if opt.K < 3 || opt.K > 32 {
		return fmt.Errorf("invalid k value: %d, valid range: [3, 32]", opt.K)
	}
	if opt.Masks < 4 {
		return fmt.Errorf("invalid numer of masks: %d, should be >=4", opt.Masks)
	}
	if opt.PrefixForCheckLC > opt.K {
		return fmt.Errorf("invalid prefix: %d, valid range: [0, k], 0 for no checking", opt.PrefixForCheckLC)
	}

	if opt.Chunks < 1 || opt.Chunks > 512 {
		return fmt.Errorf("invalid chunks: %d, valid range: [1, 512]", opt.Chunks)
	}

	if opt.GenomeBatchSize < 1 || opt.GenomeBatchSize > 1<<17 {
		return fmt.Errorf("invalid genome batch size: %d, valid range: [1, 131072]", opt.GenomeBatchSize)
	}

	if opt.MaxOpenFiles < 2 {
		return fmt.Errorf("invalid max open files: %d, should be >= 2", opt.MaxOpenFiles)
	}

	return nil
}

func BuildIndex(outdir string, infiles []string, opt *IndexBuildingOptions) error {
	// check options
	// err := CheckIndexBuildingOptions(opt)
	// if err != nil {
	// 	return err
	// }

	// generate masks
	lh, err := lexichash.NewWithSeed(opt.K, opt.Masks, opt.RandSeed, opt.PrefixForCheckLC)
	if err != nil {
		return err
	}

	// tmp dir
	tmpDir := filepath.Clean(outdir) + TmpDirExt
	err = os.RemoveAll(tmpDir)
	if err != nil {
		return err
	}

	// split the files in to batches
	nFiles := len(infiles)
	nBatches := (nFiles + opt.GenomeBatchSize - 1) / opt.GenomeBatchSize
	tmpIndexes := make([]string, 0, nBatches)

	var begin, end int
	for batch := 0; batch < nBatches; batch++ {
		// files for this batch
		begin = batch * opt.GenomeBatchSize
		end = begin + opt.GenomeBatchSize
		if end > nFiles {
			end = nFiles
		}
		files := infiles[begin:end]

		// outdir for this batch
		var outdirB string
		if nBatches > 1 {
			outdirB = filepath.Join(tmpDir, fmt.Sprintf("batch_%4d", batch))
			tmpIndexes = append(tmpIndexes, outdirB)
		} else {
			outdirB = outdir
		}

		// build index for this batch
		buildAnIndex(lh, opt, outdirB, files)
	}

	if nBatches == 1 {
		return nil
	}

	// merge indexes
	mergeIndexes(lh, opt, outdir, tmpIndexes)

	return nil
}

func buildAnIndex(lh *lexichash.LexicHash, opt *IndexBuildingOptions, outdir string, files []string) {
	// process bar
	var pbs *mpb.Progress
	var bar *mpb.Bar
	var chDuration chan time.Duration
	var doneDuration chan int
	if opt.Verbose {
		pbs = mpb.New(mpb.WithWidth(40), mpb.WithOutput(os.Stderr))
		bar = pbs.AddBar(int64(len(files)),
			mpb.PrependDecorators(
				decor.Name("processed files: ", decor.WC{W: len("processed files: "), C: decor.DindentRight}),
				decor.Name("", decor.WCSyncSpaceR),
				decor.CountersNoUnit("%d / %d", decor.WCSyncWidth),
			),
			mpb.AppendDecorators(
				decor.Name("ETA: ", decor.WC{W: len("ETA: ")}),
				decor.EwmaETA(decor.ET_STYLE_GO, 10),
				decor.OnComplete(decor.Name(""), ". done"),
			),
		)

		chDuration = make(chan time.Duration, opt.NumCPUs)
		doneDuration = make(chan int)
		go func() {
			for t := range chDuration {
				bar.Increment()
				bar.EwmaIncrBy(1, t)
			}
			doneDuration <- 1
		}()
	}

	if opt.Verbose {
		close(chDuration)
		<-doneDuration
		pbs.Wait()
	}

}

func mergeIndexes(lh *lexichash.LexicHash, opt *IndexBuildingOptions, outdir string, paths []string) {

}
