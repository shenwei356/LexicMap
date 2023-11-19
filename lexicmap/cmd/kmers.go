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
	"fmt"
	"io/fs"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"time"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/lexichash"
	"github.com/shenwei356/lexichash/index"
	"github.com/shenwei356/lexichash/tree"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
)

var kmersCmd = &cobra.Command{
	Use:   "kmers",
	Short: "view k-mers captured by the masks",
	Long: `view k-mers captured by the masks

Attentions:
  1. Mask index (column mask) is 1-based.
  2. Reference indexes (column ref_idx) are 1-based.
  3. K-mer positions (column pos) are 1-based.

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

		var err error

		// ---------------------------------------------------------------

		dbDir := getFlagString(cmd, "index")
		if dbDir == "" {
			checkError(fmt.Errorf("flag -d/--index needed"))
		}
		outFile := getFlagString(cmd, "out-file")

		mask := getFlagNonNegativeInt(cmd, "mask")

		// ---------------------------------------------------------------
		// checking index

		if outputLog {
			log.Info()
			log.Infof("checking index: %s", dbDir)
		}

		// Mask file
		fileMask := filepath.Join(dbDir, index.MaskFile)
		ok, err := pathutil.Exists(fileMask)
		if err != nil || !ok {
			checkError(index.ErrInvalidIndexDir)
		}

		lh, err := lexichash.NewFromFile(fileMask)
		if err != nil {
			checkError(index.ErrInvalidIndexDir)
		}

		if mask > len(lh.Masks) {
			log.Errorf("the index has only %s masks, but %d is given", len(lh.Masks), mask)
		}

		// Trees
		dirTrees := filepath.Join(dbDir, index.TreeDir)
		ok, err = pathutil.DirExists(dirTrees)
		if err != nil || !ok {
			checkError(index.ErrInvalidIndexDir)
		}

		// ---------------------------------------------------------------
		// output file handler
		outfh, gw, w, err := outStream(outFile, strings.HasSuffix(outFile, ".gz"), opt.CompressionLevel)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		// read and output

		decoder := lexichash.MustDecoder()

		fmt.Fprintf(outfh, "mask\tkmer\tlen_v\tref_idx\tpos\tstrand\n")

		var refpos uint64
		var idIdx uint64
		var pos uint64
		var rc uint8
		if mask == 0 {
			treePaths := make([]string, 0, len(lh.Masks))
			fs.WalkDir(os.DirFS(dirTrees), ".", func(p string, d fs.DirEntry, err error) error {
				if filepath.Ext(p) == index.TreeFileExt {
					treePaths = append(treePaths, filepath.Join(dirTrees, p))
				}
				return nil
			})
			if len(treePaths) != len(lh.Masks) {
				checkError(index.ErrTreeFileMissing)
			}

			k := uint8(lh.K)

			var t *tree.Tree
			var base string
			var i int

			for _, file := range treePaths {
				// read tree from the file
				t, err = tree.NewFromFile(file)
				checkError(err)

				base = filepath.Base(file)
				i, err = strconv.Atoi(base[0 : len(base)-len(index.TreeFileExt)])
				if err != nil {
					checkError(index.ErrInvalidIndexDir)
				}

				t.Walk(func(key uint64, v []uint64) bool {
					for _, refpos = range v {
						idIdx = refpos >> 38
						pos = refpos << 26 >> 28
						rc = uint8(refpos & 1)
						fmt.Fprintf(outfh, "%d\t%s\t%d\t%d\t%d\t%c\n",
							i+1, decoder(key, k), len(v), idIdx+1, pos+1, lexichash.Strands[rc])
					}
					return false
				})
			}

			return
		}

		// path of tree file
		idStr := fmt.Sprintf("%04d", mask-1) // convert to 0-based
		subDir := idStr[len(idStr)-2:]
		file := filepath.Join(dbDir, index.TreeDir, subDir, idStr+index.TreeFileExt)

		// read tree from the file
		t, err := tree.NewFromFile(file)
		checkError(err)

		k := uint8(t.K())
		t.Walk(func(key uint64, v []uint64) bool {
			for _, refpos = range v {
				idIdx = refpos >> 38
				pos = refpos << 26 >> 28
				rc = uint8(refpos & 1)
				fmt.Fprintf(outfh, "%d\t%s\t%d\t%d\t%d\t%c\n",
					mask, decoder(key, k), len(v), idIdx, pos+1, lexichash.Strands[rc])
			}
			return false
		})

	},
}

func init() {
	utilsCmd.AddCommand(kmersCmd)

	kmersCmd.Flags().StringP("index", "d", "",
		formatFlagUsage(`Index directory created by "lexicmap index".`))

	kmersCmd.Flags().StringP("out-file", "o", "-",
		formatFlagUsage(`Out file, supports and recommends a ".gz" suffix ("-" for stdout).`))

	kmersCmd.Flags().IntP("mask", "m", 1,
		formatFlagUsage(`View k-mers captured by Xth mask. (0 for all)`))

	kmersCmd.SetUsageTemplate(usageTemplate("-d <index path> -m <mask index> [-o out.tsv.gz]"))
}
