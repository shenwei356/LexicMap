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
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/lexichash"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
)

var masksCmd = &cobra.Command{
	Use:   "masks",
	Short: "View masks of the index or generate some new",
	Long: `View masks of the index or generate some new

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

		outFile := getFlagString(cmd, "out-file")

		k := getFlagPositiveInt(cmd, "kmer")
		if k < minK || k > 32 {
			checkError(fmt.Errorf("the value of flag -k/--kmer should be in range of [%d, 32]", minK))
		}

		nMasks := getFlagPositiveInt(cmd, "masks")
		lcPrefix := getFlagNonNegativeInt(cmd, "prefix")
		seed := getFlagPositiveInt(cmd, "seed")

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

		// ---------------------------------------------------------------

		var lh *lexichash.LexicHash

		if dbDir != "" { // from the index
			if outputLog {
				log.Info()
				log.Infof("checking index: %s", dbDir)
			}

			// Mask file
			fileMask := filepath.Join(dbDir, FileMasks)
			ok, err := pathutil.Exists(fileMask)
			if err != nil || !ok {
				checkError(fmt.Errorf("mask file not found: %s", fileMask))
			}

			lh, err = lexichash.NewFromFile(fileMask)
			if err != nil {
				checkError(fmt.Errorf("%s", err))
			}

			if outputLog {
				log.Infof("  checking passed")
				log.Infof("reading masks...")
			}
		} else { // re generate
			if outputLog {
				log.Infof("generating new mask...")
			}
			lh, err = lexichash.NewWithSeed(k, nMasks, int64(seed), lcPrefix)
			checkError(err)
		}

		decoder := lexichash.MustDecoder()
		_k := uint8(lh.K)

		for i, code := range lh.Masks {
			fmt.Fprintf(outfh, "%d\t%s\n", i+1, decoder(code, _k))
		}

	},
}

func init() {
	utilsCmd.AddCommand(masksCmd)

	masksCmd.Flags().StringP("index", "d", "",
		formatFlagUsage(`Index directory created by "lexicprof index".`))

	masksCmd.Flags().StringP("out-file", "o", "-",
		formatFlagUsage(`Out file, supports and recommends a ".gz" suffix ("-" for stdout).`))

	masksCmd.Flags().IntP("kmer", "k", 31,
		formatFlagUsage(`Maximum k-mer size. K needs to be <= 32.`))

	masksCmd.Flags().IntP("masks", "m", 1000,
		formatFlagUsage(`Number of masks.`))

	masksCmd.Flags().IntP("seed", "s", 1,
		formatFlagUsage(`The seed for generating random masks.`))

	masksCmd.Flags().IntP("prefix", "p", 15,
		formatFlagUsage(`Length of mask k-mer prefix for checking low-complexity (0 for no checking).`))

	masksCmd.SetUsageTemplate(usageTemplate("{ -d <index path> | [-k <k>] [-n <masks>] [-s <seed>] } [-o out.tsv.gz]"))
}
