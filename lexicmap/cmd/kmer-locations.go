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
	"github.com/shenwei356/lexichash/index"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
)

var kmerLocations = &cobra.Command{
	Use:   "kmer-locations",
	Short: "view locations of k-mers captured by the masks for each reference sequence",
	Long: `view locations of k-mers captured by the masks for each reference sequence

Attentions:
  1. To get the index of reference sequence, please check the file IDs.txt.
     The index of one sequence is simply the line number.
  2. K-mer positions (column pos) are 1-based.

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

		refIdx := getFlagNonNegativeInt(cmd, "ref-idx")
		outFile := getFlagString(cmd, "out-file")

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

		file := filepath.Join(dbDir, index.MaskLocationFile)
		ok, err := pathutil.Exists(file)
		checkError(err)

		var idx *index.Index
		var kmerLocations [][]uint64
		var ids [][]byte

		if !ok {
			if outputLog {
				log.Info()
				log.Infof("loading index: %s", dbDir)
			}

			idx, err = index.NewFromPath(dbDir, opt.NumCPUs)
			checkError(err)

			if outputLog {
				log.Infof("index loaded in %s", time.Since(timeStart))
				log.Info()
			}

			timeStart1 := time.Now()
			if outputLog {
				log.Infof("extracting k-mer locations ...")
			}
			idx.ExtractKmerLocations()
			if outputLog {
				log.Infof("finished extracting k-mer locations in %s", time.Since(timeStart1))
			}

			timeStart1 = time.Now()
			if outputLog {
				log.Infof("saving k-mer locations into the binary file...")
			}
			checkError(idx.WriteKmerLocations())
			if outputLog {
				log.Infof("finished saving k-mer locations binary file in %s", time.Since(timeStart1))
			}

			kmerLocations = idx.KmerLocations
			ids = idx.IDs
		} else {
			if outputLog {
				log.Infof("reading from the mask location binary file...")
			}
			kmerLocations, err = index.ReadKmerLocationsFromFile(file)
			checkError(err)

			ids, err = index.ReadIDlistFromFile(filepath.Join(dbDir, index.IDListFile))
			checkError(err)

			if outputLog {
				log.Infof("finished reading k-mer locations in %s", time.Since(timeStart))
			}
		}

		if refIdx >= len(kmerLocations) {
			log.Warningf("the value of -i/--ref-idx %d is larger than the number of reference sequences (%d)", refIdx, len(kmerLocations))
		}

		fmt.Fprintf(outfh, "ref\tpos\tstrand\n")
		var refpos uint64
		var pos uint64
		var rc uint8

		if refIdx == 0 {
			for i, locs := range kmerLocations {
				for _, refpos = range locs {
					pos = refpos >> 2
					rc = uint8(refpos & 1)
					fmt.Fprintf(outfh, "%s\t%d\t%c\n", ids[i], pos, Strands[rc])
				}
			}
		} else {
			i := refIdx - 1
			for _, refpos = range kmerLocations[i] {
				pos = refpos >> 2
				rc = uint8(refpos & 1)
				fmt.Fprintf(outfh, "%s\t%d\t%c\n", ids[i], pos, Strands[rc])
			}
		}

	},
}

func init() {
	utilsCmd.AddCommand(kmerLocations)

	kmerLocations.Flags().StringP("index", "d", "",
		formatFlagUsage(`Index directory created by "lexicmap index".`))

	kmerLocations.Flags().StringP("out-file", "o", "-",
		formatFlagUsage(`Out file, supports and recommends a ".gz" suffix ("-" for stdout).`))

	kmerLocations.Flags().IntP("ref-idx", "i", 1,
		formatFlagUsage(`View locations of k-mers for Xth reference sequence. (0 for all)`))

	kmerLocations.SetUsageTemplate(usageTemplate("-d <index path> [-i <refseq index>] [-o out.tsv.gz]"))
}
