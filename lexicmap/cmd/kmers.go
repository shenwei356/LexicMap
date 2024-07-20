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
	"io"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/lexichash"
	"github.com/spf13/cobra"
	"github.com/vbauerster/mpb/v8"
	"github.com/vbauerster/mpb/v8/decor"
)

var kmersCmd = &cobra.Command{
	Use:   "kmers",
	Short: "View k-mers captured by the masks",
	Long: `View k-mers captured by the masks

Attention:
  1. Mask index (column mask) is 1-based.
  2. Prefix means the length of shared prefix between a k-mer and the mask.
  3. K-mer positions (column pos) are 1-based.
     For reference genomes with multiple sequences, the sequences were
     concatenated to a single sequence with intervals of N's.

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
		fileMask := filepath.Join(dbDir, FileMasks)
		lh, err := lexichash.NewFromFile(fileMask)
		if err != nil {
			checkError(err)
		}

		if mask > len(lh.Masks) {
			checkError(fmt.Errorf("the index has only %d masks, but %d is given", len(lh.Masks), mask))
		}

		// info file
		fileInfo := filepath.Join(dbDir, FileInfo)
		info, err := readIndexInfo(fileInfo)
		if err != nil {
			checkError(fmt.Errorf("failed to read info file: %s", err))
		}

		if outputLog {
			log.Infof("  checking passed")
		}

		// genomes.map file for mapping index to genome id
		m, err := readGenomeMapIdx2Name(filepath.Join(dbDir, FileGenomeIndex))
		if err != nil {
			checkError(fmt.Errorf("failed to read genomes index mapping file: %s", err))
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

		fmt.Fprintf(outfh, "mask\tkmer\tprefix\tnumber\tref\tpos\tstrand\treversed\n")

		// ---------------------------------------------------------------

		// process bar
		var pbs *mpb.Progress
		var bar *mpb.Bar
		var chDuration chan time.Duration
		var doneDuration chan int
		var showProgressBar bool

		var masks []int
		if mask == 0 {
			masks = make([]int, len(lh.Masks))
			for i := range masks {
				masks[i] = i + 1
			}

			if opt.Verbose {
				showProgressBar = true

				pbs = mpb.New(mpb.WithWidth(40), mpb.WithOutput(os.Stderr))
				bar = pbs.AddBar(int64(len(lh.Masks)),
					mpb.PrependDecorators(
						decor.Name("processed masks: ", decor.WC{W: len("processed masks: "), C: decor.DindentRight}),
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
						bar.EwmaIncrBy(1, t)
					}
					doneDuration <- 1
				}()
			}
		} else {
			masks = []int{mask}
		}

		var chunkSize, chunk, iMask int
		var fileSeeds string

		var startTime time.Time
		k8 := uint8(lh.K)
		buf := make([]byte, 64)
		buf8 := make([]uint8, 8)
		var ctrlByte byte
		var first bool     // the first kmer has a different way to comput the value
		var lastPair bool  // check if this is the last pair
		var hasKmer2 bool  // check if there's a kmer2
		var _offset uint64 // offset of kmer
		var nBytes int
		var nReaded, nDecoded int
		var v1, v2 uint64
		var kmer1, kmer2 uint64
		var lenVal1, lenVal2 uint64
		var j uint64
		var v, batchIDAndRefID uint64
		var pos, rc int
		var maskCode uint64
		var rvFlag int

		// compute the chunk
		chunkSize = (len(lh.Masks) + info.Chunks - 1) / info.Chunks

		for _, mask = range masks {
			startTime = time.Now()

			chunk = (mask - 1) / chunkSize
			iMask = (mask - 1) % chunkSize

			fileSeeds = filepath.Join(dbDir, DirSeeds, chunkFile(chunk))

			// kv-data index file
			k, _, indexes, err := kv.ReadKVIndex(filepath.Clean(fileSeeds) + kv.KVIndexFileExt)
			if err != nil {
				checkError(fmt.Errorf("failed to read kv-data index file: %s", err))
			}

			if len(indexes[iMask]) == 0 { // no k-mers
				if showProgressBar {
					chDuration <- time.Duration(float64(time.Since(startTime)))
				}
				continue
			}

			// kv-data file
			r, err := os.Open(fileSeeds)
			if err != nil {
				checkError(fmt.Errorf("failed to read kv-data file: %s", err))
			}

			_, err = r.Seek(int64(indexes[iMask][1]), 0)
			if err != nil {
				checkError(fmt.Errorf("failed to seed kv-data file: %s", err))
			}

			maskCode = lh.Masks[mask-1]

			_offset = 0
			for {
				// read the control byte
				_, err = io.ReadFull(r, buf[:1])
				if err != nil {
					checkError(err)
				}
				ctrlByte = buf[0]

				lastPair = ctrlByte&128 > 0 // 1<<7
				hasKmer2 = ctrlByte&64 == 0 // 1<<6

				ctrlByte &= 63

				// parse the control byte
				nBytes = util.CtrlByte2ByteLengthsUint64(ctrlByte)

				// read encoded bytes
				nReaded, err = io.ReadFull(r, buf[:nBytes])
				if err != nil {
					checkError(err)
				}
				if nReaded < nBytes {
					checkError(kv.ErrBrokenFile)
				}

				v1, v2, nDecoded = util.Uint64s(ctrlByte, buf[:nBytes])
				if nDecoded == 0 {
					checkError(kv.ErrBrokenFile)
				}

				if first {
					kmer1 = indexes[iMask][0] // from the index
					first = false
				} else {
					kmer1 = v1 + _offset
				}
				kmer2 = kmer1 + v2
				_offset = kmer2

				// ------------------ lengths of values -------------------

				// read the control byte
				_, err = io.ReadFull(r, buf[:1])
				if err != nil {
					checkError(err)
				}
				ctrlByte = buf[0]

				// parse the control byte
				nBytes = util.CtrlByte2ByteLengthsUint64(ctrlByte)

				// read encoded bytes
				nReaded, err = io.ReadFull(r, buf[:nBytes])
				if err != nil {
					checkError(err)
				}
				if nReaded < nBytes {
					checkError(kv.ErrBrokenFile)
				}

				lenVal1, lenVal2, nDecoded = util.Uint64s(ctrlByte, buf[:nBytes])
				if nDecoded == 0 {
					checkError(kv.ErrBrokenFile)
				}

				// ------------------ values -------------------

				for j = 0; j < lenVal1; j++ {
					nReaded, err = io.ReadFull(r, buf8)
					if err != nil {
						checkError(err)
					}
					if nReaded < 8 {
						checkError(kv.ErrBrokenFile)
					}

					v = be.Uint64(buf8)
					// pos, rc = int(v<<34>>35), int(v&1)
					// batchIDAndRefID = v >> 30
					pos, rc, rvFlag = int(v<<BITS_IDX>>BITS_NONE_POS), int(v>>BITS_REVERSE&BITS_STRAND), int(v&BITS_REVERSE)
					batchIDAndRefID = v >> BITS_NONE_IDX
					fmt.Fprintf(outfh, "%d\t%s\t%d\t%d\t%s\t%d\t%c\t%s\n",
						mask, decoder(kmer1, k), util.MustKmerLongestPrefix(kmer1, maskCode, k8, k8),
						lenVal1, m[batchIDAndRefID], pos+1, lexichash.Strands[rc], reversedStr[rvFlag])
				}

				if lastPair && !hasKmer2 {
					break
				}

				for j = 0; j < lenVal2; j++ {
					nReaded, err = io.ReadFull(r, buf8)
					if err != nil {
						checkError(err)
					}
					if nReaded < 8 {
						checkError(kv.ErrBrokenFile)
					}

					v = be.Uint64(buf8)
					// pos, rc = int(v<<34>>35), int(v&1)
					// batchIDAndRefID = v >> 30
					pos, rc, rvFlag = int(v<<BITS_IDX>>BITS_NONE_POS), int(v>>BITS_REVERSE&BITS_STRAND), int(v&BITS_REVERSE)
					batchIDAndRefID = v >> BITS_NONE_IDX
					fmt.Fprintf(outfh, "%d\t%s\t%d\t%d\t%s\t%d\t%c\t%s\n",
						mask, decoder(kmer2, k), util.MustKmerLongestPrefix(kmer2, maskCode, k8, k8),
						lenVal2, m[batchIDAndRefID], pos+1, lexichash.Strands[rc], reversedStr[rvFlag])
				}

				if lastPair {
					break
				}

			}

			r.Close()

			if showProgressBar {
				chDuration <- time.Duration(float64(time.Since(startTime)))
			}
		}

		if showProgressBar {
			close(chDuration)
			<-doneDuration
			pbs.Wait()
		}

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

	kmersCmd.SetUsageTemplate(usageTemplate("-d <index path> [-m <mask index>] [-o out.tsv.gz]"))
}

var reversedStr = []string{"no", "yes"}
