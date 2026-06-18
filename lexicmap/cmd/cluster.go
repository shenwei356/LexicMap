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

var clusterCmd = &cobra.Command{
	Use:   "cluster",
	Short: "Cluster genomes in the index",
	Long: `Cluster genomes in the index


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

		// info file
		fileInfo := filepath.Join(dbDir, FileInfo)
		info, err := readIndexInfo(fileInfo)
		if err != nil {
			checkError(fmt.Errorf("failed to read info file: %s", err))
		}

		if outputLog {
			log.Infof("  checking passed")
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

		// ---------------------------------------------------------------

		// process bar
		var pbs *mpb.Progress
		var bar *mpb.Bar
		var chDuration chan time.Duration
		var doneDuration chan int
		var showProgressBar bool

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

		// ----------

		var chunkSize, chunk, iMask int
		var fileSeeds string

		var startTime time.Time
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
		_ = batchIDAndRefID

		// compute the chunk
		chunkSize = (len(lh.Masks) + info.Chunks - 1) / info.Chunks

		for mask, maskCode := range lh.Masks {
			_ = maskCode

			startTime = time.Now()

			chunk = mask / chunkSize
			iMask = mask % chunkSize

			fileSeeds = filepath.Join(dbDir, DirSeeds, chunkFile(chunk))

			// kv-data index file
			_, _, indexes, _, _, config1, err := kv.ReadKVIndex(filepath.Clean(fileSeeds) + kv.KVIndexFileExt)
			if err != nil {
				checkError(fmt.Errorf("failed to read kv-data index file: %s", err))
			}

			use3BytesForSeedPos := config1&kv.MaskUse3BytesForSeedPos > 0
			if !use3BytesForSeedPos {
				checkError(fmt.Errorf("index with genome batch number > 512 is not supported"))
			}

			bytesPos := 8
			fUint64 := be.Uint64
			if use3BytesForSeedPos {
				bytesPos = 7
				fUint64 = kv.Uint64ThreeBytes
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

			_, err = r.Seek(int64(indexes[iMask][1])>>1, 0)
			if err != nil {
				checkError(fmt.Errorf("failed to seed kv-data file: %s", err))
			}

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
					nReaded, err = io.ReadFull(r, buf8[:bytesPos])
					if err != nil {
						checkError(err)
					}
					if nReaded < bytesPos {
						checkError(kv.ErrBrokenFile)
					}

					v = fUint64(buf8)
					if v&MASK_REVERSE == 1 {
						continue
					}
					// for genome batch number <= 512, only 9 bits is used to encode batch ID.
					// so 9 + 17 = 26, which can be stored with an uint32 number.
					batchIDAndRefID = (v >> BITS_NONE_IDX) & 4294967295

				}

				if lastPair && !hasKmer2 {
					break
				}

				for j = 0; j < lenVal2; j++ {
					nReaded, err = io.ReadFull(r, buf8[:bytesPos])
					if err != nil {
						checkError(err)
					}
					if nReaded < bytesPos {
						checkError(kv.ErrBrokenFile)
					}

					v = fUint64(buf8)
					if v&MASK_REVERSE == 1 {
						continue
					}
					batchIDAndRefID = v >> BITS_NONE_IDX

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
	RootCmd.AddCommand(clusterCmd)

	clusterCmd.Flags().StringP("index", "d", "",
		formatFlagUsage(`Index directory created by "lexicmap index".`))

	clusterCmd.Flags().StringP("out-file", "o", "-",
		formatFlagUsage(`Out file, supports and recommends a ".gz" suffix ("-" for stdout).`))

	clusterCmd.SetUsageTemplate(usageTemplate("-d <index path> [-o out.tsv.gz]"))

	clusterCmd.Flags().IntP("min-prefix", "p", 21,
		formatFlagUsage(`Minimum prefix length between k-mers captured by a mask.`))
}
