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
	"sort"
	"strings"
	"sync"
	"time"

	"github.com/cznic/sortutil"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/lexichash"
	"github.com/shenwei356/lexichash/index"
	"github.com/spf13/cobra"
)

var mapCmd = &cobra.Command{
	Use:   "map",
	Short: "map sequences against an index",
	Long: `map sequences against an index

Attentions:
  1. Input format should be (gzipped) FASTA or FASTQ from files or stdin.
  2. The positions are 1-based.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		seq.ValidateSeq = false

		var fhLog *os.File
		if opt.Log2File {
			fhLog = addLog(opt.LogFile, opt.Verbose)
		}

		verbose := opt.Verbose
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
		minSubLen := getFlagPositiveInt(cmd, "min-subs")
		if minSubLen > 32 {
			checkError(fmt.Errorf("the value of flag -m/--min-subs should be <= 32"))
		}

		// ---------------------------------------------------------------

		if outputLog {
			log.Infof("LexicMap v%s", VERSION)
			log.Info("  https://github.com/shenwei356/LexicMap")
			log.Info()
		}

		// ---------------------------------------------------------------
		// input files

		if outputLog {
			log.Info("checking input files ...")
		}

		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)

		if outputLog {
			if len(files) == 1 && isStdin(files[0]) {
				log.Info("  no files given, reading from stdin")
			} else {
				log.Infof("  %d input file(s) given", len(files))
			}
		}

		outFileClean := filepath.Clean(outFile)
		for _, file := range files {
			if !isStdin(file) && filepath.Clean(file) == outFileClean {
				checkError(fmt.Errorf("out file should not be one of the input file"))
			}
		}

		// ---------------------------------------------------------------
		// loading index

		if outputLog {
			log.Info()
			log.Infof("loading index: %s", dbDir)
		}

		idx, err := index.NewFromPath(dbDir, opt.NumCPUs)
		checkError(err)

		if outputLog {
			log.Infof("index loaded in %s", time.Since(timeStart))
			log.Info()
		}

		if minSubLen > idx.K() {
			checkError(fmt.Errorf("the value of flag -m/--min-subs (%d) should be <= K (%d)", minSubLen, idx.K()))
		}

		if outputLog {
			log.Info("searching ...")
		}

		// ---------------------------------------------------------------
		// mapping

		timeStart1 := time.Now()

		outfh, gw, w, err := outStream(outFile, strings.HasSuffix(outFile, ".gz"), opt.CompressionLevel)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()

		var total, matched uint64
		var speed float64 // k reads/second

		fmt.Fprintf(outfh, "query\ttarget\tqstart\tqend\tqstrand\ttstart\ttend\ttstrand\tlen\tmatch\n")

		type Result struct {
			id      uint64
			queryID []byte
			result  *[]*index.SearchResult
		}

		decoder := lexichash.MustDecoder()

		printResult := func(queryID []byte, sr *[]*index.SearchResult) {
			total++
			if verbose {
				if (total < 4096 && total&63 == 0) || total&4095 == 0 {
					speed = float64(total) / 1000000 / time.Since(timeStart1).Minutes()
					fmt.Fprintf(os.Stderr, "processed queries: %d, speed: %.3f million queries per minute\r", total, speed)
				}
			}

			if sr == nil {
				return
			}
			matched++

			for _, r := range *sr {
				for _, v := range *r.Subs {
					fmt.Fprintf(outfh, "%s\t%s\t%d\t%d\t%c\t%d\t%d\t%c\t%d\t%s\n",
						queryID, idx.IDs[r.IdIdx],
						v.QBegin+1, v.QEnd, Strands[v.QRC],
						v.TBegin+1, v.TEnd, Strands[v.TRC],
						v.QK, decoder(v.QCode, v.QK))
				}
			}
			idx.RecycleSearchResult(sr)
		}

		// outputter
		ch := make(chan Result, opt.NumCPUs)
		done := make(chan int)
		go func() {
			var id uint64 = 1 // for keepping order
			buf := make(map[uint64]Result, 128)

			var r, r2 Result
			var ok bool

			for r := range ch {
				if id == r.id {
					printResult(r.queryID, r.result)
					id++
					continue
				}
				buf[r.id] = r

				if r2, ok = buf[id]; ok {
					printResult(r2.queryID, r2.result)
					delete(buf, r2.id)
					id++
				}
			}
			if len(buf) > 0 {
				ids := make(sortutil.Uint64Slice, len(buf))
				i := 0
				for id := range buf {
					ids[i] = id
					i++
				}
				sort.Sort(ids)

				for _, id := range ids {
					r = buf[id]
					printResult(r.queryID, r.result)
				}

			}
			done <- 1
		}()

		var wg sync.WaitGroup
		token := make(chan int, opt.NumCPUs)
		var id uint64

		var record *fastx.Record
		var fastxReader *fastx.Reader

		for _, file := range files {
			fastxReader, err = fastx.NewReader(nil, file, "")
			checkError(err)
			for {
				record, err = fastxReader.Read()
				if err != nil {
					if err == io.EOF {
						break
					}
					checkError(err)
					break
				}

				token <- 1
				wg.Add(1)
				id++
				go func(id uint64, record *fastx.Record) {
					defer func() {
						<-token
						wg.Done()
					}()

					sr, err := idx.Search(record.Seq.Seq, uint8(minSubLen))
					if err != nil {
						checkError(err)
					}

					ch <- Result{id: id, queryID: record.ID, result: sr}
				}(id, record.Clone())
			}
		}
		wg.Wait()
		close(ch)
		<-done

		if outputLog {
			fmt.Fprintf(os.Stderr, "\n")

			speed = float64(total) / 1000000 / time.Since(timeStart1).Minutes()
			log.Infof("")
			log.Infof("processed queries: %d, speed: %.3f million queries per minute\n", total, speed)
			log.Infof("%.4f%% (%d/%d) queries matched", float64(matched)/float64(total)*100, matched, total)
			log.Infof("done searching")
			if outFile != "-" {
				log.Infof("search results saved to: %s", outFile)
			}

		}

	},
}

func init() {
	RootCmd.AddCommand(mapCmd)

	mapCmd.Flags().StringP("index", "d", "",
		formatFlagUsage(`Index directory created by "lexicmap index".`))

	mapCmd.Flags().StringP("out-file", "o", "-",
		formatFlagUsage(`Out file, supports and recommends a ".gz" suffix ("-" for stdout).`))

	mapCmd.Flags().IntP("min-subs", "m", 15,
		formatFlagUsage(`Minimum length of shared substrings`))

	mapCmd.SetUsageTemplate(usageTemplate("-d <index path> [read.fq.gz ...] [-o read.tsv.gz]"))
}

// Strands could be used to output strand for a reverse complement flag
var Strands = [2]byte{'+', '-'}
