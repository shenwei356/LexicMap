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
	"sync"
	"time"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/lexichash/index"
	"github.com/shenwei356/lexichash/index/align"
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
		minPrefix := getFlagPositiveInt(cmd, "min-prefix")
		if minPrefix > 32 {
			checkError(fmt.Errorf("the value of flag -m/--min-prefix should be <= 32"))
		}
		minSinglePrefix := getFlagPositiveInt(cmd, "min-single-prefix")
		if minSinglePrefix > 32 {
			checkError(fmt.Errorf("the value of flag -M/--min-single-prefix should be <= 32"))
		}
		if minSinglePrefix < minPrefix {
			checkError(fmt.Errorf("the value of flag -M/--min-single-prefix should be >= that of -m/--min-prefix "))
		}
		maxGap := getFlagNonNegativeInt(cmd, "max-gap")
		topn := getFlagNonNegativeInt(cmd, "top-n")

		index.Threads = opt.NumCPUs

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
		defer func() {
			checkError(idx.Close())
		}()

		if outputLog {
			log.Infof("index loaded in %s", time.Since(timeStart))
			log.Info()
		}

		if minPrefix > idx.K() {
			checkError(fmt.Errorf("the value of flag -m/--min-prefix (%d) should be <= K (%d)", minPrefix, idx.K()))
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

		fmt.Fprintf(outfh, "query\tqlen\ttargets\ttarget\tchain\tident\ttlen\tqstart\tqend\ttstart\ttend\tlen\n")

		printResult := func(q *Query) {
			total++
			if q.result == nil { // seqs shorter than K or queries without matches.
				poolQuery.Put(q)
				return
			}

			if verbose {
				if (total < 4096 && total&63 == 0) || total&4095 == 0 {
					speed = float64(total) / 1000000 / time.Since(timeStart1).Minutes()
					fmt.Fprintf(os.Stderr, "processed queries: %d, speed: %.3f million queries per minute\r", total, speed)
				}
			}

			matched++

			queryID := q.seqID
			targets := len(*q.result)
			var chain *[]int
			var v *index.SubstrPair
			var c, i int
			var subs *[]*index.SubstrPair
			var aligns *[]*align.AlignResult
			var ar *align.AlignResult
			var ident float64
			for _, r := range *q.result {
				subs = r.Subs
				aligns = r.AlignResults
				// fmt.Printf("%p: %d, %p: %d\n", r.Chains, len(*r.Chains), r.AlignResults, len(*r.AlignResults))
				for c, chain = range *r.Chains {
					ar = (*aligns)[c]
					ident = float64(ar.Matches) / float64(ar.Len) * 100
					for _, i = range *chain {
						v = (*subs)[i]
						fmt.Fprintf(outfh, "%s\t%d\t%d\t%s\t%d\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\n",
							queryID, len(q.seq), targets, idx.IDs[r.IdIdx],
							c+1, ident, idx.RefSeqInfos[r.IdIdx].Len,
							v.QBegin+1, v.QBegin+v.Len,
							v.TBegin+1, v.TBegin+v.Len,
							v.Len)
					}
				}
				outfh.Flush()

			}
			idx.RecycleSearchResults(q.result)

			poolQuery.Put(q)
		}

		// outputter
		ch := make(chan *Query, opt.NumCPUs)
		done := make(chan int)
		go func() {

			for r := range ch {
				printResult(r)
			}

			done <- 1
		}()

		var wg sync.WaitGroup
		tokens := make(chan int, opt.NumCPUs)

		var record *fastx.Record
		K := idx.K()

		idx.SetSearchingOptions(&index.SearchOptions{
			MinPrefix:       uint8(minPrefix),
			MinSinglePrefix: uint8(minSinglePrefix),
			TopN:            topn,

			MaxGap: float64(maxGap),
		})
		idx.SetAlignOptions(&align.AlignOptions{
			MatchScore:    1,
			MisMatchScore: -1,
			GapScore:      -1,
		})

		for _, file := range files {
			fastxReader, err := fastx.NewReader(nil, file, "")
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

				query := poolQuery.Get().(*Query)
				query.Reset()

				if len(record.Seq.Seq) < K {
					query.result = nil
					ch <- query
					return
				}

				tokens <- 1
				wg.Add(1)

				query.seqID = append(query.seqID, record.ID...)
				query.seq = append(query.seq, record.Seq.Seq...)

				go func(query *Query) {
					defer func() {
						<-tokens
						wg.Done()
					}()

					var err error
					query.result, err = idx.Search(query.seq)
					if err != nil {
						checkError(err)
					}

					ch <- query
				}(query)
			}
			fastxReader.Close()
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

	// searching

	mapCmd.Flags().IntP("min-prefix", "m", 15,
		formatFlagUsage(`Minimum length of shared substrings`))

	mapCmd.Flags().IntP("min-single-prefix", "M", 20,
		formatFlagUsage(`Minimum length of shared substrings if there's only one pair`))

	mapCmd.Flags().IntP("top-n", "n", 10,
		formatFlagUsage(`Keep top n matches for a query`))

	mapCmd.Flags().IntP("max-gap", "g", 5000,
		formatFlagUsage(`max gap`))

	mapCmd.SetUsageTemplate(usageTemplate("-d <index path> [read.fq.gz ...] [-o read.tsv.gz]"))
}

// Strands could be used to output strand for a reverse complement flag
var Strands = [2]byte{'+', '-'}

type Query struct {
	seqID  []byte
	seq    []byte
	result *[]*index.SearchResult
}

func (q *Query) Reset() {
	q.seqID = q.seqID[:0]
	q.seq = q.seq[:0]
	q.result = nil
}

var poolQuery = &sync.Pool{New: func() interface{} {
	return &Query{
		seqID: make([]byte, 0, 128),
		seq:   make([]byte, 0, 10<<10),
	}
}}
