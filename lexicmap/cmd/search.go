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

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/spf13/cobra"
)

var mapCmd = &cobra.Command{
	Use:   "search",
	Short: "Search sequences against an index",
	Long: `Search sequences against an index

Attention:
  1. Input format should be (gzipped) FASTA or FASTQ from files or stdin.
  2. The positions in output are 1-based.

Output format:
  Tab-delimited format with 18 columns:

    1.  query,    Query ID.
    2.  qlen,     Query length.
    3.  qstart,   Query start position.
    4.  qend,     Query end position.
    5.  refs,     The number of matched reference genomes.
    6.  ref,      Reference genome ID.
    7.  seqid,    Target sequence ID.
    8.  qcovGnm,  Query coverage (percentage): $(aligned bases in the genome)/$qlen.
    9.  hit,      Nth hit in the genome.
    10. qcovHit   Query coverage (percentage): $(aligned bases in a hit)/$qlen.
    11. cmlen,    Matched bases in current hit, a hit might have >=1 segments.
    12. smlen,    Matched bases in current match/seqgment.
    13. pident,   Percentage of base identity.
    14. tlen,     Target sequence length.
    15. tstart,   Target start position.
    16. tend,     Target end position.
    17. str,      Strand of matched sequence.
    18. seeds,    Number of seeds.

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
			checkError(fmt.Errorf("the value of flag -m/--min-prefix (%d) should be <= 32", minPrefix))
		}
		maxMismatch := getFlagInt(cmd, "max-mismatch")
		minSinglePrefix := getFlagPositiveInt(cmd, "min-single-prefix")
		if minSinglePrefix > 32 {
			checkError(fmt.Errorf("the value of flag -M/--min-single-prefix (%d) should be <= 32", minSinglePrefix))
		}
		if minSinglePrefix < minPrefix {
			checkError(fmt.Errorf("the value of flag -M/--min-single-prefix (%d) should be >= that of -m/--min-prefix (%d)", minSinglePrefix, minPrefix))
		}
		maxGap := getFlagPositiveInt(cmd, "max-gap")
		maxDist := getFlagPositiveInt(cmd, "max-dist")
		extLen := getFlagNonNegativeInt(cmd, "ext-len")
		if extLen < 1000 {
			checkError(fmt.Errorf("the value of flag --ext-len should be >= 1000"))
		}
		topn := getFlagNonNegativeInt(cmd, "top-n")
		inMemorySearch := getFlagBool(cmd, "load-whole-seeds")

		minAlignLen := getFlagPositiveInt(cmd, "min-match-len")
		if minAlignLen < minSinglePrefix {
			checkError(fmt.Errorf("the value of flag -l/--min-match-len (%d) should be >= that of -M/--min-single-prefix (%d)", minAlignLen, minSinglePrefix))
		}

		minQcovGenome := getFlagNonNegativeFloat64(cmd, "min-qcov-in-genome")
		if minQcovGenome > 100 {
			checkError(fmt.Errorf("the value of flag -f/min-qcov-in-genome (%f) should be in range of [0, 100]", minQcovGenome))
		} else if minQcovGenome < 1 {
			log.Warningf("the value of flag -Q/min-qcov-in-genome is percentage in a range of [0, 100], you set: %f", minQcovGenome)
		}
		minIdent := getFlagNonNegativeFloat64(cmd, "min-match-identity")
		if minIdent > 100 {
			checkError(fmt.Errorf("the value of flag -i/min-match-identity (%f) should be in range of [0, 100]", minIdent))
		} else if minIdent < 1 {
			log.Warningf("the value of flag -i/min-match-identity is percentage in a range of [0, 100], you set: %f", minIdent)
		}
		minQcovChain := getFlagNonNegativeFloat64(cmd, "min-qcov-in-hit")
		if minQcovChain > 100 {
			checkError(fmt.Errorf("the value of flag -q/min-qcov-in-hit (%f) should be in range of [0, 100]", minIdent))
		}

		maxOpenFiles := getFlagPositiveInt(cmd, "max-open-files")

		// ---------------------------------------------------------------

		if outputLog {
			log.Infof("LexicProf v%s", VERSION)
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

		sopt := &IndexSearchingOptions{
			NumCPUs:      opt.NumCPUs,
			Verbose:      opt.Verbose,
			Log2File:     opt.Log2File,
			MaxOpenFiles: maxOpenFiles,

			MinPrefix:       uint8(minPrefix),
			MaxMismatch:     maxMismatch,
			MinSinglePrefix: uint8(minSinglePrefix),
			TopN:            topn,
			InMemorySearch:  inMemorySearch,

			MaxGap:      float64(maxGap),
			MaxDistance: float64(maxDist),

			ExtendLength: extLen,

			MinQueryAlignedFractionInAGenome: minQcovGenome,
		}

		idx, err := NewIndexSearcher(dbDir, sopt)
		checkError(err)
		defer func() {
			checkError(idx.Close())
		}()

		if outputLog {
			log.Infof("index loaded in %s", time.Since(timeStart))
			log.Info()
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

		fmt.Fprintf(outfh, "query\tqlen\tqstart\tqend\trefs\tref\tseqid\tqcovGnm\thit\tqcovHit\tcmlen\tsmlen\tpident\ttlen\ttstart\ttend\tstr\tseeds\n")

		results := make([]*SearchResult, 0, topn)
		printResult := func(q *Query) {
			total++
			if q.result == nil { // seqs shorter than K or queries without matches.
				poolQuery.Put(q)
				return
			}

			if verbose {
				if (total < 128 && total&7 == 0) || total&127 == 0 {
					speed = float64(total) / time.Since(timeStart1).Minutes()
					fmt.Fprintf(os.Stderr, "processed queries: %d, speed: %.3f queries per minute\r", total, speed)
				}
			}

			queryID := q.seqID
			// var c int
			// var v *index.SubstrPair
			// var i int
			// var subs *[]*index.SubstrPair
			var sd *SimilarityDetail
			var cr *SeqComparatorResult
			var c *Chain2Result
			var targets int

			results = results[:0]
			for _, r := range *q.result {
				if r.SimilarityDetails == nil {
					continue
				}

				if len(*r.SimilarityDetails) > 0 {
					results = append(results, r)
				}
			}
			sort.Slice(results, func(i, j int) bool {
				return (*results[i].SimilarityDetails)[0].SimilarityScore > (*results[j].SimilarityDetails)[0].SimilarityScore
			})
			targets = len(results)

			if targets > 0 {
				matched++
			}

			var strand byte
			var j int
			for _, r := range results { // each genome
				if r.SimilarityDetails == nil {
					continue
				}

				j = 1
				for _, sd = range *r.SimilarityDetails { // each chain
					cr = sd.Similarity

					for _, c = range *cr.Chains { // each match
						if c.TBegin < 0 || c.TEnd < 0 { // the extend part belongs to another contig
							continue
						}

						if sd.RC {
							strand = '-'
						} else {
							strand = '+'
						}
						fmt.Fprintf(outfh, "%s\t%d\t%d\t%d\t%d\t%s\t%s\t%.3f\t%d\t%.3f\t%d\t%d\t%.3f\t%d\t%d\t%d\t%c\t%d\n",
							queryID, len(q.seq),
							c.QBegin+1, c.QEnd+1,
							targets, r.ID,
							sd.SeqID, r.AlignedFraction, j, cr.AlignedFraction, cr.AlignedBases, c.AlignedBases, cr.PIdentity,
							sd.SeqLen,
							c.TBegin+1, c.TEnd+1, strand,
							sd.NSeeds,
						)
					}
					j++
				}
				outfh.Flush()
			}
			idx.RecycleSearchResults(q.result)

			poolQuery.Put(q)
			outfh.Flush()
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
		K := idx.k

		idx.SetSeqCompareOptions(&SeqComparatorOptions{
			K:         uint8(K),
			MinPrefix: 11, // can not be too small, or there will be a large number of anchors.

			Chaining2Options: Chaining2Options{
				// should be relative small
				MaxGap: 32,
				// better be larger than MinPrefix
				MinScore: minSinglePrefix,
				// can not be < k
				MaxDistance: 50,
				// can not be two small
				Band: 20,
			},

			MinIdentity:        minIdent,
			MinSegmentLength:   minAlignLen,
			MinAlignedFraction: minQcovChain,
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

			speed = float64(total) / time.Since(timeStart1).Minutes()
			log.Infof("")
			log.Infof("processed queries: %d, speed: %.3f queries per minute\n", total, speed)
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

	mapCmd.Flags().IntP("max-open-files", "", 512,
		formatFlagUsage(`Maximum opened files.`))

	// seed searching

	mapCmd.Flags().IntP("min-prefix", "p", 15,
		formatFlagUsage(`Minimum length of shared substrings (seeds).`))

	mapCmd.Flags().IntP("min-single-prefix", "P", 20,
		formatFlagUsage(`Minimum length of shared substrings if there's only one pair.`))

	mapCmd.Flags().IntP("max-mismatch", "m", -1,
		formatFlagUsage(`Minimum mismatch between non-prefix regions of shared substrings.`))

	mapCmd.Flags().IntP("max-gap", "g", 2000,
		formatFlagUsage(`Max gap in seed chaining.`))
	mapCmd.Flags().IntP("max-dist", "", 10000,
		formatFlagUsage(`Max distance in seed chaining.`))
	mapCmd.Flags().IntP("ext-len", "", 2000,
		formatFlagUsage(`Extend length of upstream and downstream of seed region, for extracting query and target sequences for alignment`))

	mapCmd.Flags().IntP("top-n", "n", 500,
		formatFlagUsage(`Keep top N matches for a query.`))

	mapCmd.Flags().BoolP("load-whole-seeds", "w", false,
		formatFlagUsage(`Load the whole seed data into memory for faster search.`))

	// sequence similarity

	mapCmd.Flags().IntP("min-match-len", "l", 50,
		formatFlagUsage(`Minimum matched sequence segment length`))
	mapCmd.Flags().Float64P("min-match-identity", "i", 70,
		formatFlagUsage(`Minimum base identity (percentage) between query and matched sequence segment.`))

	mapCmd.Flags().Float64P("min-qcov-in-hit", "q", 0,
		formatFlagUsage(`Minimum query coverage (percentage) in a hit with >=1 segments/matches).`))

	mapCmd.Flags().Float64P("min-qcov-in-genome", "Q", 50,
		formatFlagUsage(`Minimum query coverage (percentage) in a genome.`))

	mapCmd.SetUsageTemplate(usageTemplate("-d <index path> [query.fasta.gz ...] [-o query.tsv.gz]"))

}

// Strands could be used to output strand for a reverse complement flag
var Strands = [2]byte{'+', '-'}

// Query is an object for each query sequence, it also contains the query result.
type Query struct {
	seqID  []byte
	seq    []byte
	result *[]*SearchResult
}

// Reset reset the data for next round of using
func (q *Query) Reset() {
	q.seqID = q.seqID[:0]
	q.seq = q.seq[:0]
	q.result = nil
}

var poolQuery = &sync.Pool{New: func() interface{} {
	return &Query{
		seqID: make([]byte, 0, 128),     // the id should be not too long
		seq:   make([]byte, 0, 100<<10), // initialize with 100K
	}
}}
