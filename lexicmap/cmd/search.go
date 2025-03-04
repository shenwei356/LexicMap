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
	"bytes"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"runtime"
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
  1. Input should be (gzipped) FASTA or FASTQ records from files or stdin.
  2. For multiple queries, the order of queries might be different from the input.

Tips:
  1. When using -a/--all, the search result would be formatted to Blast-style format
     with 'lexicmap utils 2blast'. And the search speed would be slightly slowed down.
  2. Alignment result filtering is performed in the final phase, so stricter filtering criteria,
     including -q/--min-qcov-per-hsp, -Q/--min-qcov-per-genome, and -i/--align-min-match-pident,
     do not significantly accelerate the search speed. Hence, you can search with default
     parameters and then filter the result with tools like awk or csvtk.

Alignment result relationship:

  Query
  ├── Subject genome
      ├── Subject sequence
          ├── High-Scoring segment Pair (HSP)

  Here, the defination of HSP is similar with that in BLAST. Actually there are small gaps in HSPs.

  > A High-scoring Segment Pair (HSP) is a local alignment with no gaps that achieves one of the
  > highest alignment scores in a given search. https://www.ncbi.nlm.nih.gov/books/NBK62051/

Output format:
  Tab-delimited format with 19+ columns, with 1-based positions.

    1.  query,    Query sequence ID.
    2.  qlen,     Query sequence length.
    3.  hits,     Number of subject genomes.
    4.  sgenome,  Subject genome ID.
    5.  sseqid,   Subject sequence ID.
    6.  qcovGnm,  Query coverage (percentage) per genome: $(aligned bases in the genome)/$qlen.
    7.  hsp,      Nth HSP in the genome. (just for improving readability)
    8.  qcovHSP   Query coverage (percentage) per HSP: $(aligned bases in a HSP)/$qlen.
    9.  alenHSP,  Aligned length in the current HSP.
    10. pident,   Percentage of identical matches in the current HSP.
    11. gaps,     Gaps in the current HSP.
    12. qstart,   Start of alignment in query sequence.
    13. qend,     End of alignment in query sequence.
    14. sstart,   Start of alignment in subject sequence.
    15. send,     End of alignment in subject sequence.
    16. sstr,     Subject strand.
    17. slen,     Subject sequence length.
    18. evalue,   E value.
    19. bitscore, bit score.
    20. cigar,    CIGAR string of the alignment.                      (optional with -a/--all)
    21. qseq,     Aligned part of query sequence.                     (optional with -a/--all)
    22. sseq,     Aligned part of subject sequence.                   (optional with -a/--all)
    23. align,    Alignment text ("|" and " ") between qseq and sseq. (optional with -a/--all)

Result ordering:
  1. Within each subject genome, alignments (HSP) are sorted by qcovHSP*pident.
  2. Results of multiple subject genomes are sorted by qcovHSP*pident of the best alignment.

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
		minPrefix := getFlagPositiveInt(cmd, "seed-min-prefix")
		if minPrefix > 32 || minPrefix < 5 {
			checkError(fmt.Errorf("the value of flag -p/--seed-min-prefix (%d) should be in the range of [5, 32]", minPrefix))
		}
		moreColumns := getFlagBool(cmd, "all")

		// maxMismatch := getFlagInt(cmd, "seed-max-mismatch")
		minSinglePrefix := getFlagPositiveInt(cmd, "seed-min-single-prefix")
		if minSinglePrefix > 32 {
			checkError(fmt.Errorf("the value of flag -P/--seed-min-single-prefix (%d) should be <= 32", minSinglePrefix))
		}
		if minSinglePrefix < minPrefix {
			checkError(fmt.Errorf("the value of flag -P/--seed-min-single-prefix (%d) should be >= that of -p/--seed-min-prefix (%d)", minSinglePrefix, minPrefix))
		}

		// minMatches := getFlagPositiveInt(cmd, "seed-min-matches")
		// if minMatches > 32 {
		// 	checkError(fmt.Errorf("the value of flag -m/--seed-min-matches (%d) should be <= 32", minMatches))
		// }
		// if minMatches < minPrefix {
		// 	checkError(fmt.Errorf("the value of flag -m/--seed-min-matches (%d) should be >= that of -P/--seed-min-single-prefix (%d)", minMatches, minSinglePrefix))
		// }

		maxGap := getFlagPositiveInt(cmd, "seed-max-gap")
		maxDist := getFlagPositiveInt(cmd, "seed-max-dist")
		extLen := getFlagNonNegativeInt(cmd, "align-ext-len")
		// if extLen < 1000 {
		// 	checkError(fmt.Errorf("the value of flag --align-ext-len should be >= 1000"))
		// }
		topn := getFlagNonNegativeInt(cmd, "top-n-genomes")
		inMemorySearch := getFlagBool(cmd, "load-whole-seeds")

		onlyPseudoAlign := getFlagBool(cmd, "pseudo-align")

		minAlignLen := getFlagPositiveInt(cmd, "align-min-match-len")
		if minAlignLen < minSinglePrefix {
			checkError(fmt.Errorf("the value of flag -l/--align-min-match-len (%d) should be >= that of -M/--seed-min-single-prefix (%d)", minAlignLen, minSinglePrefix))
		}
		maxAlignMaxGap := getFlagPositiveInt(cmd, "align-max-gap")
		// maxAlignMismatch := getFlagPositiveInt(cmd, "align-max-kmer-dist")
		alignBand := getFlagPositiveInt(cmd, "align-band")
		if alignBand < maxAlignMaxGap {
			checkError(fmt.Errorf("the value of flag --align-band should not be smaller thant the value of --align-max-gap"))
		}

		minQcovGenome := getFlagNonNegativeFloat64(cmd, "min-qcov-per-genome")
		if minQcovGenome > 100 {
			checkError(fmt.Errorf("the value of flag -Q/--min-qcov-per-genome (%f) should be in range of [0, 100]", minQcovGenome))
		}
		// } else if minQcovGenome < 1 {
		// 	log.Warningf("the value of flag -Q/--min-qcov-per-genome is percentage in a range of [0, 100], you set: %f", minQcovGenome)
		// }
		minIdent := getFlagNonNegativeFloat64(cmd, "align-min-match-pident")
		if minIdent < 60 || minIdent > 100 {
			checkError(fmt.Errorf("the value of flag -i/--align-min-match-pident (%f) should be in range of [60, 100]", minIdent))
		}

		// } else if minIdent < 1 {
		// 	log.Warningf("the value of flag -i/--align-min-match-pident is percentage in a range of [0, 100], you set: %f", minIdent)
		// }
		minQcovChain := getFlagNonNegativeFloat64(cmd, "min-qcov-per-hsp")
		if minQcovChain > 100 {
			checkError(fmt.Errorf("the value of flag -q/--min-qcov-per-hsp (%f) should be in range of [0, 100]", minIdent))
		}

		maxOpenFiles := getFlagPositiveInt(cmd, "max-open-files")

		// ---------------------------------------------------------------

		if outputLog {
			log.Infof("LexicMap v%s (%s)", VERSION, COMMIT)
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
			if len(files) == 1 {
				if isStdin(files[0]) {
					log.Info("  no files given, reading from stdin")
				} else {
					log.Infof("  %d input file given: %s", len(files), files[0])
				}
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

		maxQueryConcurrency := getFlagNonNegativeInt(cmd, "max-query-conc")
		if maxQueryConcurrency == 0 {
			maxQueryConcurrency = runtime.NumCPU()
		}

		// maxSeedingConcurrency := getFlagNonNegativeInt(cmd, "max-seed-conc")
		// if maxSeedingConcurrency == 0 {
		// 	maxSeedingConcurrency = runtime.NumCPU()
		// }

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

			// MaxSeedingConcurrency: maxSeedingConcurrency,

			MinPrefix: uint8(minPrefix),
			// MaxMismatch:     maxMismatch,
			MinSinglePrefix: uint8(minSinglePrefix),
			// MinMatchedBases: uint8(minMatches),
			TopN:           topn,
			InMemorySearch: inMemorySearch,

			MaxGap:      float64(maxGap),
			MaxDistance: float64(maxDist),

			ExtendLength:  extLen,
			ExtendLength2: 50,

			MinQueryAlignedFractionInAGenome: minQcovGenome,

			MoreAccurateAlignment: !onlyPseudoAlign,

			OutputSeq: moreColumns,

			Debug: getFlagBool(cmd, "debug"),
		}

		// read info file to get the contig interval size
		// fileInfo := filepath.Join(dbDir, FileInfo)
		// info, err := readIndexInfo(fileInfo)
		// if err != nil {
		// 	checkError(fmt.Errorf("failed to read index info file: %s", err))
		// }

		// if extLen > info.ContigInterval {
		// 	log.Infof("the value of flag --align-ext-len (%d) is adjusted to contig interval length in database (%d)", extLen, info.ContigInterval)
		// 	sopt.ExtendLength = info.ContigInterval
		// }
		// if maxDist > info.ContigInterval {
		// 	log.Infof("the value of flag --seed-max-dist (%d) is adjusted to contig interval length in database (%d)", maxDist, info.ContigInterval)
		// 	sopt.MaxDistance = float64(info.ContigInterval)
		// }

		idx, err := NewIndexSearcher(dbDir, sopt)
		checkError(err)

		if outputLog {
			log.Infof("index loaded in %s", time.Since(timeStart))
			log.Info()
		}

		if outputLog {
			log.Infof("searching with %d threads...", opt.NumCPUs)
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

		// fmt.Fprintf(outfh, "query\tqlen\tqstart\tqend\thits\tsgenome\tsseqid\tqcovGnm\thsp\tqcovHSP\talenHSP\talenSeg\tpident\tslen\tsstart\tsend\tsstr\tseeds\n")
		fmt.Fprintf(outfh, "query\tqlen\thits\tsgenome\tsseqid\tqcovGnm\thsp\tqcovHSP\talenHSP\tpident\tgaps\tqstart\tqend\tsstart\tsend\tsstr\tslen\tevalue\tbitscore")
		if moreColumns {
			fmt.Fprintf(outfh, "\tcigar\tqseq\tsseq\talign")
		}
		fmt.Fprintln(outfh)

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
			var targets = len(*q.result)
			matched++

			var strand byte
			var j int
			for _, r := range *q.result { // each genome
				j = 1
				for _, sd = range *r.SimilarityDetails { // each chain
					cr = sd.Similarity

					if sd.RC {
						strand = '-'
					} else {
						strand = '+'
					}

					for _, c = range *cr.Chains { // each match
						if c == nil {
							continue
						}

						if sd.RC {
							strand = '-'
						} else {
							strand = '+'
						}

						fmt.Fprintf(outfh, "%s\t%d\t%d\t%s\t%s\t%.3f\t%d\t%.3f\t%d\t%.3f\t%d\t%d\t%d\t%d\t%d\t%c\t%d\t%.2e\t%d",
							queryID, len(q.seq),
							targets, r.ID, sd.SeqID, r.AlignedFraction,
							j, c.AlignedFraction, c.AlignedLength, c.PIdent, c.Gaps,
							c.QBegin+1, c.QEnd+1,
							c.TBegin+1, c.TEnd+1,
							strand, sd.SeqLen,
							c.Evalue, c.BitScore,
						)
						if moreColumns {
							if onlyPseudoAlign {
								fmt.Fprintf(outfh, "\t%s\t%s\t%s\t%s", c.CIGAR, q.seq[c.QBegin:c.QEnd+1], c.TSeq, c.Alignment)
							} else {
								fmt.Fprintf(outfh, "\t%s\t%s\t%s\t%s", c.CIGAR, c.QSeq, c.TSeq, c.Alignment)
							}
						}

						fmt.Fprintln(outfh)

						j++
					}
				}
				outfh.Flush()
			}
			idx.RecycleSearchResults(q.result)

			poolQuery.Put(q)
			outfh.Flush()
		}

		// outputter
		ch := make(chan *Query, maxQueryConcurrency)
		done := make(chan int)
		go func() {

			for r := range ch {
				printResult(r)
			}

			done <- 1
		}()

		var wg sync.WaitGroup
		tokens := make(chan int, maxQueryConcurrency)

		var record *fastx.Record
		K := idx.k

		idx.SetSeqCompareOptions(&SeqComparatorOptions{
			K:         uint8(K),
			MinPrefix: 11, // can not be too small, or there will be a large number of anchors.

			Chaining2Options: Chaining2Options{
				// should be relative small
				MaxGap: maxAlignMaxGap,
				// better be larger than MinPrefix
				MinScore:    int(float64(minAlignLen) * minIdent / 100),
				MinAlignLen: minAlignLen,
				MinIdentity: minIdent,
				// can not be < k
				// MaxDistance: maxAlignMismatch,
				// can not be two small
				BandBase:  alignBand,
				BandCount: int(alignBand / 2),
			},

			MinAlignedFraction: minQcovChain,
			MinIdentity:        minIdent,
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
					continue
				}

				tokens <- 1
				wg.Add(1)

				query.seqID = append(query.seqID, record.ID...)
				query.seq = append(query.seq, bytes.ToUpper(record.Seq.Seq)...)

				go func(query *Query) {
					defer func() {
						<-tokens
						wg.Done()
					}()

					var err error
					query.result, err = idx.Search(query)
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

		checkError(idx.Close())
	},
}

func init() {
	RootCmd.AddCommand(mapCmd)

	mapCmd.Flags().StringP("index", "d", "",
		formatFlagUsage(`Index directory created by "lexicmap index".`))

	mapCmd.Flags().StringP("out-file", "o", "-",
		formatFlagUsage(`Out file, supports a ".gz" suffix ("-" for stdout).`))

	mapCmd.Flags().IntP("max-open-files", "", 1024,
		formatFlagUsage(`Maximum opened files. It mainly affects candidate subsequence extraction. Increase this value if you have hundreds of genome batches or have multiple queries, and do not forgot to set a bigger "ulimit -n" in shell if the value is > 1024.`))

	mapCmd.Flags().BoolP("all", "a", false,
		formatFlagUsage(`Output more columns, e.g., matched sequences. Use this if you want to output blast-style format with "lexicmap utils 2blast".`))

	mapCmd.Flags().IntP("max-query-conc", "J", 12,
		formatFlagUsage(`Maximum number of concurrent queries. Bigger values do not improve the batch searching speed and consume much memory.`))

	// mapCmd.Flags().IntP("max-seed-conc", "S", 8,
	// 	formatFlagUsage(`Maximum number of concurrent seed file matching. Bigger values improve seed matching speed in SSD.`))

	// seed searching

	mapCmd.Flags().IntP("seed-min-prefix", "p", 15,
		formatFlagUsage(`Minimum (prefix) length of matched seeds.`))

	mapCmd.Flags().IntP("seed-min-single-prefix", "P", 17,
		formatFlagUsage(`Minimum (prefix) length of matched seeds if there's only one pair of seeds matched.`))

	// mapCmd.Flags().IntP("seed-min-matches", "m", 20,
	// 	formatFlagUsage(`Minimum matched bases in the only one pair of seeds.`))

	// mapCmd.Flags().IntP("seed-max-mismatch", "m", -1,
	// 	formatFlagUsage(`Maximum mismatch between non-prefix regions of shared substrings.`))

	mapCmd.Flags().IntP("seed-max-gap", "", 200,
		formatFlagUsage(`Max gap in seed chaining.`))
	mapCmd.Flags().IntP("seed-max-dist", "", 1000,
		formatFlagUsage(`Max distance between seeds in seed chaining. It should be <= contig interval length in database.`))

	mapCmd.Flags().IntP("top-n-genomes", "n", 0,
		formatFlagUsage(`Keep top N genome matches for a query (0 for all) in chaining phase. Value 1 is not recommended as the best chaining result does not always bring the best alignment, so it better be >= 10.`))

	mapCmd.Flags().BoolP("load-whole-seeds", "w", false,
		formatFlagUsage(`Load the whole seed data into memory for faster search.`))

	// pseudo alignment
	mapCmd.Flags().BoolP("pseudo-align", "", false,
		formatFlagUsage(`Only perform pseudo alignment, alignment metrics, including qcovGnm, qcovSHP and pident, will be less accurate.`))

	mapCmd.Flags().IntP("align-ext-len", "", 1000,
		formatFlagUsage(`Extend length of upstream and downstream of seed regions, for extracting query and target sequences for alignment. It should be <= contig interval length in database.`))

	mapCmd.Flags().IntP("align-max-gap", "", 20,
		formatFlagUsage(`Maximum gap in a HSP segment.`))
	// mapCmd.Flags().IntP("align-max-kmer-dist", "", 100,
	// 	formatFlagUsage(`Maximum distance of (>=11bp) k-mers in a HSP segment.`))
	mapCmd.Flags().IntP("align-band", "", 100,
		formatFlagUsage(`Band size in backtracking the score matrix (pseduo alignment phase).`))
	mapCmd.Flags().IntP("align-min-match-len", "l", 50,
		formatFlagUsage(`Minimum aligned length in a HSP segment.`))

	// general filtering thresholds

	mapCmd.Flags().Float64P("align-min-match-pident", "i", 70,
		formatFlagUsage(`Minimum base identity (percentage) in a HSP segment.`))

	mapCmd.Flags().Float64P("min-qcov-per-hsp", "q", 0,
		formatFlagUsage(`Minimum query coverage (percentage) per HSP.`))

	mapCmd.Flags().Float64P("min-qcov-per-genome", "Q", 0,
		formatFlagUsage(`Minimum query coverage (percentage) per genome.`))

	mapCmd.Flags().BoolP("debug", "", false,
		formatFlagUsage(`Print debug information.`))

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
