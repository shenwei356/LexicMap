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
	"bufio"
	"bytes"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
)

var mapCmd = &cobra.Command{
	Use:   "search",
	Short: "Search sequences against an index",
	Long: `Search sequences against an index

Attention:
  1. Input should be (gzipped) FASTA or FASTQ records from files or stdin.
  2. For multiple queries, the order of queries in output might be different from the input.

Tips:
  1. When using -a/--all, the search result would be formatted to Blast-style format
     with 'lexicmap utils 2blast'. And the search speed would be slightly slowed down.
  2. Alignment result filtering is performed in the final phase, so stricter filtering criteria,
     including -q/--min-qcov-per-hsp, -Q/--min-qcov-per-genome, and -i/--align-min-match-pident,
     do not significantly accelerate the search speed. Hence, you can search with default
     parameters and then filter the result with tools like awk or csvtk.
  3. Users can limit search by TaxId(s) via -t/--taxids or --taxid-file.
     Only genomes with descendant TaxIds of the specific ones or themselves are searched,
     in a similar way with BLAST+ 2.15.0 or later versions.
     Negative values are allowed as a black list.

     For example, searching non-Escherichia (561) genera of Enterobacteriaceae (543) family with
     -t 543,-561.

     Users only need to provide NCBI-format taxdump files (-T/--taxdump, can also create from
     any taxonomy data with TaxonKit https://bioinf.shenwei.me/taxonkit/usage/#create-taxdump )
     and a genome-ID-to-TaxId mapping file (-G/--genome2taxid).
     There's no need to rebuild the index.

Alignment result relationship:

  Query
  ├── Subject genome
      ├── Subject sequence
          ├── HSP cluster (a cluster of neighboring HSPs)
              ├── High-Scoring segment Pair (HSP)

  Here, the defination of HSP is similar with that in BLAST. Actually there are small gaps in HSPs.

  > A High-scoring Segment Pair (HSP) is a local alignment with no gaps that achieves one of the
  > highest alignment scores in a given search. https://www.ncbi.nlm.nih.gov/books/NBK62051/

Output format:
  Tab-delimited format with 20+ columns, with 1-based positions.

    1.  query,    Query sequence ID.
    2.  qlen,     Query sequence length.
    3.  hits,     Number of subject genomes.
    4.  sgenome,  Subject genome ID.
    5.  sseqid,   Subject sequence ID.
    6.  qcovGnm,  Query coverage (percentage) per genome: $(aligned bases in the genome)/$qlen.
    7.  cls,      Nth HSP cluster in the genome. (just for improving readability)
                  It's useful to show if multiple adjacent HSPs are collinear.
    8.  hsp,      Nth HSP in the genome.         (just for improving readability)
    9.  qcovHSP   Query coverage (percentage) per HSP: $(aligned bases in a HSP)/$qlen.
    10. alenHSP,  Aligned length in the current HSP.
    11. pident,   Percentage of identical matches in the current HSP.
    12. gaps,     Gaps in the current HSP.
    13. qstart,   Start of alignment in query sequence.
    14. qend,     End of alignment in query sequence.
    15. sstart,   Start of alignment in subject sequence.
    16. send,     End of alignment in subject sequence.
    17. sstr,     Subject strand.
    18. slen,     Subject sequence length.
    19. evalue,   Expect value.
    20. bitscore, Bit score.
    21. cigar,    CIGAR string of the alignment.                      (optional with -a/--all)
    22. qseq,     Aligned part of query sequence.                     (optional with -a/--all)
    23. sseq,     Aligned part of subject sequence.                   (optional with -a/--all)
    24. align,    Alignment text ("|" and " ") between qseq and sseq. (optional with -a/--all)

Result ordering:
  For a HSP cluster, SimilarityScore = max(bitscore*pident)
  1. Within each HSP cluster, HSPs are sorted by sstart.
  2. Within each subject genome, HSP clusters are sorted in descending order by SimilarityScore.
  3. Results of multiple subject genomes are sorted by the highest SimilarityScore of HSP clusters.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		seq.ValidateSeq = false

		outFile := getFlagString(cmd, "out-file")

		var fhLog *os.File
		if opt.Log2File {
			ro, err := filepath.Abs(outFile)
			if err != nil {
				checkError(fmt.Errorf("failed to check output file: %s", err))
			}
			rl, err := filepath.Abs(opt.LogFile)
			if err != nil {
				checkError(fmt.Errorf("failed to check log file: %s", err))
			}
			if ro == rl {
				checkError(fmt.Errorf("output file and log file should not be the same: %s", outFile))
			}
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
		maxEvalue := getFlagNonNegativeFloat64(cmd, "max-evalue")

		// } else if minIdent < 1 {
		// 	log.Warningf("the value of flag -i/--align-min-match-pident is percentage in a range of [0, 100], you set: %f", minIdent)
		// }
		minQcovChain := getFlagNonNegativeFloat64(cmd, "min-qcov-per-hsp")
		if minQcovChain > 100 {
			checkError(fmt.Errorf("the value of flag -q/--min-qcov-per-hsp (%f) should be in range of [0, 100]", minIdent))
		}

		maxOpenFiles := getFlagPositiveInt(cmd, "max-open-files")

		// taxonomy
		taxdumpDir := getFlagString(cmd, "taxdump")
		genome2taxidFile := getFlagString(cmd, "genome2taxid")
		taxidsStr := getFlagStringSlice(cmd, "taxids")
		keepGenomesWithoutTaxId := getFlagBool(cmd, "keep-genomes-without-taxid")
		var taxids, negativeTaxids []uint32
		taxidFile := getFlagString(cmd, "taxid-file")

		var m, negativeM map[uint32]interface{}
		var ok bool
		var v uint32

		if len(taxidsStr) > 0 {
			if !(taxdumpDir != "" && genome2taxidFile != "") {
				checkError(fmt.Errorf("flags -T/--taxdump and -G/--genome2taxid are need if -t/--taxids is given"))
			}
			m = make(map[uint32]interface{}, len(taxidsStr))
			negativeM = make(map[uint32]interface{}, len(taxidsStr))
			taxids = make([]uint32, 0, len(taxidsStr))
			negativeTaxids = make([]uint32, 0, len(taxidsStr))

			var val int64
			for _, tmp := range taxidsStr {
				val, err = strconv.ParseInt(tmp, 10, 32)
				if err != nil {
					checkError(fmt.Errorf("invalid TaxId: %s", tmp))
				}

				if val > 0 {
					v = uint32(val)
					if _, ok = m[v]; !ok {
						taxids = append(taxids, v)
						m[v] = struct{}{}
					}
				} else if val < 0 {
					v = uint32(-val)
					if _, ok = negativeM[v]; !ok {
						negativeTaxids = append(negativeTaxids, v)
						negativeM[v] = struct{}{}
					}
				}
			}
		}
		if taxidFile != "" {
			if m == nil {
				m = make(map[uint32]interface{}, len(taxidsStr))
			}

			fh, err := xopen.Ropen(taxidFile)
			if err != nil {
				checkError(fmt.Errorf("failed to read taxid file: %s", taxidFile))
			}

			scanner := bufio.NewScanner(fh)
			var line string
			var val int64
			for scanner.Scan() {
				line = strings.TrimSpace(strings.TrimRight(scanner.Text(), "\r\n"))
				if line == "" {
					continue
				}

				val, err = strconv.ParseInt(line, 10, 32)
				if err != nil {
					checkError(fmt.Errorf("invalid TaxId: %s", line))
				}

				if val > 0 {
					v = uint32(val)
					if _, ok = m[v]; !ok {
						taxids = append(taxids, v)
						m[v] = struct{}{}
					}
				} else if val < 0 {
					v = uint32(-val)
					if _, ok = negativeM[v]; !ok {
						negativeTaxids = append(negativeTaxids, v)
						negativeM[v] = struct{}{}
					}
				}
			}
			if err = scanner.Err(); err != nil {
				checkError(fmt.Errorf("failed to read taxid file: %s", taxidFile))
			}
		}
		// } else if taxdumpDir != "" {
		// 	checkError(fmt.Errorf("the flag -T/--taxdump is given, but -t/--taxids is not"))
		// } else if genome2taxidFile != "" {
		// 	checkError(fmt.Errorf("the flag -G/--genome2taxid is given, but -t/--taxids is not"))
		// }

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

		_gcInterval := getFlagNonNegativeInt(cmd, "gc-interval")
		gcInterval := uint64(_gcInterval)
		if gcInterval > 0 {
			gcInterval = uint64(roundup32(uint32(gcInterval)))
			if gcInterval == 1 {
				gcInterval = 2
			}
		}

		gc := gcInterval > 0

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
			MaxEvalue:                        maxEvalue,

			OutputSeq: moreColumns,

			Debug: getFlagBool(cmd, "debug"),

			TaxdumpDir:              taxdumpDir,
			Genome2TaxIdFile:        genome2taxidFile,
			TaxIds:                  taxids,
			NegativeTaxIds:          negativeTaxids,
			KeepGenomesWithoutTaxId: keepGenomesWithoutTaxId,
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
			if gc {
				log.Infof("  maximum number of concurrent queries: %d, force garbage collection for every %d queries", maxQueryConcurrency, gcInterval)
			}
			if len(taxids)+len(negativeTaxids) > 0 {
				log.Infof("  filtering genomes by %d TaxIds and %d negative TaxIds", len(taxids), len(negativeTaxids))
			}
		}

		// ---------------------------------------------------------------
		// mapping

		id2name := idx.BatchGenomeIndex2GenomeID

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

		fmt.Fprintf(outfh, "query\tqlen\thits\tsgenome\tsseqid\tqcovGnm\tcls\thsp\tqcovHSP\talenHSP\tpident\tgaps\tqstart\tqend\tsstart\tsend\tsstr\tslen\tevalue\tbitscore")
		if moreColumns {
			fmt.Fprintf(outfh, "\tcigar\tqseq\tsseq\talign")
		}
		fmt.Fprintln(outfh)

		gcIntervalMinus1 := gcInterval - 1

		printResult := func(q *Query) {
			total++
			if q.result == nil { // seqs shorter than K or queries without matches.
				poolQuery.Put(q)

				if gc && total&gcIntervalMinus1 == 0 {
					runtime.GC()
				}
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
			var _c, j int
			for _, r := range *q.result { // each genome
				_c = 1
				j = 1
				for _, sd = range *r.SimilarityDetails { // each chain
					cr = sd.Similarity

					// if sd.RC {
					// 	strand = '-'
					// } else {
					// 	strand = '+'
					// }

					for _, c = range *cr.Chains { // each match
						if c == nil {
							continue
						}

						if sd.RC {
							strand = '-'
						} else {
							strand = '+'
						}

						fmt.Fprintf(outfh, "%s\t%d\t%d\t%s\t%s\t%.3f\t%d\t%d\t%.3f\t%d\t%.3f\t%d\t%d\t%d\t%d\t%d\t%c\t%d\t%.2e\t%d",
							queryID, len(q.seq),
							// targets, r.ID, sd.SeqID, r.AlignedFraction,
							targets, id2name[r.BatchGenomeIndex], sd.SeqID, r.AlignedFraction,
							_c,
							j, c.AlignedFraction, c.AlignedLength, c.PIdent, c.Gaps,
							c.QBegin+1, c.QEnd+1,
							c.TBegin+1, c.TEnd+1,
							strand, sd.SeqLen,
							c.Evalue, c.BitScore,
						)
						if moreColumns {
							fmt.Fprintf(outfh, "\t%s\t%s\t%s\t%s", c.CIGAR, c.QSeq, c.TSeq, c.Alignment)
						}

						fmt.Fprintln(outfh)

						j++
					}
					_c++
				}
				outfh.Flush()
			}
			idx.RecycleSearchResults(q.result)

			poolQuery.Put(q)
			outfh.Flush()

			if gc && total&gcIntervalMinus1 == 0 {
				runtime.GC()
			}
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
			K:         uint8(31),
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

	mapCmd.Flags().IntP("max-query-conc", "J", 8,
		formatFlagUsage(`Maximum number of concurrent queries. Bigger values do not improve the batch searching speed and consume much memory.`))

	mapCmd.Flags().IntP("gc-interval", "", 64,
		formatFlagUsage(`Force garbage collection every N queries (0 for disable). The value can't be too small.`))

	// mapCmd.Flags().IntP("max-seed-conc", "S", 8,
	// 	formatFlagUsage(`Maximum number of concurrent seed file matching. Bigger values improve seed matching speed in SSD.`))

	// seed searching

	mapCmd.Flags().IntP("seed-min-prefix", "p", 15,
		formatFlagUsage(`Minimum (prefix/suffix) length of matched seeds (anchors).`))

	mapCmd.Flags().IntP("seed-min-single-prefix", "P", 17,
		formatFlagUsage(`Minimum (prefix/suffix) length of matched seeds (anchors) if there's only one pair of seeds matched.`))

	// mapCmd.Flags().IntP("seed-min-matches", "m", 20,
	// 	formatFlagUsage(`Minimum matched bases in the only one pair of seeds.`))

	// mapCmd.Flags().IntP("seed-max-mismatch", "m", -1,
	// 	formatFlagUsage(`Maximum mismatch between non-prefix regions of shared substrings.`))

	mapCmd.Flags().IntP("seed-max-gap", "", 50,
		formatFlagUsage(`Minimum gap in seed chaining.`))
	mapCmd.Flags().IntP("seed-max-dist", "", 1000,
		formatFlagUsage(`Minimum distance between seeds in seed chaining. It should be <= contig interval length in database.`))

	mapCmd.Flags().IntP("top-n-genomes", "n", 0,
		formatFlagUsage(`Keep top N genome matches for a query (0 for all) in chaining phase. Value 1 is not recommended as the best chaining result does not always bring the best alignment, so it better be >= 100. (default 0)`))

	mapCmd.Flags().BoolP("load-whole-seeds", "w", false,
		formatFlagUsage(`Load the whole seed data into memory for faster seed matching. It will consume a lot of RAM.`))

	// pseudo alignment
	mapCmd.Flags().IntP("align-ext-len", "", 1000,
		formatFlagUsage(`Extend length of upstream and downstream of seed regions, for extracting query and target sequences for alignment. It should be <= contig interval length in database.`))

	mapCmd.Flags().IntP("align-max-gap", "", 20,
		formatFlagUsage(`Maximum gap in a HSP segment.`))
	// mapCmd.Flags().IntP("align-max-kmer-dist", "", 100,
	// 	formatFlagUsage(`Maximum distance of (>=11bp) k-mers in a HSP segment.`))
	mapCmd.Flags().IntP("align-band", "", 100,
		formatFlagUsage(`Band size in backtracking the score matrix (pseudo alignment phase).`))
	mapCmd.Flags().IntP("align-min-match-len", "l", 50,
		formatFlagUsage(`Minimum aligned length in a HSP segment.`))

	// general filtering thresholds

	mapCmd.Flags().Float64P("align-min-match-pident", "i", 70,
		formatFlagUsage(`Minimum base identity (percentage) in a HSP segment.`))

	mapCmd.Flags().Float64P("min-qcov-per-hsp", "q", 0,
		formatFlagUsage(`Minimum query coverage (percentage) per HSP.`))

	mapCmd.Flags().Float64P("min-qcov-per-genome", "Q", 0,
		formatFlagUsage(`Minimum query coverage (percentage) per genome.`))

	mapCmd.Flags().Float64P("max-evalue", "e", 10,
		formatFlagUsage(`Maximum evalue of a HSP segment.`))

	mapCmd.Flags().BoolP("debug", "", false,
		formatFlagUsage(`Print debug information, including a progress bar. (recommended when searching with one query).`))

	mapCmd.SetUsageTemplate(usageTemplate("-d <index path> [query.fasta.gz ...] [-o query.tsv.gz]"))

	// filter by taxids

	mapCmd.Flags().StringP("taxdump", "T", "",
		formatFlagUsage(`Directory containing taxdump files (nodes.dmp, names.dmp, etc.), needed for filtering results with TaxIds. For other non-NCBI taxonomy data, please use 'taxonkit create-taxdump' to create taxdump files.`))
	mapCmd.Flags().StringP("genome2taxid", "G", "",
		formatFlagUsage(`Two-column tabular file for mapping genome ID to TaxId, needed for filtering results with TaxIds. Genome IDs in the index can be exported via "lexicmap utils genomes -d db.lmi/ | csvtk cut -t -f 1 | csvtk uniq -Ut"`))
	mapCmd.Flags().BoolP("keep-genomes-without-taxid", "k", false,
		formatFlagUsage(`Keep genome hits without TaxId, i.e., those without TaxId in the --genome2taxid file.`))
	mapCmd.Flags().StringSliceP("taxids", "t", []string{},
		formatFlagUsage(`TaxIds(s) for filtering results, where the taxids are equal to or are the children of the given taxids. Negative values are allowed as a black list.`))
	mapCmd.Flags().StringP("taxid-file", "", "",
		formatFlagUsage(`TaxIds from a file for filtering results, where the taxids are equal to or are the children of the given taxids. Negative values are allowed as a black list.`))

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
