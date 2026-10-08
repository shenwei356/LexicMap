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
	"bytes"
	"cmp"
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"runtime"
	"slices"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/pkg/errors"
	"github.com/shenwei356/bio/seq"
	"github.com/spf13/cobra"
)

var gsearchCmd = &cobra.Command{
	Use:   "search",
	Short: "Search genomes against an index, with ANI and AF computed",
	Long: `Search genomes against an index, with ANI and AF computed

Algorithm:
  1. Genome screening: candidate genomes are screened by the total length of shared seeds
     (k-mers longer than the value of -p/--seed-min-prefix) between the query genome and
     the genomes in the index, with lexichash masking.
  2. Alignment:
     a) (default): cut the query genome into non-overlapping fragments and align them to
        the candidate genomes with an approach similar to that used in 'lexicmap search',
        with steps of seed matching, chaining, pseudoalignment, and base-level alignment.
     b) (OrthoANI mode): cut both the query genome and candidate genomes into non-overlapping
        fragments and only orthologous fragment pairs are used for calculating the ANI and
        AF values, which is similar to the algorithm of OrthoANI.

Attention:
  1. Input should be (gzipped) FASTA records from files or stdin, with one genome per file.
     Experimental feature: Long-read FASTQ files are also accepted. For these inputs,
     --adjust-ani-by-quality can be used to estimate ANI corrected for sequencing errors:
         e = mean_i(10^(-Q_i/10)) and ANIadj = (ANI - e/3) / (1 - 4e/3),
     where Q_i are the Phred scores in aligned query regions and ANI is expressed as a
     fraction. The correction assumes independent, symmetric substitution errors.
  2. One or more input files are accepted, via positional parameters
     and/or a file list via the flag -X/--infile-list.
  3. For multiple queries, the order of queries in output might be different from the input.

Taxonomic operations:
  1. Taxonomy data, including NCBI-format taxdump files (-T/--taxdump) and a genome-ID-to-TaxId
     mapping file (-G/--genome2taxid), are needed for filtering genomes by TaxId(s) and/or
     showing taxonomic names.
     
     Taxdump files can be createted from any taxonomy data with TaxonKit, see
     https://bioinf.shenwei.me/taxonkit/usage/#create-taxdump 

     If -T and -G are not provided, it will try to find them in the index directory,
     in the subdirectory "taxdump" and the file "taxdump/taxid.map", respectively.

  2. Users can limit search by TaxId(s) via -t/--taxids or --taxid-file.
     Only genomes with descendant TaxIds of the specific ones or themselves are searched,
     in a similar way with BLAST+ 2.15.0 or later versions.
     Negative values are allowed as a black list.

     For example, searching non-Escherichia (561) genera of Enterobacteriaceae (543) family with
     -t 543,-561.

  3. Flag --show-genome-name and --show-species-name can be used to show the taxonomic/species
     name of the subject genome in the output.

Output format:
  Tab-delimited format with 10 columns. --adjust-ani-by-quality appends
  ANIadj as the last column.

    1.  query,    Query genome ID.
    2.  subject,  Subject genome ID.
    3.  ANI,      Average nucleotide identity.
    4.  qAF,      Align fraction of the query genome (sum_aligned_length / sum_fragment_length).
    5.  sAF,      Align fraction of the subject genome.
                  Default mode : sum_aligned_length / sum_fragment_length)
                  OrthoANI mode: sum_aligned_length / genome size)
    6.  qctgs,    Number of contigs in the query genome.
    7.  qsize,    Size of the query genome.
    8.  sctgs,    Number of contigs in the subject genome.
    9.  ssize,    Size of the subject genome.
    10. sname,    Taxonomic name of the subject name.
    11. ANIadj,   Quality-adjusted ANI (optional, experimental feature).
 
`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		seq.ValidateSeq = false

		outFile := getFlagString(cmd, "out-file")
		defer setupCommandLog(opt, outFile)()
		verbose := opt.Verbose
		outputLog := opt.Verbose || opt.Log2File
		timeStart := time.Now()

		var err error

		// ---------------------------------------------------------------

		dbDir := getFlagString(cmd, "index")
		if dbDir == "" {
			checkError(fmt.Errorf("flag -d/--index needed"))
		}

		nMasks := getFlagNonNegativeInt(cmd, "masks")
		if !(nMasks == 0 || (isPowerOf4(nMasks) && nMasks >= 256)) {
			checkError(fmt.Errorf("the value of -m/--masks should be 0 (for all masks in the index) or power of 4 (needs to be >= 256, e.g., 256, 1024, 4096, 16384)"))
		}

		reRefNameStr := getFlagString(cmd, "ref-name-regexp")
		var reRefName *regexp.Regexp
		if reRefNameStr != "" {
			if !regexp.MustCompile(`\(.+\)`).MatchString(reRefNameStr) {
				checkError(fmt.Errorf(`value of --ref-name-regexp must contains "(" and ")" to capture the ref name from file name`))
			}
			if !reIgnoreCase.MatchString(reRefNameStr) {
				reRefNameStr = reIgnoreCaseStr + reRefNameStr
			}

			reRefName, err = regexp.Compile(reRefNameStr)
			if err != nil {
				checkError(errors.Wrapf(err, "failed to parse regular expression for matching sequence header: %s", reRefName))
			}
		}

		windows := getFlagPositiveInt(cmd, "windows")
		fragSize := getFlagPositiveInt(cmd, "frag-size")
		if fragSize < 100 {
			checkError(fmt.Errorf("the value of flag --frag-size should be >= 100"))
		}

		minFragLen := getFlagPositiveInt(cmd, "min-frag-size")
		if minFragLen < 100 {
			checkError(fmt.Errorf("the value of flag --min-frag-size should be >= 100"))
		}
		minAF := getFlagNonNegativeFloat64(cmd, "min-af") / 100
		minANI := getFlagNonNegativeFloat64(cmd, "min-ani") / 100
		adjustANIByQuality := getFlagBool(cmd, "adjust-ani-by-quality")

		onlyGenomeScreening := getFlagBool(cmd, "only-genome-screening")
		extra := getFlagBool(cmd, "extra")
		if adjustANIByQuality && onlyGenomeScreening {
			checkError(fmt.Errorf("--adjust-ani-by-quality cannot be used with --only-genome-screening"))
		}

		minPrefix := getFlagPositiveInt(cmd, "seed-min-prefix")
		if minPrefix > 32 || minPrefix < 5 {
			checkError(fmt.Errorf("the value of flag -p/--seed-min-prefix (%d) should be in the range of [5, 32]", minPrefix))
		}

		maxSubjectGenomeSize := getFlagNonNegativeInt(cmd, "max-subject-genome-size") * 1000 * 1000

		orthoANI := getFlagBool(cmd, "OrthoANI")
		samplingScale := getFlagPositiveInt(cmd, "kmer-scale")
		if samplingScale != 2 && samplingScale != 4 && samplingScale != 8 {
			checkError(fmt.Errorf("the value of flag --kmer-scale (%d) should be one of 2, 4, or 8", samplingScale))
		}
		gsa3SamplingScale = samplingScale

		// minSinglePrefix := getFlagPositiveInt(cmd, "seed-min-single-prefix")
		// if minSinglePrefix > 32 {
		// 	checkError(fmt.Errorf("the value of flag -P/--seed-min-single-prefix (%d) should be <= 32", minSinglePrefix))
		// }
		// if minSinglePrefix < minPrefix {
		// 	checkError(fmt.Errorf("the value of flag -P/--seed-min-single-prefix (%d) should be >= that of -p/--seed-min-prefix (%d)", minSinglePrefix, minPrefix))
		// }
		minSinglePrefix := minPrefix // not used in this command

		// maxGap := getFlagPositiveInt(cmd, "seed-max-gap")
		// maxDist := getFlagPositiveInt(cmd, "seed-max-dist")
		maxGap := 1  // not used in this command
		maxDist := 1 // not used in this command

		// extLen := getFlagNonNegativeInt(cmd, "align-ext-len")
		extLen := 1 // not used in this command

		topn := getFlagNonNegativeInt(cmd, "top-n-genomes")
		topNChains := getFlagNonNegativeInt(cmd, "top-n-chains")
		// topNChains := 5 // not used in this command

		inMemorySearch := getFlagBool(cmd, "load-whole-seeds")
		minAlignLen := getFlagPositiveInt(cmd, "align-min-match-len")
		if minAlignLen < 20 {
			checkError(fmt.Errorf("the value of flag -l/--align-min-match-len (%d) should be >= 20", minAlignLen))
		}
		maxAlignMaxGap := getFlagPositiveInt(cmd, "align-max-gap")
		alignBand := getFlagPositiveInt(cmd, "align-band")
		if alignBand < maxAlignMaxGap {
			checkError(fmt.Errorf("the value of flag --align-band should not be smaller thant the value of --align-max-gap"))
		}

		// minQcovGenome := getFlagNonNegativeFloat64(cmd, "min-qcov-per-genome")
		// if minQcovGenome > 100 {
		// 	checkError(fmt.Errorf("the value of flag -Q/--min-qcov-per-genome (%f) should be in range of [0, 100]", minQcovGenome))
		// }
		minQcovGenome := 0.0 // not used in this command

		minIdent := getFlagNonNegativeFloat64(cmd, "align-min-match-pident")
		if minIdent < 60 || minIdent > 100 {
			checkError(fmt.Errorf("the value of flag -i/--align-min-match-pident (%f) should be in range of [60, 100]", minIdent))
		}
		maxEvalue := getFlagNonNegativeFloat64(cmd, "max-evalue")

		minQcovChain := getFlagNonNegativeFloat64(cmd, "min-qcov-per-hsp")
		if minQcovChain > 100 {
			checkError(fmt.Errorf("the value of flag -q/--min-qcov-per-hsp (%f) should be in range of [0, 100]", minIdent))
		}
		if orthoANI {
			minQcovChain /= 2
			minFragLen = fragSize
			log.Warningf("When using OrthoANI mode, the value of -q/--min-qcov-per-hsp is halved (%.2f%%) and the value of --min-frag-size is set with the value of --frag-size (%d)", minQcovChain, fragSize)
		} else {
			// maxDesert := getFlagPositiveInt(cmd, "seed-max-desert")
			// seedInDesertDist := getFlagPositiveInt(cmd, "seed-in-desert-dist")
			// if seedInDesertDist > maxDesert/2 {
			// 	checkError(fmt.Errorf("value of --seed-in-desert-dist should be smaller than 0.5 * --seed-max-desert"))
			// }

			// gsa3DesertMaxLen = maxDesert
			// gsa3DesertExpectedSeedDist = seedInDesertDist
		}

		maxOpenFiles := getFlagPositiveInt(cmd, "max-open-files")

		// taxonomy
		taxdumpDir := getFlagString(cmd, "taxdump")
		if !cmd.Flags().Lookup("taxdump").Changed {
			taxdumpDir = filepath.Join(dbDir, "taxdump")
		}
		genome2taxidFile := getFlagString(cmd, "genome2taxid")
		if !cmd.Flags().Lookup("genome2taxid").Changed {
			genome2taxidFile = filepath.Join(taxdumpDir, "taxid.map")
		}
		taxidsStr := getFlagStringSlice(cmd, "taxids")
		taxidFile := getFlagString(cmd, "taxid-file")
		keepGenomesWithoutTaxId := getFlagBool(cmd, "keep-genomes-without-taxid")

		taxids, negativeTaxids := parseTaxids(taxdumpDir, genome2taxidFile, taxidsStr, taxidFile)

		showTaxName := getFlagBool(cmd, "show-genome-name")
		showSpeciesName := getFlagBool(cmd, "show-species-name")
		if (showTaxName || showSpeciesName) && !(taxdumpDir != "" && genome2taxidFile != "") {
			checkError(fmt.Errorf("flags -T/--taxdump and -G/--genome2taxid are needed for --show-genome-name"))
		}

		// ---------------------------------------------------------------

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
			maxQueryConcurrency = opt.NumCPUs
		}
		// maxSeedSearchingConcurrency := getFlagPositiveInt(cmd, "max-seed-matching-conc")
		maxSeedSearchingConcurrency := maxQueryConcurrency / 2
		if maxSeedSearchingConcurrency < 2 {
			maxSeedSearchingConcurrency = 2
		}

		var threadsPerQuery int
		if len(files) < maxQueryConcurrency {
			threadsPerQuery = int(float64(opt.NumCPUs) * 1.2 / float64(len(files)))
		} else {
			threadsPerQuery = int(float64(opt.NumCPUs) * 1.2 / float64(maxQueryConcurrency))
		}
		if threadsPerQuery < 1 {
			threadsPerQuery = 1
		} else if threadsPerQuery > opt.NumCPUs {
			threadsPerQuery = opt.NumCPUs
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

			MaxSeedSearchingConcurrency: maxSeedSearchingConcurrency,

			MinPrefix:       uint8(minPrefix),
			MinSinglePrefix: uint8(minSinglePrefix),
			TopN:            topn,
			TopNChains:      topNChains,
			InMemorySearch:  inMemorySearch,

			MaxGap:      float64(maxGap),
			MaxDistance: float64(maxDist),

			ExtendLength:  extLen,
			ExtendLength2: 50,

			MinQueryAlignedFractionInAGenome: minQcovGenome,
			MaxEvalue:                        maxEvalue,

			OutputSeq: false,

			Debug: getFlagBool(cmd, "debug"),

			TaxdumpDir:              taxdumpDir,
			Genome2TaxIdFile:        genome2taxidFile,
			TaxIds:                  taxids,
			NegativeTaxIds:          negativeTaxids,
			KeepGenomesWithoutTaxId: keepGenomesWithoutTaxId,
			LoadTaxName:             showTaxName,
			LoadTaxRank:             showSpeciesName,

			MaxSubjectGenomeSize: maxSubjectGenomeSize,
			SearchMaskCount:      nMasks,
		}

		idx, err := NewIndexSearcher(dbDir, sopt)
		checkError(err)

		_nMasks := len(idx.lh.Masks)
		if nMasks > 0 {
			_nMasks = nMasks
		}

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

				HeuristicKmerPidentThreshold: 0,
			},

			MinAlignedFraction: minQcovChain,
			MinIdentity:        minIdent,
		})

		kf := 11
		minSharedKmers := MinSharedKmersThresholdExact(fragSize, uint8(kf), uint32(samplingScale), 0.80, 0.99)
		idx.SetFragmentCompareOptions(&FragmentComparatorOptions{
			K:              uint8(kf),
			MinSharedKmers: max(3, minSharedKmers),
			Scaled:         uint32(samplingScale),
			TopNFragments:  topNChains,
		})

		if outputLog {
			log.Infof("index loaded in %s", time.Since(timeStart))
			log.Info()
		}

		if outputLog {
			log.Infof("searching with %d threads...", opt.NumCPUs)
			if nMasks > 0 {
				log.Infof("  screen genomes with seed data of %d masks", _nMasks)
			} else {
				log.Infof("  screen genomes with seed data of all the %d masks", _nMasks)
			}
			if sopt.TopN > 0 {
				log.Infof("    keep the top %d genomes", sopt.TopN)
			}
			// if sopt.TopNChains > 0 {
			// 	log.Infof("  keep the top %d chains", sopt.TopNChains)
			// }
			if orthoANI {
				log.Infof("  minimum shared k-mers between genome fragments: %d", minSharedKmers)
				log.Infof("  minimum query coverage per HSP: %.2f%% (the value of -q/--min-qcov-per-hsp is halved)", minQcovChain)
			} else {
				log.Infof("  minimum query coverage per HSP: %.2f%%", minQcovChain)
			}
			log.Infof("  minimum base identity in a HSP segment: %.2f%%", minIdent)
			log.Infof("  maximum evalue: %.2e", maxEvalue)

			if gc {
				log.Infof("  maximum number of concurrent queries: %d, force garbage collection for every %d queries", maxQueryConcurrency, gcInterval)
				log.Infof("  threads per query: %d", threadsPerQuery)
			}
			if len(taxids)+len(negativeTaxids) > 0 {
				log.Infof("  filtering genomes by %d TaxIds and %d negative TaxIds", len(taxids), len(negativeTaxids))
			}

		}

		// ---------------------------------------------------------------
		// searching

		// -------  output handler -------

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

		if !onlyGenomeScreening {
			if adjustANIByQuality {
				fmt.Fprintf(outfh, "query\tsubject\tANI\tqAF\tsAF\tqcontigs\tqsize\tscontigs\tssize\tsname\tANIadj\n")
			} else {
				fmt.Fprintf(outfh, "query\tsubject\tANI\tqAF\tsAF\tqcontigs\tqsize\tscontigs\tssize\tsname\n")
			}
		} else {
			// fmt.Fprintf(outfh, "query\tsubject\tminPrefix\tfracMasks\tnMasks\tnKmers\tnBases\tavgLen\tnBestBases\tavgBestLen")
			fmt.Fprintf(outfh, "query\tsubject\tminPrefix\tfracMasks\tnMasks\tsumPrefix\tavgPrefix\tsname")
			if extra {
				fmt.Fprintf(outfh, "\tprefixes")
			}
			fmt.Fprintf(outfh, "\n")
		}

		// -------  output function -------

		gcIntervalMinus1 := gcInterval - 1
		id2name := idx.BatchGenomeIndex2GenomeID
		tax := idx.Taxonomy

		matches := make([]uint8, 0, _nMasks)
		var matchesS bytes.Buffer

		matches2strslice := func(sep byte) {
			n := len(matches)
			if n == 0 {
				return
			}

			slices.SortFunc(matches, func(a, b uint8) int {
				return cmp.Compare(b, a)
			})

			matchesS.Reset()
			matchesS.WriteString(strconv.Itoa(int(matches[0])))

			if n == 1 {
				return
			}

			for _, v := range matches[1:] {
				matchesS.WriteByte(sep)
				matchesS.WriteString(strconv.Itoa(int(v)))
			}
		}
		recycleQuery := func(q *GQuery) {
			if q.screenDetails != nil {
				idx.RecycleGSearchScreenDetailResults(q.screenDetails)
				q.screenDetails = nil
			}
			RecycleGQuery(q)
		}

		printResult := func(q *GQuery) {
			total++
			if !((q.result != nil && len(*q.result) != 0) ||
				(q.screenDetails != nil && len(*q.screenDetails) != 0)) { // seqs shorter than K or queries without matches.
				recycleQuery(q)

				if gc && total&gcIntervalMinus1 == 0 {
					runtime.GC()
				}
				return
			}

			matched++
			if verbose {
				if (total < 128 && total&7 == 0) || total&127 == 0 {
					speed = float64(total) / time.Since(timeStart1).Minutes()
					fmt.Fprintf(os.Stderr, "processed queries: %d, speed: %.3f queries per minute\r", total, speed)
				}
			}

			var vSgenome string
			var taxid uint32

			if !onlyGenomeScreening {
				for _, gr := range *q.result {
					if showSpeciesName {
						vSgenome = tax.Name(idx.genomeIdx2TaxId[gr.BatchGenomeIndex])
						for _, taxid = range tax.LineageTaxIds(idx.genomeIdx2TaxId[gr.BatchGenomeIndex]) {
							if tax.Rank(taxid) == "species" {
								vSgenome = tax.Name(taxid)
								break
							}
						}
					} else if showTaxName {
						vSgenome = tax.Name(idx.genomeIdx2TaxId[gr.BatchGenomeIndex])
					} else {
						vSgenome = ""
					}
					if adjustANIByQuality {
						fmt.Fprintf(outfh, "%s\t%s\t%.3f\t%.3f\t%.3f\t%d\t%d\t%d\t%d\t%s\t%.3f\n",
							q.id, id2name[gr.BatchGenomeIndex], gr.ANI*100, gr.AFq*100, gr.AFs*100,
							len(q.seqs), q.genomeSize, gr.NumSeqs, gr.GenomeSize, vSgenome, gr.ANIAdjusted*100)
					} else {
						fmt.Fprintf(outfh, "%s\t%s\t%.3f\t%.3f\t%.3f\t%d\t%d\t%d\t%d\t%s\n",
							q.id, id2name[gr.BatchGenomeIndex], gr.ANI*100, gr.AFq*100, gr.AFs*100,
							len(q.seqs), q.genomeSize, gr.NumSeqs, gr.GenomeSize, vSgenome)
					}
				}
			} else {
				var hitMasks uint64
				var v uint8
				for _, gr := range *q.screenDetails {
					if showSpeciesName {
						vSgenome = tax.Name(idx.genomeIdx2TaxId[gr.BatchGenomeIndex[0]])
						for _, taxid = range tax.LineageTaxIds(idx.genomeIdx2TaxId[gr.BatchGenomeIndex[0]]) {
							if tax.Rank(taxid) == "species" {
								vSgenome = tax.Name(taxid)
								break
							}
						}
					} else if showTaxName {
						vSgenome = tax.Name(idx.genomeIdx2TaxId[gr.BatchGenomeIndex[0]])
					} else {
						vSgenome = ""
					}

					hitMasks = 0
					// _stats.Reset()
					if extra {
						matches = matches[:0]
					}
					for _, v = range gr.LongestMatches {
						if v > 0 {
							hitMasks++

							if extra {
								matches = append(matches, v)
							}
						}
					}

					if extra {
						matches2strslice(',')

						fmt.Fprintf(outfh, "%s\t%s\t%d\t%.4f\t%d\t%d\t%.2f\t%s\t%s\n",
							q.id, id2name[gr.BatchGenomeIndex[0]],
							minPrefix, float64(hitMasks)/float64(_nMasks), hitMasks,
							gr.SumPrefix, float64(gr.SumPrefix)/float64(hitMasks), vSgenome,
							matchesS.Bytes(),
						)
					} else {
						fmt.Fprintf(outfh, "%s\t%s\t%d\t%.4f\t%d\t%d\t%.2f\t%s\n",
							q.id, id2name[gr.BatchGenomeIndex[0]],
							minPrefix, float64(hitMasks)/float64(_nMasks), hitMasks,
							gr.SumPrefix, float64(gr.SumPrefix)/float64(hitMasks), vSgenome,
						)
					}
				}
			}

			recycleQuery(q)
			outfh.Flush()

			if gc && total&gcIntervalMinus1 == 0 {
				runtime.GC()
			}
		}

		ch := make(chan *GQuery, maxQueryConcurrency)
		done := make(chan int)
		go func() {
			for r := range ch {
				printResult(r)
			}

			done <- 1
		}()

		// -------  input  -------

		var wg sync.WaitGroup
		tokens := make(chan int, maxQueryConcurrency)

		gr := NewGenomeReader(idx.k, reRefName)

		for _, file := range files {
			tokens <- 1
			wg.Add(1)

			go func(file string) {
				defer func() {
					<-tokens
					wg.Done()
				}()

				// 1. read all sequences of the query genome
				var query *GQuery
				var readErr error
				if adjustANIByQuality {
					query, readErr = gr.ReadWithQual(file, true, idx.softMasking)
				} else {
					query, readErr = gr.Read(file, true, idx.softMasking)
				}
				checkError(readErr)
				if query == nil { // no valid sequence
					log.Warningf("no valid sequences in %s, skipped", file)
					return
				}

				// fmt.Printf("seqs: %d, len: %d\n", len(query.seqs), len(query.bigSeq))
				if query.genomeSize < fragSize {
					log.Warningf("query genome %s is smaller than fragment size (%d), skipped", query.id, fragSize)
					ch <- query
					return
				}

				// 2. search possible genome matches
				genomeIds, rs, screenErr := idx.GSearchScreen(query, windows, nil)
				checkError(screenErr)

				if onlyGenomeScreening {
					query.screenDetails = rs
					ch <- query
					return
				} else {
					idx.RecycleGSearchScreenDetailResults(rs)
				}

				if genomeIds != nil {
					// 3. search fragments for the query
					// err = idx.GSearchAlign(query, fragSize, minFragLen, genomeIds, minAF, minANI, maxQueryConcurrency, gcInterval)
					var alignErr error
					if orthoANI {
						alignErr = idx.GSearchAlignOrthoANI(query, fragSize, minFragLen, genomeIds, minAF, minANI, threadsPerQuery)
					} else {
						alignErr = idx.GSearchAlign3Sampled(query, fragSize, minFragLen, genomeIds, minAF, minANI, threadsPerQuery)
					}

					checkError(alignErr)

					// clear up

					idx.RecycleGSearchScreenResult(genomeIds)
				}

				ch <- query
			}(file)
		}
		wg.Wait()
		close(ch)
		<-done

		// -------  final log  -------

		if verbose {
			fmt.Fprintf(os.Stderr, "\n")
		}
		if outputLog {
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
	genomeCmd.AddCommand(gsearchCmd)

	// general flags

	gsearchCmd.Flags().StringP("index", "d", "",
		formatFlagUsage(`Index directory created by "lexicmap index".`))

	gsearchCmd.Flags().IntP("masks", "m", 1024,
		formatFlagUsage(`Only use seed data of N masks for genome prefiltering/screening. It should be 0 (for all masks in the index) or power of 4 (needs to be >= 256, e.g., 256, 1024, 4096, 16384).`))

	gsearchCmd.Flags().StringP("out-file", "o", "-",
		formatFlagUsage(`Out file, supports a ".gz" suffix ("-" for stdout).`))

	gsearchCmd.Flags().IntP("max-open-files", "", 1024,
		formatFlagUsage(`Maximum opened files. It mainly affects candidate genome extraction. Increase this value if you have hundreds of genome batches or have multiple queries, and do not forgot to set a bigger "ulimit -n" in shell if the value is > 1024.`))

	gsearchCmd.Flags().IntP("max-query-conc", "J", 8,
		formatFlagUsage(`Maximum number of concurrent queries.`))

	gsearchCmd.Flags().IntP("gc-interval", "", 4,
		formatFlagUsage(`Force garbage collection every N queries (0 for disable). The value can't be too small.`))

	// genome filtering/screening

	gsearchCmd.Flags().StringP("ref-name-regexp", "", `(?i)(.+)\.(f[aq](st[aq])?|fna)(\.gz|\.xz|\.zst|\.bz2)?$`,
		formatFlagUsage(`Regular expression (must contains "(" and ")") for extracting the reference name from the input filename. Attention: use double quotation marks for patterns containing commas, e.g., -p '"A{2,}"'.`))

	gsearchCmd.Flags().IntP("seed-min-prefix", "p", 21,
		formatFlagUsage(`Minimum prefix length of matched seeds in the genome filtering phase.`))

	gsearchCmd.Flags().IntP("windows", "", 1,
		formatFlagUsage(`The number of windows in lexichash masking, for genome screening (1 is sufficient).`))
	gsearchCmd.Flags().IntP("frag-size", "", 1020,
		formatFlagUsage(`The size of non-overlap fragments cut for ANI computation.`))
	gsearchCmd.Flags().IntP("min-frag-size", "", 100,
		formatFlagUsage(`The minimum length of fragments in the end of a sequence during cutting fragments.`))

	gsearchCmd.Flags().IntP("top-n-genomes", "n", 10,
		formatFlagUsage(`Keep the top N genome matches for a query (0 for all). Screening includes all candidates tied with the Nth score; ANI ranking then returns at most N matches.`))

	gsearchCmd.Flags().IntP("max-subject-genome-size", "", 20,
		formatFlagUsage(`Maximum size of subject genomes to be considered (in MB).`))

	// alignment

	// gsearchCmd.Flags().IntP("seed-min-single-prefix", "P", 17,
	// 	formatFlagUsage(`Minimum (prefix/suffix) length of matched seeds (anchors) if there's only one pair of seeds matched.`))

	// gsearchCmd.Flags().IntP("seed-max-gap", "", 50,
	// 	formatFlagUsage(`Minimum gap in seed chaining.`))
	// gsearchCmd.Flags().IntP("seed-max-dist", "", 1000,
	// 	formatFlagUsage(`Minimum distance between seeds in seed chaining. It should be <= contig interval length in database.`))

	gsearchCmd.Flags().IntP("top-n-chains", "N", 5,
		formatFlagUsage(`Keep the top N chains in a genome for the query (0 for all) in the chaining phase. Value 1 is not recommended as the best chaining result does not always bring the best alignment.`))

	gsearchCmd.Flags().BoolP("load-whole-seeds", "w", false,
		formatFlagUsage(`Load the whole seed data into memory for faster seed matching. It will consume a lot of RAM.`))

	// pseudo alignment
	// gsearchCmd.Flags().IntP("align-ext-len", "", 1000,
	// 	formatFlagUsage(`Extend length of upstream and downstream of seed regions, for extracting query and target sequences for alignment. It should be <= contig interval length in database.`))

	// gsearchCmd.Flags().IntP("seed-max-desert", "", 50,
	// 	formatFlagUsage(`Maximum length of sketching deserts, or maximum seed distance. Deserts with seed distance larger than this value will be filled by choosing k-mers roughly every --seed-in-desert-dist bases.`))
	// gsearchCmd.Flags().IntP("seed-in-desert-dist", "", 25,
	// 	formatFlagUsage(`Distance of k-mers to fill deserts.`))

	gsearchCmd.Flags().IntP("align-max-gap", "", 100,
		formatFlagUsage(`Maximum gap in a HSP segment.`))
	gsearchCmd.Flags().IntP("align-band", "", 100,
		formatFlagUsage(`Band size in backtracking the score matrix (pseudo alignment phase).`))
	gsearchCmd.Flags().IntP("align-min-match-len", "l", 30,
		formatFlagUsage(`Minimum aligned length in a HSP segment.`))

	// general filtering thresholds

	gsearchCmd.Flags().Float64P("align-min-match-pident", "i", 70,
		formatFlagUsage(`Minimum base identity (percentage) in a HSP segment.`))

	gsearchCmd.Flags().Float64P("min-qcov-per-hsp", "q", 30,
		formatFlagUsage(`Minimum query coverage (percentage) per HSP.`))

	// gsearchCmd.Flags().Float64P("min-qcov-per-genome", "Q", 0,
	// 	formatFlagUsage(`Minimum query coverage (percentage) per genome.`))

	gsearchCmd.Flags().Float64P("max-evalue", "e", 1e-15,
		formatFlagUsage(`Maximum evalue of a HSP segment.`))

	gsearchCmd.Flags().BoolP("debug", "", false,
		formatFlagUsage(`Print debug information, including a progress bar. (recommended when searching with one query).`))

	gsearchCmd.SetUsageTemplate(usageTemplate("-d <index path> [query.fasta[.gz] ...] [-o result.tsv[.gz]]"))

	// filter by taxids

	gsearchCmd.Flags().StringP("taxdump", "T", "",
		formatFlagUsage(`Directory containing taxdump files (nodes.dmp, names.dmp, etc.), needed for filtering results with TaxIds. For other non-NCBI taxonomy data, please use 'taxonkit create-taxdump' to create taxdump files. (default: <index path>/taxdump)`))
	gsearchCmd.Flags().StringP("genome2taxid", "G", "",
		formatFlagUsage(`Two-column tabular file for mapping genome ID to TaxId, needed for filtering results with TaxIds. Genome IDs in the index can be exported via "lexicmap utils genomes -d db.lmi/ | csvtk cut -t -f 1 | csvtk uniq -Ut". (default: <taxdump path>/taxid.map)`))
	gsearchCmd.Flags().BoolP("keep-genomes-without-taxid", "k", false,
		formatFlagUsage(`Keep genome hits without TaxId, i.e., those without TaxId in the --genome2taxid file.`))
	gsearchCmd.Flags().StringSliceP("taxids", "t", []string{},
		formatFlagUsage(`TaxIds(s) for filtering results, where the taxids are equal to or are the children of the given taxids. Negative values are allowed as a black list.`))
	gsearchCmd.Flags().StringP("taxid-file", "", "",
		formatFlagUsage(`TaxIds from a file for filtering results, where the taxids are equal to or are the children of the given taxids. Negative values are allowed as a black list.`))

	// ani-af related filtering

	gsearchCmd.Flags().Float64P("min-af", "F", 15.0,
		formatFlagUsage(`Only output results where one genome has aligned fraction > than this value (percentage).`))

	gsearchCmd.Flags().Float64P("min-ani", "I", 70,
		formatFlagUsage(`Only output results where one genome has ANI > than this value (percentage).`))

	gsearchCmd.Flags().BoolP("adjust-ani-by-quality", "", false,
		formatFlagUsage(`(Experimental feature) Estimate ANI corrected for Phred+33 sequencing errors in the aligned query regions and add the ANIAdjusted column. FASTQ input is required. The raw ANI is still used for filtering and sorting.`))

	// OrthoANI
	gsearchCmd.Flags().BoolP("OrthoANI", "O", false,
		formatFlagUsage(`Compute OrthoANI using reciprocal best hit of fragment pairs. Type 'lexicmap genome search --help' for details.`))

	gsearchCmd.Flags().IntP("kmer-scale", "", 4,
		formatFlagUsage(`Using 1/scale of k-mers for seeding (default mode) or fragment comparison (OrthoANI mode). Available values: 2, 4, 8.`))

	// only for research
	gsearchCmd.Flags().BoolP("only-genome-screening", "S", false,
		formatFlagUsage(`Only perform genome screening, no ANI computation.`))
	gsearchCmd.Flags().BoolP("extra", "", false,
		formatFlagUsage(`Show extra columns for -S/--only-genome-screening.`))

	// extra info
	gsearchCmd.Flags().BoolP("show-genome-name", "g", false,
		formatFlagUsage(`Show the taxonomic name of subject genome in field 'sname'. Flags -T/--taxdump and -G/--genome2taxid are needed.`))
	gsearchCmd.Flags().BoolP("show-species-name", "s", false,
		formatFlagUsage(`Show the species name of subject genome in field 'sname'. Flags -T/--taxdump and -G/--genome2taxid are needed.`))
}
