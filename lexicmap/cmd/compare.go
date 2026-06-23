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
	"os"
	"path/filepath"
	"regexp"
	"runtime"
	"strings"
	"time"

	"github.com/pkg/errors"
	"github.com/shenwei356/bio/seq"
	"github.com/spf13/cobra"
)

var compareCmd = &cobra.Command{
	Use:   "compare",
	Short: "Compare genome pairs",
	Long: `Compare genome pairs

Input:
  - 
  - 
 
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

		fragSize := getFlagPositiveInt(cmd, "frag-size")
		if fragSize < 100 {
			checkError(fmt.Errorf("the value of flag --frag-size should be >= 100"))
		}

		minFragLen := getFlagPositiveInt(cmd, "min-frag-size")
		if minFragLen < 100 {
			checkError(fmt.Errorf("the value of flag --min-frag-size should be >= 100"))
		}
		// minAF := getFlagNonNegativeFloat64(cmd, "min-af") / 100
		// minAF := 0

		minPrefix := getFlagPositiveInt(cmd, "seed-min-prefix")
		if minPrefix > 32 || minPrefix < 5 {
			checkError(fmt.Errorf("the value of flag -p/--seed-min-prefix (%d) should be in the range of [5, 32]", minPrefix))
		}

		// maxSubjectGenomeSize := getFlagNonNegativeInt(cmd, "max-subject-genome-size") * 1000 * 1000
		maxSubjectGenomeSize := 1 << 30

		samplingScale := getFlagPositiveInt(cmd, "kmer-scale")
		if samplingScale != 2 && samplingScale != 4 && samplingScale != 8 {
			checkError(fmt.Errorf("the value of flag --kmer-scale (%d) should be one of 2, 4, or 8", samplingScale))
		}
		gsa3SamplingScale = samplingScale

		minSinglePrefix := minPrefix // not used in this command

		maxGap := 1  // not used in this command
		maxDist := 1 // not used in this command
		extLen := 1  // not used in this command

		topn := 10
		topNChains := 5 // not used in this command

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

		maxOpenFiles := getFlagPositiveInt(cmd, "max-open-files")

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
		var threadsPerQuery int
		if len(files) < maxQueryConcurrency {
			threadsPerQuery = opt.NumCPUs / len(files)
		} else {
			threadsPerQuery = opt.NumCPUs / maxQueryConcurrency
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

			// TaxdumpDir:              taxdumpDir,
			// Genome2TaxIdFile:        genome2taxidFile,
			// TaxIds:                  taxids,
			// NegativeTaxIds:          negativeTaxids,
			// KeepGenomesWithoutTaxId: keepGenomesWithoutTaxId,

			MaxSubjectGenomeSize: maxSubjectGenomeSize,
		}

		idx, err := NewIndexSearcher(dbDir, sopt)
		checkError(err)

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
			log.Infof("  minimum query coverage per HSP: %.2f%%", minQcovChain)
			log.Infof("  minimum base identity in a HSP segment: %.2f%%", minIdent)
			log.Infof("  maximum evalue: %.2e", maxEvalue)

			if gc {
				log.Infof("  maximum number of concurrent comparing: %d, force garbage collection for every %d genome pairs", maxQueryConcurrency, gcInterval)
				log.Infof("  threads per query: %d", threadsPerQuery)
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

		var total uint64
		var speed float64 // k reads/second

		// -------  output function -------

		// gcIntervalMinus1 := gcInterval - 1

		// ------- input -------
		// for _, file := range files {

		// }

		// -------  final log  -------

		if verbose {
			fmt.Fprintf(os.Stderr, "\n")
		}
		if outputLog {
			speed = float64(total) / time.Since(timeStart1).Minutes()
			log.Infof("")
			log.Infof("processed genome pairs: %d, speed: %.3f pairs per minute\n", total, speed)
			log.Infof("done comparing")
			if outFile != "-" {
				log.Infof("search results saved to: %s", outFile)
			}
		}

		checkError(idx.Close())
	},
}

func init() {
	genomeCmd.AddCommand(compareCmd)

	// general flags

	compareCmd.Flags().StringP("index", "d", "",
		formatFlagUsage(`Index directory created by "lexicmap index".`))

	compareCmd.Flags().StringP("out-file", "o", "-",
		formatFlagUsage(`Out file, supports a ".gz" suffix ("-" for stdout).`))

	compareCmd.Flags().IntP("max-open-files", "", 1024,
		formatFlagUsage(`Maximum opened files. It mainly affects candidate genome extraction. Increase this value if you have hundreds of genome batches or have multiple queries, and do not forgot to set a bigger "ulimit -n" in shell if the value is > 1024.`))

	compareCmd.Flags().IntP("max-query-conc", "J", 4,
		formatFlagUsage(`Maximum number of concurrent queries.`))

	compareCmd.Flags().IntP("gc-interval", "", 4,
		formatFlagUsage(`Force garbage collection every N queries (0 for disable). The value can't be too small.`))

	// genome filtering/screening

	compareCmd.Flags().StringP("ref-name-regexp", "", `(?i)(.+)\.(f[aq](st[aq])?|fna)(\.gz|\.xz|\.zst|\.bz2)?$`,
		formatFlagUsage(`Regular expression (must contains "(" and ")") for extracting the reference name from the input filename. Attention: use double quotation marks for patterns containing commas, e.g., -p '"A{2,}"'.`))

	compareCmd.Flags().IntP("seed-min-prefix", "p", 21,
		formatFlagUsage(`Minimum prefix length of matched seeds in the genome filtering phase.`))

	compareCmd.Flags().IntP("frag-size", "", 1020,
		formatFlagUsage(`The size of non-overlap fragments cut for ANI computation.`))
	compareCmd.Flags().IntP("min-frag-size", "", 100,
		formatFlagUsage(`The minimum length of fragments in the end of a sequence during cutting fragments.`))

	// alignment

	// compareCmd.Flags().IntP("top-n-chains", "N", 5,
	// 	formatFlagUsage(`Keep the top N chains in a genome for the query (0 for all) in the chaining phase. Value 1 is not recommended as the best chaining result does not always bring the best alignment.`))

	compareCmd.Flags().BoolP("load-whole-seeds", "w", false,
		formatFlagUsage(`Load the whole seed data into memory for faster seed matching. It will consume a lot of RAM.`))

	compareCmd.Flags().IntP("align-max-gap", "", 100,
		formatFlagUsage(`Maximum gap in a HSP segment.`))
	compareCmd.Flags().IntP("align-band", "", 100,
		formatFlagUsage(`Band size in backtracking the score matrix (pseudo alignment phase).`))
	compareCmd.Flags().IntP("align-min-match-len", "l", 30,
		formatFlagUsage(`Minimum aligned length in a HSP segment.`))

	// general filtering thresholds

	compareCmd.Flags().Float64P("align-min-match-pident", "i", 70,
		formatFlagUsage(`Minimum base identity (percentage) in a HSP segment.`))

	compareCmd.Flags().Float64P("min-qcov-per-hsp", "q", 30,
		formatFlagUsage(`Minimum query coverage (percentage) per HSP.`))

	compareCmd.Flags().Float64P("max-evalue", "e", 1e-15,
		formatFlagUsage(`Maximum evalue of a HSP segment.`))

	compareCmd.Flags().BoolP("debug", "", false,
		formatFlagUsage(`Print debug information, including a progress bar. (recommended when searching with one query).`))

	compareCmd.SetUsageTemplate(usageTemplate("-d <index path> [query.fasta[.gz] ...] [-o result.tsv[.gz]]"))

	// ani-af related filtering

	compareCmd.Flags().IntP("kmer-scale", "", 4,
		formatFlagUsage(`Using 1/scale of k-mers for seeding (default mode) or fragment comparison (OrthoANI mode). Available values: 2, 4, 8.`))

}
