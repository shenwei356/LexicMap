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
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"runtime"
	"strings"
	"sync"
	"time"

	"github.com/pkg/errors"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
	"github.com/vbauerster/mpb/v8"
	"github.com/vbauerster/mpb/v8/decor"
	"gonum.org/v1/gonum/stat/combin"
)

var compareCmd = &cobra.Command{
	Use:   "compare",
	Short: "Compare genome pairs and compute ANI and AF",
	Long: `Compare genome pairs and compute ANI and AF

Input:
  - Option 1:
    Two or more FASTA files. All combinations of 2 genomes will be compared.
  - Option 2: 
    Tab-delimited file(s), with genome IDs in the first two columns.
    The file can be the output of 'lexicmap genome pair', and the flag '-H' is needed
    to skip the header line.

Output format:
  Tab-delimited format with 11 columns.

    1.  genome1,  Genome 1.
    2.  genome2,  Genome 2.
    3.  tANI,     Total Average nucleotide identity, calculated by dividing the total matched
                  bases in pairwise alignments (genome 1 vs genome 2 and genome 2 vs genome 1)
                  by the total genome sizes of two genomes.
    4.  ANI1,     Average nucleotide identity when aligning genome 1 to genome 2.
    5.  ANI2,     Average nucleotide identity when aligning genome 2 to genome 1.
    6.  AF1,      Align fraction of genome 1.
    7.  AF2,      Align fraction of genome 2.
    8.  contigs1, Number of contigs in genome 1.
    9.  size1,    Size of the genome 1.
    10. contigs2, Number of contigs in genome 2.
    11. size2,    Size of the genome 1.

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
		// if dbDir == "" {
		// 	checkError(fmt.Errorf("flag -d/--index needed"))
		// }

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
		minAF := getFlagNonNegativeFloat64(cmd, "min-af") / 100

		// minPrefix := getFlagPositiveInt(cmd, "seed-min-prefix")
		// if minPrefix > 32 || minPrefix < 5 {
		// 	checkError(fmt.Errorf("the value of flag -p/--seed-min-prefix (%d) should be in the range of [5, 32]", minPrefix))
		// }
		minPrefix := 21 // not used in this command

		maxSubjectGenomeSize := getFlagNonNegativeInt(cmd, "max-genome-size") * 1000 * 1000

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

		// inMemorySearch := getFlagBool(cmd, "load-whole-seeds")
		inMemorySearch := false

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

		hasHeaderLine := getFlagBool(cmd, "skip-header-line")

		debug := getFlagBool(cmd, "debug")

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

		loadIndex := dbDir != ""

		if !loadIndex && len(files) < 2 {
			checkError(fmt.Errorf("if no LexicMap index given via '-d', >=2 input sequence files are required"))
		}

		if outputLog {
			log.Info()
			if loadIndex {
				log.Infof("loading index: %s", dbDir)
			}
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

			NoIndex: !loadIndex,
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

		if outputLog && loadIndex {
			log.Infof("index loaded in %s", time.Since(timeStart))
			log.Info()
		}

		if outputLog {
			log.Infof("compare with %d threads...", opt.NumCPUs)
			log.Infof("  minimum query coverage per HSP: %.2f%%", minQcovChain)
			log.Infof("  minimum base identity in a HSP segment: %.2f%%", minIdent)
			log.Infof("  maximum evalue: %.2e", maxEvalue)

			if gc {
				log.Infof("  force garbage collection for every %d genome pairs", gcInterval)
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

		gcIntervalMinus1 := gcInterval - 1

		fmt.Fprintf(outfh, "genome1\tgenome2\ttANI\tANI1\tANI2\tAF1\tAF2\tcontigs1\tsize1\tcontigs2\tsize2\n")

		printResult := func(q *GPair) {
			total++

			if verbose {
				if (total < 15 && total&3 == 0) || total&15 == 0 {
					speed = float64(total) / time.Since(timeStart1).Minutes()
					fmt.Fprintf(os.Stderr, "processed genome pairs: %d, speed: %.3f pairs per minute\r", total, speed)
				}
			}

			// ------------------------------

			g1, g2 := q.g1, q.g2
			gr1, gr2 := (*g1.result)[0], (*g2.result)[0]

			fmt.Fprintf(outfh, "%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%d\t%d\n",
				g1.id, g2.id,
				float64(gr1.AlignedMatches+gr2.AlignedMatches)/float64(g1.genomeSize+g2.genomeSize)*100,
				gr1.ANI*100, gr2.ANI*100,
				gr1.AFq*100, gr2.AFq*100,
				len(g1.seqs), g1.genomeSize,
				len(g2.seqs), g2.genomeSize,
			)

			// ------------------------------

			outfh.Flush()
			RecycleGPair(q)

			if gc && total&gcIntervalMinus1 == 0 {
				runtime.GC()
			}
		}

		ch := make(chan *GPair, opt.NumCPUs)
		done := make(chan int)
		go func() {
			for r := range ch {
				printResult(r)
			}

			done <- 1
		}()

		// ------- input -------

		var wg sync.WaitGroup
		tokens := make(chan int, opt.NumCPUs/2) // cause each pair use 2 threads
		fcpus := float64(idx.opt.NumCPUs / 2)

		var gname2idx map[string]*[]uint64

		var pairs []string
		if loadIndex {
			// genomes.map file for mapping index to genome id
			gname2idx, err = readGenomeMapName2Idx(filepath.Join(dbDir, FileGenomeIndex))
			if err != nil {
				checkError(fmt.Errorf("failed to read genomes index mapping file: %s", err))
			}

			for _, file := range files {
				pairs1, err := readPairs(file, hasHeaderLine)
				if err != nil {
					checkError(fmt.Errorf("failed to parse input, two-column tab-delimited input needed: %s\n", err))
				}
				if pairs == nil {
					pairs = pairs1
				} else {
					pairs = append(pairs, pairs1...)
				}
			}
			if len(pairs)>>1 < 1 {
				checkError(fmt.Errorf("no valid genome pairs given from %d pair list file(s)", len(files)))
			}
		} else {
			combs := combin.Combinations(len(files), 2)
			pairs = make([]string, 0, len(combs)<<1)
			for _, pair := range combs {
				pairs = append(pairs, files[pair[0]])
				pairs = append(pairs, files[pair[1]])
			}
		}

		nPairs := len(pairs) >> 1

		if outputLog {
			log.Info()
			log.Infof("%d genome pairs loaded from %d file(s)", nPairs, len(files))
		}

		// -----------------------------------------------------------
		// process bar
		var pbs *mpb.Progress
		var bar *mpb.Bar
		var chDuration chan time.Duration
		var doneDuration chan int
		if debug {
			pbs = mpb.New(mpb.WithWidth(40), mpb.WithOutput(os.Stderr))
			bar = pbs.AddBar(int64(nPairs),
				mpb.PrependDecorators(
					decor.Name("compared genome pairs: ", decor.WC{W: len("compared genome pairs: "), C: decor.DindentRight}),
					decor.Name("", decor.WCSyncSpaceR),
					decor.CountersNoUnit("%d / %d", decor.WCSyncWidth),
				),
				mpb.AppendDecorators(
					decor.Name("ETA: ", decor.WC{W: len("ETA: ")}),
					decor.EwmaETA(decor.ET_STYLE_GO, 1024),
					decor.OnComplete(decor.Name(""), ". done"),
				),
			)

			chDuration = make(chan time.Duration, idx.opt.NumCPUs)
			doneDuration = make(chan int)
			go func() {
				for t := range chDuration {
					bar.EwmaIncrBy(1, t)
				}
				doneDuration <- 1
			}()
		}

		// -----------------------------------------------------------

		for i := 0; i < len(pairs); i += 2 {
			wg.Add(1)
			tokens <- 1

			go func(genome1, genome2 string) {
				timeStart := time.Now()
				defer func() {
					wg.Done()
					<-tokens

					if debug {
						chDuration <- time.Duration(float64(time.Since(timeStart)) / fcpus)
					}
				}()

				q := poolGPair.Get().(*GPair)

				var _wg sync.WaitGroup

				// read genome sequences
				_wg.Add(2)
				go func() {
					if loadIndex {
						batchIDAndRefIDs, ok := gname2idx[genome1]
						if !ok {
							checkError(fmt.Errorf("reference name not found: %s, you might need use '-H' to skip the header line", genome1))
						}
						q.g1, err = idx.ReadGenome(batchIDAndRefIDs)
						checkError(err)
						q.g1.id = append(q.g1.id, []byte(genome1)...)
					} else {
						q.g1, err = ReadGenomeFromFile(genome1, reRefName)
						checkError(err)
					}

					_wg.Done()
				}()
				go func() {
					if loadIndex {
						batchIDAndRefIDs, ok := gname2idx[genome2]
						if !ok {
							checkError(fmt.Errorf("reference name not found: %s, you might need use '-H' to skip the header line", genome2))
						}
						q.g2, err = idx.ReadGenome(batchIDAndRefIDs)
						checkError(err)
						q.g2.id = append(q.g2.id, []byte(genome2)...)
					} else {
						q.g2, err = ReadGenomeFromFile(genome2, reRefName)
						checkError(err)
					}

					_wg.Done()
				}()
				_wg.Wait()

				// compre genomes
				_wg.Add(2)
				go func() {
					err = idx.CompareTwoGenomes(q.g1, q.g2, fragSize, minAlignLen, minAF)
					if err != nil {
						checkError(fmt.Errorf("compare %s to %s: %s", q.g1.id, q.g2.id, err))
					}

					_wg.Done()
				}()
				go func() {
					err = idx.CompareTwoGenomes(q.g2, q.g1, fragSize, minAlignLen, minAF)
					if err != nil {
						checkError(fmt.Errorf("compare %s to %s: %s", q.g2.id, q.g1.id, err))
					}

					_wg.Done()
				}()
				_wg.Wait()

				ch <- q

			}(pairs[i], pairs[i+1])
		}

		wg.Wait()

		// -----------------------------------------------------------
		if debug {
			close(chDuration)
			<-doneDuration
			pbs.Wait()
		}

		close(ch)
		<-done

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
				log.Infof("results saved to: %s", outFile)
			}
		}

		checkError(idx.Close())
	},
}

func init() {
	genomeCmd.AddCommand(compareCmd)

	compareCmd.Flags().BoolP("skip-header-line", "H", false,
		formatFlagUsage(`Skip the header line in the input file.`))

	// general flags

	compareCmd.Flags().StringP("index", "d", "",
		formatFlagUsage(`Index directory created by "lexicmap index".`))

	compareCmd.Flags().StringP("out-file", "o", "-",
		formatFlagUsage(`Out file, supports a ".gz" suffix ("-" for stdout).`))

	compareCmd.Flags().IntP("max-open-files", "", 1024,
		formatFlagUsage(`Maximum opened files. It mainly affects candidate genome extraction. Increase this value if you have hundreds of genome batches or have multiple queries, and do not forgot to set a bigger "ulimit -n" in shell if the value is > 1024.`))

	compareCmd.Flags().IntP("gc-interval", "", 128,
		formatFlagUsage(`Force garbage collection every N queries (0 for disable). The value can't be too small.`))

	// genome filtering/screening

	compareCmd.Flags().StringP("ref-name-regexp", "", `(?i)(.+)\.(f[aq](st[aq])?|fna)(\.gz|\.xz|\.zst|\.bz2)?$`,
		formatFlagUsage(`Regular expression (must contains "(" and ")") for extracting the reference name from the input filename. Attention: use double quotation marks for patterns containing commas, e.g., -p '"A{2,}"'.`))

	// compareCmd.Flags().IntP("seed-min-prefix", "p", 21,
	// 	formatFlagUsage(`Minimum prefix length of matched seeds in the genome filtering phase.`))

	compareCmd.Flags().IntP("frag-size", "", 1020,
		formatFlagUsage(`The size of non-overlap fragments cut for ANI computation.`))
	compareCmd.Flags().IntP("min-frag-size", "", 100,
		formatFlagUsage(`The minimum length of fragments in the end of a sequence during cutting fragments.`))

	compareCmd.Flags().IntP("max-genome-size", "", 20,
		formatFlagUsage(`Maximum size of genomes to be considered (in MB) when reading genome from an index.`))

	// alignment
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

	compareCmd.SetUsageTemplate(usageTemplate(""))

	// ani-af related filtering

	compareCmd.Flags().Float64P("min-af", "", 15.0,
		formatFlagUsage(`Only output results where one genome has aligned fraction > than this value (percentage).`))

	compareCmd.Flags().IntP("kmer-scale", "", 4,
		formatFlagUsage(`Using 1/scale of k-mers for seeding (default mode) or fragment comparison (OrthoANI mode). Available values: 2, 4, 8.`))

}

func readPairs(file string, hasHeaderLine bool) ([]string, error) {
	fh, err := xopen.Ropen(file)
	if err != nil {
		return nil, err
	}

	list := make([]string, 0, 4096)

	items := make([]string, 2)
	scanner := bufio.NewScanner(fh)
	var line string
	var i int
	for scanner.Scan() {
		line = strings.TrimRight(scanner.Text(), "\r\n")
		if line == "" {
			continue
		}

		if hasHeaderLine {
			hasHeaderLine = false
			continue
		}

		stringSplitNByByte(line, '\t', 2, &items)
		if len(items) < 2 {
			continue
		}

		if i = strings.IndexByte(items[1], '\t'); i > 0 { // a\tb\tc
			list = append(list, items[0])
			list = append(list, items[1][:i])
		} else if i == 0 { // a\t\tc
			continue
		} else { // a\tb
			list = append(list, items[0])
			list = append(list, items[1])
		}
	}
	if err = scanner.Err(); err != nil {
		return nil, err
	}

	return list, fh.Close()
}

type GPair struct {
	g1, g2 *GQuery
}

var poolGPair = &sync.Pool{
	New: func() interface{} {
		return &GPair{}
	},
}

func (q *GPair) Reset() {
	q.g1 = nil
	q.g2 = nil
}

func RecycleGPair(q *GPair) {
	RecycleGQuery(q.g1)
	RecycleGQuery(q.g2)
	q.Reset()
	poolGPair.Put(q)
}
