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
	"container/heap"
	"fmt"
	"io"
	"math"
	"math/rand"
	"os"
	"path/filepath"
	"regexp"
	"sort"
	"strings"
	"sync"
	"time"

	"github.com/pkg/errors"
	"github.com/rdleal/intervalst/interval"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/kmers"
	"github.com/shenwei356/lexichash"
	"github.com/shenwei356/lexichash/iterator"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
	"github.com/twotwotwo/sorts"
	"github.com/twotwotwo/sorts/sortutil"
	"github.com/vbauerster/mpb/v8"
	"github.com/vbauerster/mpb/v8/decor"
)

var geneMasksCmd = &cobra.Command{
	Use:   "gen-masks",
	Short: "Generate masks from the top N largest genomes",
	Long: `Generate masks from the top N largest genomes

How:

    |ATTATAACGCCACGGGGAGCCGCGGGGTTTC One k-bp mask
    |--------========_______________
    |
    |-------- Prefixes for covering all possible P-bp DNA.
    |         The length is the largest number for 4^P <= #masks
    |
    |--------======== Extend prefixes, chosen from the most frequent extended prefixes
    |                 of which the prefix-derived k-mers do not overlap masked k-mers.
    |
    |                _______________ Randomly generated bases

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

		// ---------------------------------------------------------------
		// input

		var err error

		inDir := getFlagString(cmd, "in-dir")
		skipFileCheck := getFlagBool(cmd, "skip-file-check")

		readFromDir := inDir != ""
		if readFromDir {
			var isDir bool
			isDir, err = pathutil.IsDir(inDir)
			if err != nil {
				checkError(errors.Wrapf(err, "checking -I/--in-dir"))
			}
			if !isDir {
				checkError(fmt.Errorf("value of -I/--in-dir should be a directory: %s", inDir))
			}
		}

		reFileStr := getFlagString(cmd, "file-regexp")
		var reFile *regexp.Regexp
		if reFileStr != "" {
			if !reIgnoreCase.MatchString(reFileStr) {
				reFileStr = reIgnoreCaseStr + reFileStr
			}
			reFile, err = regexp.Compile(reFileStr)
			checkError(errors.Wrapf(err, "failed to parse regular expression for matching file: %s", reFileStr))
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

		reSeqNameStrs := getFlagStringSlice(cmd, "seq-name-filter")
		reSeqNames := make([]*regexp.Regexp, 0, len(reSeqNameStrs))
		for _, kw := range reSeqNameStrs {
			if !reIgnoreCase.MatchString(kw) {
				kw = reIgnoreCaseStr + kw
			}
			re, err := regexp.Compile(kw)
			if err != nil {
				checkError(errors.Wrapf(err, "failed to parse regular expression for matching sequence header: %s", kw))
			}
			reSeqNames = append(reSeqNames, re)
		}

		// ---------------------------------------------------------------
		// basic flags

		maxGenomeSize := getFlagPositiveInt(cmd, "max-genome")
		fileBigGenomes := getFlagString(cmd, "big-genomes")
		topN := getFlagPositiveInt(cmd, "top-n")
		_prefix := getFlagPositiveInt(cmd, "prefix-ext")

		k := getFlagPositiveInt(cmd, "kmer")
		if k < minK || k > 32 {
			checkError(fmt.Errorf("the value of flag -k/--kmer should be in range of [%d, 32]", minK))
		}

		nMasks := getFlagPositiveInt(cmd, "masks")
		seed := getFlagPositiveInt(cmd, "rand-seed")

		outFile := getFlagString(cmd, "out-file")
		seedPosFile := getFlagString(cmd, "seed-pos")
		if seedPosFile == "-" {
			checkError(fmt.Errorf(`stdin ("-") not supported for the flag --seed-pos`))
		}

		// contigInterval := getFlagPositiveInt(cmd, "contig-interval")
		// if contigInterval <= k-1 {
		// 	checkError(fmt.Errorf("the value of --contig-interval should be >= k-1"))
		// }

		bopt := &IndexBuildingOptions{
			// general
			NumCPUs:      opt.NumCPUs,
			Verbose:      opt.Verbose,
			Log2File:     opt.Log2File,
			Force:        false,
			MaxOpenFiles: 512,

			// skip extremely large genomes
			MaxGenomeSize: maxGenomeSize,
			BigGenomeFile: fileBigGenomes,

			// LexicHash
			K:        k,
			Masks:    nMasks,
			RandSeed: int64(seed),

			// generate masks
			TopN:      topN,
			PrefixExt: _prefix,

			// genome
			ReRefName:    reRefName,
			ReSeqExclude: reSeqNames,

			ContigInterval: 1000,
		}

		// ---------------------------------------------------------------
		// checking input files

		if opt.Verbose || opt.Log2File {
			log.Infof("LexicMap v%s", VERSION)
			log.Info("  https://github.com/shenwei356/LexicMap")
			log.Info()

			log.Info("checking input files ...")
		}

		var files []string
		if readFromDir {
			files, err = getFileListFromDir(inDir, reFile, opt.NumCPUs)
			if err != nil {
				checkError(errors.Wrapf(err, "walking dir: %s", inDir))
			}
			if len(files) == 0 {
				log.Warningf("  no files matching regular expression: %s", reFileStr)
			}
		} else {
			files = getFileListFromArgsAndFile(cmd, args, !skipFileCheck, "infile-list", !skipFileCheck)
			if opt.Verbose || opt.Log2File {
				if len(files) == 1 && isStdin(files[0]) {
					log.Info("  no files given, reading from stdin")
				}
			}
		}
		if len(files) < 1 {
			checkError(fmt.Errorf("FASTA/Q files needed"))
		} else if opt.Verbose || opt.Log2File {
			log.Infof("  %d input file(s) given", len(files))
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

		masks, skippedFiles, err := GenerateMasks(files, bopt, seedPosFile)
		checkError(err)

		for _, mask := range masks {
			fmt.Fprintf(outfh, "%s\n", kmers.MustDecode(mask, bopt.K))
		}

		if fileBigGenomes != "" {
			outfh2, err := os.Create(fileBigGenomes)
			if err != nil {
				checkError(fmt.Errorf("failed to write file: %s", fileBigGenomes))
				return
			}
			for _, file := range skippedFiles {
				fmt.Fprintf(outfh2, "%s\n", file)
			}
			outfh2.Close()
		}
	},
}

func init() {
	utilsCmd.AddCommand(geneMasksCmd)

	// -----------------------------  input  -----------------------------

	geneMasksCmd.Flags().StringP("in-dir", "I", "",
		formatFlagUsage(`Directory containing FASTA/Q files. Directory symlinks are followed.`))

	geneMasksCmd.Flags().StringP("file-regexp", "r", `\.(f[aq](st[aq])?|fna)(.gz)?$`,
		formatFlagUsage(`Regular expression for matching sequence files in -I/--in-dir, case ignored.`))

	geneMasksCmd.Flags().StringP("ref-name-regexp", "N", `(?i)(.+)\.(f[aq](st[aq])?|fna)(.gz)?$`,
		formatFlagUsage(`Regular expression (must contains "(" and ")") for extracting the reference name from the filename.`))

	geneMasksCmd.Flags().StringSliceP("seq-name-filter", "B", []string{},
		formatFlagUsage(`List of regular expressions for filtering out sequences by header/name, case ignored.`))

	geneMasksCmd.Flags().BoolP("skip-file-check", "S", false,
		formatFlagUsage(`Skip input file checking when given files or a file list.`))

	geneMasksCmd.Flags().IntP("max-genome", "g", 15000000,
		formatFlagUsage(`Maximum genome size. Extremely large genomes (non-isolate assemblies) will be skipped.`))

	// geneMasksCmd.Flags().IntP("contig-interval", "", 1000,
	// 	formatFlagUsage(`Length of interval (N's) between contigs in a genome (>=K).`))

	// -----------------------------  output   -----------------------------

	geneMasksCmd.Flags().StringP("out-file", "o", "-",
		formatFlagUsage(`Out file of generated masks. The ".gz" suffix is not recommended. ("-" for stdout).`))

	geneMasksCmd.Flags().StringP("big-genomes", "G", "",
		formatFlagUsage(`Out file of skipped files with genomes >= -G/--max-genome`))

	geneMasksCmd.Flags().StringP("seed-pos", "D", "",
		formatFlagUsage(`Out file of seed postions and distances, supports and recommends a ".gz" suffix.`))

	// -----------------------------  lexichash   -----------------------------

	geneMasksCmd.Flags().IntP("kmer", "k", 31,
		formatFlagUsage(`Maximum k-mer size. K needs to be <= 32.`))

	geneMasksCmd.Flags().IntP("masks", "m", 40000,
		formatFlagUsage(`Number of masks.`))

	geneMasksCmd.Flags().IntP("rand-seed", "s", 1,
		formatFlagUsage(`Rand seed for generating random masks.`))

	// -----------------------------  generate mask   -----------------------------

	geneMasksCmd.Flags().IntP("top-n", "n", 20,
		formatFlagUsage(`Select the top N largest genomes for generating masks.`))

	geneMasksCmd.Flags().IntP("prefix-ext", "P", 8,
		formatFlagUsage(`Extension length of prefixes, higher values -> smaller maximum seed distances.`))

	geneMasksCmd.SetUsageTemplate(usageTemplate("[-k <k>] [-n <masks>] [-n <top-n>] [-D <seeds.tsv.gz>] [-o masks.txt] { -I <seqs dir> | -X <file list>}"))
}

// GenerateMasks generates Masks from the top N largest genomes.
func GenerateMasks(files []string, opt *IndexBuildingOptions, outFile string) ([]uint64, []string, error) {
	maxGenomeSize := opt.MaxGenomeSize
	topN := opt.TopN
	_lenPrefix := opt.PrefixExt

	var outfh *bufio.Writer
	if outFile != "" { // output seed locations and distances
		var gw io.WriteCloser
		var w *os.File
		var err error
		outfh, gw, w, err = outStream(outFile, strings.HasSuffix(outFile, ".gz"), 5)
		checkError(err)
		defer func() {
			outfh.Flush()
			if gw != nil {
				gw.Close()
			}
			w.Close()
		}()
	}

	var timeStart time.Time
	if opt.Verbose || opt.Log2File {
		timeStart = time.Now()
		log.Info()
		log.Infof("generating masks from the top %d out of %d genomes...", topN, len(files))
		defer func() {
			log.Info()
			log.Infof("  finished generating masks in: %s", time.Since(timeStart))
		}()
	}

	// --------------------------------------------------------------------
	// find top n files that contain the largest genome.

	var filesTop []string
	file2gsizes := &File2GSizes{}

	var maxGenome int

	// find top n genome files with biggest genome size
	if opt.Verbose || opt.Log2File {
		log.Info()
		log.Infof("  checking genomes sizes of %d files...", len(files))
	}

	// process bar
	var pbs *mpb.Progress
	var bar *mpb.Bar
	var chDuration chan time.Duration
	var doneDuration chan int
	if opt.Verbose {
		pbs = mpb.New(mpb.WithWidth(40), mpb.WithOutput(os.Stderr))
		bar = pbs.AddBar(int64(len(files)),
			mpb.PrependDecorators(
				decor.Name("processed files: ", decor.WC{W: len("processed files: "), C: decor.DindentRight}),
				decor.Name("", decor.WCSyncSpaceR),
				decor.CountersNoUnit("%d / %d", decor.WCSyncWidth),
			),
			mpb.AppendDecorators(
				decor.Name("ETA: ", decor.WC{W: len("ETA: ")}),
				decor.EwmaETA(decor.ET_STYLE_GO, 20),
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

	// genome size
	chGS := make(chan File2GSize, opt.NumCPUs)
	doneGS := make(chan int)
	heap.Init(file2gsizes)
	go func() {
		for gs := range chGS {
			heap.Push(file2gsizes, File2GSize{File: gs.File, Size: gs.Size})
			if len(*file2gsizes) == topN {
				heap.Init(file2gsizes)
			} else if len(*file2gsizes) > topN {
				heap.Pop(file2gsizes)
			}
		}
		doneGS <- 1
	}()

	// skipped extreme big genomes
	chSkippedFiles := make(chan string, 8)
	doneSkip := make(chan int)
	skippedFiles := make([]string, 0, 8)
	go func() {
		for file := range chSkippedFiles {
			skippedFiles = append(skippedFiles, file)
		}
		doneSkip <- 1
	}()

	// read seq file and count bases
	var wg sync.WaitGroup                 // ensure all jobs done
	tokens := make(chan int, opt.NumCPUs) // control the max concurrency number
	k := opt.K
	filterNames := len(opt.ReSeqExclude) > 0
	threadsFloat := float64(opt.NumCPUs)
	for _, file := range files {
		tokens <- 1
		wg.Add(1)

		go func(file string) {
			defer func() {
				wg.Done()
				<-tokens
			}()
			startTime := time.Now()

			fastxReader, err := fastx.NewReader(nil, file, "")
			if err != nil {
				checkError(fmt.Errorf("failed to read seq file: %s", err))
			}
			defer fastxReader.Close()

			var record *fastx.Record

			var ignoreSeq bool
			var re *regexp.Regexp

			var genomeSize int
			var i int = 0
			for {
				record, err = fastxReader.Read()
				if err != nil {
					if err == io.EOF {
						break
					}
					checkError(fmt.Errorf("read seq %d in %s: %s", i, file, err))
					break
				}

				// filter out sequences shorter than k
				if len(record.Seq.Seq) < k {
					continue
				}

				// filter out sequences with names in the blast list
				if filterNames {
					ignoreSeq = false
					for _, re = range opt.ReSeqExclude {
						if re.Match(record.Name) {
							ignoreSeq = true
							break
						}
					}
					if ignoreSeq {
						continue
					}
				}

				genomeSize += len(record.Seq.Seq)
			}

			if genomeSize == 0 {
				log.Warningf("skipping %s: no valid sequences", file)
				log.Info()
				if opt.Verbose {
					chDuration <- time.Microsecond // important, or the progress bar will get hung
				}
				return
			}
			if genomeSize > maxGenomeSize {
				if opt.Verbose || opt.Log2File {
					log.Warningf("skipping big genome (%d bp): %s", genomeSize, file)
					if !opt.Log2File {
						log.Info()
					}
				}
				chSkippedFiles <- file
				if opt.Verbose {
					chDuration <- time.Duration(float64(time.Since(startTime)) / threadsFloat)
				}
				return
			}

			chGS <- File2GSize{File: file, Size: genomeSize}

			if opt.Verbose {
				chDuration <- time.Duration(float64(time.Since(startTime)) / threadsFloat)
			}
		}(file)
	}

	wg.Wait()
	close(chGS)
	<-doneGS
	close(chSkippedFiles)
	<-doneSkip

	// process bar
	if opt.Verbose {
		close(chDuration)
		<-doneDuration
		pbs.Wait()
	}

	// topN might be > the number of available files
	if topN > len(files)-len(skippedFiles) {
		topN = len(files) - len(skippedFiles)
	}

	filesTop = make([]string, 0, topN)
	sort.Sort(*file2gsizes) // make sure it's sorted
	maxGenome = (*file2gsizes)[len(*file2gsizes)-1].Size

	if opt.Verbose || opt.Log2File {
		log.Infof("    %d genomes longer than %d are filtered out", len(skippedFiles), maxGenomeSize)
		log.Infof("    genome size range in the top %d files: [%d, %d]",
			topN, (*file2gsizes)[0].Size, maxGenome)
	}

	for _, gs := range *file2gsizes {
		filesTop = append(filesTop, gs.File)
	}

	// --------------------------------------------------------------------
	// count k-mers from the top n files

	if opt.Verbose || opt.Log2File {
		log.Info()
		log.Infof("  collecting k-mers from %d files...", len(filesTop))
	}

	lenPrefix := 1
	for 1<<(lenPrefix<<1) <= opt.Masks {
		lenPrefix++
	}
	lenPrefix-- // 7 for 20,000 masks

	nPrefix := int(math.Pow(4, float64(lenPrefix)))

	// prefix (list, 16384) -> refs (list, #refs) -> kmers (map, might >10k) -> location (list, small)
	data := make([][]map[uint64]*[]int32, nPrefix)
	for i := range data {
		data[i] = make([]map[uint64]*[]int32, topN)
		for j := range data[i] {
			data[i][j] = make(map[uint64]*[]int32, 8)
		}
	}

	var m map[uint64]*[]int32
	var prefix int
	k8 := uint8(opt.K)
	p8 := uint8(lenPrefix)
	genomeIDs := make([]string, topN)

	nnn := bytes.Repeat([]byte{'N'}, opt.ContigInterval)
	reRefName := opt.ReRefName
	extractRefName := reRefName != nil
	_seq := make([]byte, 0, 10<<20)
	var iter *iterator.Iterator

	var i int

	var record *fastx.Record

	var ignoreSeq bool
	var re *regexp.Regexp
	var baseFile string

	var kmer, kmerRC uint64
	var ok bool
	var j int
	var locs *[]int32
	// store k-mer-location pairs in a list, sort them and handle them.
	// sorted k-mers have the same prefix, so the k-mer maps
	// (map[uint64]*[]int32) would be next to each other.
	// This increases the data locality and improves the speed of later processes
	// which need frequently iterate these maps.
	_kmers := make([][2]uint64, 40<<20)
	var _kmers2 [][2]uint64
	var kmer2loc [2]uint64
	var loc int32
	var iK, lenKmers int

	for iG, file := range filesTop {
		_seq = _seq[:0]

		baseFile = filepath.Base(file)

		fastxReader, err := fastx.NewReader(nil, file, "")
		if err != nil {
			checkError(fmt.Errorf("failed to read seq file: %s", err))
		}

		// --------------------------------
		// concatenate contigs

		i = 0
		for {
			record, err = fastxReader.Read()
			if err != nil {
				if err == io.EOF {
					break
				}
				checkError(fmt.Errorf("read seq %d in %s: %s", i, file, err))
				break
			}

			// filter out sequences shorter than k
			if len(record.Seq.Seq) < k {
				continue
			}

			// filter out sequences with names in the blast list
			if filterNames {
				ignoreSeq = false
				for _, re = range opt.ReSeqExclude {
					if re.Match(record.Name) {
						ignoreSeq = true
						break
					}
				}
				if ignoreSeq {
					continue
				}
			}

			if i > 0 { // add N's between two contigs
				_seq = append(_seq, nnn...)
			}
			_seq = append(_seq, record.Seq.Seq...)

			i++
		}

		if len(_seq) == 0 {
			fastxReader.Close()

			log.Warningf("skipping %s: no valid sequences", file)
			log.Info()
			continue
		}

		// genome id
		var genomeID string
		if extractRefName {
			if reRefName.MatchString(baseFile) {
				genomeID = reRefName.FindAllStringSubmatch(baseFile, 1)[0][1]
			} else {
				genomeID, _, _ = filepathTrimExtension(baseFile, nil)
			}
		} else {
			genomeID, _, _ = filepathTrimExtension(baseFile, nil)
		}
		genomeIDs[iG] = genomeID

		// --------------------------------
		// generate k-mers

		iter, err = iterator.NewKmerIterator(_seq, k)
		if err != nil {
			checkError(fmt.Errorf("count kmer for %s: %s", file, err))
		}

		lenKmers = len(_kmers)
		iK = 0
		for {
			kmer, kmerRC, ok, _ = iter.NextKmer()
			if !ok {
				break
			}
			if kmer == 0 { // all bases are A's or N's.
				continue
			}

			j = iter.Index()

			if iK < lenKmers { //
				_kmers[iK][0] = kmer
				_kmers[iK][1] = uint64(j)
			} else {
				_kmers = append(_kmers, [2]uint64{kmer, uint64(j)})
			}
			iK++

			if iK < lenKmers {
				_kmers[iK][0] = kmerRC
				_kmers[iK][1] = uint64(j)
			} else {
				_kmers = append(_kmers, [2]uint64{kmerRC, uint64(j)})
			}
			iK++
		}

		fastxReader.Close()

		_kmers2 = _kmers[:iK] // only used data
		// sort.Slice(_kmers2, func(i, j int) bool { return _kmers2[i][0] < _kmers2[j][0] })
		sorts.Quicksort(Kmer2Locs(_kmers2))

		for _, kmer2loc = range _kmers2 {
			kmer = kmer2loc[0]
			loc = int32(kmer2loc[1])

			prefix = int(util.KmerPrefix(kmer, k8, p8))
			m = data[prefix][iG]

			if locs, ok = m[kmer]; !ok {
				tmp := []int32{loc}
				m[kmer] = &tmp
			} else {
				*locs = append(*locs, loc)
			}
			// fmt.Printf("%s, %d, %d\n", kmers.Decode(uint64(prefix), lenPrefix), iG, len(data[prefix][iG]))
		}

		if opt.Verbose {
			fmt.Fprintf(os.Stderr, "\rprocessed files: %d/%d", iG+1, topN)
		}
	}
	if opt.Verbose {
		fmt.Fprintln(os.Stderr)
	}

	// --------------------------------------------------------------------
	// generate masks

	if opt.Verbose || opt.Log2File {
		log.Info()
		log.Infof("  generating masks...")
	}

	// count prefixes' k-mers and sort in descending order.
	counts := make([][2]int, nPrefix)
	var dRefs []map[uint64]*[]int32 // kmer data of all refs
	var n int
	for i, dRefs = range data {
		n = 0
		for _, m = range dRefs {
			n += len(m)
		}
		counts[i] = [2]int{i, n} // prefix, count
	}
	sort.Slice(counts, func(i, j int) bool { return counts[i][1] < counts[j][1] })

	// ATTATAACGCCACGGGGAGCCGCGGGGTTTC
	// ------- prefix
	//        -------- _prefix
	//                ---------------- will be generated randomly

	// frequency table for bases behind the prefix.

	_nPrefix := int(math.Pow(4, float64(_lenPrefix)))
	var _prefix int
	_p8 := uint8(_lenPrefix)
	lenPrefix8 := uint8(lenPrefix)

	freqs := make([]map[int]interface{}, _nPrefix) // list of _prefix
	for i := range freqs {
		freqs[i] = make(map[int]interface{}, topN)
	}
	_counts := make([][2]int, _nPrefix)
	for i = range _counts {
		_counts[i] = [2]int{}
	}
	var _count [2]int
	_sortFunc := func(i, j int) bool { return _counts[i][1] > _counts[j][1] }
	var _m map[int]interface{}

	r := rand.New(rand.NewSource(opt.RandSeed))
	var _mask uint64 = 1<<(uint64(k-lenPrefix-_lenPrefix)<<1) - 1
	shiftP := uint64(k-lenPrefix) << 1
	_shiftP := uint64(k-lenPrefix-_lenPrefix) << 1
	var prefix64 uint64
	var mask uint64
	var hash, minHash, minKmer uint64
	var minLocs *[]int32
	// var loc int32

	// for intervals
	flank := int32(maxGenome / nPrefix)
	cmpFn := func(x, y int32) int { return int(x - y) }
	itrees := make([]*interval.SearchTree[uint64, int32], topN)
	for i := 0; i < topN; i++ {
		itrees[i] = interval.NewSearchTree[uint64, int32](cmpFn)
	}
	var overlap bool
	k32 := int32(k)

	locations := make([][]int32, topN) // k-mer locations in all genomes
	for i := 0; i < topN; i++ {
		locations[i] = make([]int32, 0, opt.Masks)
	}

	// prefix -> map of k-mers
	if opt.Verbose || opt.Log2File {
		log.Infof("    generating %d masks covering all %d-bp prefixes...", nPrefix, lenPrefix)
	}

	// masks: prefix -> mask
	masks := make([]map[uint64]interface{}, nPrefix)

	// for deduplicating: prefix -> _prefix
	masks2 := make([]map[int]interface{}, nPrefix)

	for j, count := range counts { // for all prefix
		if opt.Verbose {
			fmt.Fprintf(os.Stderr, "\rprocessed prefixes: %d/%d", j+1, nPrefix)
		}

		// a prefix not existing in any genome.
		// maybe it's just because we don't have enough genomes.
		// in a test of top 5 genomes with 7-bp prefix, it did not happen.
		if count[1] == 0 {
			prefix = count[0]
			prefix64 = uint64(prefix)
			_prefix = r.Intn(_nPrefix)
			// generate one random mask with a prefix of "prefix"+"_prefix"
			// mask = prefix + _ prefix + random
			mask = prefix64<<shiftP | uint64(_prefix)<<_shiftP | util.Hash64(r.Uint64())&_mask

			// record it: prefix -> mask
			masks[prefix] = map[uint64]interface{}{mask: struct{}{}}
			masks2[prefix] = map[int]interface{}{_prefix: struct{}{}}

			continue
		}

		prefix = count[0]

		prefix64 = uint64(prefix)

		// -----------------------------------------------------------------
		// extend the prefix according to existing k-mers and generate a mask

		// clear the freqs table
		for i = range freqs {
			clear(freqs[i])
		}
		// count how many genomes have k-mers with this _prefix
		for i, m = range data[prefix] { // i is the genome idx
			for kmer, locs = range m {
				// check if this kmer overlap with existing intervals
				overlap = false
				for _, loc = range *locs {
					if _, ok = itrees[i].AnyIntersection(loc, loc+k32); ok {
						overlap = true
					}
				}
				if overlap {
					continue
				}

				// extract _prefix and count it
				_prefix = int(util.KmerPrefix(util.KmerSuffix(kmer, k8, lenPrefix8), k8-p8, _p8))

				freqs[_prefix][i] = struct{}{}
			}
		}

		// clear the count table
		for _prefix, _m = range freqs {
			if len(_m) == 0 { // a _prefix not existing in any genome, reset the count
				_counts[_prefix][0] = 0
				_counts[_prefix][1] = 0
				continue
			}
			_counts[_prefix][0] = _prefix
			_counts[_prefix][1] = len(_m)
		}
		sort.Slice(_counts, _sortFunc) // sort in descending order of frequencies

		if _counts[0][1] == 0 { // there's no k-mers available, just generate a random one
			_prefix = r.Intn(_nPrefix)
		} else { // randomly choose one of the most frequent _prefixes.
			_prefix = -1
			for _, _count = range _counts {
				if _, ok = masks2[prefix][_count[0]]; !ok { // make sure that prefix+_prefix is unique.
					_prefix = _count[0]
					break
				}
			}
			if _prefix == -1 { // have no choice
				_prefix = _counts[0][0]
			}
		}

		// generate one random mask with a prefix of "prefix"+"_prefix"
		// mask = prefix + _ prefix + random
		mask = prefix64<<shiftP | uint64(_prefix)<<_shiftP | util.Hash64(r.Uint64())&_mask

		// record it: prefix -> mask
		masks[prefix] = map[uint64]interface{}{mask: struct{}{}}
		masks2[prefix] = map[int]interface{}{_prefix: struct{}{}}

		// -----------------------------------------------------------------
		// capture the most similar k-mer in each genome

		for i, m = range data[prefix] { // each genome, i is the genome idx
			if len(m) == 0 { // this genome does not have k-mers starting with the prefix
				continue
			}

			minHash = math.MaxUint64
			for kmer, locs = range m { // each kmer
				hash = mask ^ kmer // lexichash

				if hash < minHash { // hash == minHash would not happen, because they are saved in a map
					minKmer = kmer
					minLocs = locs
					minHash = hash
				}
			}
			for _, loc = range *minLocs {
				// add the region of this k-mer to the interval tree of the genome
				itrees[i].Insert(loc-flank, loc+flank, minKmer)

				// store locations for further report
				locations[i] = append(locations[i], loc)
			}
		}
	}
	if opt.Verbose {
		fmt.Fprintln(os.Stderr)
	}

	// generate left (Masks - nPrefix) masks
	if opt.Masks-nPrefix > 0 {
		leftMasks := opt.Masks - nPrefix // left masks to generate
		if opt.Verbose || opt.Log2File {
			log.Infof("    generating left %d masks...", leftMasks)
		}
		// shuffle all prefixes
		prefixes := make([]int, 0, nPrefix) // only use non-low-complexity prefix
		for i = 0; i < nPrefix; i++ {
			if lexichash.IsLowComplexity(uint64(i), lenPrefix) {
				continue
			}
			prefixes = append(prefixes, i)
		}
		nPrefix2 := len(prefixes)
		r.Shuffle(nPrefix2, func(i, j int) { prefixes[i], prefixes[j] = prefixes[j], prefixes[i] })

		rounds := leftMasks/nPrefix2 + 1 // the round of using the prefixes
		var last int                     // the last element
		j = 0
		for round := 0; round < rounds; round++ {
			if round < rounds-1 { // for previous rounds, all the prefixes are used
				last = nPrefix2
			} else { // for the last round, only a part of prefixes is used
				last = leftMasks % nPrefix2
			}

			for _, prefix = range prefixes[:last] {
				if opt.Verbose {
					j++
					fmt.Fprintf(os.Stderr, "\rprocessed prefixes: %d/%d", j, leftMasks)
				}

				prefix64 = uint64(prefix)

				// -----------------------------------------------------------------
				// extend the prefix according to existing k-mers and generate a mask

				// clear the freqs table
				for i = range freqs {
					clear(freqs[i])
				}
				// count how many genomes have k-mers with this _prefix
				for i, m = range data[prefix] { // i is the genome idx
					for kmer, locs = range m {
						// check if this kmer overlap with existing intervals
						overlap = false
						for _, loc = range *locs {
							if _, ok = itrees[i].AnyIntersection(loc, loc+k32); ok {
								overlap = true
							}
						}
						if overlap {
							continue
						}

						// extract _prefix and count it
						_prefix = int(util.KmerPrefix(util.KmerSuffix(kmer, k8, lenPrefix8), k8-p8, _p8))

						freqs[_prefix][i] = struct{}{}
					}
				}

				// clear the count table
				for _prefix, _m = range freqs {
					if len(_m) == 0 { // a _prefix not existing in any genome, reset the count
						_counts[_prefix][0] = 0
						_counts[_prefix][1] = 0
						continue
					}
					_counts[_prefix][0] = _prefix
					_counts[_prefix][1] = len(_m)
				}
				sort.Slice(_counts, _sortFunc) // sort in descending order of frequencies

				if _counts[0][1] == 0 { // there's no k-mers available, just generate a random one
					_prefix = r.Intn(_nPrefix)
				} else { // randomly choose one of the most frequent _prefixes.
					_prefix = -1
					for _, _count = range _counts {
						if _, ok = masks2[prefix][_count[0]]; !ok { // make sure that prefix+_prefix is unique.
							_prefix = _count[0]
							break
						}
					}
					if _prefix == -1 { // have no choice
						_prefix = _counts[0][0]
					}
				}

				// generate one random mask with a prefix of "prefix"+"_prefix"
				// mask = prefix + _ prefix + random
				mask = prefix64<<shiftP | uint64(_prefix)<<_shiftP | util.Hash64(r.Uint64())&_mask

				// record it
				// it's different here !!!!!!!!!!!!1
				// masks[prefix] = map[uint64]interface{}{mask: struct{}{}}
				masks[prefix][mask] = struct{}{}
				masks2[prefix][_prefix] = struct{}{}

				// -----------------------------------------------------------------
				// capture the most similar k-mer in each genome

				for i, m = range data[prefix] { // each genome, i is the genome idx
					if len(m) == 0 { // this genome does not have k-mers starting with the prefix
						continue
					}
					minHash = math.MaxUint64
					for kmer, locs = range m { // each kmer
						hash = mask ^ kmer // lexichash

						if hash < minHash { // hash == minHash would not happen, because they are saved in a map
							minKmer = kmer
							minLocs = locs
							minHash = hash
						}
					}
					for _, loc = range *minLocs {
						// add the region of this k-mer to the interval tree of the genome
						itrees[i].Insert(loc-flank, loc+flank, minKmer)

						// store locations for further report
						locations[i] = append(locations[i], loc)
					}
				}
			}
		}
		if opt.Verbose {
			fmt.Fprintln(os.Stderr)
		}
	}

	// check the max distance in each genome
	if opt.Verbose || opt.Log2File {
		log.Info()
		log.Infof("  maximum distance between seeds:")
	}

	writeLocs := outFile != ""

	var _loc, dist, maxDist int32
	if writeLocs {
		fmt.Fprintf(outfh, "ref\tpos\tdist\n")
	}
	for i, locs2 := range locations {
		sortutil.Int32s(locs2)
		_loc = 0
		maxDist = 0
		for j, loc = range locs2 {
			dist = loc - _loc

			if writeLocs {
				fmt.Fprintf(outfh, "%s\t%d\t%d\n", genomeIDs[i], loc, dist)
			}

			if dist > maxDist {
				maxDist = dist
			}
			_loc = loc
		}

		if opt.Verbose {
			log.Infof("    %s: %d\n", genomeIDs[i], maxDist)
		}
	}
	if writeLocs {
		if opt.Verbose || opt.Log2File {
			log.Info()
			log.Infof("  seed locations and distances of the %d genomes are saved to %s, ", topN, outFile)
			log.Infof("  you can plot the histogram of seed distances:")
			log.Infof("    csvtk grep -t -f ref -p %s %s | csvtk plot hist -t -f dist -o hist.png", genomeIDs[i], outFile)
		}
	}

	// collect masks to return
	_masks := make([]uint64, 0, opt.Masks)
	for _, m := range masks {
		for kmer = range m {
			_masks = append(_masks, kmer)
		}
	}
	sort.Slice(_masks, func(i, j int) bool { return _masks[i] < _masks[j] })

	return _masks, skippedFiles, nil
}

type File2GSize struct {
	Size int
	File string
} // An IntHeap is a min-heap of ints.
type File2GSizes []File2GSize

func (h File2GSizes) Len() int           { return len(h) }
func (h File2GSizes) Less(i, j int) bool { return h[i].Size < h[j].Size }
func (h File2GSizes) Swap(i, j int)      { h[i], h[j] = h[j], h[i] }

func (h *File2GSizes) Push(x any) {
	// Push and Pop use pointer receivers because they modify the slice's length,
	// not just its contents.
	*h = append(*h, x.(File2GSize))
}

func (h *File2GSizes) Pop() any {
	old := *h
	n := len(old)
	x := old[n-1]
	*h = old[0 : n-1]
	return x
}

type Kmer2Locs [][2]uint64

func (h Kmer2Locs) Len() int           { return len(h) }
func (h Kmer2Locs) Less(i, j int) bool { return h[i][0] < h[j][0] }
func (h Kmer2Locs) Swap(i, j int)      { h[i], h[j] = h[j], h[i] }
