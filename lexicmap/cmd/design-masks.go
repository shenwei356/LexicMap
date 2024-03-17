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
	"container/heap"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"regexp"
	"sort"
	"strings"
	"sync"
	"time"

	"github.com/pkg/errors"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/lexichash/iterator"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
	"github.com/vbauerster/mpb/v8"
	"github.com/vbauerster/mpb/v8/decor"
)

var designMasksCmd = &cobra.Command{
	Use:   "design-masks",
	Short: "Design-masks from genomes",
	Long: `Design-masks from genomes

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		seq.ValidateSeq = false

		// ------------------- input -------------------------

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

		// ------------------------------------

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

		topN := getFlagPositiveInt(cmd, "top-n")
		// ---------------------------------------------------------------

		outFile := getFlagString(cmd, "out-file")

		k := getFlagPositiveInt(cmd, "kmer")
		if k < minK || k > 32 {
			checkError(fmt.Errorf("the value of flag -k/--kmer should be in range of [%d, 32]", minK))
		}

		nMasks := getFlagPositiveInt(cmd, "masks")
		lcPrefix := getFlagNonNegativeInt(cmd, "prefix")
		seed := getFlagPositiveInt(cmd, "seed")

		bopt := &IndexBuildingOptions{
			// general
			NumCPUs:      opt.NumCPUs,
			Verbose:      opt.Verbose,
			Log2File:     opt.Log2File,
			Force:        false,
			MaxOpenFiles: 512,

			// LexicHash
			K:                k,
			Masks:            nMasks,
			RandSeed:         int64(seed),
			PrefixForCheckLC: lcPrefix,

			// genome
			ReRefName:    reRefName,
			ReSeqExclude: reSeqNames,
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

		DesignMasks(files, bopt, topN)
	},
}

func init() {
	utilsCmd.AddCommand(designMasksCmd)

	// -----------------------------  input  -----------------------------

	designMasksCmd.Flags().StringP("in-dir", "I", "",
		formatFlagUsage(`Directory containing FASTA/Q files. Directory symlinks are followed.`))

	designMasksCmd.Flags().StringP("file-regexp", "r", `\.(f[aq](st[aq])?|fna)(.gz)?$`,
		formatFlagUsage(`Regular expression for matching sequence files in -I/--in-dir, case ignored.`))

	designMasksCmd.Flags().StringP("ref-name-regexp", "N", `(?i)(.+)\.(f[aq](st[aq])?|fna)(.gz)?$`,
		formatFlagUsage(`Regular expression (must contains "(" and ")") for extracting the reference name from the filename.`))

	designMasksCmd.Flags().StringSliceP("seq-name-filter", "B", []string{},
		formatFlagUsage(`List of regular expressions for filtering out sequences by header/name, case ignored.`))

	designMasksCmd.Flags().BoolP("skip-file-check", "S", false,
		formatFlagUsage(`skip input file checking when given files or a file list.`))

	// -----------------------------  output   -----------------------------

	designMasksCmd.Flags().StringP("out-file", "o", "-",
		formatFlagUsage(`Out file, supports and recommends a ".gz" suffix ("-" for stdout).`))

	designMasksCmd.Flags().IntP("top-n", "n", 100,
		formatFlagUsage(`Choose top n fils with the biggest genome size`))

	// -----------------------------  lexichash   -----------------------------

	designMasksCmd.Flags().IntP("kmer", "k", 31,
		formatFlagUsage(`Maximum k-mer size. K needs to be <= 32.`))

	designMasksCmd.Flags().IntP("masks", "m", 20480,
		formatFlagUsage(`Number of masks.`))

	designMasksCmd.Flags().IntP("seed", "s", 1,
		formatFlagUsage(`The seed for generating random masks.`))

	designMasksCmd.Flags().IntP("prefix", "p", 15,
		formatFlagUsage(`Length of mask k-mer prefix for checking low-complexity (0 for no checking).`))

	designMasksCmd.SetUsageTemplate(usageTemplate("{ -d <index path> | [-k <k>] [-n <masks>] [-s <seed>] } [-o out.tsv.gz]"))
}

func DesignMasks(files []string, opt *IndexBuildingOptions, topN int) error {
	// if topN < 100 {
	// 	topN = 100
	// }
	var timeStart time.Time
	if opt.Verbose || opt.Log2File {
		timeStart = time.Now()
		defer func() {
			log.Infof("  finished in: %s", time.Since(timeStart))
		}()
	}

	// --------------------------------------------------------------------
	// find top n files that contain the largest genome.

	var filesTop []string
	file2gsizes := &File2GSizes{}
	if len(files) > topN {
		// find topn genome files with biggest genome size
		if opt.Verbose || opt.Log2File {
			log.Infof("check genomes sizes of %d files...", len(files))
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
					decor.EwmaETA(decor.ET_STYLE_GO, 3),
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

		// receiver

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

		// do it
		var wg sync.WaitGroup                 // ensure all jobs done
		tokens := make(chan int, opt.NumCPUs) // control the max concurrency number
		k := opt.K
		filterNames := len(opt.ReSeqExclude) > 0
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

				chGS <- File2GSize{File: file, Size: genomeSize}

				chDuration <- time.Since(startTime)
			}(file)
		}

		wg.Wait()
		close(chGS)
		<-doneGS

		// process bar
		if opt.Verbose {
			close(chDuration)
			<-doneDuration
			pbs.Wait()
		}

		filesTop = make([]string, 0, topN)
		sort.Sort(*file2gsizes) // make sure it's sorted
		if opt.Verbose || opt.Log2File {
			log.Infof("  genome size range in the top %d files: [%d, %d]",
				topN, (*file2gsizes)[0].Size, (*file2gsizes)[len(*file2gsizes)-1].Size)
		}

		for _, gs := range *file2gsizes {
			filesTop = append(filesTop, gs.File)
		}
	} else {
		filesTop = files
	}

	// --------------------------------------------------------------------
	// anlysis
	if opt.Verbose || opt.Log2File {
		log.Info()
		log.Infof("analyzing %d files...", len(filesTop))
	}

	// -------------------------------------------------
	// count k-mers from the top n files
	// process bar

	var pbs *mpb.Progress
	var bar *mpb.Bar
	var chDuration chan time.Duration
	var doneDuration chan int
	if opt.Verbose {
		pbs = mpb.New(mpb.WithWidth(40), mpb.WithOutput(os.Stderr))
		bar = pbs.AddBar(int64(len(filesTop)),
			mpb.PrependDecorators(
				decor.Name("processed files: ", decor.WC{W: len("processed files: "), C: decor.DindentRight}),
				decor.Name("", decor.WCSyncSpaceR),
				decor.CountersNoUnit("%d / %d", decor.WCSyncWidth),
			),
			mpb.AppendDecorators(
				decor.Name("ETA: ", decor.WC{W: len("ETA: ")}),
				decor.EwmaETA(decor.ET_STYLE_GO, 3),
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

	// 2. Collect k-mers

	// prefix (list, 16384) -> refs (list, #refs) -> kmers (map, might >10k) -> location (list, small)
	// var data [][]map[uint64][]uint32
	nPrefix := 1
	for 1<<(nPrefix<<1) <= opt.Masks {
		nPrefix++
	}
	nPrefix-- // 7 for 20,000 files

	poolGenomes := &sync.Pool{New: func() interface{} {
		g := &Genome{
			ID:    make([]byte, 0, 128),
			Seq:   make([]byte, 0, 10<<20),
			Kmers: make(map[uint64]*[]uint32, 10<<20),
		}
		return g
	}}

	genomes := make(chan *Genome, opt.NumCPUs)
	done := make(chan int)

	go func() {
		threadsFloat := float64(opt.NumCPUs) // just avoid repeated type conversion
		for refseq := range genomes {        // each genome
			fmt.Printf("refseq: %s, size: %d bp, kmers: %d\n",
				refseq.ID, refseq.GenomeSize, len(refseq.Kmers))

			poolGenomes.Put(refseq)

			if opt.Verbose {
				chDuration <- time.Duration(float64(time.Since(refseq.StartTime)) / threadsFloat)
			}
		}
		done <- 1
	}()

	// --------------------------

	// 1. Compute k-mers from the topn files
	k := opt.K
	nnn := bytes.Repeat([]byte{'N'}, k-1)
	reRefName := opt.ReRefName
	extractRefName := reRefName != nil
	filterNames := len(opt.ReSeqExclude) > 0

	var wg sync.WaitGroup                 // ensure all jobs done
	tokens := make(chan int, opt.NumCPUs) // control the max concurrency number

	for _, file := range filesTop {
		tokens <- 1
		wg.Add(1)

		go func(file string) {
			defer func() {
				wg.Done()
				<-tokens
			}()
			startTime := time.Now()

			// --------------------------------
			// read sequence

			fastxReader, err := fastx.NewReader(nil, file, "")
			if err != nil {
				checkError(fmt.Errorf("failed to read seq file: %s", err))
			}
			defer fastxReader.Close()

			var record *fastx.Record

			var ignoreSeq bool
			var re *regexp.Regexp
			var baseFile = filepath.Base(file)

			// object for storing the genome data
			refseq := poolGenomes.Get().(*Genome)
			refseq.Reset()

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

				if i > 0 { // add N's between two contigs
					refseq.Seq = append(refseq.Seq, nnn...)
				}
				refseq.Seq = append(refseq.Seq, record.Seq.Seq...)
				refseq.GenomeSize += len(record.Seq.Seq)

				i++
			}

			if len(refseq.Seq) == 0 {
				genome.PoolGenome.Put(refseq) // important
				log.Warningf("skipping %s: no valid sequences", file)
				log.Info()
				if opt.Verbose {
					chDuration <- time.Microsecond // important, or the progress bar will get hung
				}
				return
			}

			var genomeID string // genome id
			if extractRefName {
				if reRefName.MatchString(baseFile) {
					genomeID = reRefName.FindAllStringSubmatch(baseFile, 1)[0][1]
				} else {
					genomeID, _, _ = filepathTrimExtension(baseFile, nil)
				}
			} else {
				genomeID, _, _ = filepathTrimExtension(baseFile, nil)
			}

			refseq.ID = []byte(genomeID)
			refseq.StartTime = startTime

			// --------------------------------

			iter, err := iterator.NewKmerIterator(refseq.Seq, k)
			if err != nil {
				checkError(fmt.Errorf("count kmer for %s: %s", file, err))
			}
			var kmer, kmerRC uint64
			var ok bool
			var j int
			kmers := refseq.Kmers
			var locs *[]uint32
			for {
				kmer, kmerRC, ok, _ = iter.NextKmer()
				if !ok {
					break
				}
				if kmer == 0 { // all bases are A's or N's.
					continue
				}

				j = iter.Index()

				if locs, ok = kmers[kmer]; !ok {
					tmp := []uint32{uint32(j)}
					kmers[kmer] = &tmp
				} else {
					*locs = append(*locs, uint32(j))
				}

				if locs, ok = kmers[kmerRC]; !ok {
					tmp := []uint32{uint32(j)}
					kmers[kmerRC] = &tmp
				} else {
					*locs = append(*locs, uint32(j))
				}
			}

			genomes <- refseq

		}(file)
	}

	wg.Wait() // all infiles are parsed
	close(genomes)
	<-done // all k-mer data are collected

	// process bar
	if opt.Verbose {
		close(chDuration)
		<-doneDuration
		pbs.Wait()
	}

	return nil
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

type Genome struct {
	GenomeSize int
	ID         []byte
	Seq        []byte
	Kmers      map[uint64]*[]uint32
	StartTime  time.Time
}

func (g *Genome) Reset() {
	g.ID = g.ID[:0]
	g.Seq = g.Seq[:0]
	clear(g.Kmers)

}
