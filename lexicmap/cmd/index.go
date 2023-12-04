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
	"regexp"
	"sync"
	"time"

	"github.com/pkg/errors"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/lexichash/index"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
	"github.com/vbauerster/mpb/v8"
	"github.com/vbauerster/mpb/v8/decor"
)

var indexCmd = &cobra.Command{
	Use:   "index",
	Short: "Generate index from FASTA/Q sequences",
	Long: `Generate index from FASTA/Q sequences

Input:
  1. Input plain or gzipped FASTA/Q files can be given via positional
     arguments or the flag -X/--infile-list with the list of input files,
  2. Or a directory containing sequence files via the flag -I/--in-dir,
     with multiple-level sub-directories allowed. A regular expression
     for matching sequencing files is available via the flag -r/--file-regexp.
 *3. For taxonomic profiling, the sequences of each reference genome should be
     saved in a separate file, with the reference identifier in the file name.

  Attention:
    You may rename the sequence files for convenience because the 
  sequence/genome identifier in the index and search results would be:
    1). For the default mode (computing k-mers for the whole file):
          the basename of file with common FASTA/Q file extension removed,
          captured via the flag -N/--ref-name-regexp.  
    2). For computing k-mers for each sequence:
          the sequence identifier.

Attentions:
  1. Unwanted sequences like plasmid can be filtered out by
     the name via regular expressions (-B/--seq-name-filter).
  2. By default, LexicMap index each file as a whole genome,
     you can also use --by-seq to compute for every sequence,
     where sequence IDs in all input files are better to be distinct.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		seq.ValidateSeq = false

		var fhLog *os.File
		if opt.Log2File {
			fhLog = addLog(opt.LogFile, opt.Verbose)
		}
		timeStart := time.Now()
		defer func() {
			if opt.Verbose || opt.Log2File {
				log.Info()
				log.Infof("elapsed time: %s", time.Since(timeStart))
				log.Info()
			}
			if opt.Log2File {
				fhLog.Close()
			}
		}()

		// ---------------------------------------------------------------
		// basic flags

		k := getFlagPositiveInt(cmd, "kmer")
		if k < 4 || k > 32 {
			checkError(fmt.Errorf("the value of flag -k/--kmer should be in range of [4, 32]"))
		}

		nMasks := getFlagPositiveInt(cmd, "masks")
		seed := getFlagPositiveInt(cmd, "seed")

		outDir := getFlagString(cmd, "out-dir")
		force := getFlagBool(cmd, "force")
		bySeq := getFlagBool(cmd, "by-seq")
		skipFileCheck := getFlagBool(cmd, "skip-file-check")

		if outDir == "" {
			checkError(fmt.Errorf("flag -O/--out-dir is needed"))
		}

		var err error

		inDir := getFlagString(cmd, "in-dir")

		if filepath.Clean(inDir) == filepath.Clean(outDir) {
			checkError(fmt.Errorf("intput and output paths should not be the same: %s", outDir))
		}

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
		var extractRefName bool
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
			extractRefName = true
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
		filterNames := len(reSeqNames) > 0

		// ---------------------------------------------------------------
		// out dir

		outputDir := outDir != ""
		if outputDir {
			makeOutDir(outDir, force)
		}

		// ---------------------------------------------------------------
		// input files

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
		// log

		if opt.Verbose || opt.Log2File {
			log.Info()
			log.Infof("-------------------- [main parameters] --------------------")
			log.Info()
			log.Info("input and output:")
			log.Infof("  input directory: %s", inDir)
			log.Infof("    regular expression of input files: %s", reFileStr)
			log.Infof("    *regular expression for extracting reference name from file name: %s", reRefNameStr)
			log.Infof("    *regular expressions for filtering out sequences: %s", reSeqNameStrs)
			log.Infof("  output directory: %s", outDir)
			log.Info()
			log.Infof("k-mer size: %d", k)
			log.Infof("number of masks: %d", nMasks)
			log.Infof("rand seed: %d", seed)
			log.Info()

			if bySeq {
				log.Infof("index each sequence: %v", bySeq)
				log.Info()
			}

			log.Infof("-------------------- [main parameters] --------------------")
			log.Info()
			log.Infof("computing ...")
		}

		// ---------------------------------------------------------------

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
					decor.EwmaETA(decor.ET_STYLE_GO, 10),
					decor.OnComplete(decor.Name(""), ". done"),
				),
			)

			chDuration = make(chan time.Duration, opt.NumCPUs)
			doneDuration = make(chan int)
			go func() {
				for t := range chDuration {
					bar.Increment()
					bar.EwmaIncrBy(1, t)
				}
				doneDuration <- 1
			}()
		}

		// index
		idx, err := index.NewIndexWithSeed(k, nMasks, int64(seed))
		checkError(err)
		// BatchInsert is faster than Insert()
		input, done := idx.BatchInsert()

		// wait group
		var wg sync.WaitGroup                 // ensure all jobs done
		tokens := make(chan int, opt.NumCPUs) // control the max concurrency number
		threadsFloat := float64(opt.NumCPUs)  // just avoid repeated type conversion

		nnn := bytes.Repeat([]byte{'N'}, k-1)

		for _, file := range files {
			tokens <- 1
			wg.Add(1)

			go func(file string) {
				defer func() {
					wg.Done()
					<-tokens
				}()
				startTime := time.Now()

				var record *fastx.Record
				var fastxReader *fastx.Reader

				fastxReader, err = fastx.NewReader(nil, file, "")
				checkError(err)

				if !bySeq {
					var ignoreSeq bool
					var re *regexp.Regexp
					var baseFile = filepath.Base(file)

					refseq := index.PoolRefSeq.Get().(*index.RefSeq)
					refseq.Reset()

					var i int = 0
					for {
						record, err = fastxReader.Read()
						if err != nil {
							if err == io.EOF {
								break
							}
							checkError(err)
							break
						}

						// filter out sequences shorter than k
						if len(record.Seq.Seq) < k {
							continue
						}

						// filter out sequences with names in the blast list
						if filterNames {
							ignoreSeq = false
							for _, re = range reSeqNames {
								if re.Match(record.Name) {
									ignoreSeq = true
									break
								}
							}
							if ignoreSeq {
								continue
							}
						}
						if i > 0 {
							refseq.Seq = append(refseq.Seq, nnn...)
						}
						refseq.Seq = append(refseq.Seq, record.Seq.Seq...)
						refseq.SeqSizes = append(refseq.SeqSizes, len(record.Seq.Seq))
						refseq.RefSeqSize += len(record.Seq.Seq)

						i++
					}

					if len(refseq.Seq) == 0 {
						index.PoolRefSeq.Put(refseq)
						log.Warningf("skipping %s: no valid sequences", file)
						log.Info()
						return
					}

					var seqID string
					if extractRefName {
						if reRefName.MatchString(baseFile) {
							seqID = reRefName.FindAllStringSubmatch(baseFile, 1)[0][1]
						} else {
							seqID, _ = filepathTrimExtension(baseFile)
						}
					} else {
						seqID, _ = filepathTrimExtension(baseFile)
					}

					refseq.ID = []byte(seqID)

					input <- refseq

					if opt.Verbose || opt.Log2File {
						chDuration <- time.Duration(float64(time.Since(startTime)) / threadsFloat)
					}

					return
				}

				for {
					record, err = fastxReader.Read()
					if err != nil {
						if err == io.EOF {
							break
						}
						checkError(err)
						break
					}

					if len(record.Seq.Seq) < k {
						continue
					}

					refseq := index.PoolRefSeq.Get().(*index.RefSeq)
					refseq.Reset()

					refseq.ID = append(refseq.ID, record.ID...)
					refseq.Seq = append(refseq.Seq, record.Seq.Seq...)
					refseq.SeqSizes = append(refseq.SeqSizes, len(record.Seq.Seq))
					refseq.RefSeqSize = len(record.Seq.Seq)

					input <- refseq
				}

				if opt.Verbose || opt.Log2File {
					chDuration <- time.Duration(float64(time.Since(startTime)) / threadsFloat)
				}

			}(file)
		}

		wg.Wait()
		close(input) // wait BatchInsert
		<-done       // wait BatchInsert

		if opt.Verbose {
			close(chDuration)
			<-doneDuration
			pbs.Wait()
		}

		if opt.Verbose || opt.Log2File {
			log.Infof("finished building LexicMap index in %s from %d files with %d masks",
				time.Since(timeStart), len(files), nMasks)
			log.Info()
			log.Infof("writing to directory: %s ...", outDir)
		}

		timeStart2 := time.Now()
		err = idx.WriteToPath(outDir, force, opt.NumCPUs)
		checkError(err)

		if opt.Verbose || opt.Log2File {
			log.Infof("finished writing to disk in %s", time.Since(timeStart2))
		}
	},
}

func init() {
	RootCmd.AddCommand(indexCmd)

	// -----------------------------  input  -----------------------------

	indexCmd.Flags().StringP("in-dir", "I", "",
		formatFlagUsage(`Directory containing FASTA/Q files. Directory symlinks are followed.`))

	indexCmd.Flags().StringP("file-regexp", "r", `\.(f[aq](st[aq])?|fna)(.gz)?$`,
		formatFlagUsage(`Regular expression for matching sequence files in -I/--in-dir, case ignored.`))

	indexCmd.Flags().BoolP("by-seq", "", false,
		formatFlagUsage(`Compute k-mers (sketches) for each sequence, instead of the whole file.`))

	indexCmd.Flags().StringP("ref-name-regexp", "N", `(?i)(.+)\.(f[aq](st[aq])?|fna)(.gz)?$`,
		formatFlagUsage(`Regular expression (must contains "(" and ")") for extracting reference name from filename.`))

	indexCmd.Flags().StringSliceP("seq-name-filter", "B", []string{},
		formatFlagUsage(`List of regular expressions for filtering out sequences by header/name, case ignored.`))

	indexCmd.Flags().BoolP("skip-file-check", "S", false,
		formatFlagUsage(`skip input file checking when given files or a file list`))

	// -----------------------------  output  -----------------------------

	indexCmd.Flags().StringP("out-dir", "O", "",
		formatFlagUsage(`Output directory.`))

	indexCmd.Flags().BoolP("force", "", false,
		formatFlagUsage(`Overwrite existed output directory.`))

	// -----------------------------  main parameters   -----------------------------

	indexCmd.Flags().IntP("kmer", "k", 31,
		formatFlagUsage(`Maximum k-mer size. K needs to be <= 32.`))

	indexCmd.Flags().IntP("masks", "n", 1024,
		formatFlagUsage(`Number of masks.`))

	indexCmd.Flags().IntP("seed", "s", 1,
		formatFlagUsage(`The seed for generating random masks.`))

	indexCmd.SetUsageTemplate(usageTemplate("[-k <k>] [-n <masks>] [-s <seed>] {[-I <seqs dir>] | <seq files> | -X <file list>} -O <out dir>"))
}

var reIgnoreCaseStr = "(?i)"
var reIgnoreCase = regexp.MustCompile(`\(\?i\)`)
