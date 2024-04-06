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
	"os"
	"path/filepath"
	"regexp"
	"runtime"
	"time"

	"github.com/pkg/errors"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
)

var indexCmd = &cobra.Command{
	Use:   "index",
	Short: "Generate an index from FASTA/Q sequences",
	Long: `Generate an index from FASTA/Q sequences

Input:
 *1. Sequences of each reference genome should be saved in a separate FASTA/Q file, with the reference
     identifier in the file name.
  2. Input plain or gzipped FASTA/Q files can be given via positional arguments or the flag -X/--infile-list
     with the list of input files.
  3. Input can also be a directory containing sequence files via the flag -I/--in-dir, with multiple-level
     sub-directories allowed. A regular expression for matching sequencing files is available via the flag
     -r/--file-regexp.
  4. Some none-isolate assemblies might have extremely large genomes (e.g., GCA_000765055.1, >150 mb).
     The flag -g/--max-genome is used to skip these input files, and the file list would be write to a file
     (-G/--big-genomes).

  Attentions:
   *1) ► You can rename the sequence files for convenience because the genome identifiers in the index and
       search result would be: the basenames of files with common FASTA/Q file extensions removed, which
       are extracted via the flag -N/--ref-name-regexp.
       ► The extracted genome identifiers better be distinct/unique, which will be shown in search results
       and are used to extract subsequences in the command "lexicmap utils subseq".
    2) ► Unwanted sequences like plasmids can be filtered out by the name via regular expressions
       (-B/--seq-name-filter).

Important parameters:
  --- LexicHash mask generation ---
  0. -M/--mask-file,   ► File with custom masks, which could be genrated by "lexicmap utils gen-masks".
                       This flag oversides -k/--kmer, -m/--masks, -s/--rand-seed, -n/--top-n, and -P/--prefix-ext.
  1. -k/--kmer,        ► K-mer size (maximum: 32, default: 31).
                       ► Bigger values improve the search specificity and do not increase the index size.
  2. -m/--masks,       ► Number of masks (default: 20480).
                       ► Bigger values improve the search sensitivity and increase the index size.
  3. -n/--top-n,       ► The top N largest genomes for generating mask (default 20).
                       ► Bigger values increase the indexing memory occupation (1.0~1.5 GB per genome for 20k masks).
  4. -P/--prefix-ext,  ► Extension length of prefixes.
                       ► Bigger values improve the search sensitivity by decreasing the maximum seed distances.

  --- seeds (k-mer-value data) ---
  1. -c/--chunks,      ► Number of seed file chunks (maximum: 128, default: #CPUs).
                       ► Bigger values accelerate the search speed at the cost of a high disk reading load.
                       The maximum number should not exceed the maximum number of open files set by the
                       operating systems.
  2. -p/--partitions,  ► Number of partitions for indexing each seed file (default: 1024).
                       ► Bigger values slightly improve the search speed at the cost of slightly higher
                       memory occupation.
  3. --max-open-files, ► Maximum number of open files (default: 512).
                       ► It's only used in merging indexes of multiple genome batches.

  --- genome data ---
  1. -b/--batch-size,  ► Maximum number of genomes in each batch (maximum: 131072,default: 10000).
                       ► If the number of input files exceeds this number, input files are into multiple
                       batches and indexes are built for all batches. Next, seed files are merged into a
                       big one, while genome data files are kept unchanged and collected.
                       ► Bigger values increase indexing memory occupation, while improving the search speed.

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
		if k < minK || k > 32 {
			checkError(fmt.Errorf("the value of flag -k/--kmer should be in range of [%d, 32]", minK))
		}

		nMasks := getFlagPositiveInt(cmd, "masks")
		// lcPrefix := getFlagNonNegativeInt(cmd, "prefix")
		seed := getFlagPositiveInt(cmd, "rand-seed")
		maskFile := getFlagString(cmd, "mask-file")
		chunks := getFlagPositiveInt(cmd, "chunks")
		partitions := getFlagPositiveInt(cmd, "partitions")
		batchSize := getFlagPositiveInt(cmd, "batch-size")
		maxOpenFiles := getFlagPositiveInt(cmd, "max-open-files")

		maxGenomeSize := getFlagNonNegativeInt(cmd, "max-genome")
		fileBigGenomes := getFlagString(cmd, "big-genomes")

		topN := getFlagPositiveInt(cmd, "top-n")
		prefixExt := getFlagPositiveInt(cmd, "prefix-ext")

		outDir := getFlagString(cmd, "out-dir")
		force := getFlagBool(cmd, "force")

		if outDir == "" {
			checkError(fmt.Errorf("flag -O/--out-dir is needed"))
		}

		var err error

		inDir := getFlagString(cmd, "in-dir")
		// skipFileCheck := getFlagBool(cmd, "skip-file-check")

		outDir = filepath.Clean(outDir)

		if filepath.Clean(inDir) == outDir {
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
		// options for building index
		bopt := &IndexBuildingOptions{
			// general
			NumCPUs:      opt.NumCPUs,
			Verbose:      opt.Verbose,
			Log2File:     opt.Log2File,
			Force:        force,
			MaxOpenFiles: maxOpenFiles,

			// skip extremely large genomes
			MaxGenomeSize: maxGenomeSize,
			BigGenomeFile: fileBigGenomes,

			// LexicHash
			MaskFile: maskFile,
			K:        k,
			Masks:    nMasks,
			RandSeed: int64(seed),
			// PrefixForCheckLC: lcPrefix,

			// generate masks
			TopN:      topN,
			PrefixExt: prefixExt,

			// k-mer-value data
			Chunks:     chunks,
			Partitions: partitions,

			// genome batches
			GenomeBatchSize: batchSize,

			// genome
			ReRefName:    reRefName,
			ReSeqExclude: reSeqNames,

			SaveSeedPositions: getFlagBool(cmd, "save-seed-pos"),
		}
		err = CheckIndexBuildingOptions(bopt)
		checkError(err)

		// ---------------------------------------------------------------
		// out dir

		outputDir := outDir != ""
		if outputDir {
			makeOutDir(outDir, force, "out-dir")
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
			files = getFileListFromArgsAndFile(cmd, args, false, "infile-list", false)
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
			log.Infof("--------------------- [ main parameters ] ---------------------")
			log.Info()
			log.Info("input and output:")
			log.Infof("  input directory: %s", inDir)
			log.Infof("    regular expression of input files: %s", reFileStr)
			log.Infof("    *regular expression for extracting reference name from file name: %s", reRefNameStr)
			log.Infof("    *regular expressions for filtering out sequences: %s", reSeqNameStrs)
			log.Infof("  max genome size: %d", maxGenomeSize)
			log.Infof("  output directory: %s", outDir)
			if fileBigGenomes != "" {
				log.Infof("  output file of skipped genomes: %s", fileBigGenomes)
			}
			log.Info()
			if maskFile != "" {
				log.Infof("  custom mask file: %s", maskFile)
			} else {
				log.Infof("k-mer size: %d", k)
				log.Infof("number of masks: %d", nMasks)
				log.Infof("rand seed: %d", seed)
				log.Infof("top N genomes for generating mask: %d", topN)
				log.Infof("prefix extension length: %d", prefixExt)
			}
			log.Info()
			log.Infof("seeds data chunks: %d", chunks)
			log.Infof("seeds data indexing partitions: %d", partitions)
			log.Info()
			log.Infof("genome batch size: %d", batchSize)
			log.Info()
		}

		// ---------------------------------------------------------------

		// index
		err = BuildIndex(outDir, files, bopt)
		if err != nil {
			checkError(fmt.Errorf("failed to create a new index: %s", err))
		}

		if opt.Verbose || opt.Log2File {
			log.Info()
			log.Infof("finished building LexicMap index from %d files with %d masks in %s",
				len(files), nMasks, time.Since(timeStart))
			log.Infof("LexicMap index saved: %s", outDir)
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

	indexCmd.Flags().StringP("ref-name-regexp", "N", `(?i)(.+)\.(f[aq](st[aq])?|fna)(.gz)?$`,
		formatFlagUsage(`Regular expression (must contains "(" and ")") for extracting the reference name from the filename.`))

	indexCmd.Flags().StringSliceP("seq-name-filter", "B", []string{},
		formatFlagUsage(`List of regular expressions for filtering out sequences by header/name, case ignored.`))

	// indexCmd.Flags().BoolP("skip-file-check", "S", false,
	// 	formatFlagUsage(`Skip input file checking when given files or a file list.`))

	indexCmd.Flags().IntP("max-genome", "g", 15000000,
		formatFlagUsage(`Maximum genome size. Extremely large genomes (non-isolate assemblies) will be skipped.`))

	// -----------------------------  output  -----------------------------

	indexCmd.Flags().StringP("out-dir", "O", "",
		formatFlagUsage(`Output directory.`))

	indexCmd.Flags().StringP("big-genomes", "G", "",
		formatFlagUsage(`Out file of skipped files with genomes >= -G/--max-genome`))

	indexCmd.Flags().BoolP("force", "", false,
		formatFlagUsage(`Overwrite existing output directory.`))

	// -----------------------------  lexichash masks   -----------------------------

	indexCmd.Flags().IntP("kmer", "k", 31,
		formatFlagUsage(`Maximum k-mer size. K needs to be <= 32.`))

	indexCmd.Flags().IntP("masks", "m", 20000,
		formatFlagUsage(`Number of masks.`))

	indexCmd.Flags().IntP("rand-seed", "s", 1,
		formatFlagUsage(`Rand seed for generating random masks.`))

	indexCmd.Flags().StringP("mask-file", "M", "",
		formatFlagUsage(`File of custom masks. This flag oversides -k/--kmer, -m/--masks, -s/--rand-seed, --prefix, -n/--top-n, and -P/--prefix-ext.`))

	// ------  generate masks randomly

	// indexCmd.Flags().IntP("prefix", "", 15,
	// 	formatFlagUsage(`Length of mask k-mer prefix for checking low-complexity (0 for no checking).`))

	// ------  generate mask from the top N biggest genomes

	indexCmd.Flags().IntP("top-n", "n", 20,
		formatFlagUsage(`Select the top N largest genomes for generating masks.`))

	indexCmd.Flags().IntP("prefix-ext", "P", 8,
		formatFlagUsage(`Extension length of prefixes, higher values -> smaller maximum seed distances.`))

	// -----------------------------  kmer-value data   -----------------------------

	defaultChunks := runtime.NumCPU()
	if defaultChunks > 128 {
		defaultChunks = 128
	}
	indexCmd.Flags().IntP("chunks", "c", defaultChunks,
		formatFlagUsage(`Number of chunks for storing seeds (k-mer-value data) files.`))
	indexCmd.Flags().IntP("partitions", "p", 1024,
		formatFlagUsage(`Number of partitions for indexing seeds (k-mer-value data) files.`))
	indexCmd.Flags().IntP("max-open-files", "", 512,
		formatFlagUsage(`Maximum opened files, used in merging indexes.`))

	indexCmd.Flags().BoolP("save-seed-pos", "", false,
		formatFlagUsage(`Save seed positions, which can be inspect with "lexicmap utils seed-positions".`))

	// -----------------------------  genome batches   -----------------------------

	indexCmd.Flags().IntP("batch-size", "b", 10000,
		formatFlagUsage(`Maximum number of genomes in each batch (maximum value: 131072)`))

	// ----------------------------------------------------------

	indexCmd.SetUsageTemplate(usageTemplate("[-k <k>] [-m <masks>] { -I <seqs dir> | -X <file list>} -O <out dir>"))
}

var reIgnoreCaseStr = "(?i)"
var reIgnoreCase = regexp.MustCompile(`\(\?i\)`)
