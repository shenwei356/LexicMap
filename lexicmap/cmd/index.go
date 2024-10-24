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
 *1. Sequences of each reference genome should be saved in separate FASTA/Q files, with reference identifiers
     in the file names.
  2. Input plain or gzip/xz/zstd/bzip2 compressed FASTA/Q files can be given via positional arguments or
     the flag -X/--infile-list with a list of input files.
     Flag -S/--skip-file-check is optional for skipping file checking if you trust the file list.
  3. Input can also be a directory containing sequence files via the flag -I/--in-dir, with multiple-level
     sub-directories allowed. A regular expression for matching sequencing files is available via the flag
     -r/--file-regexp.
  4. Some non-isolate assemblies might have extremely large genomes (e.g., GCA_000765055.1, >150 mb).
     The flag -g/--max-genome is used to skip these input files, and the file list would be written to a file
     (-G/--big-genomes).
     Changes since v0.5.0: 
	   - Genomes with any single contig larger than the threshold will be skipped as before.
       - However, fragmented (with many contigs) genomes with the total bases larger than the threshold will
         be split into chunks and alignments from these chunks will be merged in "lexicmap search".
     You need to increase the value for indexing fungi genomes.
  5. Maximum genome size: 268,435,456.
     More precisely: $total_bases + ($num_contigs - 1) * 1000 <= 268,435,456, as we concatenate contigs with
     1000-bp intervals of N’s to reduce the sequence scale to index.
  6. A flag -l/--min-seq-len can filter out sequences shorter than the threshold (default is the k value).

  Attention:
   *1) ► You can rename the sequence files for convenience, e.g., GCF_000017205.1.fa.gz, because the genome
       identifiers in the index and search result would be: the basenames of files with common FASTA/Q file
       extensions removed, which are extracted via the flag -N/--ref-name-regexp.
       ► The extracted genome identifiers better be distinct, which will be shown in search results
       and are used to extract subsequences in the command "lexicmap utils subseq".
    2) ► Unwanted sequences like plasmids can be filtered out by content in FASTA/Q header via regular
       expressions (-B/--seq-name-filter).
    3) All degenerate bases are converted to their lexicographic first bases. E.g., N is converted to A.
        code  bases    saved
        A     A        A
        C     C        C
        G     G        G
        T/U   T        T

        M     A/C      A
        R     A/G      A
        W     A/T      A
        S     C/G      C
        Y     C/T      C
        K     G/T      G

        V     A/C/G    A
        H     A/C/T    A
        D     A/G/T    A
        B     C/G/T    C

        N     A/C/G/T  A

Important parameters:

  --- Genome data ---
 *1. -b/--batch-size,       ► Maximum number of genomes in each batch (maximum: 131072, default: 5000).
                            ► If the number of input files exceeds this number, input files are split into multiple
                            batches and indexes are built for all batches. In the end, seed files are merged, while
                            genome data files are kept unchanged and collected.
                            ■ Bigger values increase indexing memory occupation and increase batch searching speed,
                            while single query searching speed is not affected.

  --- LexicHash mask generation ---
  0. -M/--mask-file,        ► File with custom masks, which could be exported from an existing index or newly
                            generated by "lexicmap utils masks".
                            This flag oversides -k/--kmer, -m/--masks, -s/--rand-seed, etc.
 *1. -k/--kmer,             ► K-mer size (maximum: 32, default: 31).
                            ■ Bigger values improve the search specificity and do not increase the index size.
 *2. -m/--masks,            ► Number of LexicHash masks (default: 40000).
                            ■ Bigger values improve the search sensitivity, increase the index size, and slow down
                            the search speed.

  --- Seeds data (k-mer-value data) ---
 *1. --seed-max-desert      ► Maximum length of distances between seeds (default: 200).
                            The default value of 200 guarantees queries >=200 bp would match at least one seed.
                            ► Large regions with no seeds are called sketching deserts. Deserts with seed distance
                            larger than this value will be filled by choosing k-mers roughly every
                            --seed-in-desert-dist (50 by default) bases.
                            ■ Big values decrease the search sensitivity for distant targets, speed up the indexing
                            speed, decrease the indexing memory occupation and decrease the index size. While the
                            alignment speed is almost not affected.
  2. -c/--chunks,           ► Number of seed file chunks (maximum: 128, default: value of -j/--threads).
                            ► Bigger values accelerate the search speed at the cost of a high disk reading load.
                            The maximum number should not exceed the maximum number of open files set by the
                            operating systems.
                            ► Make sure the value of '-j/--threads' in 'lexicmap search' is >= this value.
 *3. -J/--seed-data-threads ► Number of threads for writing seed data and merging seed chunks from all batches
                            (maximum: -c/--chunks, default: 8).
                            ■ The actual value is min(--seed-data-threads, max(1, --max-open-files/($batches_1_round + 2))),
                            where $batches_1_round = min(int($input_files / --batch-size), --max-open-files).
                            ■ Bigger values increase indexing speed at the cost of slightly higher memory occupation.
  4. --partitions,          ► Number of partitions for indexing each seed file (default: 4096).
                            ► Bigger values bring a little higher memory occupation.
                            ► After indexing, "lexicmap utils reindex-seeds" can be used to reindex the seeds data
                            with another value of this flag.
 *5. --max-open-files,      ► Maximum number of open files (default: 1024).
                            ► It's only used in merging indexes of multiple genome batches. If there are >100 batches,
                            ($input_files / --batch-size), please increase this value and set a bigger "ulimit -n" in shell.

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
		minSeqLen := getFlagInt(cmd, "min-seq-len")
		if minSeqLen < k {
			minSeqLen = -1
			// checkError(fmt.Errorf("the value (%d) of flag -l/--min-seq-len should be >= k (%d)", minSeqLen, k))
		}
		if minSeqLen <= 0 {
			minSeqLen = k
		}
		minSeqLen = max(minSeqLen, k)

		nMasks := getFlagPositiveInt(cmd, "masks")
		seed := getFlagPositiveInt(cmd, "rand-seed")
		maskFile := getFlagString(cmd, "mask-file")

		chunks := opt.NumCPUs // the default value is equal to -j/--threads
		_chunks := getFlagPositiveInt(cmd, "chunks")
		if chunks != _chunks && cmd.Flags().Lookup("chunks").Changed {
			chunks = _chunks
		}
		if chunks > 128 {
			chunks = 128
		}

		mergeThreads := getFlagPositiveInt(cmd, "seed-data-threads")
		if mergeThreads > chunks {
			mergeThreads = chunks
		}
		partitions := getFlagPositiveInt(cmd, "partitions")
		batchSize := getFlagPositiveInt(cmd, "batch-size")
		maxOpenFiles := getFlagPositiveInt(cmd, "max-open-files")

		maxGenomeSize := getFlagNonNegativeInt(cmd, "max-genome")
		if maxGenomeSize > MAX_GENOME_SIZE {
			checkError(fmt.Errorf("value of -g/--max-genome (%d) should not be greater than the maximum supported genome size (%d)", maxGenomeSize, MAX_GENOME_SIZE))
		}
		fileBigGenomes := getFlagString(cmd, "big-genomes")

		// minPrefix := getFlagNonNegativeInt(cmd, "seed-min-prefix")
		// if minPrefix > 32 || minPrefix < 5 {
		// 	checkError(fmt.Errorf("the value of flag -p/--seed-min-prefix (%d) should be in the range of [5, 32]", minPrefix))
		// }
		maxDesert := getFlagPositiveInt(cmd, "seed-max-desert")
		seedInDesertDist := getFlagPositiveInt(cmd, "seed-in-desert-dist")
		if seedInDesertDist > maxDesert/2 {
			checkError(fmt.Errorf("value of --seed-in-desert-dist should be smaller than 0.5 * --seed-max-desert"))
		}
		noDesertFilling := getFlagBool(cmd, "no-desert-filling")

		// topN := getFlagPositiveInt(cmd, "top-n")
		// prefixExt := getFlagPositiveInt(cmd, "prefix-ext")

		outDir := getFlagString(cmd, "out-dir")
		force := getFlagBool(cmd, "force")

		if outDir == "" {
			checkError(fmt.Errorf("flag -O/--out-dir is needed"))
		}

		var err error

		inDir := getFlagString(cmd, "in-dir")
		skipFileCheck := getFlagBool(cmd, "skip-file-check")

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

		contigInterval := getFlagPositiveInt(cmd, "contig-interval")
		if contigInterval < maxDesert {
			checkError(fmt.Errorf("the value of --contig-interval (%d) should be >= -D/--seed-max-desert (%d)", contigInterval, maxDesert))
		}

		// refNameStr := getFlagString(cmd, "ref-name-info")
		// var name2info map[string]string

		// ---------------------------------------------------------------
		// options for building index
		bopt := &IndexBuildingOptions{
			// general
			NumCPUs:      opt.NumCPUs,
			Verbose:      opt.Verbose,
			Log2File:     opt.Log2File,
			Force:        force,
			MaxOpenFiles: maxOpenFiles,
			MergeThreads: mergeThreads,

			MinSeqLen: minSeqLen,

			// skip extremely large genomes
			MaxGenomeSize: maxGenomeSize,
			BigGenomeFile: fileBigGenomes,

			// LexicHash
			MaskFile: maskFile,
			K:        k,
			Masks:    nMasks,
			RandSeed: int64(seed),

			// randomly generating
			// Prefix: minPrefix,

			// filling sketching deserts
			DisableDesertFilling:   noDesertFilling,      // disable desert filling (just for analysis index)
			DesertMaxLen:           uint32(maxDesert),    // maxi length of sketching deserts
			DesertExpectedSeedDist: seedInDesertDist,     // expected distance between seeds
			DesertSeedPosRange:     seedInDesertDist / 2, // the upstream and down stream region for adding a seeds

			// generate masks
			// TopN:      topN,
			// PrefixExt: prefixExt,

			// k-mer-value data
			Chunks:     chunks,
			Partitions: partitions,

			// genome batches
			GenomeBatchSize: batchSize,

			// genome
			ReRefName:    reRefName,
			ReSeqExclude: reSeqNames,

			ContigInterval: contigInterval,

			SaveSeedPositions: getFlagBool(cmd, "save-seed-pos"),

			Debug: getFlagBool(cmd, "debug"),
		}
		err = CheckIndexBuildingOptions(bopt)
		checkError(err)

		// ---------------------------------------------------------------
		// out dir

		outputDir := outDir != ""
		if outputDir {
			makeOutDir(outDir, force, "out-dir", opt.Verbose || opt.Log2File)
		}

		// ---------------------------------------------------------------
		// input files

		if opt.Verbose || opt.Log2File {
			log.Infof("LexicMap v%s (%s)", VERSION, COMMIT)
			log.Info("  https://github.com/shenwei356/LexicMap")
			log.Info()

		}

		// if refNameStr != "" {
		// 	name2info, err = readKVs(refNameStr, false)
		// 	checkError(err)
		// 	if opt.Verbose || opt.Log2File {
		// 		log.Infof("%d reference name information records loaded", len(name2info))
		// 	}
		// }

		if opt.Verbose || opt.Log2File {
			log.Info("checking input files ...")
		}

		var files []string
		if readFromDir {
			if opt.Verbose || opt.Log2File {
				log.Infof("  scanning files from directory: %s", inDir)
			}
			files, err = getFileListFromDir(inDir, reFile, opt.NumCPUs)
			if err != nil {
				checkError(errors.Wrapf(err, "walking dir: %s", inDir))
			}
			if len(files) == 0 {
				log.Warningf("  no files matching regular expression: %s", reFileStr)
			}
		} else {
			if opt.Verbose || opt.Log2File {
				log.Info("  checking files from command-line argument or/and file list ...")
			}
			files = getFileListFromArgsAndFile(cmd, args, !skipFileCheck, "infile-list", !skipFileCheck)
			if opt.Verbose || opt.Log2File {
				if len(files) == 1 && isStdin(files[0]) {
					log.Info("  no files given, reading from stdin")
				}
			}
		}
		if len(files) < 1 {
			checkError(fmt.Errorf("FASTA/Q files needed"))
		} else if len(files) > 1<<BITS_IDX { // 1<< 34
			checkError(fmt.Errorf("at most %d files supported, given: %d", 1<<BITS_IDX, len(files)))
		} else if opt.Verbose || opt.Log2File {
			log.Infof("  %d input file(s) given", len(files))
		}

		// sort files according to taxonomic information
		// if len(name2info) > 0 {
		// 	if opt.Verbose || opt.Log2File {
		// 		log.Info("sorting input files according to reference name information...")
		// 	}
		// 	file2info := make([][2]string, len(files))

		// 	var baseFile, genomeID string
		// 	for i, file := range files {
		// 		baseFile = filepath.Base(file)
		// 		if reRefName.MatchString(baseFile) {
		// 			genomeID = reRefName.FindAllStringSubmatch(baseFile, 1)[0][1]
		// 		} else {
		// 			genomeID, _, _ = filepathTrimExtension(baseFile, nil)
		// 		}

		// 		file2info[i] = [2]string{file, name2info[genomeID]}
		// 	}
		// 	sort.Slice(file2info, func(i, j int) bool {
		// 		a, b := file2info[i][1], file2info[j][1]
		// 		if a == b {
		// 			return strings.Compare(file2info[i][0], file2info[j][0]) < 0
		// 		}
		// 		return strings.Compare(a, b) < 0
		// 	})
		// 	for i := range file2info {
		// 		files[i] = file2info[i][0]
		// 		// fmt.Printf("%s, %s\n", files[i], file2info[i][1])
		// 	}
		// 	if opt.Verbose || opt.Log2File {
		// 		log.Info("  input files sorted")
		// 	}
		// }

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
			log.Infof("  min sequence length: %d", minSeqLen)
			log.Infof("  max genome size: %d", maxGenomeSize)
			log.Infof("  output directory: %s", outDir)
			if fileBigGenomes != "" {
				log.Infof("  output file of skipped genomes: %s", fileBigGenomes)
			}
			log.Info()

			log.Info("mask generation:")
			if maskFile != "" {
				log.Infof("  custom mask file: %s", maskFile)
			} else {
				log.Infof("  k-mer size: %d", k)
				log.Infof("  number of masks: %d", nMasks)
				log.Infof("  rand seed: %d", seed)
				// log.Infof("  prefix length for checking low-complexity in mask generation: %d", minPrefix)
			}

			log.Info()
			log.Info("seed data:")
			if noDesertFilling {
				log.Infof("  disable desert filling: %v", noDesertFilling)
			} else {
				log.Infof("  maximum sketching desert length: %d", maxDesert)
				log.Infof("  distance of k-mers to fill deserts: %d", seedInDesertDist)
			}
			log.Infof("  seeds data chunks: %d", chunks)
			log.Infof("  seeds data indexing partitions: %d", partitions)
			log.Info()
			log.Info("general:")
			log.Infof("  genome batch size: %d", batchSize)
			log.Infof("  batch merge threads: %d", mergeThreads)
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
		formatFlagUsage(`Input directory containing FASTA/Q files. Directory and file symlinks are followed.`))

	indexCmd.Flags().StringP("file-regexp", "r", `\.(f[aq](st[aq])?|fna)(\.gz|\.xz|\.zst|\.bz2)?$`,
		formatFlagUsage(`Regular expression for matching sequence files in -I/--in-dir, case ignored. Attention: use double quotation marks for patterns containing commas, e.g., -p '"A{2,}"'.`))

	indexCmd.Flags().StringP("ref-name-regexp", "N", `(?i)(.+)\.(f[aq](st[aq])?|fna)(\.gz|\.xz|\.zst|\.bz2)?$`,
		formatFlagUsage(`Regular expression (must contains "(" and ")") for extracting the reference name from the filename. Attention: use double quotation marks for patterns containing commas, e.g., -p '"A{2,}"'.`))

	indexCmd.Flags().StringSliceP("seq-name-filter", "B", []string{},
		formatFlagUsage(`List of regular expressions for filtering out sequences by contents in FASTA/Q header/name, case ignored.`))

	indexCmd.Flags().BoolP("skip-file-check", "S", false,
		formatFlagUsage(`Skip input file checking when given files or a file list.`))

	indexCmd.Flags().IntP("min-seq-len", "l", -1,
		formatFlagUsage(`Maximum sequence length to index. The value would be k for values <= 0.`))

	indexCmd.Flags().IntP("max-genome", "g", 15000000,
		formatFlagUsage(fmt.Sprintf(`Maximum genome size. Genomes with any single contig larger than the threshold will be skipped, while fragmented (with many contigs) genomes larger than the threshold will be split into chunks and alignments from these chunks will be merged in "lexicmap search". The value needs to be smaller than the maximum supported genome size: %d.`, MAX_GENOME_SIZE)))

	// indexCmd.Flags().StringP("ref-name-info", "", ``,
	// 	formatFlagUsage(`A two-column tab-delimted file for mapping reference names (extracted by --ref-name-regexp) to taxonomic information such as species names. It helps to reduce memory usage.`))

	// -----------------------------  output  -----------------------------

	indexCmd.Flags().StringP("out-dir", "O", "",
		formatFlagUsage(`Output LexicMap index directory.`))

	indexCmd.Flags().StringP("big-genomes", "G", "",
		formatFlagUsage(`Out file of skipped files with $total_bases + ($num_contigs - 1) * $contig_interval >= -g/--max-genome. The second column is one of the skip types: no_valid_seqs, too_large_genome, too_many_seqs.`))

	indexCmd.Flags().BoolP("force", "", false,
		formatFlagUsage(`Overwrite existing output directory.`))

	// -----------------------------  lexichash masks   -----------------------------

	indexCmd.Flags().IntP("kmer", "k", 31,
		formatFlagUsage(`Maximum k-mer size. K needs to be <= 32.`))

	indexCmd.Flags().IntP("masks", "m", 40000,
		formatFlagUsage(`Number of LexicHash masks.`))

	indexCmd.Flags().IntP("rand-seed", "s", 1,
		formatFlagUsage(`Rand seed for generating random masks.`))

	indexCmd.Flags().StringP("mask-file", "M", "",
		formatFlagUsage(`File of custom masks. This flag oversides -k/--kmer, -m/--masks, -s/--rand-seed etc.`))
	// formatFlagUsage(`File of custom masks. This flag oversides -k/--kmer, -m/--masks, -s/--rand-seed, -p/--seed-min-prefix, etc.`))

	// ------  generate masks randomly

	indexCmd.Flags().BoolP("no-desert-filling", "", false,
		formatFlagUsage(`Disable sketching desert filling (only for debug).`))
	// indexCmd.Flags().IntP("seed-min-prefix", "p", 15,
	// 	formatFlagUsage(`Minimum length of shared substrings (anchors) in searching. Here, this value is used to remove low-complexity masks and choose k-mers to fill sketching deserts.`))
	indexCmd.Flags().IntP("seed-max-desert", "D", 200,
		formatFlagUsage(`Maximum length of sketching deserts, or maximum seed distance. Deserts with seed distance larger than this value will be filled by choosing k-mers roughly every --seed-in-desert-dist bases.`))
	indexCmd.Flags().IntP("seed-in-desert-dist", "d", 50,
		formatFlagUsage(`Distance of k-mers to fill deserts.`))

	// ------  generate mask from the top N biggest genomes

	// indexCmd.Flags().IntP("top-n", "n", 20,
	// 	formatFlagUsage(`Select the top N largest genomes for generating masks.`))

	// indexCmd.Flags().IntP("prefix-ext", "P", 8,
	// 	formatFlagUsage(`Extension length of prefixes, higher values -> smaller maximum seed distances.`))

	// -----------------------------  kmer-value data   -----------------------------

	defaultChunks = runtime.NumCPU()
	if defaultChunks > 128 {
		defaultChunks = 128
	}
	indexCmd.Flags().IntP("chunks", "c", defaultChunks,
		formatFlagUsage(`Number of chunks for storing seeds (k-mer-value data) files. Max: 128. Default: the value of -j/--threads.`))
	indexCmd.Flags().IntP("partitions", "", 4096,
		formatFlagUsage(`Number of partitions for indexing seeds (k-mer-value data) files. The value needs to be the power of 4.`))
	indexCmd.Flags().IntP("max-open-files", "", 1024,
		formatFlagUsage(`Maximum opened files, used in merging indexes. If there are >100 batches, please increase this value and set a bigger "ulimit -n" in shell.`))

	indexCmd.Flags().BoolP("save-seed-pos", "", false,
		formatFlagUsage(`Save seed positions, which can be inspected with "lexicmap utils seed-pos".`))

	// -----------------------------  genome batches   -----------------------------

	indexCmd.Flags().IntP("batch-size", "b", 5000,
		formatFlagUsage(fmt.Sprintf(`Maximum number of genomes in each batch (maximum value: %d)`, 1<<BITS_GENOME_IDX)))

	indexCmd.Flags().IntP("seed-data-threads", "J", 8,
		formatFlagUsage(`Number of threads for writing seed data and merging seed chunks from all batches, the value should be in range of [1, -c/--chunks]. If there are >100 batches, please also increase the value of --max-open-files and set a bigger "ulimit -n" in shell.`))

	indexCmd.Flags().IntP("contig-interval", "", 1000,
		formatFlagUsage(`Length of interval (N's) between contigs in a genome.`))

	// ----------------------------------------------------------

	indexCmd.Flags().BoolP("debug", "", false,
		formatFlagUsage(`Print debug information.`))

	indexCmd.SetUsageTemplate(usageTemplate("[-k <k>] [-m <masks>] { -I <seqs dir> | -X <file list>} -O <out dir>"))
}

var defaultChunks int
var reIgnoreCaseStr = "(?i)"
var reIgnoreCase = regexp.MustCompile(`\(\?i\)`)
