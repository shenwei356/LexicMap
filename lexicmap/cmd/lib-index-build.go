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
	"encoding/binary"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"regexp"
	"slices"
	"sync"
	"time"

	"github.com/pelletier/go-toml/v2"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/seedposition"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/kmers"
	"github.com/shenwei356/lexichash"
	"github.com/shenwei356/lexichash/iterator"
	"github.com/twotwotwo/sorts/sortutil"
	"github.com/vbauerster/mpb/v8"
	"github.com/vbauerster/mpb/v8/decor"

	itree "github.com/rdleal/intervalst/interval"
)

var be = binary.BigEndian

// MainVersion is use for checking compatibility
var MainVersion uint8 = 3

// MinorVersion is less important
var MinorVersion uint8 = 1

// ExtTmpDir is the path extension for temporary files
const ExtTmpDir = ".tmp"

// FileMasks is the file for storing lexichash mask
const FileMasks = "masks.bin"

// DirSeeds is the directory of k-mer-value data files
const DirSeeds = "seeds"

// ExtSeeds is file extention of  k-mer-value data files
const ExtSeeds = ".bin"

// DirGenomes is the directory of genomes datas
const DirGenomes = "genomes"

// FileGenomes is the name of each genome file
const FileGenomes = "genomes.bin"

// FileSeedPositions is the name of seed position file
const FileSeedPositions = "seed_positions.bin"

// FileInfo is the summary file
const FileInfo = "info.toml"

// FileGenomeIndex maps genome id to genome batch id and index in the batch
const FileGenomeIndex = "genomes.map.bin"

// FileGenomeChunks store lists of batch+genome index of genome chunks
const FileGenomeChunks = "genomes.chunks.bin"

// batchDir returns the direcotry name of a genome batch
func batchDir(batch int) string {
	return fmt.Sprintf("batch_%04d", batch)
}

// chunkFile returns the file name of a k-mer-value file
func chunkFile(chunk int) string {
	return fmt.Sprintf("chunk_%03d%s", chunk, ExtSeeds)
}

// IndexBuildingOptions contains all options for building an LexicMap index.
type IndexBuildingOptions struct {
	// general
	NumCPUs      int
	Verbose      bool // show log
	Log2File     bool // log file
	Force        bool // force overwrite existed index
	MaxOpenFiles int  // maximum opened files, used in merging indexes
	MergeThreads int  // Maximum Concurrent Merge Jobs

	MinSeqLen int // minimum sequence length, should be >= k

	// skipping extremely large genome
	MaxGenomeSize int    // Maximum genome size. Extremely large genomes (non-isolate assemblies) will be skipped
	BigGenomeFile string // Out file of skipped files with genomes

	// LexicHash
	MaskFile string // file of custom masks

	K        int   // k-mer size
	Masks    int   // number of masks
	RandSeed int64 // random seed

	// generate mask randomly
	// Prefix int // length of prefix for checking low-complexity and choosing k-mers to fill deserts

	// filling sketching deserts
	DisableDesertFilling   bool   // disable desert filling (just for analysis index)
	DesertMaxLen           uint32 // maxi length of sketching deserts
	DesertExpectedSeedDist int    // expected distance between seeds
	DesertSeedPosRange     int    // the upstream and down stream region for adding a seeds

	// generate mask from the top N biggest genomes
	TopN      int // Select the the top N largest genomes for generating masks
	PrefixExt int // Extension length of prefixes

	// k-mer-value data

	Chunks     int // the number of chunks for storing k-mer data
	Partitions int // the number of partitions for indexing k-mer data

	// genome batches

	GenomeBatchSize int // the maximum number of genomes of a batch

	// genome

	ReRefName    *regexp.Regexp   // for extracting genome id from the file name
	ReSeqExclude []*regexp.Regexp // for excluding sequences according to name pattern

	ContigInterval int // the length of N's between contigs

	SaveSeedPositions bool

	Debug bool
}

// CheckIndexBuildingOptions checks some important options
func CheckIndexBuildingOptions(opt *IndexBuildingOptions) error {
	if opt.K < minK || opt.K > 32 {
		return fmt.Errorf("invalid k value: %d, valid range: [%d, 32]", opt.K, minK)
	}
	if opt.Masks < 64 {
		return fmt.Errorf("invalid numer of masks: %d, should be >=64", opt.Masks)
	}
	// if opt.Prefix > opt.K {
	// 	return fmt.Errorf("invalid prefix: %d, valid range: [0, k], 0 for no checking", opt.Prefix)
	// }

	if opt.Chunks < 1 || opt.Chunks > 128 {
		return fmt.Errorf("invalid chunks: %d, valid range: [1, 128]", opt.Chunks)
	}

	if opt.Chunks > opt.Masks {
		return fmt.Errorf("invalid chunks: %d, should be <= masks (%d)", opt.Chunks, opt.Masks)
	}

	if opt.Partitions < 1 {
		return fmt.Errorf("invalid numer of partitions in indexing k-mer data: %d, should be >=1", opt.Partitions)
	}

	if opt.GenomeBatchSize < 1 || opt.GenomeBatchSize > 1<<BITS_BATCH_IDX {
		return fmt.Errorf("invalid genome batch size: %d, valid range: [1, %d]", opt.GenomeBatchSize, 1<<BITS_BATCH_IDX)
	}

	// ------------------------

	if opt.NumCPUs < 1 {
		return fmt.Errorf("invalid number of CPUs: %d, should be >= 1", opt.NumCPUs)
	}
	if opt.MergeThreads < 1 {
		return fmt.Errorf("invalid number of merge threads: %d, should be >= 1", opt.MergeThreads)
	}
	openFiles := opt.Chunks + 2
	if opt.MaxOpenFiles < openFiles {
		return fmt.Errorf("invalid max open files: %d, should be >= %d", opt.MaxOpenFiles, openFiles)
	}

	return nil
}

// BuildIndex builds index from a list of input files
func BuildIndex(outdir string, infiles []string, opt *IndexBuildingOptions) error {
	// they are already checked.
	//
	// check options
	// err := CheckIndexBuildingOptions(opt)
	// if err != nil {
	// 	return err
	// }

	if opt.Verbose || opt.Log2File {
		log.Info()
		log.Infof("--------------------- [ generating masks ] ---------------------")
	}

	// generate masks
	var lh *lexichash.LexicHash
	var err error

	if opt.MaskFile != "" {
		if opt.Verbose || opt.Log2File {
			log.Info()
			log.Infof("reading masks from file: %s", opt.MaskFile)
		}
		lh, err = lexichash.NewFromTextFile(opt.MaskFile)
		checkError(err)
		if len(lh.Masks) < 64 {
			return fmt.Errorf("invalid numer of masks: %d, should be >=64", opt.Masks)
		}
		opt.K = lh.K
	} else {
		lh, err = lexichash.NewWithSeed(opt.K, opt.Masks, opt.RandSeed, 0)
		if err != nil {
			return err
		}
	}

	// ----------------------------------
	// mask prefix length
	maskPrefix := 1
	for 1<<(maskPrefix<<1) <= len(lh.Masks) {
		maskPrefix++
	}
	maskPrefix--
	if maskPrefix < 1 {
		maskPrefix = 1
	}

	anchorPrefix := 0
	partitions := opt.Partitions
	for partitions > 0 {
		partitions >>= 2
		anchorPrefix++
	}
	anchorPrefix--
	if anchorPrefix < 1 {
		anchorPrefix = 1
	}

	// output failed genome
	outputBigGenomes := opt.BigGenomeFile != ""
	var outfhBG *os.File
	var chBG chan string
	var doneBG chan int
	var nBG int
	if outputBigGenomes {
		outfhBG, err = os.Create(opt.BigGenomeFile)
		if err != nil {
			checkError(fmt.Errorf("failed to write file: %s", opt.BigGenomeFile))
		}

		chBG = make(chan string, opt.NumCPUs)
		doneBG = make(chan int)

		go func() {
			for r := range chBG {
				nBG++
				outfhBG.WriteString(r)
			}

			doneBG <- 1
		}()
	}

	// create a lookup table for faster masking
	lenPrefix := 1
	for 1<<(lenPrefix<<1) <= len(lh.Masks) {
		lenPrefix++
	}
	lenPrefix--
	err = lh.IndexMasks(lenPrefix)
	if err != nil {
		checkError(fmt.Errorf("indexing masks: %s", err))
	}
	err = lh.IndexMasksWithDistinctPrefixes(lenPrefix + 1)
	if err != nil {
		checkError(fmt.Errorf("indexing masks for distinct prefixes: %s", err))
	}

	// save mask later

	if opt.Verbose || opt.Log2File {
		log.Info()
		log.Infof("--------------------- [ building index ] ---------------------")
	}

	datas := make([]*map[uint64]*[]uint64, opt.Masks)
	for i := 0; i < opt.Masks; i++ {
		m := kv.PoolKmerData.Get().(*map[uint64]*[]uint64)
		datas[i] = m
	}

	// split the files in to batches
	nFiles := len(infiles)
	nBatches := (nFiles + opt.GenomeBatchSize - 1) / opt.GenomeBatchSize
	tmpIndexes := make([]string, 0, nBatches)

	// tmp dir
	tmpDir := filepath.Clean(outdir) + ExtTmpDir
	err = os.RemoveAll(tmpDir)
	if err != nil {
		return err
	}
	if nBatches > 1 { // only used for > 1 batches
		err = os.MkdirAll(tmpDir, 0755)
		if err != nil {
			checkError(fmt.Errorf("failed to create dir: %s", err))
		}
	} else if nBatches > 1<<BITS_BATCH_IDX { // 1<<17
		checkError(fmt.Errorf("at most %d batches supported. current: %d", 1<<BITS_BATCH_IDX, nBatches))
	}

	var begin, end int
	var kvChunks int
	var hasSomeGenomes bool
	for batch := 0; batch < nBatches; batch++ {
		// files for this batch
		begin = batch * opt.GenomeBatchSize
		end = begin + opt.GenomeBatchSize
		if end > nFiles {
			end = nFiles
		}
		files := infiles[begin:end]

		// outdir for this batch
		var outdirB string
		if nBatches > 1 {
			outdirB = filepath.Join(tmpDir, batchDir(batch))
		} else {
			outdirB = outdir
		}

		// build index for this batch
		kvChunks, hasSomeGenomes = buildAnIndex(lh, uint8(maskPrefix), uint8(anchorPrefix), opt, &datas, outdirB, files, batch, nBatches, outputBigGenomes, chBG)

		if nBatches > 1 && hasSomeGenomes { // only merge indexes with valid genomes
			tmpIndexes = append(tmpIndexes, outdirB)
		}
	}

	if outputBigGenomes {
		close(chBG)
		<-doneBG
		outfhBG.Close()
		if opt.Verbose || opt.Log2File {
			log.Infof("  finished saving %d skipped genome files: %s", nBG, opt.BigGenomeFile)
		}
	}

	for _, data := range datas {
		kv.PoolKmerData.Put(data)
	}

	if nBatches == 1 {
		return nil
	}

	// merge indexes
	if opt.Verbose || opt.Log2File {
		log.Info()
		log.Infof("merging %d indexes...", len(tmpIndexes))
	}
	err = mergeIndexes(lh, uint8(maskPrefix), uint8(anchorPrefix), opt, kvChunks, outdir, tmpIndexes, tmpDir, 1)
	if err != nil {
		return fmt.Errorf("failed to merge indexes: %s", err)
	}

	// clean tmp dir
	err = os.RemoveAll(tmpDir)
	if err != nil {
		checkError(fmt.Errorf("failed to remove tmp directory: %s", err))
	}

	return err
}

// ----------------------------------

// BITS_BATCH_IDX is the number of bits to store the genome batch index.
const BITS_BATCH_IDX = 17

// BITS_GENOME_IDX is the number of bits to store the genome index.
const BITS_GENOME_IDX = 17

// MASK_GENOME_IDX is the mask of genome index.
const MASK_GENOME_IDX = (1 << BITS_GENOME_IDX) - 1

// BITS_POSITION is the number of bits to store the k-mer position/coordinate.
const BITS_POSITION = 28

// MAX_GENOME_SIZE is the maximum genome size, 268435456
const MAX_GENOME_SIZE = 1 << BITS_POSITION

// BITS_STRAND is the flag to indicate if the k-mer is from the reverse complement strand.
const BITS_STRAND = 1

// BITS_SUFFIX_IDX is the flag to indicate if the k-mer is reversed.
const BITS_REVERSE = 1

// MASK_REVERSE is the mask of reversed flag
const MASK_REVERSE = 1

// BITS_IDX is the number of bits to strore batch index and genome index.
const BITS_IDX = BITS_BATCH_IDX + BITS_BATCH_IDX

// BITS_NONE_POS is the number of bits except for position
const BITS_NONE_POS = 64 - BITS_POSITION

// BITS_NONE_IDX is the number of bits to store data except for batch index and genome index.
const BITS_NONE_IDX = 64 - BITS_BATCH_IDX - BITS_GENOME_IDX

// MASK_NONE_IDX is the mask of non-index data
const MASK_NONE_IDX = (1 << BITS_NONE_IDX) - 1

// BITS_FLAGS is the number of bits to store two bits
const BITS_FLAGS = BITS_STRAND + BITS_REVERSE

// BITS_IDX_FLAGS is the sum of BITS_IDX and BITS_FLAGS
const BITS_IDX_FLAGS = BITS_IDX + BITS_FLAGS

// ----------------------------------

// NO_VALID_SEQS means there are no valid sequences in a genome file.
const NO_VALID_SEQS = "no_valid_seqs"

// TOO_LARGE_GENOME means the genome is too big to index.
const TOO_LARGE_GENOME = "too_large_genome"

// TOO_MANY_SEQS means there are too many sequences, as we require: $total_bases + ($num_contigs - 1) * $interval_size <= 268,435,456
const TOO_MANY_SEQS = "too_many_seqs"

// build an index for the files of one batch
func buildAnIndex(lh *lexichash.LexicHash, maskPrefix uint8, anchorPrefix uint8, opt *IndexBuildingOptions,
	datas *[]*map[uint64]*[]uint64,
	outdir string, files []string, batch int, nbatches int, outputBigGenomes bool, chBG chan string) (int, bool) {

	debug := opt.Debug

	var timeStart time.Time
	if opt.Verbose || opt.Log2File {
		timeStart = time.Now()
		log.Info()
		log.Infof("  ------------------------[ batch %d/%d ]------------------------", batch+1, nbatches)
		log.Infof("  building index for batch %d with %d files...", batch+1, len(files))
		defer func() {
			log.Infof("  finished building index for batch %d in: %s", batch+1, time.Since(timeStart))
		}()
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
				decor.EwmaETA(decor.ET_STYLE_GO, 10),
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

	// -------------------------------------------------------------------
	// dir structure

	err := os.MkdirAll(outdir, 0755)
	if err != nil {
		checkError(fmt.Errorf("failed to create dir: %s", err))
	}

	// masks
	fileMask := filepath.Join(outdir, FileMasks)
	_, err = lh.WriteToFile(fileMask)
	if err != nil {
		checkError(fmt.Errorf("failed to write masks: %s", err))
	}

	// genomes
	dirGenomes := filepath.Join(outdir, DirGenomes, batchDir(batch))
	err = os.MkdirAll(dirGenomes, 0755)
	if err != nil {
		checkError(fmt.Errorf("failed to create dir: %s", err))
	}

	// seeds
	dirSeeds := filepath.Join(outdir, DirSeeds)
	err = os.MkdirAll(dirSeeds, 0755)
	if err != nil {
		checkError(fmt.Errorf("failed to create dir: %s", err))
	}

	// -------------------------------------------------------------------

	// --------------------------------
	// 2) collect k-mers data & write genomes to file
	for _, data := range *datas { // reset all maps
		clear(*data)
	}

	threadsFloat := float64(opt.NumCPUs) // just avoid repeated type conversion

	genomes := make(chan *genome.Genome, opt.NumCPUs)
	genomesW := make(chan *genome.Genome, opt.NumCPUs)
	done := make(chan int)

	// genome writer
	fileGenomes := filepath.Join(dirGenomes, FileGenomes)
	gw, err := genome.NewWriter(fileGenomes, uint32(batch))
	if err != nil {
		checkError(fmt.Errorf("failed to write genome file: %s", err))
	}
	doneGW := make(chan int)

	// seed positions
	var fileSeedLoc string
	var locw *seedposition.Writer
	if opt.SaveSeedPositions {
		fileSeedLoc = filepath.Join(dirGenomes, FileSeedPositions)
		locw, err = seedposition.NewWriter(fileSeedLoc, uint32(batch))
		if err != nil {
			checkError(fmt.Errorf("failed to write seed position file: %s", err))
		}
	}

	// 2.2) write genomes to file
	var nFiles int // the total number of indexed files
	go func() {

		for refseq := range genomesW { // each genome, one by one
			nFiles++

			// write the genome to file
			err = gw.Write(refseq)
			if err != nil {
				checkError(fmt.Errorf("failed to write genome: %s", err))
			}

			// --------------------------------
			// seed positions
			if opt.SaveSeedPositions {
				err = locw.Write(*refseq.Locs)
				if err != nil {
					checkError(fmt.Errorf("failed to write seed position: %s", err))
				}
			}

			if debug {
				log.Debugf("batch: %d, file #%d, genome data saved: %s", batch, refseq.GenomeIdx, refseq.ID)
				log.Info()
			}

			// send signal of genome being written
			refseq.Done <- 1
		}
		doneGW <- 1
	}()

	// genome id -> a list of batch+ref-index
	mGenomeChunks := make(map[int]*[]uint64, 1024)

	// 2.1) collect k-mer data
	go func() {
		var wg sync.WaitGroup
		threads := opt.NumCPUs
		tokens := make(chan int, threads)
		nMasks := opt.Masks
		chunkSize := (nMasks + threads - 1) / threads
		var j, begin, end int

		// genome-index mapping file
		fileGenomeIndex := filepath.Join(outdir, FileGenomeIndex)
		fhGI, err := os.Create(fileGenomeIndex)
		if err != nil {
			checkError(fmt.Errorf("%s", err))
		}
		bw := bufio.NewWriter(fhGI)

		var batchIDAndRefID, batchIDAndRefIDShift, refIdx uint64 // genome number
		buf := make([]byte, 8)

		lockers := make([]sync.Mutex, len(lh.Masks)) // to avoid concurrent data writes

		var li *[]uint64
		var ok bool
		for refseq := range genomes { // each genome, one by one
			genomesW <- refseq // send to save to file, asynchronously writing

			_kmers := refseq.Kmers
			loces := refseq.Locses
			extraKmers := refseq.ExtraKmers

			// refseq id -> this
			batchIDAndRefID = uint64(batch)<<BITS_GENOME_IDX | (refIdx & MASK_GENOME_IDX)
			be.PutUint16(buf[:2], uint16(len(refseq.ID)))
			bw.Write(buf[:2])
			bw.Write(refseq.ID)
			be.PutUint64(buf, batchIDAndRefID)
			bw.Write(buf)

			// genome id -> a list of batch+ref-index
			if li, ok = mGenomeChunks[refseq.GenomeIdx]; !ok {
				tmp := make([]uint64, 0, 8)
				li = &tmp
				mGenomeChunks[refseq.GenomeIdx] = &tmp
			}
			*li = append(*li, batchIDAndRefID)

			batchIDAndRefIDShift = batchIDAndRefID << BITS_NONE_IDX

			// save k-mer data into all masks by chunks
			for j = 0; j < threads; j++ { // each chunk for storing kmer-value data
				begin = j * chunkSize
				end = begin + chunkSize
				if end > nMasks {
					end = nMasks
				}

				wg.Add(1)
				tokens <- 1
				go func(begin, end int) { // a chunk of masks
					var kmer uint64
					var loc int
					var value uint64
					var ok bool
					var values *[]uint64
					var knl *[]uint64
					var _j, _end int

					var data *map[uint64]*[]uint64

					// normal k-mers
					for i := begin; i < end; i++ {
						data = (*datas)[i] // the map to save into

						if len((*loces)[i]) > 0 { // locations from from MaskKnownPrefixes might be empty.
							kmer = (*_kmers)[i] // captured k-mer by the mask
							if values, ok = (*data)[kmer]; !ok {
								// values = &[]uint64{}
								tmp := make([]uint64, 0, len((*loces)[i]))
								values = &tmp
								(*data)[kmer] = values
							}
							for _, loc = range (*loces)[i] { // position information of the captured k-mer
								//  batch idx: 17 bits
								//  ref idx:   17 bits
								//  pos:       28 bits
								//  strand:     1 bits
								//  reverse:    1 bits
								// here, the position from Mask() already contains the strand information.
								// value = uint64(batch)<<47 | ((refIdx & 131071) << 30) |
								// value = batchIDAndRefIDShift | (uint64(loc) & 1073741823)
								value = batchIDAndRefIDShift | ((uint64(loc) << BITS_REVERSE) & MASK_NONE_IDX)
								// fmt.Printf("%s, batch: %d, refIdx: %d, value: %064b\n", refseq.ID, batch, refIdx, value)

								*values = append(*values, value)
							}
						}

						// -----------------------------
						// extra k-mers

						if extraKmers == nil {
							continue
						}

						knl = (*extraKmers)[i]
						if knl == nil {
							continue
						}
						_end = len(*knl) - 2
						for _j = 0; _j <= _end; _j += 2 {
							kmer = (*knl)[_j]
							if values, ok = (*data)[kmer]; !ok {
								values = &[]uint64{}
								(*data)[kmer] = values
							}

							// value = batchIDAndRefIDShift | ((*knl)[_j+1] & 1073741823)
							value = batchIDAndRefIDShift | (((*knl)[_j+1] << BITS_REVERSE) & MASK_NONE_IDX)
							*values = append(*values, value)
						}
					}

					wg.Done()
					<-tokens
				}(begin, end)
			}

			wg.Wait() // wait all mask chunks

			// ---------------------------------------------------------
			// another round for reversed k-mers
			for j = 0; j < threads; j++ { // each chunk for storing kmer-value data
				begin = j * chunkSize
				end = begin + chunkSize
				if end > nMasks {
					end = nMasks
				}

				wg.Add(1)
				tokens <- 1
				go func(begin, end int) { // a chunk of masks
					var kmer uint64
					var loc int
					var value uint64
					var ok bool
					var values *[]uint64
					var knl *[]uint64
					var _j, _end int

					var iMasks *[]int
					var j int
					var minj int
					var mask, h, minh uint64
					k := lh.K

					var data *map[uint64]*[]uint64

					// reversed k-mers
					for i := begin; i < end; i++ {

						if len((*loces)[i]) > 0 { // locations from from MaskKnownPrefixes might be empty.
							kmer = (*_kmers)[i] // captured k-mer by the mask

							// ----------------------------------------------------
							// save reversed k-mer to another mask. (suffix index)
							kmer = kmers.MustReverse(kmer, k) // reverse the k-mer
							iMasks = lh.MaskKmer(kmer)
							minh = math.MaxUint64
							for _, j = range *iMasks {
								mask = lh.Masks[j]
								h = mask ^ kmer
								if h < minh {
									minj, minh = j, h
								}
							}
							// mask = lh.Masks[minj]
							// fmt.Printf("%s -> %s, %d\n", kmers.MustDecode(kmer, k), kmers.MustDecode(mask, k), len(*iMasks))

							lockers[minj].Lock()
							data = (*datas)[minj] // the map to save into
							if values, ok = (*data)[kmer]; !ok {
								// values = &[]uint64{}
								tmp := make([]uint64, 0, len((*loces)[i]))
								values = &tmp
								(*data)[kmer] = values
							}
							for _, loc = range (*loces)[i] {
								value = batchIDAndRefIDShift | ((uint64(loc)<<BITS_REVERSE | MASK_REVERSE) & MASK_NONE_IDX)
								*values = append(*values, value)
							}
							lockers[minj].Unlock()

							lh.RecycleMaskKmerResult(iMasks)
						}

						// -----------------------------
						// extra k-mers

						if extraKmers == nil {
							continue
						}

						knl = (*extraKmers)[i]
						if knl == nil {
							continue
						}
						_end = len(*knl) - 2
						for _j = 0; _j <= _end; _j += 2 {
							kmer = (*knl)[_j]

							// ----------------------------------------------------
							// save reversed k-mer to another mask. (suffix index)
							kmer = kmers.MustReverse(kmer, k) // reverse the k-mer
							iMasks = lh.MaskKmer(kmer)
							minh = math.MaxUint64
							for _, j = range *iMasks {
								mask = lh.Masks[j]
								h = mask ^ kmer
								if h < minh {
									minj, minh = j, h
								}
							}
							// mask = lh.Masks[minj]
							// fmt.Printf("%s -> %s, %d\n", kmers.MustDecode(kmer, k), kmers.MustDecode(mask, k), len(*iMasks))
							lockers[minj].Lock()
							data = (*datas)[minj] // the map to save into
							if values, ok = (*data)[kmer]; !ok {
								values = &[]uint64{}
								(*data)[kmer] = values
							}
							value = batchIDAndRefIDShift | (((*knl)[_j+1]<<BITS_REVERSE | MASK_REVERSE) & MASK_NONE_IDX)
							*values = append(*values, value)
							lockers[minj].Unlock()

							lh.RecycleMaskKmerResult(iMasks)
						}
					}

					wg.Done()
					<-tokens
				}(begin, end)
			}

			wg.Wait() // wait all mask chunks

			if debug {
				log.Debugf("batch: %d, file #%d, seeds collected: %s", batch, refseq.GenomeIdx, refseq.ID)
				log.Info()
			}

			// wait the genome data being written
			<-refseq.Done

			// recycle the genome
			if refseq.Kmers != nil {
				lh.RecycleMaskResult(refseq.Kmers, refseq.Locses)
			}
			genome.RecycleGenome(refseq)

			if opt.Verbose && !refseq.StartTime.IsZero() {
				chDuration <- time.Duration(float64(time.Since(refseq.StartTime)) / threadsFloat)
			}

			refIdx++
		}

		// genome index
		bw.Flush()
		checkError(fhGI.Close())

		// genome data
		close(genomesW)
		done <- 1
	}()

	// --------------------------------
	// 1) parsing input genome files & mask & pack sequences
	nnn := bytes.Repeat([]byte{'A'}, opt.ContigInterval)
	reRefName := opt.ReRefName
	extractRefName := reRefName != nil
	filterNames := len(opt.ReSeqExclude) > 0

	reGaps := regexp.MustCompile(fmt.Sprintf(`[Nn]{%d,}`, 5))

	var wg sync.WaitGroup                 // ensure all jobs done
	tokens := make(chan int, opt.NumCPUs) // control the max concurrency number

	// 1.2) mask & pack sequences

	var wgMask sync.WaitGroup                 // ensure all jobs done
	tokensMask := make(chan int, opt.NumCPUs) // control the max concurrency number

	genomesMask := make(chan *genome.Genome, opt.NumCPUs)
	doneMask := make(chan int)

	go func() {
		for refseq := range genomesMask { // each genome
			tokensMask <- 1
			wgMask.Add(1)

			go func(refseq *genome.Genome) {
				defer func() {
					if debug {
						log.Debugf("batch: %d, file #%d, finished computing seeds: %s", batch, refseq.GenomeIdx, refseq.ID)
						log.Info()
					}
					wgMask.Done()
					<-tokensMask
				}()

				if debug {
					log.Debugf("batch: %d, file #%d, start to compute seeds: %s", batch, refseq.GenomeIdx, refseq.ID)
					log.Info()
				}

				// --------------------------------
				// mask with lexichash

				k := lh.K

				// skip regions around junctions of two sequences.

				// interval tree of 1000-bp interval regions, used in filling deserts
				_itree := itree.NewSearchTree[uint8, int](cmpFn)

				// because lh.Mask accepts a list, while when skipRegions is nil, *skipRegions is illegal.
				var _skipRegions [][2]int
				var skipRegions *[][2]int
				var interval int = opt.ContigInterval
				if len(refseq.SeqSizes) > 1 {
					skipRegions = poolSkipRegions.Get().(*[][2]int)
					*skipRegions = (*skipRegions)[:0]
					var n int // len of concatenated seqs
					for i, s := range refseq.SeqSizes {
						if i > 0 {
							// 0-based region. [n, n+interval-1]
							*skipRegions = append(*skipRegions, [2]int{n, n + interval - 1})

							_itree.Insert(n-k+1, n+interval-1, 1)

							n += interval
						}
						n += s
					}
					_skipRegions = *skipRegions
				}

				// skip gap regions (N's)
				gaps := reGaps.FindAllSubmatchIndex(refseq.Seq, -1)
				if gaps != nil {
					if _skipRegions == nil {
						skipRegions = poolSkipRegions.Get().(*[][2]int)
						*skipRegions = (*skipRegions)[:0]
					}

					for _, gap := range gaps {
						*skipRegions = append(*skipRegions, [2]int{gap[0], gap[1] - 1})

						_itree.Insert(gap[0]-k+1, gap[1]-1, 1)
					}

					// sort.Slice(*skipRegions, func(i, j int) bool {
					// 	return (*skipRegions)[i][0] < (*skipRegions)[j][0]
					// })
					slices.SortFunc(*skipRegions, func(a, b [2]int) int {
						return a[0] - b[0]
					})

					_skipRegions = *skipRegions
				}
				//

				var _kmers *[]uint64
				var locses *[][]int

				// if len(lh.Masks) > 1024 && len(refseq.Seq) > 1048576 {
				// 	_kmers, locses, err = lh.MaskLongSeqs(refseq.Seq, _skipRegions)
				// } else {
				// 	_kmers, locses, err = lh.Mask(refseq.Seq, _skipRegions)
				// }
				// _kmers, locses, err = lh.MaskKnownPrefixes(refseq.Seq, _skipRegions)
				_kmers, locses, err = lh.MaskKnownDistinctPrefixes(refseq.Seq, _skipRegions, true)
				if err != nil {
					panic(err)
				}

				// remove low-complexity k-mers
				k8 := uint8(lh.K)
				ccc := util.Ns(0b01, k8)
				ggg := util.Ns(0b10, k8)
				ttt := (uint64(1) << (k << 1)) - 1
				for i, kmer := range *_kmers {
					if kmer == 0 || kmer == ccc || kmer == ggg || kmer == ttt ||
						util.IsLowComplexityDust(kmer, k8) {
						// fmt.Printf("low-complexity k-mer #%d: %s\n", i, lexichash.MustDecode(kmer, k8))
						(*_kmers)[i] = 0
						(*locses)[i] = (*locses)[i][:0]
						continue
					}
				}

				refseq.Kmers = _kmers
				refseq.Locses = locses

				// --------------------------------
				// bit-packed sequences
				refseq.TwoBit = genome.Seq2TwoBit(refseq.Seq)

				// --------------------------------
				// locations

				if refseq.Locs == nil { // existing positions, will be updated with new k-mers
					tmp := make([]uint32, 0, 40<<10)
					refseq.Locs = &tmp
				}
				locs := refseq.Locs
				*locs = (*locs)[:0]

				var loc int
				for _, _locs := range *locses { // insert positions
					for _, loc = range _locs {
						*locs = append(*locs, uint32(loc&4294967295)) // only posision | strand flag
					}
				}
				sortutil.Uint32s(*locs)

				if refseq.ExtraKmers == nil { // extra k-mers
					tmp := make([]*[]uint64, opt.Masks)
					refseq.ExtraKmers = &tmp
				}
				extraKmers := refseq.ExtraKmers
				for _i, k2l := range *extraKmers { // reset
					if k2l != nil {
						poolKmerAndLocs.Put(k2l)
						(*extraKmers)[_i] = nil
					}
				}

				// ----------------------------------------------------------------
				// fill sketching deserts

				var extraLocs *[]int
				var loc2maskidx *[]int
				var loc2maskidxRC *[]int
				var kmerList *[]uint64

				if !opt.DisableDesertFilling {
					var pos2str, pos, pre, d uint32
					maxDesert := opt.DesertMaxLen
					seedDist := opt.DesertExpectedSeedDist
					seedPosR := opt.DesertSeedPosRange

					lenSeq := len(refseq.Seq)
					var start, end int
					var _kmers2 *[]uint64
					var _locses2 *[][]int
					var _locs []int
					var _i int
					// kmer2maskidx := poolKmer2MaskIdx.Get().(*map[uint64]int)
					loc2maskidx = poolLoc2MaskIdx.Get().(*[]int)
					loc2maskidxRC = poolLoc2MaskIdx.Get().(*[]int)
					kmerList = poolKmerKmerRC.Get().(*[]uint64)

					var iter *iterator.Iterator
					var kmer, kmerRC, kmerPos uint64
					var ok bool
					var _j, posOfPre, posOfCur, _start, _end int

					var _im int
					var knl *[]uint64

					var inIntervalRegion bool

					// solved by adding gap regions into the skip regions and the interval tree.
					//
					// it's hard to skip k-mers at the edge of a gap
					//       ACCAAAAAAAA      AAAAAAGCCAGA
					// ----------AAAAAAAAAAAAAAAAAAA-------
					// ----------TTTTTTTTTTTTTTTTTTT-------
					//      TTTTTTAC         ACGTTTTTTT
					// so we just detect these k-mers by it's suffix
					// var lenSuffix uint8 = 5
					// var lenPrefix uint8 = 3                          // that means the firt n base can't be As, just in case.
					// tttPrefix := uint64((1 << (lenPrefix << 1)) - 1) // for k-mer on the negative strand
					// tttSuffix := uint64((1 << (lenSuffix << 1)) - 1) // for k-mer on the negative strand

					extraLocs = poolInts.Get().(*[]int) // extra positions
					*extraLocs = (*extraLocs)[:0]

					pre = 0
					var iD int
					for _, pos2str = range *locs {
						pos = pos2str >> 1
						d = pos - pre

						if d < maxDesert { // small distance, cool
							pre = pos
							continue
						}

						// there's a really big gap in it, it might be the interval between contigs or a assembly gap
						if float64(lengthAAs(refseq.Seq[pre:pos]))/float64(d) >= 0.7 {
							pre = pos
							continue
						}

						// range of desert region +- 1000 bp,
						// as we don't want other kmers with the same prefix exist around.
						start = int(pre) - 1000 // start position in the sequence
						posOfPre = 1000         // the location of previous seed in the list
						if start < 0 {
							posOfPre += start
							start = 0
						}
						end = int(pos) + 1000 + k // end position in the sequence
						if end > lenSeq {
							end = lenSeq
						}

						iD++
						// fmt.Printf("desert %d: %d-%d, len: %d, region: %d-%d, list size: %d\n",
						// 	iD, pre, pos, d, start, end, end-start+1)

						posOfCur = posOfPre + int(d) // their distance keeps the same

						// .                       desert
						// .                 -------------------
						// 0123      start                             end
						// ----o-----[-o-----o------------------o---o--]----o-------
						//           0123    posOfPre           posOfCur

						// fmt.Printf("  posOfPre: %d, posOfCur: %d\n", posOfPre, posOfCur)

						// iterate k-mers
						iter, err = iterator.NewKmerIterator(refseq.Seq[start:end], k)
						if err != nil {
							checkError(err)
						}

						*kmerList = (*kmerList)[:0]
						for {
							kmer, kmerRC, ok, _ = iter.NextKmer()
							if !ok {
								break
							}
							*kmerList = append(*kmerList, kmer)
							*kmerList = append(*kmerList, kmerRC)
						}

						// masks this region, just treat it as a query sequence
						// _kmers2, _locses2, _ = lh.MaskKnownPrefixes(refseq.Seq[start:end], nil)
						// here, checkShorterPrefix can be false, as we do not need all probes to capture there k-mers,
						// we only need a few.
						_kmers2, _locses2, _ = lh.MaskKnownDistinctPrefixes(refseq.Seq[start:end], nil, false)

						// // remove low-complexity k-mers
						// // k8 := uint8(lh.K)
						// for i, kmer := range *_kmers2 {
						// 	if kmer == ttt || (kmer != 0 && util.IsLowComplexity(kmer, k8)) {
						// 		// fmt.Printf("low-complexity k-mer #%d: %s\n", i, lexichash.MustDecode(kmer, k8))
						// 		(*_kmers2)[i] = 0
						// 		continue
						// 	}
						// }

						// clear(*kmer2maskidx)
						// for _i, kmer = range *_kmers2 {
						// 	// mulitple masks probably capture more than one k-mer in such a short sequence,
						// 	// we just record the last mask.
						// 	(*kmer2maskidx)[kmer] = _i
						// }

						*loc2maskidx = (*loc2maskidx)[:0]
						*loc2maskidxRC = (*loc2maskidxRC)[:0]
						for _i = start; _i < end; _i++ {
							*loc2maskidx = append(*loc2maskidx, -1)
							*loc2maskidxRC = append(*loc2maskidxRC, -1)
						}
						for _i, _locs = range *_locses2 {
							for _, loc = range _locs {
								if loc&1 == 0 {
									(*loc2maskidx)[loc>>1] = _i
								} else {
									(*loc2maskidxRC)[loc>>1] = _i
								}
							}
						}

						lh.RecycleMaskResult(_kmers2, _locses2)

						// start from the previous seed
						_j = posOfPre + seedDist
						for {
							if _j >= posOfCur {
								break
							}

							// upstream scan range: _j-seedPosR, _j
							_start = _j + 1 // start of downstream scan range
							_end = _j - seedPosR

							// fmt.Printf("  check %d, <-end: %d, ->start: %d\n", _j, _end, _start)

							ok = false
							for ; _j > _end; _j-- {
								// fmt.Printf("    test u %d\n", _j)

								if _, inIntervalRegion = _itree.AnyIntersection(start+_j, start+_j); inIntervalRegion {
									continue
								}

								// strand +
								kmer = (*kmerList)[_j<<1]
								// if kmer != 0 &&
								// 	!util.MustKmerHasSuffix(kmer, 0, k8, lenSuffix) &&
								// 	!util.MustKmerHasPrefix(kmer, 0, k8, lenPrefix) {
								// if kmer != 0 {
								if kmer != 0 && kmer != ccc && kmer != ggg && kmer != ttt &&
									!util.IsLowComplexityDust(kmer, k8) {
									// if _im, ok = (*kmer2maskidx)[kmer]; ok {
									// 	kmerPos = uint64(start+_j) << 1
									// 	break
									// }
									_im = (*loc2maskidx)[_j]
									if _im >= 0 {
										kmerPos = uint64(start+_j) << 1
										ok = true
										break
									}
								}

								// strand -
								kmer = (*kmerList)[(_j<<1)+1]
								// if kmer != 0 &&
								// 	!util.MustKmerHasSuffix(kmer, tttSuffix, k8, lenSuffix) &&
								// 	!util.MustKmerHasPrefix(kmer, tttPrefix, k8, lenPrefix) {
								// if kmer != 0 {
								if kmer != 0 && kmer != ccc && kmer != ggg && kmer != ttt &&
									!util.IsLowComplexityDust(kmer, k8) {
									// if _im, ok = (*kmer2maskidx)[kmer]; ok {
									// 	kmerPos = uint64(start+_j)<<1 | 1
									// 	break
									// }
									_im = (*loc2maskidxRC)[_j]
									if _im >= 0 {
										kmerPos = uint64(start+_j)<<1 | 1
										ok = true
										break
									}
								}
							}
							if ok {
								// fmt.Printf("    uadd: %s at %d (%d)\n", kmers.MustDecode(kmer, k), _j, start+_j)

								knl = (*extraKmers)[_im]
								if knl == nil {
									knl = poolKmerAndLocs.Get().(*[]uint64)
									*knl = (*knl)[:0]
									(*extraKmers)[_im] = knl
								}
								*knl = append(*knl, kmer)
								*knl = append(*knl, kmerPos)

								// fmt.Printf("  ADD to mask %d with %s, from %d\n", _im+1, lexichash.MustDecode(kmer, k8), (kmerPos>>1)+1)

								*extraLocs = append(*extraLocs, int(kmerPos))

								_j += seedDist
								continue
							}

							if _start >= posOfCur {
								break
							}

							// downstream scan range: _j+1, _j+seedpoR
							_end = _start + seedPosR
							if _end >= posOfCur {
								_end = posOfCur - 1
							}
							for _j = _start; _j < _end; _j++ {
								// fmt.Printf("    test d %d\n", _j)

								if _, inIntervalRegion = _itree.AnyIntersection(start+_j, start+_j); inIntervalRegion {
									continue
								}

								// strand +
								kmer = (*kmerList)[_j<<1]
								// if kmer != 0 &&
								// 	!util.MustKmerHasSuffix(kmer, 0, k8, lenSuffix) &&
								// 	!util.MustKmerHasPrefix(kmer, 0, k8, lenPrefix) {
								// if kmer != 0 {
								if kmer != 0 && kmer != ccc && kmer != ggg && kmer != ttt &&
									!util.IsLowComplexityDust(kmer, k8) {
									// if _im, ok = (*kmer2maskidx)[kmer]; ok {
									// 	kmerPos = uint64(start+_j) << 1
									// 	break
									// }
									_im = (*loc2maskidx)[_j]
									if _im >= 0 {
										kmerPos = uint64(start+_j) << 1
										ok = true
										break
									}
								}

								// strand -
								kmer = (*kmerList)[(_j<<1)+1]
								// if kmer != 0 &&
								// 	!util.MustKmerHasSuffix(kmer, tttSuffix, k8, lenSuffix) &&
								// 	!util.MustKmerHasPrefix(kmer, tttPrefix, k8, lenPrefix) {
								// if kmer != 0 {
								if kmer != 0 && kmer != ccc && kmer != ggg && kmer != ttt &&
									!util.IsLowComplexityDust(kmer, k8) {
									// if _im, ok = (*kmer2maskidx)[kmer]; ok {
									// 	kmerPos = uint64(start+_j)<<1 | 1
									// 	break
									// }
									_im = (*loc2maskidxRC)[_j]
									if _im >= 0 {
										kmerPos = uint64(start+_j)<<1 | 1
										ok = true
										break
									}
								}
							}
							if ok {
								// fmt.Printf("    uadd: %s at %d (%d)\n", kmers.MustDecode(kmer, k), _j, start+_j)
								knl = (*extraKmers)[_im]
								if knl == nil {
									knl = poolKmerAndLocs.Get().(*[]uint64)
									*knl = (*knl)[:0]
									(*extraKmers)[_im] = knl
								}
								*knl = append(*knl, kmer)
								*knl = append(*knl, kmerPos)

								// fmt.Printf("  ADD to mask %d with %s, from %d\n", _im+1, lexichash.MustDecode(kmer, k8), (kmerPos>>1)+1)

								*extraLocs = append(*extraLocs, int(kmerPos))

								_j += seedDist
								continue
							}

							// it might fail to fill current region of the desert.
							//   1. there's a gap here.
							//   2. it's the interval region between two contigs.

							// fmt.Printf("desert %d: %d-%d, len: %d, region: %d-%d, list size: %d\n",
							// 	iD, pre, pos, d, start, end, end-start+1)

							_j += seedDist
						}

						pre = pos
					}

				}

				if opt.SaveSeedPositions {
					if !opt.DisableDesertFilling {
						// add extra locs
						for _, loc = range *extraLocs {
							*locs = append(*locs, uint32(loc&4294967295))
						}
						sortutil.Uint32s(*locs)
					}

					// add an extra flag so we can skip these seed pairs accrossing interval regions.

					var nRegions int
					checkRegion := skipRegions != nil
					var _i, _j, _ri, _rs int
					var _loc uint32
					if checkRegion {
						nRegions = len(*skipRegions)

						_ri = 0
						_rs = (*skipRegions)[_ri][1] // end position of a interval region
					}

					if !checkRegion {
						for _i, _loc = range *locs {
							(*locs)[_i] = _loc << 1
						}
					} else {
						for _i, _loc = range *locs {
							_j = int(_loc >> 1)

							// fmt.Printf("checkregion: %v, _j:%d, _rs: %d, nRegions: %d\n", checkRegion, _j, _rs, nRegions)

							if checkRegion && _j >= _rs { // this is the first pos after an interval region
								(*locs)[_i] = _loc<<1 | 1 // add a flag

								// fmt.Printf("  the first pos: %d after a region: %d-%d\n", _j, _rs, (*skipRegions)[_ri][1])

								_ri++

								// some short contigs might do not have seeds
								for _ri+1 < nRegions && _j > (*skipRegions)[_ri+1][1] {
									_rs += (*skipRegions)[_ri+1][1]
									_ri++
									if _ri == nRegions-1 {
										checkRegion = false
										break
									}
								}

								if _ri == nRegions { // this is already the last one
									checkRegion = false
								} else {
									_rs = (*skipRegions)[_ri][1]
								}
							} else {
								(*locs)[_i] = _loc << 1
								// fmt.Printf("  %d do not have the flag\n", _j)
							}
						}
					}
				}

				if !opt.DisableDesertFilling {
					poolInts.Put(extraLocs)
					poolKmerKmerRC.Put(kmerList)
					// poolKmer2MaskIdx.Put(kmer2maskidx)
					poolLoc2MaskIdx.Put(loc2maskidx)
					poolLoc2MaskIdx.Put(loc2maskidxRC)
				}

				// recycle
				if skipRegions != nil {
					poolSkipRegions.Put(skipRegions)
				}

				// send to collect data
				genomes <- refseq

			}(refseq)
		}

		doneMask <- 1
	}()

	// 1.1) parsing input genome files
	iFileBase := batch * opt.GenomeBatchSize
	for iFile, file := range files {
		tokens <- 1
		wg.Add(1)

		go func(file string, iFile int) {
			defer func() {
				wg.Done()
				<-tokens
			}()
			startTime := time.Now()

			if debug {
				log.Debugf("batch: %d, file #%d: %s, begin to parse file", batch, iFile, file)
				log.Info()
			}

			k := lh.K
			// k8 := uint8(lh.K)

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

			maxGenomeSize := opt.MaxGenomeSize

			// object for storing the genome data
			refseq := genome.PoolGenome.Get().(*genome.Genome)
			refseq.Reset()

			refseq.GenomeIdx = iFile

			refseq.StartTime = startTime

			minSeqLen := opt.MinSeqLen
			if minSeqLen <= 0 {
				minSeqLen = k
			}
			minSeqLen = max(minSeqLen, k)

			var i int = 0
			var chunks int = 1
			var seqSize int
			for {
				record, err = fastxReader.Read()
				if err != nil {
					if err == io.EOF {
						break
					}
					checkError(fmt.Errorf("read seq %d in %s: %s", i, file, err))
					break
				}

				// filter out sequences shorter than k or minSeqLen
				if len(record.Seq.Seq) < minSeqLen {
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

				// --------------------------------------------------------

				// check the length of the concatenated sequence, and decide whether to split the genome into multiple chunks.
				//   1. each chunk only saves corresponding sequence meta data: id, length
				//   2. other genome information is the same to the whole genome
				//
				// further processing:
				//   1. use an extra file to store the batch+genome ID list of these split genome chunks.
				//      How: use a map: map[genome id][]batch+ref-index
				//      File format: #chunks, batch+ref indexes.
				//      This format can be simply be concatenated in the merging step
				//
				//   2. after alignment, check if there are genome chunks from the same reference genomes
				//      If yes, merge them, and recompute qcovGnm.
				//        use a map, map[id1][id2]interface{}, map[id2][id1]interface{}

				seqSize = len(refseq.Seq) + len(record.Seq.Seq)
				if seqSize > maxGenomeSize {
					if len(record.Seq.Seq) > maxGenomeSize { // the current sequence is larger than maxGenomeSize
						if outputBigGenomes {
							chBG <- file + "\t" + TOO_LARGE_GENOME + "\n"
						}
						genome.PoolGenome.Put(refseq) // important
						if opt.Verbose || opt.Log2File {
							log.Warningf("  skipping a big genome with a sequence of %d bp: %s", len(record.Seq.Seq), file)
							if !opt.Log2File {
								log.Info()
							}
						}
						if opt.Verbose {
							chDuration <- time.Microsecond // important, or the progress bar will get hung
						}
						return
					}

					refseq.NumSeqs = i

					// ---------------

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

					if debug {
						log.Debugf("batch: %d, file #%d: %s, send split genome: %s", batch, iFile, file, genomeID)
						log.Info()
					}

					// send to mask
					genomesMask <- refseq

					chunks++

					// ---------------

					refseq = genome.PoolGenome.Get().(*genome.Genome)
					refseq.Reset()

					refseq.GenomeIdx = iFile

					// set a zero time, and won't send the time to the progress bar
					refseq.StartTime = time.Time{}

					// ---------------

					i = 0
				}

				// --------------------------------------------------------

				if i > 0 { // add N's between two contigs
					refseq.Seq = append(refseq.Seq, nnn...)
					refseq.Len += len(nnn)
				}
				refseq.Seq = append(refseq.Seq, record.Seq.Seq...)
				refseq.Len += len(record.Seq.Seq)

				// sizes of all contigs
				refseq.SeqSizes = append(refseq.SeqSizes, len(record.Seq.Seq))
				// ids of all contigs
				seqid := []byte(string(record.ID))
				refseq.SeqIDs = append(refseq.SeqIDs, &seqid)
				refseq.GenomeSize += len(record.Seq.Seq)

				i++
			}
			refseq.NumSeqs = i

			if len(refseq.Seq) == 0 {
				if outputBigGenomes {
					chBG <- file + "\t" + NO_VALID_SEQS + "\n"
				}
				genome.PoolGenome.Put(refseq) // important
				if opt.Verbose || opt.Log2File {
					log.Warningf("  skipping %s: no valid sequences", file)
					log.Info()
				}
				if opt.Verbose {
					chDuration <- time.Microsecond // important, or the progress bar will get hung
				}
				return
			}

			// if refseq.GenomeSize > opt.MaxGenomeSize {
			// 	if outputBigGenomes {
			// 		chBG <- file + "\t" + TOO_LARGE_GENOME + "\n"
			// 	}
			// 	genome.PoolGenome.Put(refseq) // important
			// 	if opt.Verbose || opt.Log2File {
			// 		log.Warningf("skipping a big genome (%d bp, %d sequences): %s", refseq.GenomeSize, refseq.NumSeqs, file)
			// 		if !opt.Log2File {
			// 			log.Info()
			// 		}
			// 	}
			// 	if opt.Verbose {
			// 		chDuration <- time.Microsecond // important, or the progress bar will get hung
			// 	}
			// 	return
			// }

			// if len(refseq.Seq) > MAX_GENOME_SIZE {
			// 	if outputBigGenomes {
			// 		chBG <- file + "\t" + TOO_MANY_SEQS + "\n"
			// 	}
			// 	genome.PoolGenome.Put(refseq) // important
			// 	if opt.Verbose || opt.Log2File {
			// 		log.Warningf("skipping a genome with lots of sequences (%d bp, %d sequences): %s", refseq.GenomeSize, refseq.NumSeqs, file)
			// 		if !opt.Log2File {
			// 			log.Info()
			// 		}
			// 	}
			// 	if opt.Verbose {
			// 		chDuration <- time.Microsecond // important, or the progress bar will get hung
			// 	}
			// 	return

			// 	// can split the genome (with thousands of contigs) into several chunks.
			// }

			if chunks > 1 && (opt.Verbose || opt.Log2File) {
				log.Warningf("  splitting a big genome into %d chunks: %s", chunks, file)
				if !opt.Log2File {
					log.Info()
				}
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

			if debug {
				if chunks > 1 {
					log.Debugf("batch: %d, file #%d: %s, send split genome: %s", batch, iFile, file, genomeID)
				} else {
					log.Debugf("batch: %d, file #%d: %s, send genome: %s", batch, iFile, file, genomeID)
				}
				log.Info()
			}

			// send to mask
			genomesMask <- refseq

		}(file, iFileBase+iFile)
	}

	wg.Wait() // all infiles are parsed

	close(genomesMask)
	<-doneMask // all genomes are masked and sent to collect k-mers

	wgMask.Wait() // all genomes are masked

	close(genomes)
	<-done // all k-mer data are collected

	<-doneGW // all genome data are saved
	checkError(gw.Close())
	if opt.SaveSeedPositions {
		checkError(locw.Close())
	}

	// genome chunk lists
	fileGenomeChunks := filepath.Join(outdir, FileGenomeChunks)
	fhGC, err := os.Create(fileGenomeChunks)
	if err != nil {
		checkError(fmt.Errorf("%s", err))
	}
	bw := bufio.NewWriter(fhGC)
	buf := make([]byte, 8)
	var batchIDAndRefID uint64
	for _, li := range mGenomeChunks {
		if len(*li) <= 1 {
			continue
		}

		be.PutUint64(buf, uint64(len(*li)))
		bw.Write(buf)
		for _, batchIDAndRefID = range *li {
			be.PutUint64(buf, batchIDAndRefID)
			bw.Write(buf)
		}
	}
	bw.Flush()
	checkError(fhGC.Close())

	// process bar
	if opt.Verbose {
		close(chDuration)
		<-doneDuration
		pbs.Wait()
	}

	// --------------------------------
	// 4) Summary file
	doneInfo := make(chan int)
	go func() {
		// index summary
		info := &IndexInfo{
			MainVersion:  MainVersion,
			MinorVersion: MinorVersion,

			K:        uint8(lh.K),
			Masks:    len(lh.Masks),
			RandSeed: lh.Seed,

			MaxDesert:        int(opt.DesertMaxLen),
			SeedDistInDesert: int(opt.DesertExpectedSeedDist),

			Chunks:     opt.Chunks,
			Partitions: opt.Partitions,

			InputGenomes:    len(mGenomeChunks), // original genome number. TODO
			Genomes:         nFiles,
			GenomeBatchSize: nFiles, // just for this batch
			GenomeBatches:   1,      // just for this batch
			ContigInterval:  opt.ContigInterval,
		}
		err = writeIndexInfo(filepath.Join(outdir, FileInfo), info)
		if err != nil {
			checkError(fmt.Errorf("failed to write index summary: %s", err))
		}

		doneInfo <- 1
	}()

	// --------------------------------
	// 3) write k-mer data to files

	timeStart2 := time.Now()
	if opt.Verbose || opt.Log2File {
		log.Infof("  writing seeds...")
	}
	chunks := opt.Chunks
	nMasks := opt.Masks
	chunkSize := (nMasks + chunks - 1) / opt.Chunks
	var j, begin, end int
	k8 := uint8(lh.K)
	tokens = make(chan int, opt.MergeThreads)

	for j = 0; j < chunks; j++ { // each chunk
		begin = j * chunkSize
		end = begin + chunkSize
		if end > nMasks {
			end = nMasks
		}
		// fmt.Printf("chunk %d, masks: [%d, %d)\n", j, begin, end)
		wg.Add(1)
		tokens <- 1
		go func(chunk, begin, end int) { // a chunk of masks
			file := filepath.Join(dirSeeds, chunkFile(chunk))

			// for m, data := range (*datas)[begin:end] {
			// 	for key, values := range data {
			// 		for _, v := range *values {
			// 			fmt.Printf("mask: %d, key: %s, value: %064b\n",
			// 				m, kmers.MustDecode(key, lh.K), v)
			// 		}
			// 	}
			// }

			_, err := kv.WriteKVData(k8, begin, (*datas)[begin:end], file, uint8(maskPrefix), uint8(anchorPrefix))
			if err != nil {
				checkError(fmt.Errorf("failed to write seeds data: %s", err))
			}

			// if opt.Verbose || opt.Log2File {
			// 	log.Infof("    seeds file of chunk %d saved", chunk)
			// }

			wg.Done()
			<-tokens
		}(j, begin, end)
	}
	wg.Wait() // all k-mer-value data are saved.
	if opt.Verbose || opt.Log2File {
		log.Infof("  finished writing seeds in %s", time.Since(timeStart2))
	}

	<-doneInfo // info file

	// rare case: no valid genome in a batch: nFiles == 0
	return chunks, nFiles > 0
}

// IndexInfo contains summary of the index
type IndexInfo struct {
	MainVersion      uint8 `toml:"main-version" comment:"Index format"`
	MinorVersion     uint8 `toml:"minor-version"`
	K                uint8 `toml:"max-K" comment:"LexicHash"`
	Masks            int   `toml:"masks"`
	RandSeed         int64 `toml:"rand-seed"`
	MaxDesert        int   `toml:"max-seed-dist" comment:"Seed distance"`
	SeedDistInDesert int   `toml:"seed-dist-in-desert"`
	Chunks           int   `toml:"chunks" comment:"Seeds (k-mer-value data) files"`
	Partitions       int   `toml:"index-partitions"`
	InputGenomes     int   `toml:"input-genomes" comment:"Input genomes"`
	Genomes          int   `toml:"genomes" comment:"Genome data. \n'genomes' might be larger than 'input-genomes', as some big fragmented genomes are split into multiple chunks.\nIn this case, 'genome-batch-size' is not accurate, being variable in different batches."`
	GenomeBatchSize  int   `toml:"genome-batch-size"`
	GenomeBatches    int   `toml:"genome-batches"`
	ContigInterval   int   `toml:"contig-interval"`
}

// writeIndexInfo writes summary of one index
func writeIndexInfo(file string, info *IndexInfo) error {
	fh, err := os.Create(file)
	if err != nil {
		return err
	}

	data, err := toml.Marshal(info)
	if err != nil {
		log.Fatalf("error: %v", err)
	}

	fh.Write(data)

	return fh.Close()
}

// readIndexInfo reads summary frm a file
func readIndexInfo(file string) (*IndexInfo, error) {
	data, err := os.ReadFile(file)
	if err != nil {
		return nil, err
	}

	v := &IndexInfo{}
	err = toml.Unmarshal(data, v)
	return v, err
}

var poolSkipRegions = &sync.Pool{New: func() interface{} {
	tmp := make([][2]int, 0, 128)
	return &tmp
}}

// readGenomeMapIdx2Name reads genome-index mapping file
func readGenomeMapIdx2Name(file string) (map[uint64][]byte, error) {
	fh, err := os.Open(file)
	if err != nil {
		return nil, err
	}
	defer fh.Close()

	r := bufio.NewReader(fh)
	m := make(map[uint64][]byte, 1024)

	buf := make([]byte, 8)
	var n, lenID int
	var batchIDAndRefID uint64
	for {
		n, err = io.ReadFull(r, buf[:2])
		if err != nil {
			if err == io.EOF {
				break
			}
			return nil, err
		}
		if n < 2 {
			return nil, fmt.Errorf("broken genome map file")
		}
		lenID = int(be.Uint16(buf[:2]))
		id := make([]byte, lenID)

		n, err = io.ReadFull(r, id)
		if err != nil {
			return nil, err
		}
		if n < lenID {
			return nil, fmt.Errorf("broken genome map file")
		}

		n, err = io.ReadFull(r, buf)
		if err != nil {
			return nil, err
		}
		if n < 8 {
			return nil, fmt.Errorf("broken genome map file")
		}

		batchIDAndRefID = be.Uint64(buf)

		m[batchIDAndRefID] = id
	}
	return m, nil
}

// readGenomeMap reads genome-index mapping file
func readGenomeMapName2Idx(file string) (map[string]*[]uint64, error) {
	fh, err := os.Open(file)
	if err != nil {
		return nil, err
	}
	defer fh.Close()

	r := bufio.NewReader(fh)
	m := make(map[string]*[]uint64, 1024)

	buf := make([]byte, 8)
	var n, lenID int
	var batchIDAndRefID uint64
	var li *[]uint64
	var ok bool
	for {
		n, err = io.ReadFull(r, buf[:2])
		if err != nil {
			if err == io.EOF {
				break
			}
			return nil, err
		}
		if n < 2 {
			return nil, fmt.Errorf("broken genome map file")
		}
		lenID = int(be.Uint16(buf[:2]))
		id := make([]byte, lenID)

		n, err = io.ReadFull(r, id)
		if err != nil {
			return nil, err
		}
		if n < lenID {
			return nil, fmt.Errorf("broken genome map file")
		}

		n, err = io.ReadFull(r, buf)
		if err != nil {
			return nil, err
		}
		if n < 8 {
			return nil, fmt.Errorf("broken genome map file")
		}

		batchIDAndRefID = be.Uint64(buf)

		if li, ok = m[string(id)]; !ok {
			tmp := make([]uint64, 0, 1)
			li = &tmp
			m[string(id)] = &tmp
		}
		*li = append(*li, batchIDAndRefID)
	}
	return m, nil
}

// readGenomeList reads genome-index mapping file and return the list of genomes
func readGenomeList(file string) ([]string, error) {
	fh, err := os.Open(file)
	if err != nil {
		return nil, err
	}
	defer fh.Close()

	r := bufio.NewReader(fh)
	m := make([]string, 0, 1024)

	buf := make([]byte, 8)
	var n, lenID int
	for {
		n, err = io.ReadFull(r, buf[:2])
		if err != nil {
			if err == io.EOF {
				break
			}
			return nil, err
		}
		if n < 2 {
			return nil, fmt.Errorf("broken genome map file")
		}
		lenID = int(be.Uint16(buf[:2]))
		id := make([]byte, lenID)

		n, err = io.ReadFull(r, id)
		if err != nil {
			return nil, err
		}
		if n < lenID {
			return nil, fmt.Errorf("broken genome map file")
		}

		n, err = io.ReadFull(r, buf)
		if err != nil {
			return nil, err
		}
		if n < 8 {
			return nil, fmt.Errorf("broken genome map file")
		}

		m = append(m, string(id))
	}
	return m, nil
}

// readGenomeChunksMapBig2Small reads the genome chunkfile and return a map
// with bigger batch+ref index mapping to a smaller one.
func readGenomeChunksMapBig2Small(file string) (map[uint64]map[uint64]interface{}, error) {
	fh, err := os.Open(file)
	if err != nil {
		if os.IsNotExist(err) { // no file
			return nil, nil
		}
		return nil, err
	}
	defer fh.Close()

	r := bufio.NewReader(fh)
	data := make(map[uint64]map[uint64]interface{}, 1024)

	buf := make([]byte, 8)
	var n, chunks, i, j int
	var a, b uint64

	list := make([]uint64, 0, 1024)
	for {
		n, err = io.ReadFull(r, buf)
		if err != nil {
			if err == io.EOF {
				break
			}
			return nil, err
		}
		if n < 8 {
			return nil, fmt.Errorf("broken genome chunk file")
		}

		chunks = int(be.Uint64(buf))

		list = list[:0]

		for i = 0; i < chunks; i++ {
			n, err = io.ReadFull(r, buf)
			if err != nil {
				return nil, err
			}
			if n < 8 {
				return nil, fmt.Errorf("broken genome chunk file")
			}

			a = be.Uint64(buf)

			m := make(map[uint64]interface{}, max(8, chunks))

			if i > 0 {
				for j = 0; j < i; j++ {
					b = list[j]

					// data[a][b] = struct{}{} // a > b
					m[b] = struct{}{} // a > b, and no need to sort
				}
			}

			data[a] = m

			list = append(list, a)
		}
	}
	return data, nil
}

// readGenomeChunksMapBig2Small reads the genome chunkfile and return a map
// with all batch+ref indexes of the same genome mapping to the same list.
func readGenomeChunksLists(file string) ([][]uint64, error) {
	fh, err := os.Open(file)
	if err != nil {
		if os.IsNotExist(err) { // no file
			return nil, nil
		}
		return nil, err
	}
	defer fh.Close()

	r := bufio.NewReader(fh)

	data := make([][]uint64, 0, 1024)

	buf := make([]byte, 8)
	var n, chunks, i int
	var idx uint64

	for {
		n, err = io.ReadFull(r, buf)
		if err != nil {
			if err == io.EOF {
				break
			}
			return nil, err
		}
		if n < 8 {
			return nil, fmt.Errorf("broken genome chunk file")
		}

		chunks = int(be.Uint64(buf))

		m := make([]uint64, 0, max(8, chunks))

		for i = 0; i < chunks; i++ {
			n, err = io.ReadFull(r, buf)
			if err != nil {
				return nil, err
			}
			if n < 8 {
				return nil, fmt.Errorf("broken genome chunk file")
			}

			idx = be.Uint64(buf)

			m = append(m, idx)

		}

		data = append(data, m)
	}

	return data, nil
}

// readGenomeChunksMap reads the genome chunkfile and return a map
// with batch+ref index as the key.
func readGenomeChunksMap(file string) (map[uint64]interface{}, error) {
	fh, err := os.Open(file)
	if err != nil {
		if os.IsNotExist(err) { // no file
			return nil, nil
		}
		return nil, err
	}
	defer fh.Close()

	r := bufio.NewReader(fh)
	data := make(map[uint64]interface{}, 1024)

	buf := make([]byte, 8)
	var n, chunks, i int
	var a uint64

	for {
		n, err = io.ReadFull(r, buf)
		if err != nil {
			if err == io.EOF {
				break
			}
			return nil, err
		}
		if n < 8 {
			return nil, fmt.Errorf("broken genome chunk file")
		}

		chunks = int(be.Uint64(buf))

		for i = 0; i < chunks; i++ {
			n, err = io.ReadFull(r, buf)
			if err != nil {
				return nil, err
			}
			if n < 8 {
				return nil, fmt.Errorf("broken genome chunk file")
			}

			a = be.Uint64(buf)

			data[a] = struct{}{}
		}
	}
	return data, nil
}

var poolPrefix2Kmers = &sync.Pool{New: func() interface{} {
	tmp := make([]*[4]uint64, 0, 1024)
	return &tmp
}}

// prefix1, kmer, prefixRC, kmerRC
var poolPrefix2Kmer = &sync.Pool{New: func() interface{} {
	return &[4]uint64{}
}}

var poolPrefxCounter = &sync.Pool{New: func() interface{} {
	tmp := make(map[uint64]uint32)
	return &tmp
}}

var poolKmer2MaskIdx = &sync.Pool{New: func() interface{} {
	tmp := make(map[uint64]int, 1024)
	return &tmp
}}

var poolLoc2MaskIdx = &sync.Pool{New: func() interface{} {
	tmp := make([]int, 1024)
	return &tmp
}}

// kmer, kmerRC
var poolKmerKmerRC = &sync.Pool{New: func() interface{} {
	tmp := make([]uint64, 1024)
	return &tmp
}}

var poolKmerAndLocs = &sync.Pool{New: func() interface{} {
	tmp := make([]uint64, 0, 128)
	return &tmp
}}

var poolInts = &sync.Pool{New: func() interface{} {
	tmp := make([]int, 0, 1024)
	return &tmp
}}

var cmpFn = func(x, y int) int { return int(x - y) }
