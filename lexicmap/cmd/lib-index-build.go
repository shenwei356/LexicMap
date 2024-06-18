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
	"os"
	"path/filepath"
	"regexp"
	"sync"
	"time"

	"github.com/pelletier/go-toml/v2"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/seedposition"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/lexichash"
	"github.com/shenwei356/lexichash/iterator"
	"github.com/twotwotwo/sorts/sortutil"
	"github.com/vbauerster/mpb/v8"
	"github.com/vbauerster/mpb/v8/decor"

	itree "github.com/rdleal/intervalst/interval"
)

var be = binary.BigEndian

// MainVersion is use for checking compatibility
var MainVersion uint8 = 1

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

	// skipping extremely large genome
	MaxGenomeSize int    // Maximum genome size. Extremely large genomes (non-isolate assemblies) will be skipped
	BigGenomeFile string // Out file of skipped files with genomes

	// LexicHash
	MaskFile string // file of custom masks

	K        int   // k-mer size
	Masks    int   // number of masks
	RandSeed int64 // random seed

	// generate mask randomly
	Prefix int // length of prefix for checking low-complexity and choosing k-mers to fill deserts

	// filling sketching deserts
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
}

// CheckIndexBuildingOptions checks some important options
func CheckIndexBuildingOptions(opt *IndexBuildingOptions) error {
	if opt.K < 3 || opt.K > 32 {
		return fmt.Errorf("invalid k value: %d, valid range: [3, 32]", opt.K)
	}
	if opt.Masks < 64 {
		return fmt.Errorf("invalid numer of masks: %d, should be >=64", opt.Masks)
	}
	if opt.Prefix > opt.K {
		return fmt.Errorf("invalid prefix: %d, valid range: [0, k], 0 for no checking", opt.Prefix)
	}

	if opt.Chunks < 1 || opt.Chunks > 128 {
		return fmt.Errorf("invalid chunks: %d, valid range: [1, 128]", opt.Chunks)
	}

	if opt.Chunks > opt.Masks {
		return fmt.Errorf("invalid chunks: %d, should be <= masks (%d)", opt.Chunks, opt.Masks)
	}

	if opt.Partitions < 1 {
		return fmt.Errorf("invalid numer of partitions in indexing k-mer data: %d, should be >=1", opt.Partitions)
	}

	if opt.GenomeBatchSize < 1 || opt.GenomeBatchSize > 1<<17 {
		return fmt.Errorf("invalid genome batch size: %d, valid range: [1, 131072]", opt.GenomeBatchSize)
	}

	// ------------------------

	if opt.NumCPUs < 1 {
		return fmt.Errorf("invalid number of CPUs: %d, should be >= 1", opt.NumCPUs)
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
		// } else {
		// 	masks, skippedFiles, err := GenerateMasks(infiles, opt, "")
		// 	checkError(err)

		// 	lh, err = lexichash.NewWithMasks(opt.K, masks)
		// 	checkError(err)

		// 	if opt.BigGenomeFile != "" {
		// 		outfh2, err := os.Create(opt.BigGenomeFile)
		// 		if err != nil {
		// 			checkError(fmt.Errorf("failed to write file: %s", opt.BigGenomeFile))
		// 		}
		// 		for _, file := range skippedFiles {
		// 			fmt.Fprintf(outfh2, "%s\n", file)
		// 		}
		// 		outfh2.Close()
		// 		if opt.Verbose || opt.Log2File {
		// 			log.Infof("  finished saving skipped genome files: %s", opt.BigGenomeFile)
		// 		}
		// 	}
		// }
	} else {
		lh, err = lexichash.NewWithSeed(opt.K, opt.Masks, opt.RandSeed, opt.Prefix)
		if err != nil {
			return err
		}
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
	}

	var begin, end int
	var kvChunks int
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
			tmpIndexes = append(tmpIndexes, outdirB)
		} else {
			outdirB = outdir
		}

		// build index for this batch
		kvChunks = buildAnIndex(lh, opt, &datas, outdirB, files, batch)
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
	err = mergeIndexes(lh, opt, kvChunks, outdir, tmpIndexes, tmpDir, 1)
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

// build an index for the files of one batch
func buildAnIndex(lh *lexichash.LexicHash, opt *IndexBuildingOptions,
	datas *[]*map[uint64]*[]uint64,
	outdir string, files []string, batch int) int {

	var timeStart time.Time
	if opt.Verbose || opt.Log2File {
		timeStart = time.Now()
		log.Info()
		log.Infof("  ------------------------[ batch %d ]------------------------", batch)
		log.Infof("  building index for batch %d with %d files...", batch, len(files))
		defer func() {
			log.Infof("  finished building index for batch %d in: %s", batch, time.Since(timeStart))
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

		for refseq := range genomesW { // each genome
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

			// send signal of genome being written
			refseq.Done <- 1
		}
		doneGW <- 1
	}()

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

		for refseq := range genomes { // each genome
			genomesW <- refseq // send to save to file, asynchronously writing

			_kmers := refseq.Kmers
			loces := refseq.Locses
			extraKmers := refseq.ExtraKmers

			// refseq id -> this
			batchIDAndRefID = uint64(batch)<<17 | (refIdx & 131071)
			be.PutUint16(buf[:2], uint16(len(refseq.ID)))
			bw.Write(buf[:2])
			bw.Write(refseq.ID)
			be.PutUint64(buf, batchIDAndRefID)
			bw.Write(buf)

			batchIDAndRefIDShift = batchIDAndRefID << 30

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
					for i := begin; i < end; i++ {
						data := (*datas)[i] // the map to save into

						kmer = (*_kmers)[i]       // captured k-mer by the mask
						if len((*loces)[i]) > 0 { // locations from from MaskKnownPrefixes might be empty.
							if values, ok = (*data)[kmer]; !ok {
								values = &[]uint64{}
								(*data)[kmer] = values
							}
							for _, loc = range (*loces)[i] { // position information of the captured k-mer
								//  batch idx: 17 bits
								//  ref idx:   17 bits
								//  pos:       29 bits
								//  strand:     1 bits
								// here, the position from Mask() already contains the strand information.
								// value = uint64(batch)<<47 | ((refIdx & 131071) << 30) |
								value = batchIDAndRefIDShift | (uint64(loc) & 1073741823)
								// fmt.Printf("%s, batch: %d, refIdx: %d, value: %064b\n", refseq.ID, batch, refIdx, value)

								*values = append(*values, value)
							}
						}

						// -----------------------------
						// extra k-mers
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

							value = batchIDAndRefIDShift | ((*knl)[_j+1] & 1073741823)
							*values = append(*values, value)
						}

					}
					wg.Done()
					<-tokens
				}(begin, end)
			}

			wg.Wait() // wait all mask chunks

			// wait the genome data being written
			<-refseq.Done

			// recycle the genome
			if refseq.Kmers != nil {
				lh.RecycleMaskResult(refseq.Kmers, refseq.Locses)
			}
			genome.RecycleGenome(refseq)

			if opt.Verbose {
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
	nnn := bytes.Repeat([]byte{'N'}, opt.ContigInterval)
	reRefName := opt.ReRefName
	extractRefName := reRefName != nil
	filterNames := len(opt.ReSeqExclude) > 0

	var wg sync.WaitGroup                 // ensure all jobs done
	tokens := make(chan int, opt.NumCPUs) // control the max concurrency number

	for _, file := range files {
		tokens <- 1
		wg.Add(1)

		go func(file string) {
			defer func() {
				wg.Done()
				<-tokens
			}()
			startTime := time.Now()

			k := lh.K
			k8 := uint8(lh.K)

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
			refseq := genome.PoolGenome.Get().(*genome.Genome)
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
				genome.PoolGenome.Put(refseq) // important
				if opt.Verbose || opt.Log2File {
					log.Warningf("skipping %s: no valid sequences", file)
					log.Info()
				}
				if opt.Verbose {
					chDuration <- time.Microsecond // important, or the progress bar will get hung
				}
				return
			}

			if refseq.GenomeSize > opt.MaxGenomeSize {
				genome.PoolGenome.Put(refseq) // important
				if opt.Verbose || opt.Log2File {
					log.Warningf("skipping big genome (%d bp): %s", refseq.GenomeSize, file)
					if !opt.Log2File {
						log.Info()
					}
				}
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
			// mask with lexichash

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
			var _kmers *[]uint64
			var locses *[][]int

			// if len(lh.Masks) > 1024 && len(refseq.Seq) > 1048576 {
			// 	_kmers, locses, err = lh.MaskLongSeqs(refseq.Seq, _skipRegions)
			// } else {
			// 	_kmers, locses, err = lh.Mask(refseq.Seq, _skipRegions)
			// }
			_kmers, locses, err = lh.MaskKnownPrefixes(refseq.Seq, _skipRegions)

			if err != nil {
				panic(err)
			}
			refseq.Kmers = _kmers
			refseq.Locses = locses

			// --------------------------------
			// bit-packed sequences
			refseq.TwoBit = genome.Seq2TwoBit(refseq.Seq)

			// ----------------------------------------------------------------
			// fill sketching deserts

			if refseq.Locs == nil { // existing positions, will be updated with new k-mers
				tmp := make([]uint32, 0, 40<<10)
				refseq.Locs = &tmp
			}
			locs := refseq.Locs
			*locs = (*locs)[:0]

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

			var loc int
			for _, _locs := range *locses { // insert positions
				for _, loc = range _locs {
					*locs = append(*locs, uint32(loc&4294967295)) // only posision | strand flag
				}
			}
			sortutil.Uint32s(*locs)

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
			loc2maskidx := poolLoc2MaskIdx.Get().(*[]int)
			loc2maskidxRC := poolLoc2MaskIdx.Get().(*[]int)
			kmerList := poolKmerKmerRC.Get().(*[]uint64)

			var iter *iterator.Iterator
			var kmer, kmerRC, kmerPos uint64
			var ok bool
			var _j, posOfPre, posOfCur, _start, _end int

			var _im int
			var knl *[]uint64

			var inIntervalRegion bool

			// it's hard to skip k-mers at the edge of a gap
			//       ACCAAAAAAAA      AAAAAAGCCAGA
			// ----------AAAAAAAAAAAAAAAAAAA-------
			// ----------TTTTTTTTTTTTTTTTTTT-------
			//      TTTTTTAC         ACGTTTTTTT
			// so we just detect these k-mers by it's suffix
			var lenSuffix uint8 = 5
			var lenPrefix uint8 = 3                          // that means the firt n base can't be As, just in case.
			tttPrefix := uint64((1 << (lenPrefix << 1)) - 1) // for k-mer on the negative strand
			tttSuffix := uint64((1 << (lenSuffix << 1)) - 1) // for k-mer on the negative strand

			extraLocs := poolInts.Get().(*[]int) // extra positions
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
				_kmers2, _locses2, _ = lh.MaskKnownPrefixes(refseq.Seq[start:end], nil)

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
						if kmer != 0 &&
							!util.MustKmerHasSuffix(kmer, 0, k8, lenSuffix) &&
							!util.MustKmerHasPrefix(kmer, 0, k8, lenPrefix) {
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
						if kmer != 0 &&
							!util.MustKmerHasSuffix(kmer, tttSuffix, k8, lenSuffix) &&
							!util.MustKmerHasPrefix(kmer, tttPrefix, k8, lenPrefix) {
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
						if kmer != 0 &&
							!util.MustKmerHasSuffix(kmer, 0, k8, lenSuffix) &&
							!util.MustKmerHasPrefix(kmer, 0, k8, lenPrefix) {
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
						if kmer != 0 &&
							!util.MustKmerHasSuffix(kmer, tttSuffix, k8, lenSuffix) &&
							!util.MustKmerHasPrefix(kmer, tttPrefix, k8, lenPrefix) {
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

			if opt.SaveSeedPositions {
				// add extra locs
				for _, loc = range *extraLocs {
					*locs = append(*locs, uint32(loc&4294967295))
				}
				sortutil.Uint32s(*locs)

				// add an extra flag so we can skip these seed pairs accrossing interval regions.

				var nRegions int
				checkRegion := len(refseq.SeqSizes) > 1
				var _i, _j, _ri, _rs int
				var _loc uint32
				if checkRegion {
					nRegions = len(*skipRegions)

					_ri = 0
					_rs = (*skipRegions)[_ri][1] // start  of a interval region
				}

				if !checkRegion {
					for _i, _loc = range *locs {
						(*locs)[_i] = _loc << 1
					}
				} else {
					for _i, _loc = range *locs {
						_j = int(_loc >> 1)

						if checkRegion && _j >= _rs { // this is the first pos after an interval region
							(*locs)[_i] = _loc<<1 | 1 // add a flag

							// fmt.Printf("the first pos: %d after a region: %d-%d\n", _j, _rs, (*skipRegions)[_ri][1])

							_ri++
							if _ri == nRegions { // this is already the last one
								checkRegion = false
							} else {
								_rs = (*skipRegions)[_ri][1]
							}
						} else {
							(*locs)[_i] = _loc << 1
						}
					}
				}
			}

			poolInts.Put(extraLocs)

			poolKmerKmerRC.Put(kmerList)
			// poolKmer2MaskIdx.Put(kmer2maskidx)
			poolLoc2MaskIdx.Put(loc2maskidx)
			poolLoc2MaskIdx.Put(loc2maskidxRC)

			if skipRegions != nil {
				poolSkipRegions.Put(skipRegions)
			}

			genomes <- refseq

		}(file)
	}

	wg.Wait() // all infiles are parsed
	close(genomes)
	<-done // all k-mer data are collected

	<-doneGW // all genome data are saved
	checkError(gw.Close())
	if opt.SaveSeedPositions {
		checkError(locw.Close())
	}

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
	tokens = make(chan int, 1)   // hope it reduces memory
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

			_, err := kv.WriteKVData(k8, begin, (*datas)[begin:end], file, opt.Partitions)
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

	return chunks
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
	Genomes          int   `toml:"genomes" comment:"Genome data"`
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
func readGenomeMapName2Idx(file string) (map[string]uint64, error) {
	fh, err := os.Open(file)
	if err != nil {
		return nil, err
	}
	defer fh.Close()

	r := bufio.NewReader(fh)
	m := make(map[string]uint64, 1024)

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

		m[string(id)] = batchIDAndRefID
	}
	return m, nil
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
