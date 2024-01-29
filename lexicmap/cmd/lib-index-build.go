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
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/lexichash"
	"github.com/vbauerster/mpb/v8"
	"github.com/vbauerster/mpb/v8/decor"
)

var be = binary.BigEndian

// MainVersion is use for checking compatibility
var MainVersion uint8 = 0

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

	// LexicHash

	K                int   // k-mer size
	Masks            int   // number of masks
	RandSeed         int64 // random seed
	PrefixForCheckLC int   // length of prefix for checking low-complexity

	// k-mer-value data

	Chunks     int // the number of chunks for storing k-mer data
	Partitions int // the number of partitions for indexing k-mer data

	// genome batches

	GenomeBatchSize int // the maximum number of genomes of a batch

	// genome

	ReRefName    *regexp.Regexp   // for extracting genome id from the file name
	ReSeqExclude []*regexp.Regexp // for excluding sequences according to name pattern
}

// CheckIndexBuildingOptions checks some important options
func CheckIndexBuildingOptions(opt *IndexBuildingOptions) error {
	if opt.K < 3 || opt.K > 32 {
		return fmt.Errorf("invalid k value: %d, valid range: [3, 32]", opt.K)
	}
	if opt.Masks < 4 {
		return fmt.Errorf("invalid numer of masks: %d, should be >=4", opt.Masks)
	}
	if opt.PrefixForCheckLC > opt.K {
		return fmt.Errorf("invalid prefix: %d, valid range: [0, k], 0 for no checking", opt.PrefixForCheckLC)
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

	// generate masks
	lh, err := lexichash.NewWithSeed(opt.K, opt.Masks, opt.RandSeed, opt.PrefixForCheckLC)
	if err != nil {
		return err
	}
	// save mask later

	// it's private variable because the size comes from opt.Masks
	var poolKmerDatas = &sync.Pool{New: func() interface{} {
		datas := make([]map[uint64]*[]uint64, opt.Masks)
		for i := 0; i < opt.Masks; i++ {
			datas[i] = make(map[uint64]*[]uint64, 1024)
		}
		return &datas
	}}

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
		kvChunks = buildAnIndex(lh, opt, poolKmerDatas, outdirB, files, batch)
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
	poolKmerDatas *sync.Pool,
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
	datas := poolKmerDatas.Get().(*[]map[uint64]*[]uint64)
	for _, data := range *datas { // reset all maps
		clear(data)
	}

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

	// 2.2) write genomes to file
	go func() {
		threadsFloat := float64(opt.NumCPUs) // just avoid repeated type conversion
		for refseq := range genomesW {       // each genome
			// write the genome to file
			err = gw.Write(refseq)
			if err != nil {
				checkError(fmt.Errorf("failed to write genome: %s", err))
			}

			// recycle the genome
			if refseq.Kmers != nil {
				lh.RecycleMaskResult(refseq.Kmers, refseq.Locses)
			}
			genome.RecycleGenome(refseq)

			if opt.Verbose {
				chDuration <- time.Duration(float64(time.Since(refseq.StartTime)) / threadsFloat)
			}
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
					for i := begin; i < end; i++ {
						data := (*datas)[i] // the map to save into
						kmer = (*_kmers)[i] // captured k-mer by the mask
						if values, ok = data[kmer]; !ok {
							values = &[]uint64{}
							data[kmer] = values
						}
						for _, loc = range (*loces)[i] { // location information of the captured k-mer
							//  batch idx: 17 bits
							//  ref idx:   17 bits
							//  pos:       29 bits
							//  strand:     1 bits
							// here, the location from Mask() already contains the strand information.
							// value = uint64(batch)<<47 | ((refIdx & 131071) << 30) |
							value = batchIDAndRefIDShift | (uint64(loc) & 1073741823)
							// fmt.Printf("%s, batch: %d, refIdx: %d, value: %064b\n", refseq.ID, batch, refIdx, value)

							*values = append(*values, value)
						}
					}
					wg.Done()
					<-tokens
				}(begin, end)
			}

			wg.Wait() // wait all mask chunks
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
	k := lh.K
	nnn := bytes.Repeat([]byte{'N'}, k-1)
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
				log.Warningf("skipping %s: no valid sequences", file)
				log.Info()
				if opt.Verbose {
					chDuration <- time.Duration(0) // important, or the progress bar will get hung
				}
				return
			}

			var genomeID string // genome id
			if extractRefName {
				if reRefName.MatchString(baseFile) {
					genomeID = reRefName.FindAllStringSubmatch(baseFile, 1)[0][1]
				} else {
					genomeID, _ = filepathTrimExtension(baseFile)
				}
			} else {
				genomeID, _ = filepathTrimExtension(baseFile)
			}

			refseq.ID = []byte(genomeID)
			refseq.StartTime = startTime

			// --------------------------------
			// mask with lexichash

			// skip regions around junctions of two sequences.

			// because lh.Mask accepts a list, while when skipRegions is nil, *skipRegions is illegal.
			var _skipRegions [][2]int
			var skipRegions *[][2]int
			if len(refseq.SeqSizes) > 1 {
				skipRegions = poolSkipRegions.Get().(*[][2]int)
				*skipRegions = (*skipRegions)[:0]
				var n int // len of concatenated seqs
				for i, s := range refseq.SeqSizes {
					if i > 0 {
						*skipRegions = append(*skipRegions, [2]int{n, n + k - 1})
						n += k - 1
					}
					n += s
				}
				_skipRegions = *skipRegions
			}
			kmers, locses, err := lh.Mask(refseq.Seq, _skipRegions)
			if err != nil {
				panic(err)
			}
			refseq.Kmers = kmers
			refseq.Locses = locses
			if skipRegions != nil {
				poolSkipRegions.Put(skipRegions)
			}

			// --------------------------------
			// bit-packed sequences
			refseq.TwoBit = genome.Seq2TwoBit(refseq.Seq)

			genomes <- refseq

		}(file)
	}

	wg.Wait() // all infiles are parsed
	close(genomes)
	<-done // all k-mer data are collected

	<-doneGW // all genome data are saved
	checkError(gw.Close())

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

			Chunks:     opt.Chunks,
			Partitions: opt.Partitions,

			Genomes:         len(files),
			GenomeBatchSize: len(files), // just for this batch
			GenomeBatches:   1,          // just for this batch
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

			_, err := kv.WriteKVData(uint8(k), begin, (*datas)[begin:end], file, opt.Partitions)
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

	poolKmerDatas.Put(datas)

	return chunks
}

// IndexInfo contains summary of the index
type IndexInfo struct {
	MainVersion     uint8 `toml:"main-version" comment:"Index format"`
	MinorVersion    uint8 `toml:"minor-version"`
	K               uint8 `toml:"max-K" comment:"LexicHash"`
	Masks           int   `toml:"masks"`
	RandSeed        int64 `toml:"rand-seed"`
	Chunks          int   `toml:"chunks" comment:"Seeds (k-mer-value data) files"`
	Partitions      int   `toml:"index-partitions"`
	Genomes         int   `toml:"genomes" comment:"Genome data"`
	GenomeBatchSize int   `toml:"genome-batch-size"`
	GenomeBatches   int   `toml:"genome-batches"`
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

// readGenomeMap reads genome-index mapping file
func readGenomeMap(file string) (map[uint64][]byte, error) {
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
