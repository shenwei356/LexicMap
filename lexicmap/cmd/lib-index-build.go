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
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OFTestSerializationTestSerialization ANY KIND, EXPRESS OR
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

	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/lexichash"
	"github.com/vbauerster/mpb/v8"
	"github.com/vbauerster/mpb/v8/decor"
)

const FileMasks = "masks.bin"
const DirSeeds = "seeds"
const DirGenomes = "genomes"
const FileGenomes = "genomes.bin"

// TmpDirExt is the path extension for temporary files
const TmpDirExt = ".tmp"

type IndexBuildingOptions struct {
	// general
	NumCPUs      int
	Verbose      bool // show log
	Log2File     bool
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

	ReRefName    *regexp.Regexp
	ReSeqExclude []*regexp.Regexp
}

// CheckIndexBuildingOptions checks the important options
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

	if opt.Chunks < 1 || opt.Chunks > 512 {
		return fmt.Errorf("invalid chunks: %d, valid range: [1, 512]", opt.Chunks)
	}

	if opt.Partitions < 1 {
		return fmt.Errorf("invalid numer of partitions in indexing k-mer data: %d, should be >=1", opt.Partitions)
	}

	if opt.GenomeBatchSize < 1 || opt.GenomeBatchSize > 1<<17 {
		return fmt.Errorf("invalid genome batch size: %d, valid range: [1, 131072]", opt.GenomeBatchSize)
	}

	if opt.MaxOpenFiles < 2 {
		return fmt.Errorf("invalid max open files: %d, should be >= 2", opt.MaxOpenFiles)
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
	tmpDir := filepath.Clean(outdir) + TmpDirExt
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
			outdirB = filepath.Join(tmpDir, fmt.Sprintf("batch_%4d", batch))
			tmpIndexes = append(tmpIndexes, outdirB)
		} else {
			outdirB = outdir
		}

		// build index for this batch
		buildAnIndex(lh, opt, poolKmerDatas, outdirB, files, batch)
	}

	if nBatches == 1 {
		return nil
	}

	// merge indexes
	mergeIndexes(lh, opt, outdir, tmpIndexes)

	return nil
}

// build an index for the files of one batch
func buildAnIndex(lh *lexichash.LexicHash, opt *IndexBuildingOptions,
	poolKmerDatas *sync.Pool,
	outdir string, files []string, batch int) {

	var timeStart time.Time
	if opt.Verbose || opt.Log2File {
		timeStart = time.Now()
		log.Info()
		log.Infof("  building index for batch %d with %d files...", batch, len(files))
		defer func() {
			log.Infof("  elapsed time: %s", time.Since(timeStart))
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
				bar.Increment()
				bar.EwmaIncrBy(1, t)
			}
			doneDuration <- 1
		}()
	}

	// -------------------------------------------------------------------
	// dir structure

	// masks
	fileMask := filepath.Join(outdir, FileMasks)
	_, err := lh.WriteToFile(fileMask)
	if err != nil {
		checkError(fmt.Errorf("failed to write masks: %s", err))
	}

	// genomes
	dirGenomes := filepath.Join(outdir, DirGenomes, fmt.Sprintf("batch_%04d", batch))
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

	// write genomes to file
	go func() {
		for refseq := range genomesW { // each genome
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

	// collect k-mer data
	go func() {
		posMask := uint64((1 << 30) - 1)
		var wg sync.WaitGroup
		threads := opt.NumCPUs
		tokens := make(chan int, threads)
		nMasks := opt.Masks
		chunkSize := (nMasks + threads - 1) / threads
		var j, begin, end int

		var refIdx uint32 // genome number

		for refseq := range genomes { // each genome
			genomesW <- refseq // send to save to file, asynchronously writing

			_kmers := refseq.Kmers
			loces := refseq.Locses

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
						data := (*datas)[i]
						kmer = (*_kmers)[i]
						if values, ok = data[kmer]; !ok {
							values = &[]uint64{}
							data[kmer] = values
						}
						for _, loc = range (*loces)[i] {
							//  batch idx: 17 bits
							//  ref idx:   17 bits
							//  pos:       29 bits
							//  strand:     1 bits
							// here, the location from Mask() already contains the strand information.
							value = uint64(batch)<<47 | uint64(refIdx)<<30 | (uint64(loc) & posMask)

							*values = append(*values, value)
						}
					}
					wg.Done()
					<-tokens
				}(begin, end)
			}

			wg.Wait()
			refIdx++
		}
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
				if i > 0 {
					refseq.Seq = append(refseq.Seq, nnn...)
					refseq.Len += len(nnn)
				}
				refseq.Seq = append(refseq.Seq, record.Seq.Seq...)
				refseq.Len += len(record.Seq.Seq)

				refseq.SeqSizes = append(refseq.SeqSizes, len(record.Seq.Seq))
				refseq.GenomeSize += len(record.Seq.Seq)

				i++
			}

			if len(refseq.Seq) == 0 {
				genome.PoolGenome.Put(refseq) // important
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

	// --------------------------------
	// 3) write k-mer data to files

	if opt.Verbose || opt.Log2File {
		log.Infof("  writing seeds...")
	}
	chunks := opt.Chunks
	nMasks := opt.Masks
	chunkSize := (nMasks + chunks - 1) / opt.Chunks
	var j, begin, end int

	for j = 0; j < chunks; j++ { // each chunk
		begin = j * chunkSize
		end = begin + chunkSize
		if end > nMasks {
			end = nMasks
		}

		wg.Add(1)
		tokens <- 1
		go func(chunk, begin, end int) { // a chunk of masks
			file := filepath.Join(dirSeeds, fmt.Sprintf("chunk_%03d.bin", chunk))

			// for m, data := range (*datas)[begin:end] {
			// 	for key, values := range data {
			// 		for _, v := range *values {
			// 			fmt.Printf("mask: %d, key: %s, value: %d\n",
			// 				m, kmers.MustDecode(key, lh.K), v)
			// 		}
			// 	}
			// }

			_, err := kv.WriteKVData(uint8(k), begin, (*datas)[begin:end], file, opt.Partitions)
			if err != nil {
				checkError(fmt.Errorf("failed to write seeds data: %s", err))
			}

			wg.Done()
			<-tokens
		}(j, begin, end)
	}
	wg.Wait() // all k-mer-value data are saved.

	poolKmerDatas.Put(datas)

	// -------------------------------------------------------------------
	if opt.Verbose {
		close(chDuration)
		<-doneDuration
		pbs.Wait()
	}
}

var poolSkipRegions = &sync.Pool{New: func() interface{} {
	tmp := make([][2]int, 0, 128)
	return &tmp
}}

func mergeIndexes(lh *lexichash.LexicHash, opt *IndexBuildingOptions, outdir string, paths []string) {

}
