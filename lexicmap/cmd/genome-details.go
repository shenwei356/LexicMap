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
	"bytes"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"sync"
	"time"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
)

var genomeDetailsCmd = &cobra.Command{
	Use:   "genome-details",
	Short: "Extract or view genome details in the index",
	Long: `Extract or view genome details in the index

On the first run, this command will extract genome details and save them to 'genomes.details.bin'.
If the file exists, it will be read directly and details will be printed.

Output format:
  Tab-delimited format.

    ref             genome id
    genome_size     genome size (sum of all genome chunks)
    chunks          the number of genome chunks
    chunk           nth genome chunk
    chunk_size      genome (chunk) size
    seqs            the number of sequences in the genome (chunk)
    seqsizes        comma-separated sequence sizes in the genome (chunk)    (optional with -e/--extra)
    seqids          comma-separated sequence ids in the genome (chunk)      (optional with -e/--extra)
                    only available when the genome details file is created with the -i/--save-seqids flag.

  Note that genome chunks are created when a genome is too large, and the chunk size is determined by the
  "-g/--max-genome" parameter in "lexicmap index". If a genome is not chunked, it will be treated as one chunk.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		seq.ValidateSeq = false

		// ------------------------------

		dbDir := getFlagString(cmd, "index")
		if dbDir == "" {
			checkError(fmt.Errorf("flag -d/--index needed"))
		}

		outFile := getFlagString(cmd, "out-file")

		extra := getFlagBool(cmd, "extra")

		saveSeqIDs := getFlagBool(cmd, "save-seqids")

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

		// -----------------------------------------------------
		// Extract genome details if not existed, otherwise read it directly.

		fileGenomeDetails := filepath.Join(dbDir, FileGenomeDetails)
		existed, err := pathutil.Exists(fileGenomeDetails)
		checkError(err)
		if !existed {
			if opt.Verbose {
				log.Infof("extracting genome details and saving to %s", fileGenomeDetails)
			}
			timeStart := time.Now()
			checkError(extractGenomeDetails(opt, dbDir, saveSeqIDs))
			if opt.Verbose {
				log.Infof("  elapsed time: %s", time.Since(timeStart))
				log.Info()
			}
		}

		// -----------------------------------------------------
		if opt.Verbose {
			log.Infof("reading genome details from %s", fileGenomeDetails)
		}
		timeStart1 := time.Now()
		checkError(readGenomeDetails(fileGenomeDetails, outfh, extra))
		if opt.Verbose {
			log.Infof("  elapsed time: %s", time.Since(timeStart1))
		}
	},
}

func init() {
	utilsCmd.AddCommand(genomeDetailsCmd)

	genomeDetailsCmd.Flags().StringP("index", "d", "",
		formatFlagUsage(`Index directory created by "lexicmap index".`))

	genomeDetailsCmd.Flags().StringP("out-file", "o", "-",
		formatFlagUsage(`Out file, supports the ".gz" suffix ("-" for stdout).`))

	genomeDetailsCmd.Flags().BoolP("save-seqids", "i", false,
		formatFlagUsage(`Extract and save sequence ids. This will increase the file size`))

	genomeDetailsCmd.Flags().BoolP("extra", "e", false,
		formatFlagUsage(`Show extra columns, including seqsizes and seqids.`))

	genomeDetailsCmd.SetUsageTemplate(usageTemplate(""))
}

// FileGenomeDetails store lists of genome details.
//
//	genome id:
//	    chunk 1:
//	       gsize
//	       n_seqs
//	           size 1, len_seqid1, seqid1
//	           ...
const FileGenomeDetails = "genomes.details.bin"
const FLAG_SAVE_SEQIDS = 1

func extractGenomeDetails(opt *Options, dbDir string, saveSeqIDs bool) error {

	// ---------------------------------------------------------------
	// info file
	if opt.Verbose {
		log.Infof("  reading index info file")
	}
	fileInfo := filepath.Join(dbDir, FileInfo)
	info, err := readIndexInfo(fileInfo)
	if err != nil {
		return fmt.Errorf("failed to read info file: %s", err)
	}
	if info.MainVersion != MainVersion {
		return fmt.Errorf("index main versions do not match: %d (index) != %d (tool). please re-create the index", info.MainVersion, MainVersion)
	}

	// ---------------------------------------------------------------
	// genome readers
	nReaders := 1

	if opt.Verbose {
		log.Infof("  creating reader pools for %d genome batches, each with %d reader(s)...", info.GenomeBatches, nReaders)
	}
	poolGenomeRdrs := make([]chan *genome.Reader, info.GenomeBatches)
	for i := 0; i < info.GenomeBatches; i++ {
		poolGenomeRdrs[i] = make(chan *genome.Reader, nReaders)
	}

	// parallelize it
	var wg sync.WaitGroup
	tokens := make(chan int, opt.NumCPUs)
	for i := 0; i < info.GenomeBatches; i++ {
		for j := 0; j < nReaders; j++ {
			tokens <- 1
			wg.Add(1)
			go func(i int) {
				fileGenomes := filepath.Join(dbDir, DirGenomes, batchDir(i), FileGenomes)
				rdr, err := genome.NewReader(fileGenomes)
				if err != nil {
					checkError(fmt.Errorf("failed to create genome reader: %s", err))
				}
				poolGenomeRdrs[i] <- rdr

				wg.Done()
				<-tokens
			}(i)
		}
	}
	wg.Wait()

	// ---------------------------------------------------------------
	// read genome chunks data if existed
	if opt.Verbose {
		log.Infof("  reading genome chunk data files")
	}
	genomeChunks, err := readGenomeChunksLists(filepath.Join(dbDir, FileGenomeChunks))
	if err != nil {
		checkError(fmt.Errorf("failed to read genome chunk file: %s", err))
	}
	var hasGenomeChunks bool
	var genomeChunkFlags map[uint64]*bool
	var genomeChunkLists map[uint64]*[]uint64

	if len(genomeChunks) > 0 {
		hasGenomeChunks = true

		genomeChunkFlags = make(map[uint64]*bool)
		genomeChunkLists = make(map[uint64]*[]uint64)
		for _, idxs := range genomeChunks {
			var flag bool
			for _, idx := range idxs {
				genomeChunkFlags[idx] = &flag
				genomeChunkLists[idx] = &idxs
			}
		}
	}

	// ---------------------------------------------------------------
	// writer
	fileGenomeDetails := filepath.Join(dbDir, FileGenomeDetails)
	fhw, err := os.Create(fileGenomeDetails)
	if err != nil {
		return err
	}
	bw := bufio.NewWriter(fhw)

	// ---------------------------------------------------------------
	// genomes.map file for mapping index to genome id
	fh, err := os.Open(filepath.Join(dbDir, FileGenomeIndex))
	if err != nil {
		checkError(fmt.Errorf("failed to read genome index mapping file: %s", err))
	}
	defer fh.Close()

	r := bufio.NewReader(fh)

	buf := make([]byte, 8)
	var n, lenID int
	var batchIDAndRefID, idx uint64
	var genomeBatch, genomeIdx int
	var chunked bool
	var processed *bool
	_genomes := make([]*genome.Genome, 0, 1024)
	_idxs := make([]uint64, 0, 1024)
	var g *genome.Genome
	var seqid []byte
	var size int
	var buf1 bytes.Buffer

	var flags uint64

	if saveSeqIDs {
		flags |= FLAG_SAVE_SEQIDS
	}

	// flags
	be.PutUint64(buf, flags)
	bw.Write(buf)

	for {
		n, err = io.ReadFull(r, buf[:2])
		if err != nil {
			if err == io.EOF {
				break
			}
			checkError(fmt.Errorf("failed to read genome index mapping file: %s", err))
		}
		if n < 2 {
			checkError(fmt.Errorf("broken genome map file"))
		}
		lenID = int(be.Uint16(buf[:2]))
		id := make([]byte, lenID)

		n, err = io.ReadFull(r, id)
		if err != nil {
			checkError(fmt.Errorf("broken genome map file"))
		}
		if n < lenID {
			checkError(fmt.Errorf("broken genome map file"))
		}

		n, err = io.ReadFull(r, buf)
		if err != nil {
			checkError(fmt.Errorf("broken genome map file"))
		}
		if n < 8 {
			checkError(fmt.Errorf("broken genome map file"))
		}

		batchIDAndRefID = be.Uint64(buf)

		chunked = false
		if hasGenomeChunks {
			if processed, chunked = genomeChunkFlags[batchIDAndRefID]; chunked {
				if *processed {
					continue
				}
				*processed = true

				_genomes = _genomes[:0]
				_idxs = _idxs[:0]
				for _, idx = range *genomeChunkLists[batchIDAndRefID] {
					genomeBatch = int(idx >> BITS_GENOME_IDX)
					genomeIdx = int(idx & MASK_GENOME_IDX)

					rdr := <-poolGenomeRdrs[genomeBatch]
					g, err = rdr.GenomeInfo(int(genomeIdx))
					if err != nil {
						return fmt.Errorf("failed to read genome info: %s", err)
					}
					poolGenomeRdrs[genomeBatch] <- rdr

					_genomes = append(_genomes, g)
					_idxs = append(_idxs, idx)
				}

			}
		}

		if !chunked { // not chunked, directly read genome details
			idx = batchIDAndRefID
			_genomes = _genomes[:0]
			_idxs = _idxs[:0]
			genomeBatch = int(idx >> BITS_GENOME_IDX)
			genomeIdx = int(idx & MASK_GENOME_IDX)

			rdr := <-poolGenomeRdrs[genomeBatch]
			g, err = rdr.GenomeInfo(int(genomeIdx))
			if err != nil {
				return fmt.Errorf("failed to read genome info: %s", err)
			}
			poolGenomeRdrs[genomeBatch] <- rdr

			_genomes = append(_genomes, g)
			_idxs = append(_idxs, idx)
		}

		for i, g := range _genomes {
			if i == 0 {
				// batch+ref index
				be.PutUint64(buf, _idxs[i])
				bw.Write(buf)

				// genome id from genome data
				// be.PutUint16(buf[:2], uint16(len(g.ID)))
				// bw.Write(buf[:2])
				// bw.Write(g.ID)

				// genome id from genomes.map.bin
				be.PutUint16(buf[:2], uint16(lenID))
				// number of chunks
				be.PutUint32(buf[2:6], uint32(len(_genomes)))
				bw.Write(buf[:6])
				bw.Write(id)
			}

			// genome size
			be.PutUint32(buf[:4], uint32(g.GenomeSize))
			// number of sequences
			be.PutUint32(buf[4:8], uint32(g.NumSeqs))
			bw.Write(buf)

			// The size and id of a sequence are not saved side by side,
			// so users can optionally skip reading sequence ids.
			//
			// seq sizes
			buf1.Reset()
			for i, size = range g.SeqSizes {
				// seq size
				be.PutUint32(buf[:4], uint32(size))
				bw.Write(buf[:4])

				if saveSeqIDs {
					// seq id
					seqid = *g.SeqIDs[i]
					be.PutUint16(buf[:2], uint16(len(seqid))) // length of id
					buf1.Write(buf[:2])
					buf1.Write(seqid)
				}
			}

			if saveSeqIDs {
				// seq ids data
				be.PutUint32(buf[:4], uint32(buf1.Len())) // length of seqid data
				bw.Write(buf[:4])
				bw.Write(buf1.Bytes())
			}

			genome.RecycleGenome(g) // do not forget to recycle it.
		}
	}

	bw.Flush()
	err = fhw.Close()
	if err != nil {
		return err
	}

	// ---------------------------------------------------------------
	// close genome readers
	var _err error
	for _, pool := range poolGenomeRdrs {
		wg.Add(1)
		go func(pool chan *genome.Reader) {
			close(pool)
			for rdr := range pool {
				err := rdr.Close()
				if err != nil {
					_err = err
				}
			}
			wg.Done()
		}(pool)
	}
	wg.Wait()

	return _err
}

// ErrBrokenFile means the file is not complete.
var ErrBrokenFile = errors.New("genome chunk detail data: broken file")

func readGenomeDetails(fileGenomeDetails string, outfh *bufio.Writer, extra bool) error {
	fh, err := os.Open(fileGenomeDetails)
	checkError(err)
	br := bufio.NewReader(fh)

	buf := make([]byte, 1024)
	var n int
	// var batchIDAndRefID uint64
	var l16 uint16
	var nChunks, nSeqs uint32
	genomeID := make([]byte, 0, 1024)
	var genomeSize, seqSize uint32
	var i, j, lenID int

	genomeSizes := make([]uint32, 0, 1024)
	seqSizes := make([][]uint32, 0, 1024)
	seqIDs := make([][][]byte, 0, 1024)
	var totalGenomeSize uint64
	var buf1, buf2 bytes.Buffer

	if extra {
		fmt.Fprintf(outfh, "ref\tgenome_size\tchunks\tchunk\tchunk_size\tseqs\tseqsizes\tseqids\n")
	} else {
		fmt.Fprintf(outfh, "ref\tgenome_size\tchunks\tchunk\tchunk_size\tseqs\n")
	}

	// flags
	n, err = io.ReadFull(br, buf[:8])
	if n < 8 {
		return ErrBrokenFile
	}
	flags := be.Uint64(buf[:8])
	hasSeqIDs := (flags & FLAG_SAVE_SEQIDS) != 0

	for {
		// batch+ref index (8 bytes), the length of genome id (2 bytes), the number of chunks (4 bytes)
		n, err = io.ReadFull(br, buf[:14])
		if err != nil {
			if err == io.EOF {
				break
			}
			return ErrBrokenFile
		}
		if n < 14 {
			return ErrBrokenFile
		}
		// batchIDAndRefID = be.Uint64(buf[:8]) // batch+ref index

		l16 = be.Uint16(buf[8:10])      // length of genome id
		nChunks = be.Uint32(buf[10:14]) // number of chunks

		n, err = io.ReadFull(br, buf[:l16])
		if n < int(l16) {
			return ErrBrokenFile
		}
		genomeID = append(genomeID[:0], buf[:l16]...) // genome id

		genomeSizes = genomeSizes[:0]
		seqSizes = seqSizes[:0]
		seqIDs = seqIDs[:0]

		for i = 0; i < int(nChunks); i++ {
			// genome size (4 bytes), number of sequences (4 bytes)
			n, err = io.ReadFull(br, buf[:8])
			if n < 8 {
				return ErrBrokenFile
			}
			genomeSize = be.Uint32(buf[:4]) // genome size
			genomeSizes = append(genomeSizes, genomeSize)

			nSeqs = be.Uint32(buf[4:8]) // number of sequences

			// seq sizes
			_seqSizes := make([]uint32, nSeqs)
			for j = 0; j < int(nSeqs); j++ {
				// seq size (4 bytes)
				n, err = io.ReadFull(br, buf[:4])
				if n < 4 {
					return ErrBrokenFile
				}
				_seqSizes[j] = be.Uint32(buf[:4])
			}
			seqSizes = append(seqSizes, _seqSizes)

			if hasSeqIDs {
				// length of seqid data
				n, err = io.ReadFull(br, buf[:4])
				if n < 4 {
					return ErrBrokenFile
				}

				// seq ids data
				_seqIDs := make([][]byte, nSeqs)
				for j = 0; j < int(nSeqs); j++ {
					n, err = io.ReadFull(br, buf[:2])
					if n < 2 {
						return ErrBrokenFile
					}
					lenID = int(be.Uint16(buf[:2])) // length of seqid

					seqid := make([]byte, lenID)
					n, _ = io.ReadFull(br, seqid)
					if n < lenID {
						return ErrBrokenFile
					}
					_seqIDs[j] = seqid
				}
				seqIDs = append(seqIDs, _seqIDs)
			}
		}

		totalGenomeSize = 0
		for _, size := range genomeSizes {
			totalGenomeSize += uint64(size)
		}

		if extra {
			for i = 0; i < int(nChunks); i++ {
				buf1.Reset()
				buf2.Reset()

				for j, seqSize = range seqSizes[i] {
					buf1.WriteString(fmt.Sprintf("%d", seqSize))

					if hasSeqIDs {
						buf2.WriteString(fmt.Sprintf("%s", seqIDs[i][j]))
					}

					if j < len(seqSizes[i])-1 {
						buf1.WriteByte(',')

						if hasSeqIDs {
							buf2.WriteByte(',')
						}
					}
				}

				fmt.Fprintf(outfh, "%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n",
					genomeID, totalGenomeSize, nChunks, i+1, genomeSizes[i], len(seqSizes[i]), buf1.String(), buf2.String())

			}
		} else {
			for i = 0; i < int(nChunks); i++ {
				fmt.Fprintf(outfh, "%s\t%d\t%d\t%d\t%d\t%d\n",
					genomeID, totalGenomeSize, nChunks, i+1, genomeSizes[i], len(seqSizes[i]))

			}
		}
	}

	return nil
}
