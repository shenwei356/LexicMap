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
	"time"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
	"github.com/vbauerster/mpb/v8"
	"github.com/vbauerster/mpb/v8/decor"
)

var genomeDetailsCmd = &cobra.Command{
	Use:   "genome-details",
	Short: "Extract or view genome details in the index",
	Long: `Extract or view genome details in the index

On the first run, this command will extract genome details and save them to 'genomes.details.bin'.
If the file exists, it will be read directly and details will be printed.
The cache is published only after extraction succeeds. Extraction keeps one genome batch open at a time.
If an older cache is reported as incomplete, remove genomes.details.bin and rerun.
Renaming IDs with edit-genome-ids removes this cache. The next run rebuilds it;
use -i/--save-seqids again if sequence IDs are needed.

Output format:
  Tab-delimited format.

    ref             genome id
    genome_size     genome size (sum of all genome chunks)
    chunks          the number of genome chunks
    chunk           nth genome chunk (1-based)
    cidx            the index of the genome chunk (0-based)
    gidx            the index of the genome in the chunk (0-based)
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

		outFile := getFlagString(cmd, "out-file")
		defer setupCommandLog(opt, outFile)()
		outputLog := opt.Verbose || opt.Log2File

		// ------------------------------

		dbDir := getFlagString(cmd, "index")
		if dbDir == "" {
			checkError(fmt.Errorf("flag -d/--index needed"))
		}

		extra := getFlagBool(cmd, "extra")

		saveSeqIDs := getFlagBool(cmd, "save-seqids")

		// output file handler
		outfh, gw, w, err := outStream(outFile, strings.HasSuffix(outFile, ".gz"), opt.CompressionLevel)
		checkError(err)
		defer func() {
			checkError(closeOutStream(outfh, gw, w))
		}()

		// -----------------------------------------------------
		// Extract genome details if not existed, otherwise read it directly.

		fileGenomeDetails := filepath.Join(dbDir, FileGenomeDetails)
		existed, err := pathutil.Exists(fileGenomeDetails)
		checkError(err)
		if !existed {
			if outputLog {
				log.Infof("extracting genome details and saving to %s", fileGenomeDetails)
			}
			timeStart := time.Now()
			checkError(extractGenomeDetails(opt, dbDir, saveSeqIDs))
			if outputLog {
				log.Infof("  elapsed time: %s", time.Since(timeStart))
				log.Info()
			}
		}

		// -----------------------------------------------------
		if outputLog {
			log.Infof("reading genome details from %s", fileGenomeDetails)
		}
		timeStart1 := time.Now()
		info, err := readIndexInfo(filepath.Join(dbDir, FileInfo))
		checkError(err)
		checkError(readGenomeDetails(fileGenomeDetails, outfh, extra, info.Genomes))
		if outputLog {
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

	genomeDetailsCmd.SetUsageTemplate(usageTemplate("-d <index path> [-o out.tsv.gz]"))
}

// FileGenomeDetails store lists of genome details.
//
// Format:
//
//	flag (8 bytes)
//	len(genome id, 2 bytes), n_genome_chunks (4 bytes), genome id (X bytes)
//	# genome chunk 1:
//	    batch+ref index (8 bytes)
//	    gsize (4 bytes)
//	    n_seqs (4 bytes)
//	    # sequence size of each sequence:
//	        size 1 (4 bytes)
//	        size 2 (4 bytes)
//	        ...
//	    # sequence id of each sequence (optional, the flag is in FLAG_SAVE_SEQIDS):
//	        len_seqid1 (2 bytes), seqid1 (X bytes)
//	        len_seqid1 (2 bytes), seqid1 (X bytes)
//	        ...
const FileGenomeDetails = "genomes.details.bin"
const FLAG_SAVE_SEQIDS = 1

func extractGenomeDetails(opt *Options, dbDir string, saveSeqIDs bool) error {
	outputLog := opt.Verbose || opt.Log2File

	// ---------------------------------------------------------------
	// info file
	if outputLog {
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
	// Keep only the current batch open. Mapping records are normally grouped
	// by batch; split genomes may revisit a batch, which is reopened on demand.
	var rdr *genome.Reader
	readerBatch := -1
	defer func() {
		if rdr != nil {
			_ = rdr.Close()
		}
	}()
	readGenomeInfo := func(batch, ref int) (*genome.Genome, error) {
		if rdr == nil || readerBatch != batch {
			if rdr != nil {
				closeErr := rdr.Close()
				rdr = nil
				if closeErr != nil {
					return nil, closeErr
				}
			}
			file := filepath.Join(dbDir, DirGenomes, batchDir(batch), FileGenomes)
			var openErr error
			rdr, openErr = genome.NewReader(file)
			if openErr != nil {
				return nil, fmt.Errorf("failed to open genome batch %d: %w", batch, openErr)
			}
			readerBatch = batch
		}
		return rdr.GenomeInfo(ref)
	}

	// ---------------------------------------------------------------
	// read genome chunks data if existed
	if outputLog {
		log.Infof("  reading genome chunk data files")
	}
	genomeChunks, err := readGenomeChunksLists(filepath.Join(dbDir, FileGenomeChunks))
	if err != nil {
		return fmt.Errorf("failed to read genome chunk file: %w", err)
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
	// Publish only a fully written and closed cache, on the same filesystem.
	fhw, err := os.CreateTemp(dbDir, "."+FileGenomeDetails+"-*")
	if err != nil {
		return err
	}
	defer func() {
		_ = fhw.Close()
		_ = os.Remove(fhw.Name())
	}()
	bw := bufio.NewWriterSize(fhw, IOBufferSize)

	// ---------------------------------------------------------------
	// genomes.map file for mapping index to genome id
	if outputLog {
		log.Infof("  reading genomes.map file and genome data")
	}

	fh, err := os.Open(filepath.Join(dbDir, FileGenomeIndex))
	if err != nil {
		return fmt.Errorf("failed to read genome index mapping file: %w", err)
	}
	defer fh.Close()
	// Keep shared indexes readable to the same users as the genome mapping.
	stat, err := fh.Stat()
	if err != nil {
		return err
	}
	if err = fhw.Chmod(stat.Mode().Perm()); err != nil {
		return err
	}

	r := bufio.NewReaderSize(fh, IOBufferSize)

	buf := make([]byte, 8)
	var n, lenID int
	var batchIDAndRefID, idx uint64
	var genomeBatch, genomeIdx int
	var chunked bool
	var processed *bool
	_genomes := make([]*genome.Genome, 0, 1024)
	defer func() {
		for _, g := range _genomes {
			if g != nil {
				genome.RecycleGenome(g)
			}
		}
	}()
	_idxs := make([]uint64, 0, 1024)
	var g *genome.Genome
	var seqid []byte
	var buf1 bytes.Buffer

	var flags uint64

	if saveSeqIDs {
		flags |= FLAG_SAVE_SEQIDS
	}

	// flags
	be.PutUint64(buf, flags)
	bw.Write(buf)

	// -----------------------------------------------------------
	// process bar
	var pbs *mpb.Progress
	var bar *mpb.Bar
	var chDuration chan time.Duration
	var doneDuration chan int
	if opt.Verbose {
		pbs = mpb.New(mpb.WithWidth(40), mpb.WithOutput(os.Stderr))
		bar = pbs.AddBar(int64(info.InputGenomes),
			mpb.PrependDecorators(
				decor.Name("processed genomes: ", decor.WC{W: len("processed genomes: "), C: decor.DindentRight}),
				decor.Name("", decor.WCSyncSpaceR),
				decor.CountersNoUnit("%d / %d", decor.WCSyncWidth),
			),
			mpb.AppendDecorators(progressTimeDecorators(1024)...),
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
	// Stop the progress worker on read/write failures as well as success.
	if opt.Verbose {
		defer func() {
			close(chDuration)
			<-doneDuration
			bar.Abort(true)
			pbs.Wait()
		}()
	}
	// -----------------------------------------------------------

	var timeStart time.Time
	for {
		timeStart = time.Now()

		n, err = io.ReadFull(r, buf[:2])
		if err != nil {
			if err == io.EOF {
				break
			}
			return fmt.Errorf("failed to read genome index mapping file: %w", err)
		}
		if n < 2 {
			return fmt.Errorf("broken genome map file")
		}
		lenID = int(be.Uint16(buf[:2]))
		id := make([]byte, lenID)

		n, _ = io.ReadFull(r, id)
		if n < lenID {
			return fmt.Errorf("broken genome map file")
		}

		n, _ = io.ReadFull(r, buf)
		if n < 8 {
			return fmt.Errorf("broken genome map file")
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

					g, err = readGenomeInfo(genomeBatch, genomeIdx)
					if err != nil {
						return fmt.Errorf("failed to read genome info: %s", err)
					}

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

			g, err = readGenomeInfo(genomeBatch, genomeIdx)
			if err != nil {
				return fmt.Errorf("failed to read genome info: %s", err)
			}

			_genomes = append(_genomes, g)
			_idxs = append(_idxs, idx)
		}

		for i, g := range _genomes {
			if i == 0 {
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

			// batch+ref index of each chunk
			be.PutUint64(buf, _idxs[i])
			bw.Write(buf)

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
			for i, size := range g.SeqSizes {
				// seq size
				be.PutUint32(buf[:4], uint32(size))
				bw.Write(buf[:4])

				if saveSeqIDs {
					// seq id
					seqid = *g.SeqIDs[i]
					if len(seqid) > 65535 { // truncate super-long sequence ID
						seqid = seqid[:65535]
					}
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
			_genomes[i] = nil
		}
		clear(_genomes)

		// Only progress-bar mode creates and drains this channel.
		if opt.Verbose {
			chDuration <- time.Since(timeStart)
		}
	}

	if err = bw.Flush(); err != nil {
		_ = fhw.Close()
		return err
	}
	err = fhw.Close()
	if err != nil {
		return err
	}

	if rdr != nil {
		err = rdr.Close()
		rdr = nil
		if err != nil {
			return err
		}
	}
	if err = fh.Close(); err != nil {
		return err
	}
	return os.Rename(fhw.Name(), fileGenomeDetails)
}

// ErrBrokenFile means the file is not complete.
var ErrBrokenFile = errors.New("genome chunk detail data: broken file")

// readGenomeDetails checks the chunk count so legacy caches truncated at a
// record boundary cannot silently appear complete. The binary format is unchanged.
func readGenomeDetails(fileGenomeDetails string, outfh *bufio.Writer, extra bool, expectedChunks int) error {
	fh, err := os.Open(fileGenomeDetails)
	if err != nil {
		return err
	}
	defer fh.Close()

	br := bufio.NewReaderSize(fh, IOBufferSize)

	buf := make([]byte, 1<<16) // 64K
	var n int
	var batchIDAndRefID uint64
	var l16 uint16
	var nChunks, nSeqs uint32
	var chunksRead uint64
	genomeID := make([]byte, 0, 1024)
	var genomeSize, seqSize uint32
	var i, j, lenID int

	batchIDAndRefIDs := make([]uint64, 0, 1024)
	genomeSizes := make([]uint32, 0, 1024)
	seqSizes := make([][]uint32, 0, 1024)
	seqIDs := make([][][]byte, 0, 1024)
	var totalGenomeSize uint64
	var buf1, buf2 bytes.Buffer

	if extra {
		fmt.Fprintf(outfh, "ref\tgenome_size\tchunks\tchunk\tcidx\tgidx\tchunk_size\tseqs\tseqsizes\tseqids\n")
	} else {
		fmt.Fprintf(outfh, "ref\tgenome_size\tchunks\tchunk\tcidx\tgidx\tchunk_size\tseqs\n")
	}

	// flags
	n, err = io.ReadFull(br, buf[:8])
	if n < 8 {
		return ErrBrokenFile
	}
	flags := be.Uint64(buf[:8])
	hasSeqIDs := (flags & FLAG_SAVE_SEQIDS) != 0

	for {
		// the length of genome id (2 bytes), the number of chunks (4 bytes)
		n, err = io.ReadFull(br, buf[:6])
		if err != nil {
			if err == io.EOF {
				break
			}
			return ErrBrokenFile
		}
		if n < 6 {
			return ErrBrokenFile
		}

		l16 = be.Uint16(buf[:2])      // length of genome id
		nChunks = be.Uint32(buf[2:6]) // number of chunks
		chunksRead += uint64(nChunks)
		if nChunks == 0 || chunksRead > uint64(expectedChunks) {
			return ErrBrokenFile
		}

		n, err = io.ReadFull(br, buf[:l16])
		if n < int(l16) {
			return ErrBrokenFile
		}
		genomeID = append(genomeID[:0], buf[:l16]...) // genome id

		batchIDAndRefIDs = batchIDAndRefIDs[:0]
		genomeSizes = genomeSizes[:0]
		clear(seqSizes)
		seqSizes = seqSizes[:0]
		clear(seqIDs)
		seqIDs = seqIDs[:0]

		for i = 0; i < int(nChunks); i++ {
			// batch+ref index (8 bytes), genome size (4 bytes), number of sequences (4 bytes)
			n, err = io.ReadFull(br, buf[:16])
			if n < 16 {
				return ErrBrokenFile
			}

			batchIDAndRefID = be.Uint64(buf[:8]) // batch+ref index
			batchIDAndRefIDs = append(batchIDAndRefIDs, batchIDAndRefID)
			genomeSize = be.Uint32(buf[8:12])             // genome size
			genomeSizes = append(genomeSizes, genomeSize) // genome sizes

			nSeqs = be.Uint32(buf[12:16]) // number of sequences

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

				fmt.Fprintf(outfh, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n",
					genomeID, totalGenomeSize, nChunks, i+1,
					batchIDAndRefIDs[i]>>BITS_GENOME_IDX,
					batchIDAndRefIDs[i]&MASK_GENOME_IDX,
					genomeSizes[i], len(seqSizes[i]), buf1.String(), buf2.String())

			}
		} else {
			for i = 0; i < int(nChunks); i++ {
				fmt.Fprintf(outfh, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
					genomeID, totalGenomeSize, nChunks, i+1,
					batchIDAndRefIDs[i]>>BITS_GENOME_IDX,
					batchIDAndRefIDs[i]&MASK_GENOME_IDX,
					genomeSizes[i], len(seqSizes[i]))

			}
		}
	}

	if chunksRead != uint64(expectedChunks) {
		return fmt.Errorf("%w: read %d genome chunks, expected %d; remove %s and rerun to rebuild the cache", ErrBrokenFile, chunksRead, expectedChunks, fileGenomeDetails)
	}
	return nil
}
