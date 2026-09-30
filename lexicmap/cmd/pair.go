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
	"fmt"
	"io"
	"math"
	"math/bits"
	"os"
	"path/filepath"
	"runtime"
	"slices"
	"strings"
	"sync"
	"time"

	"github.com/dustin/go-humanize"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/lexichash"
	"github.com/spf13/cobra"
	"github.com/twotwotwo/sorts"
	"github.com/vbauerster/mpb/v8"
	"github.com/vbauerster/mpb/v8/decor"
)

var pairCmd = &cobra.Command{
	Use:   "pair",
	Short: "Find similar genome pairs in the index",
	Long: `Find similar genome pairs in the index

Output format:
  Tab-delimited format with 7 columns.

    1.  genome1,    Genome 1.
    2.  genome2,    Genome 2.
    3.  minPrefix,  Minimum common prefix length between two seeds (-p/--min-prefix).
    4.  fracMasks,  Fraction of masks with seeds sharing a common prefix >= minPrefix.
    5.  nMasks,     Number   of masks with seeds sharing a common prefix >= minPrefix.
    6.  sumPrefix,  Total length of common prefixes, with only the best seed pair
                    (with the longest common prefix) considered for each mask.
    7.  avgPrefix,  Average prefix length (sumPrefix / nMasks).

Limitations:
  1. Not suitable for very large genome sets (e.g., > 1 million genomes) due to memory constraints.
  2. Genome sets with many highly similar genomes may result in a large number of genome pairs, 
     which can require substantial time and memory.
  3. Genomes stored in multiple chunks are not evaluated as a whole.

Notes on pair-table memory:
  1. For indexes with no more than --max-dense-genomes genomes, compact dense pair
     tables are used. They are faster and use less memory when many genome pairs
     match, but their size grows quadratically with the number of genomes. Set the
     value to 0 to force the sparse-map implementation.
  2. Raising --max-dense-genomes enables dense tables for a larger index and can
     greatly increase baseline memory. Lowering it can save memory for indexes in
     which relatively few genome pairs match. At the default limit of 20,000
     genomes, the global tables need about 1.14 GiB and each mask worker needs
     about 215 MiB; the default 768 MiB budget therefore permits up to three mask
     workers. These values exclude k-mer windows and final sorting results.
  3. --dense-mask-memory limits only the temporary tables used by concurrent mask
     workers; it is not a total RSS limit. One table is allocated up front for each
     permitted worker. Raising the value may increase parallelism, but speed usually
     has diminishing returns and can even regress because workers compete for memory
     bandwidth and CPU cache, while completed masks are merged into the global dense
     table by one collector goroutine. The global pair statistics, k-mer windows, and
     final sorting results are additional memory.

Notes on probabilistic pruning:
  1. -f/--min-mask-fraction is the final reporting threshold. In contrast,
     -s/--prob-threshold controls an early heuristic that discards pairs unlikely
     to reach that threshold. Higher values prune more aggressively; 0 disables
     the heuristic and computes exact pair statistics.
  2. In sparse-map mode, pruning can substantially reduce the number of map entries
     and therefore memory use. In dense-table mode, all pair slots are allocated
     up front, so pruning does not reduce the main table memory or skip pair
     comparisons within a mask. It can still save some time by avoiding updates
     to discarded pairs and by keeping the active result set smaller.
  3. The default -s 0.001 is retained for compatibility and sparse-map workloads.
     With -s 0, results are independent of mask completion order, so changing
     --dense-mask-memory affects only memory and parallelism. With -s > 0, changing
     concurrency can alter the completion order and may affect marginal pairs
     considered by the pruning heuristic.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		seq.ValidateSeq = false

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

		var err error

		// -------------------------------------------------------------------------

		dbDir := getFlagString(cmd, "index")
		if dbDir == "" {
			checkError(fmt.Errorf("flag -d/--index needed"))
		}
		outFile := getFlagString(cmd, "out-file")
		minPrefix := getFlagPositiveInt(cmd, "min-prefix")
		minMaskFraction := getFlagNonNegativeFloat64(cmd, "min-mask-fraction")
		probThreshold := getFlagNonNegativeFloat64(cmd, "prob-threshold")
		maxDenseGenomes := getFlagNonNegativeInt(cmd, "max-dense-genomes")
		if maxDenseGenomes > maxSupportedDenseGenomes {
			checkError(fmt.Errorf("value of --max-dense-genomes (%d) should not be greater than %d",
				maxDenseGenomes, maxSupportedDenseGenomes))
		}
		denseMaskMemoryMiB := getFlagPositiveInt(cmd, "dense-mask-memory")

		nMasks := getFlagNonNegativeInt(cmd, "masks")
		if !(nMasks == 0 || (isPowerOf4(nMasks) && nMasks >= 64)) {
			checkError(fmt.Errorf("the value of -m/--masks should be 0 (for all masks in the index) or power of 4 (needs to be >= 64, e.g., 64, 256, 1024, 4096, 16384)"))
		}

		// -------------------------------------------------------------------------

		if outputLog {
			log.Infof("LexicMap v%s", VERSION)
			log.Info("  https://github.com/shenwei356/LexicMap")
			log.Info()
		}

		// -------------------------------------------------------------------------
		// checking index

		if outputLog {
			log.Infof("checking index: %s", dbDir)
		}

		// Mask file
		fileMask := filepath.Join(dbDir, FileMasks)
		lh, err := lexichash.NewFromFile(fileMask)
		if err != nil {
			checkError(err)
		}

		if nMasks > len(lh.Masks) {
			checkError(fmt.Errorf("the value of -m/--mask (%d) is bigger than the number of masks in the index (%d)", nMasks, len(lh.Masks)))
		}

		// info file
		fileInfo := filepath.Join(dbDir, FileInfo)
		info, err := readIndexInfo(fileInfo)
		if err != nil {
			checkError(fmt.Errorf("failed to read info file: %s", err))
		}

		if outputLog {
			log.Infof("  checking passed")
			log.Infof("reading seed data of all masks...")
		}

		id2name, err := readGenomeMapIdx2Name(filepath.Join(dbDir, FileGenomeIndex))
		if err != nil {
			checkError(fmt.Errorf("failed to read %s: %s", filepath.Join(dbDir, FileGenomeIndex), err))
		}
		genomeCodes, genomeOrdinals, err := buildGenomeOrdinals(id2name, info.GenomeBatches)
		if err != nil {
			checkError(err)
		}

		// -------------------------------------------------------------------------
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

		// -------------------------------------------------------------------------
		// choose masks
		var maskPrefix int
		selectedMasks := make([]bool, len(lh.Masks))
		totalMasks := 0
		if nMasks == 0 {
			for i := range selectedMasks {
				selectedMasks[i] = true
			}
			totalMasks = len(selectedMasks)
		} else {
			// maskPrefix = int(math.Log2(float64(nMasks)) / 2)
			maskPrefix = bits.TrailingZeros(uint(nMasks)) / 2
			seenPrefixes := make(map[uint64]struct{}, nMasks)
			shift := uint64(lh.K-maskPrefix) << 1
			for i, mask := range lh.Masks {
				prefix := mask >> shift
				if _, ok := seenPrefixes[prefix]; !ok {
					selectedMasks[i] = true
					totalMasks++
					seenPrefixes[prefix] = struct{}{}
				}
			}
		}

		// -------------------------------------------------------------------------
		// process bar
		var pbs *mpb.Progress
		var bar *mpb.Bar
		var chDuration chan time.Duration
		var doneDuration chan int
		var showProgressBar bool

		if opt.Verbose {
			showProgressBar = true

			pbs = mpb.New(mpb.WithWidth(40), mpb.WithOutput(os.Stderr))
			bar = pbs.AddBar(int64(totalMasks),
				mpb.PrependDecorators(
					decor.Name("processed masks: ", decor.WC{W: len("processed masks: "), C: decor.DindentRight}),
					decor.Name("", decor.WCSyncSpaceR),
					decor.CountersNoUnit("%d / %d", decor.WCSyncWidth),
				),
				mpb.AppendDecorators(
					decor.Name("ETA: ", decor.WC{W: len("ETA: ")}),
					decor.EwmaETA(decor.ET_STYLE_GO, 64),
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

		fcpus := float64(opt.NumCPUs)

		// -------------------------------------------------------------------------

		nGenomes := uint32(len(genomeCodes))
		var nGenomePairs uint64
		if nGenomes > 1 {
			nGenomePairs = uint64(nGenomes) * uint64(nGenomes-1) / 2
		}
		// Dense tables need one slot for every possible unordered genome pair.
		// Keep the switch genome-based for a more intuitive CLI, but calculate
		// allocations from the exact n*(n-1)/2 pair count.
		useDensePairs := maxDenseGenomes > 0 && nGenomes <= uint32(maxDenseGenomes) &&
			nGenomePairs <= math.MaxUint32 && totalMasks <= math.MaxUint16
		var pairStats map[uint64]PairStats
		var denseStats *DensePairStats
		var densePairStarts []uint32
		if useDensePairs {
			// Remap the possibly sparse batch/ref genome codes to [0, nGenomes).
			// Sorted codes preserve the old pair-key tie-breaking order.
			denseStats = newDensePairStats(int(nGenomePairs))
			densePairStarts = make([]uint32, nGenomes)
			for g := range nGenomes {
				densePairStarts[g] = uint32(densePairRowStart(g, nGenomes))
			}
		} else {
			// Keep match count and prefix sum in one map to avoid storing and
			// hashing every pair key twice.
			pairStats = make(map[uint64]PairStats, 10240)
		}

		// Calculate threshold for minimum prefix length
		// threshold = 1 << ((k - minPrefix) * 2)
		k := int(info.K)
		kMinus32 := k - 32 // precompute to avoid repeated calculation
		threshold := uint64(1) << ((k - minPrefix) * 2)
		minPrefixU8 := uint8(minPrefix) // convert to uint8 for comparison

		requiredMatches := int(minMaskFraction * float64(totalMasks))

		if outputLog {
			log.Infof("  minimum prefix length between k-mers captured by a mask: %d", minPrefix)
			log.Infof("  total masks: %d, required matches: %d (%.1f%%)", totalMasks, requiredMatches, minMaskFraction*100)
			if useDensePairs {
				log.Infof("  using compact pair tables for %s possible genome pairs", humanize.Comma(int64(nGenomePairs)))
			}
		}

		// -------------------------------------------------------------------------
		// collect counting results
		maskWorkers := opt.NumCPUs
		var denseMaskPool chan *DenseMaskCounts
		if useDensePairs {
			// Each concurrent mask needs one byte per possible pair plus a touched
			// bitset. --dense-mask-memory limits the sum of only these temporary
			// tables; the global statistics and k-mer windows are additional.
			bytesPerMask := nGenomePairs + ((nGenomePairs + 63) >> 6 << 3)
			denseMaskMemoryBudget := uint64(denseMaskMemoryMiB) << 20
			workersByMemory := int(denseMaskMemoryBudget / max(uint64(1), bytesPerMask))
			if workersByMemory < 1 {
				workersByMemory = 1
			}
			if maskWorkers > workersByMemory {
				maskWorkers = workersByMemory
			}
			denseMaskPool = make(chan *DenseMaskCounts, maskWorkers)
			for range maskWorkers {
				denseMaskPool <- newDenseMaskCounts(int(nGenomePairs))
			}
			if outputLog {
				log.Infof("  dense mask workers: %d (--dense-mask-memory: %d MiB)", maskWorkers, denseMaskMemoryMiB)
			}
		}

		type Result struct {
			Counts    *map[uint64]uint8
			Dense     *DenseMaskCounts
			StartTime time.Time
		}

		ch := make(chan Result, maskWorkers)
		done := make(chan int)
		// A single collector owns the global pair statistics. Mask workers write
		// only to their private tables, so the hot dense-table updates need no locks.
		go func() {
			processedMasks := 0
			remaining := totalMasks
			var pruneDecisions []int8

			for result := range ch {
				processedMasks++
				remaining--

				if result.Counts == nil && result.Dense == nil { // no k-mers
					if showProgressBar {
						chDuration <- time.Duration(float64(time.Since(result.StartTime)) / fcpus)
					}
					continue
				}

				if useDensePairs {
					// Merge one mask at a time so matches and prefix sums have the same
					// semantics as the sparse map path.
					mergeDenseMaskCounts(denseStats, result.Dense, processedMasks, remaining,
						requiredMatches, totalMasks, minMaskFraction, probThreshold)
					denseMaskPool <- result.Dense
					if probThreshold > 0 && processedMasks < totalMasks && processedMasks&7 == 0 {
						// Dense pruning clears statistics and active bits for unlikely
						// pairs, but it cannot release the fixed-size dense arrays.
						pruneDecisions = pruneDensePairStats(denseStats, pruneDecisions, processedMasks,
							minMaskFraction, totalMasks, probThreshold)
					}
				} else {
					maskCounts := result.Counts
					if len(*maskCounts) > 0 {
						if probThreshold == 0 { // no pruning
							// Simply accumulate all pairs.
							for pair, prefixLen := range *maskCounts {
								stats := pairStats[pair]
								stats.matches++
								stats.sumPrefix += uint32(prefixLen)
								pairStats[pair] = stats
							}
						} else {
							// Check if a new pair can still reach the required number of
							// matching masks. The probability check is identical for all
							// pairs first observed in this mask, so compute it once.
							shouldAddNewPair := false
							if 1+remaining >= requiredMatches {
								shouldAddNewPair = shouldKeepPair(processedMasks, 1, minMaskFraction, totalMasks, probThreshold)
							}
							// Update match counts for pairs that matched in this mask.
							for pair, prefixLen := range *maskCounts {
								stats, ok := pairStats[pair]
								if !ok {
									// New pair: retain it only if it passes the shared
									// probability check above.
									if shouldAddNewPair {
										pairStats[pair] = PairStats{matches: 1, sumPrefix: uint32(prefixLen)}
									}
								} else {
									// Existing pair: accumulate this mask's best prefix.
									stats.matches++
									stats.sumPrefix += uint32(prefixLen)
									pairStats[pair] = stats
								}
							}
							// Probabilistic pruning: every eight processed masks, remove
							// active pairs that are unlikely to reach the final threshold.
							if processedMasks < totalMasks && processedMasks&7 == 0 {
								if cap(pruneDecisions) <= processedMasks {
									pruneDecisions = make([]int8, processedMasks+1)
								} else {
									pruneDecisions = pruneDecisions[:processedMasks+1]
									clear(pruneDecisions)
								}
								for pair, stats := range pairStats {
									if stats.matches <= 1 {
										continue
									}
									decision := pruneDecisions[stats.matches]
									if decision == 0 {
										if shouldKeepPair(processedMasks, int(stats.matches), minMaskFraction, totalMasks, probThreshold) {
											decision = 1
										} else {
											decision = -1
										}
										pruneDecisions[stats.matches] = decision
									}
									if decision < 0 {
										delete(pairStats, pair)
									}
								}
							}
						}
					}
					clear(*maskCounts)
					poolMaskCounts.Put(maskCounts)
				}

				if showProgressBar {
					chDuration <- time.Duration(float64(time.Since(result.StartTime)) / fcpus)
				}
				if processedMasks&63 == 0 {
					runtime.GC()
				}
			}

			done <- 1
		}()

		// -------------------------------------------------------------------------
		// read seed data files

		var wg sync.WaitGroup
		tokens := make(chan int, maskWorkers)

		for chunk := range info.Chunks {
			wg.Add(1)
			tokens <- 1

			go func(chunk int) {
				defer func() {
					wg.Done()
					<-tokens
				}()

				fileSeeds := filepath.Join(dbDir, DirSeeds, chunkFile(chunk))

				// -------------------------------
				// header

				buf8 := make([]uint8, 8)
				var config1 uint8
				var use3BytesForSeedPos bool
				var bytesPos int
				var fUint64 func([]byte) uint64

				// the header of kv-data file
				fh, err := os.Open(fileSeeds)
				if err != nil {
					checkError(err)
				}
				defer fh.Close()

				r := bufio.NewReaderSize(fh, 64<<10)

				var n int

				// check the magic number
				n, err = io.ReadFull(r, buf8)
				if n < 8 {
					checkError(ErrBrokenFile)
				}
				same := true
				for i := 0; i < 8; i++ {
					if kv.Magic[i] != buf8[i] {
						same = false
						break
					}
				}
				if !same {
					checkError(kv.ErrInvalidFileFormat)
				}

				// read version information
				n, err = io.ReadFull(r, buf8)
				if n < 8 {
					checkError(ErrBrokenFile)
				}
				// check compatibility
				if kv.MainVersion != buf8[0] {
					checkError(kv.ErrVersionMismatch)
				}

				config1 = buf8[3]

				// index of the first mask in current chunk.
				n, err = io.ReadFull(r, buf8)
				if n < 8 {
					checkError(ErrBrokenFile)
				}
				iFirstMask := int(be.Uint64(buf8))

				// mask chunk size
				n, err = io.ReadFull(r, buf8)
				if n < 8 {
					checkError(ErrBrokenFile)
				}
				nMasks := int(be.Uint64(buf8))

				use3BytesForSeedPos = config1&kv.MaskUse3BytesForSeedPos > 0
				if !use3BytesForSeedPos {
					checkError(fmt.Errorf("index with genome batch number > 512 is not supported"))
				}
				bytesPos = 8
				fUint64 = be.Uint64
				if use3BytesForSeedPos {
					bytesPos = 7
					fUint64 = kv.Uint64ThreeBytes
				}

				// kv-data index file
				indexes, err := kv.ReadKVIndexStarts(filepath.Clean(fileSeeds) + kv.KVIndexFileExt)
				if err != nil {
					checkError(fmt.Errorf("failed to read kv-data index file: %s", err))
				}

				// -------------------------------
				// data of all masks

				buf := make([]byte, 64)
				valueBuf := make([]byte, 0, 4096)
				var ctrlByte byte
				var first bool     // the first kmer has a different way to compute the value
				var lastPair bool  // check if this is the last pair
				var hasKmer2 bool  // check if there's a kmer2
				var _offset uint64 // offset of kmer
				var nBytes int
				var nReaded, nDecoded int
				var v1, v2 uint64
				var kmer1, kmer2 uint64
				var lenVal1, lenVal2 uint64
				var j uint64
				var v, batchIDAndRefID uint64
				var i, valueOffset, nValueBytes int
				var lastGenome uint32
				var hasLastGenome bool

				for iMask := 0; iMask < nMasks; iMask++ {
					maskIndex := iFirstMask + iMask
					if !selectedMasks[maskIndex] {
						continue
					}

					var maskStart time.Time
					if showProgressBar {
						maskStart = time.Now()
					}

					if indexes[iMask][1] == 0 { // no k-mers
						ch <- Result{
							Counts:    nil,
							StartTime: maskStart,
						}
						continue
					}

					// genome id list
					genomes := poolGenomes.Get().(*[]uint32)

					// Sliding window for all-to-all comparison
					window := poolKmerWindow.Get().(*KmerWindow)

					// Per-mask tracking: keep the maximum prefix for each pair.
					var maskCounts *map[uint64]uint8
					var denseCounts *DenseMaskCounts
					if useDensePairs {
						denseCounts = <-denseMaskPool
					} else {
						maskCounts = poolMaskCounts.Get().(*map[uint64]uint8)
					}

					// seek
					_, err = fh.Seek(int64(indexes[iMask][1])>>1, 0)
					if err != nil {
						checkError(fmt.Errorf("failed to seek kv-data file: %s", err))
					}

					r.Reset(fh) // use buffer

					// -------------------------------
					// read data of a mask

					_offset = 0
					first = true
					for {
						// read the control byte
						_, err = io.ReadFull(r, buf[:1])
						if err != nil {
							checkError(err)
						}
						ctrlByte = buf[0]

						lastPair = ctrlByte&128 > 0 // 1<<7
						hasKmer2 = ctrlByte&64 == 0 // 1<<6

						ctrlByte &= 63

						// parse the control byte
						nBytes = util.CtrlByte2ByteLengthsUint64(ctrlByte)

						// read encoded bytes
						nReaded, err = io.ReadFull(r, buf[:nBytes])
						if nReaded < nBytes {
							checkError(kv.ErrBrokenFile)
						}

						v1, v2, nDecoded = util.Uint64s(ctrlByte, buf[:nBytes])
						if nDecoded == 0 {
							checkError(kv.ErrBrokenFile)
						}

						if first {
							kmer1 = indexes[iMask][0] // from the index
							first = false
						} else {
							kmer1 = v1 + _offset
						}
						kmer2 = kmer1 + v2
						_offset = kmer2

						// ------------------ lengths of values -------------------

						// read the control byte
						_, err = io.ReadFull(r, buf[:1])
						if err != nil {
							checkError(err)
						}
						ctrlByte = buf[0]

						// parse the control byte
						nBytes = util.CtrlByte2ByteLengthsUint64(ctrlByte)

						// read encoded bytes
						nReaded, err = io.ReadFull(r, buf[:nBytes])
						if nReaded < nBytes {
							checkError(kv.ErrBrokenFile)
						}

						lenVal1, lenVal2, nDecoded = util.Uint64s(ctrlByte, buf[:nBytes])
						if nDecoded == 0 {
							checkError(kv.ErrBrokenFile)
						}

						// Values of the two k-mers are contiguous in the file. Read
						// them in one operation instead of one io.ReadFull call per
						// encoded value.
						nValues := lenVal1
						if hasKmer2 {
							nValues += lenVal2
						}
						nValueBytes = int(nValues) * bytesPos
						if cap(valueBuf) < nValueBytes {
							valueBuf = make([]byte, nValueBytes)
						} else {
							valueBuf = valueBuf[:nValueBytes]
						}
						if nValueBytes > 0 {
							nReaded, err = io.ReadFull(r, valueBuf)
							if nReaded < nValueBytes || err != nil {
								checkError(kv.ErrBrokenFile)
							}
						}
						valueOffset = 0

						// ------------------ values for kmer1 -------------------

						*genomes = (*genomes)[:0] // reuse slice
						hasLastGenome = false
						for j = 0; j < lenVal1; j++ {
							v = fUint64(valueBuf[valueOffset : valueOffset+bytesPos])
							valueOffset += bytesPos
							if v&MASK_REVERSE == 1 {
								continue // skip reverse complement
							}
							// Extract genome ID (batchID + refID)
							batchIDAndRefID = (v >> BITS_NONE_IDX) & 4294967295
							genomeCode := uint32(batchIDAndRefID)
							if hasLastGenome && genomeCode == lastGenome {
								continue
							}
							genome := genomeCode
							if useDensePairs {
								genome = genomeOrdinal(genomeCode, genomeOrdinals)
							}
							*genomes = append(*genomes, genome)
							lastGenome = genomeCode
							hasLastGenome = true
						}

						// Process kmer1 with sliding window
						if len(*genomes) > 0 {
							processKmerWithWindow(kmer1, genomes, window, maskCounts, denseCounts, densePairStarts,
								threshold, kMinus32, minPrefixU8)
						}

						if lastPair && !hasKmer2 {
							break
						}

						// ------------------ values for kmer2 -------------------

						*genomes = (*genomes)[:0] // reuse slice
						hasLastGenome = false
						for j = 0; j < lenVal2; j++ {
							v = fUint64(valueBuf[valueOffset : valueOffset+bytesPos])
							valueOffset += bytesPos
							if v&MASK_REVERSE == 1 {
								continue // skip reverse complement
							}
							batchIDAndRefID = (v >> BITS_NONE_IDX) & 4294967295
							genomeCode := uint32(batchIDAndRefID)
							if hasLastGenome && genomeCode == lastGenome {
								continue
							}
							genome := genomeCode
							if useDensePairs {
								genome = genomeOrdinal(genomeCode, genomeOrdinals)
							}
							*genomes = append(*genomes, genome)
							lastGenome = genomeCode
							hasLastGenome = true
						}

						// Process kmer2 with sliding window
						if len(*genomes) > 0 {
							processKmerWithWindow(kmer2, genomes, window, maskCounts, denseCounts, densePairStarts,
								threshold, kMinus32, minPrefixU8)
						}

						if lastPair {
							break
						}
					}

					// recycle objects and send result
					poolGenomes.Put(genomes)

					for i = window.head; i < len(window.records); i++ {
						window.records[i].genomes = window.records[i].genomes[:0]
						poolKmerRecord.Put(window.records[i])
					}
					clear(window.records)
					window.records = window.records[:0]
					window.head = 0
					poolKmerWindow.Put(window)

					ch <- Result{
						Counts:    maskCounts,
						Dense:     denseCounts,
						StartTime: maskStart,
					}
				}

			}(chunk)
		}

		wg.Wait()
		close(ch)
		<-done

		if showProgressBar {
			close(chDuration)
			<-doneDuration
			pbs.Wait()
		}

		// ---------------------------------------------------------------
		// Output results
		if outputLog {
			log.Info()
			log.Info("sorting and writing results...")
		}

		// Write header
		outfh.WriteString("genome1\tgenome2\tminPrefix\tfracMasks\tnMasks\tsumPrefix\tavgPrefix\n")

		var nResults int
		if useDensePairs {
			results := make([]uint32, 0, countDenseResults(denseStats, requiredMatches))
			for wordIndex, word := range denseStats.active {
				for word != 0 {
					bit := bits.TrailingZeros64(word)
					pairIndex := wordIndex<<6 + bit
					if pairIndex < len(denseStats.matches) && int(denseStats.matches[pairIndex]) >= requiredMatches {
						results = append(results, uint32(pairIndex))
					}
					word &= word - 1
				}
			}
			nResults = len(results)
			sorts.Quicksort(DensePairResults{pairs: results, stats: denseStats})

			for _, pairIndex := range results {
				g1, g2 := densePairGenomes(pairIndex, nGenomes)
				statsIndex := int(pairIndex)
				nMatchedMasks := int(denseStats.matches[statsIndex])
				sumPrefix := denseStats.sumPrefixes[statsIndex]
				fracMasks := float64(nMatchedMasks) / float64(totalMasks)
				fmt.Fprintf(outfh, "%s\t%s\t%d\t%.4f\t%d\t%d\t%.2f\n",
					id2name[uint64(genomeCodes[g1])], id2name[uint64(genomeCodes[g2])], minPrefix,
					fracMasks, nMatchedMasks, sumPrefix, float64(sumPrefix)/float64(nMatchedMasks))
			}
		} else {
			results := make([]PairResult, 0, len(pairStats))
			for pair, stats := range pairStats {
				if int(stats.matches) >= requiredMatches {
					results = append(results, PairResult{
						pair:      pair,
						nMasks:    int(stats.matches),
						sumPrefix: stats.sumPrefix,
					})
				}
			}
			nResults = len(results)
			sorts.Quicksort(PairResults(results))

			for _, result := range results {
				gid1 := result.pair >> 32
				gid2 := result.pair & 0xFFFFFFFF
				fracMasks := float64(result.nMasks) / float64(totalMasks)
				fmt.Fprintf(outfh, "%s\t%s\t%d\t%.4f\t%d\t%d\t%.2f\n",
					id2name[gid1], id2name[gid2], minPrefix, fracMasks, result.nMasks,
					result.sumPrefix, float64(result.sumPrefix)/float64(result.nMasks))
			}
		}
		if outputLog {
			log.Info()
			log.Infof("total genome pairs: %d", nResults)
		}

		if outputLog && outFile != "-" {
			log.Infof("results saved to: %s", outFile)
		}

	},
}

func init() {
	genomeCmd.AddCommand(pairCmd)

	pairCmd.Flags().StringP("index", "d", "",
		formatFlagUsage(`Index directory created by "lexicmap index".`))

	pairCmd.Flags().IntP("masks", "m", 1024,
		formatFlagUsage(`Only use seed data of N masks. It should be 0 (for all masks in the index) or power of 4 (needs to be >= 64, e.g., 64, 256, 1024, 4096, 16384).`))

	pairCmd.Flags().StringP("out-file", "o", "-",
		formatFlagUsage(`Out file, supports and recommends a ".gz" suffix ("-" for stdout).`))

	pairCmd.SetUsageTemplate(usageTemplate("-d <index path> [-o out.tsv.gz]"))

	pairCmd.Flags().IntP("min-prefix", "p", 21,
		formatFlagUsage(`Minimum prefix length between k-mers captured by a mask.`))

	pairCmd.Flags().Float64P("min-mask-fraction", "f", 0.25,
		formatFlagUsage(`Minimum fraction of masks that must match for a genome pair to be reported.`))

	pairCmd.Flags().Float64P("prob-threshold", "s", 0.001,
		formatFlagUsage(`Probability threshold for early pruning (higher = more aggressive; 0 disables pruning and computes exact pair statistics).`))

	pairCmd.Flags().Int("max-dense-genomes", defaultMaxDenseGenomes,
		formatFlagUsage(fmt.Sprintf(`Maximum number of genomes for using compact dense pair tables (maximum: %d; 0 disables them). Bigger values can greatly increase memory usage.`, maxSupportedDenseGenomes)))

	pairCmd.Flags().Int("dense-mask-memory", defaultDenseMaskMemoryMiB,
		formatFlagUsage(`Memory budget in MiB for concurrent temporary dense-mask tables (not total RSS). Bigger values permit more mask workers but may have diminishing returns because of memory-bandwidth contention and serial global merging. See the notes above.`))

}

// KmerRecord stores a k-mer code and its associated genome IDs
type KmerRecord struct {
	code    uint64
	genomes []uint32
}

// computeProbabilityUpperBound computes the upper bound of P(X >= t*S | X = k, n partitions processed)
// using the Agievich bound approximation from the Onika paper
// (https://doi.org/10.1101/2025.11.21.689685, Section 2.3.1).
// Returns true if the probability is above the threshold.
// n: number of partitions processed so far
// k: number of matches observed
// t: minimum similarity threshold (minMaskFraction)
// S: total number of partitions (masks)
// probThreshold: minimum probability threshold
func shouldKeepPair(n, k int, t float64, S int, probThreshold float64) bool {
	// If no probability threshold, keep it
	// if probThreshold <= 0.0 {
	// 	return true
	// }

	// if n == 0 || n < k {
	// 	return true
	// }

	// We want to estimate if the pair can reach t*S matches in the remaining partitions
	requiredMatches := int(t * float64(S))

	// If already reached the threshold, keep the pair
	if k >= requiredMatches {
		return true
	}

	// If impossible to reach even if all remaining partitions match
	remaining := S - n
	if k+remaining < requiredMatches {
		return false
	}

	// Estimate the probability using the binomial approximation
	// Following Onika's implementation: use log space to avoid overflow
	fn := float64(n)
	fk := float64(k)

	// Use the observed rate or threshold, whichever is higher
	p := t
	if n > 0 {
		observedRate := fk / fn
		if observedRate > p {
			p = observedRate
		}
	}

	// Clamp p to avoid log(0)
	p = math.Max(1e-12, math.Min(1.0-1e-12, p))
	q := 1.0 - p

	// Compute log probability using Agievich approximation
	diff := fk - 0.5*fn

	// Log coefficient: n*ln(2) - 0.5*ln(pi*n/2) - 2*diff^2/n + 23/(18n)
	logCoeff := fn*math.Ln2 - 0.5*math.Log(math.Pi*fn/2.0) - 2.0*diff*diff/fn + 23.0/(18.0*fn)

	// Log mass: logCoeff + k*ln(p) + (n-k)*ln(q)
	logMass := logCoeff + fk*math.Log(p) + (fn-fk)*math.Log(q)

	// Clamp to 0 if positive (probability can't exceed 1)
	if logMass > 0.0 {
		return true
	}

	// Compare in log space
	logThreshold := math.Log(probThreshold)
	return logMass >= logThreshold
}

// Pool for reusing KmerRecord objects
var poolKmerRecord = &sync.Pool{New: func() interface{} {
	return &KmerRecord{
		genomes: make([]uint32, 0, 128),
	}
},
}

var poolGenomes = &sync.Pool{New: func() interface{} {
	tmp := make([]uint32, 4096)
	return &tmp
}}

const WindowInitialSize = 1 << 18

const (
	defaultMaxDenseGenomes   = 20_000 // 199,990,000 possible unordered pairs.
	maxSupportedDenseGenomes = 92_682 // Pair indexes are stored as uint32.
	// Budget for per-mask dense tables only. It controls worker concurrency,
	// not total RSS; larger values often become memory-bandwidth or collector bound.
	defaultDenseMaskMemoryMiB = 768
)

// DenseMaskCounts stores the best prefix for one mask in a triangular pair
// table. A zero prefix means unseen. The bitset records touched slots so merging
// and resetting cost O(number of observed pairs), not O(number of all pairs).
type DenseMaskCounts struct {
	prefixes []uint8
	touched  []uint64
}

func newDenseMaskCounts(nPairs int) *DenseMaskCounts {
	return &DenseMaskCounts{
		prefixes: make([]uint8, nPairs),
		touched:  make([]uint64, (nPairs+63)>>6),
	}
}

type DensePairStats struct {
	// uint16 is sufficient because the dense path is limited to 65,535 masks.
	matches     []uint16
	sumPrefixes []uint32
	// active marks pairs currently retained by probabilistic pruning.
	active []uint64
}

func newDensePairStats(nPairs int) *DensePairStats {
	return &DensePairStats{
		matches:     make([]uint16, nPairs),
		sumPrefixes: make([]uint32, nPairs),
		active:      make([]uint64, (nPairs+63)>>6),
	}
}

func densePairIndex(g1, g2, nGenomes uint32) int {
	if g1 > g2 {
		g1, g2 = g2, g1
	}
	return int(uint64(g1)*uint64(2*nGenomes-g1-1)/2 + uint64(g2-g1-1))
}

// densePairRowStart returns the first flat index for all pairs whose smaller
// genome ordinal is g1. Rows contain (g1,g1+1), ..., (g1,nGenomes-1).
func densePairRowStart(g1, nGenomes uint32) uint64 {
	return uint64(g1) * uint64(2*nGenomes-g1-1) / 2
}

func densePairGenomes(pairIndex uint32, nGenomes uint32) (uint32, uint32) {
	// Locate the triangular row by binary search, then obtain the column from
	// the offset within that row. This runs only while writing final results.
	idx := uint64(pairIndex)
	lo, hi := uint32(0), nGenomes
	for lo+1 < hi {
		mid := lo + (hi-lo)/2
		if densePairRowStart(mid, nGenomes) <= idx {
			lo = mid
		} else {
			hi = mid
		}
	}
	return lo, lo + 1 + uint32(idx-densePairRowStart(lo, nGenomes))
}

func mergeDenseMaskCounts(stats *DensePairStats, counts *DenseMaskCounts, processedMasks, remaining,
	requiredMatches, totalMasks int, minMaskFraction, probThreshold float64) {
	// As in the sparse path, all pairs first observed in this mask share the
	// same keep/drop decision because they each have exactly one match so far.
	shouldAddNewPair := probThreshold == 0
	if probThreshold > 0 && 1+remaining >= requiredMatches {
		shouldAddNewPair = shouldKeepPair(processedMasks, 1, minMaskFraction, totalMasks, probThreshold)
	}

	// Visit only pair slots touched by the worker. Clear both the prefix and
	// bitset in place before returning this table to the worker pool.
	for wordIndex, word := range counts.touched {
		for word != 0 {
			bit := bits.TrailingZeros64(word)
			pairIndex := wordIndex<<6 + bit
			if pairIndex < len(counts.prefixes) {
				prefixLen := counts.prefixes[pairIndex]
				if stats.matches[pairIndex] != 0 || shouldAddNewPair {
					if stats.matches[pairIndex] == 0 {
						stats.active[wordIndex] |= uint64(1) << bit
					}
					stats.matches[pairIndex]++
					stats.sumPrefixes[pairIndex] += uint32(prefixLen)
				}
				counts.prefixes[pairIndex] = 0
			}
			word &= word - 1
		}
		counts.touched[wordIndex] = 0
	}
}

func pruneDensePairStats(stats *DensePairStats, decisions []int8, processedMasks int,
	minMaskFraction float64, totalMasks int, probThreshold float64) []int8 {
	// The pruning decision depends only on the current match count. Cache it
	// once per count instead of recomputing the probability for every pair.
	if cap(decisions) <= processedMasks {
		decisions = make([]int8, processedMasks+1)
	} else {
		decisions = decisions[:processedMasks+1]
		clear(decisions)
	}

	// Scan active bits rather than all possible pair slots.
	for wordIndex, word := range stats.active {
		for word != 0 {
			bit := bits.TrailingZeros64(word)
			pairIndex := wordIndex<<6 + bit
			if pairIndex >= len(stats.matches) {
				break
			}
			matches := stats.matches[pairIndex]
			if matches > 1 {
				decision := decisions[matches]
				if decision == 0 {
					if shouldKeepPair(processedMasks, int(matches), minMaskFraction, totalMasks, probThreshold) {
						decision = 1
					} else {
						decision = -1
					}
					decisions[matches] = decision
				}
				if decision < 0 {
					stats.active[wordIndex] &^= uint64(1) << bit
					stats.matches[pairIndex] = 0
					stats.sumPrefixes[pairIndex] = 0
				}
			}
			word &= word - 1
		}
	}
	return decisions
}

func countDenseResults(stats *DensePairStats, requiredMatches int) int {
	n := 0
	for wordIndex, word := range stats.active {
		for word != 0 {
			bit := bits.TrailingZeros64(word)
			pairIndex := wordIndex<<6 + bit
			if pairIndex < len(stats.matches) && int(stats.matches[pairIndex]) >= requiredMatches {
				n++
			}
			word &= word - 1
		}
	}
	return n
}

func buildGenomeOrdinals(id2name map[uint64][]byte, genomeBatches int) ([]uint32, [][]uint32, error) {
	// Encoded batch/ref IDs may be sparse. Dense pair tables require contiguous
	// ordinals, while genomeCodes provides the lossless mapping back for output.
	genomeCodes := make([]uint32, 0, len(id2name))
	maxGenomeIdx := make([]int, genomeBatches)
	for i := range maxGenomeIdx {
		maxGenomeIdx[i] = -1
	}
	for id := range id2name {
		if id > math.MaxUint32 {
			return nil, nil, fmt.Errorf("genome index out of range: %d", id)
		}
		code := uint32(id)
		batch := int(code >> BITS_GENOME_IDX)
		if batch >= genomeBatches {
			return nil, nil, fmt.Errorf("genome batch out of range: %d", batch)
		}
		genomeIdx := int(code & MASK_GENOME_IDX)
		if genomeIdx > maxGenomeIdx[batch] {
			maxGenomeIdx[batch] = genomeIdx
		}
		genomeCodes = append(genomeCodes, code)
	}
	slices.Sort(genomeCodes)

	ordinals := make([][]uint32, genomeBatches)
	for batch, maxIdx := range maxGenomeIdx {
		if maxIdx < 0 {
			continue
		}
		ordinals[batch] = make([]uint32, maxIdx+1)
		for i := range ordinals[batch] {
			ordinals[batch][i] = math.MaxUint32
		}
	}
	for ordinal, code := range genomeCodes {
		batch := int(code >> BITS_GENOME_IDX)
		genomeIdx := int(code & MASK_GENOME_IDX)
		ordinals[batch][genomeIdx] = uint32(ordinal)
	}
	return genomeCodes, ordinals, nil
}

func genomeOrdinal(code uint32, ordinals [][]uint32) uint32 {
	return ordinals[int(code>>BITS_GENOME_IDX)][int(code&MASK_GENOME_IDX)]
}

type KmerWindow struct {
	records []*KmerRecord
	head    int
}

var poolKmerWindow = &sync.Pool{New: func() interface{} {
	return &KmerWindow{records: make([]*KmerRecord, 0, WindowInitialSize)}
}}

var poolMaskCounts = &sync.Pool{New: func() interface{} {
	tmp := make(map[uint64]uint8, 4096)
	return &tmp
}}

func setDensePairPrefix(dense *DenseMaskCounts, pairStarts []uint32, g1, g2 uint32, prefixLen uint8) {
	if g1 > g2 {
		g1, g2 = g2, g1
	}
	pairIndex := int(pairStarts[g1] + g2 - g1 - 1)
	if prefixLen > dense.prefixes[pairIndex] {
		// Mark a slot only on its first update in this mask; later updates only
		// replace its maximum prefix length.
		if dense.prefixes[pairIndex] == 0 {
			dense.touched[pairIndex>>6] |= uint64(1) << (pairIndex & 63)
		}
		dense.prefixes[pairIndex] = prefixLen
	}
}

func setMapPairPrefix(counts *map[uint64]uint8, g1, g2 uint32, prefixLen uint8) {
	if g1 > g2 {
		g1, g2 = g2, g1
	}
	key := uint64(g1)<<32 | uint64(g2)
	if prefixLen > (*counts)[key] {
		(*counts)[key] = prefixLen
	}
}

// processKmerWithWindow processes a k-mer against the sliding window.
func processKmerWithWindow(currentCode uint64, currentGenomes *[]uint32, window *KmerWindow,
	counts *map[uint64]uint8, dense *DenseMaskCounts, pairStarts []uint32,
	threshold uint64, kMinus32 int, minPrefix uint8) {
	records := window.records
	head := window.head

	// Clean up window: remove k-mers that are too far away
	for head < len(records) && currentCode-records[head].code >= threshold {
		// Return KmerRecord to pool
		records[head].genomes = records[head].genomes[:0]
		poolKmerRecord.Put(records[head])
		records[head] = nil
		head++
	}

	// Compact occasionally. Keeping a head index avoids copying on every k-mer
	// while retaining the full backing slice for the next mask.
	if head >= 4096 && head >= len(records)/2 {
		n := copy(records, records[head:])
		clear(records[n:])
		records = records[:n]
		head = 0
	}

	// Compare with all k-mers in the window
	var g1, g2 uint32
	var prefixLen uint8
	for i := head; i < len(records); i++ {
		// Calculate exact prefix length using XOR and leading zeros
		// prefixLen = (bits.LeadingZeros64(kmer1^kmer2) >> 1) + kMinus32
		prefixLen = uint8((bits.LeadingZeros64(currentCode^records[i].code) >> 1) + kMinus32)

		// Skip if prefix length is less than minimum
		if prefixLen < minPrefix {
			continue
		}

		// Cartesian product of genome IDs. Select the storage path outside the
		// inner loop; this loop dominates highly similar genome collections.
		if dense != nil {
			for _, g1 = range records[i].genomes {
				for _, g2 = range *currentGenomes {
					if g1 != g2 {
						setDensePairPrefix(dense, pairStarts, g1, g2, prefixLen)
					}
				}
			}
		} else {
			for _, g1 = range records[i].genomes {
				for _, g2 = range *currentGenomes {
					if g1 != g2 {
						setMapPairPrefix(counts, g1, g2, prefixLen)
					}
				}
			}
		}
	}

	// Also handle pairs within currentGenomes (same k-mer code, prefix = k)
	// This handles cases where multiple genomes share the exact same k-mer

	n := len(*currentGenomes)
	var i, j int
	if n > 1 {
		prefixLen = uint8(kMinus32 + 32) // full k-mer length
		if dense != nil {
			for i = 0; i < n; i++ {
				g1 = (*currentGenomes)[i]
				rowStart := int(pairStarts[g1]) - int(g1) - 1
				for j = i + 1; j < n; j++ {
					g2 = (*currentGenomes)[j]
					pairIndex := rowStart + int(g2)
					if prefixLen > dense.prefixes[pairIndex] {
						if dense.prefixes[pairIndex] == 0 {
							dense.touched[pairIndex>>6] |= uint64(1) << (pairIndex & 63)
						}
						dense.prefixes[pairIndex] = prefixLen
					}
				}
			}
		} else {
			for i = 0; i < n; i++ {
				for j = i + 1; j < n; j++ {
					g1, g2 = (*currentGenomes)[i], (*currentGenomes)[j]
					if g1 != g2 {
						setMapPairPrefix(counts, g1, g2, prefixLen)
					}
				}
			}
		}
	}

	// Get a KmerRecord from pool
	record := poolKmerRecord.Get().(*KmerRecord)
	record.code = currentCode
	record.genomes = append(record.genomes, (*currentGenomes)...)

	// Add current k-mer to window
	records = append(records, record)
	window.records = records
	window.head = head
}

type PairStats struct {
	matches   uint32
	sumPrefix uint32
}

// Collect results into slice for sorting
type PairResult struct {
	pair      uint64
	nMasks    int
	sumPrefix uint32
}

type PairResults []PairResult

func (s PairResults) Len() int { return len(s) }
func (s PairResults) Less(i, j int) bool {
	if s[i].nMasks == s[j].nMasks { // 1. number of matched masks
		if s[i].sumPrefix == s[j].sumPrefix { // 2. total matched bases
			return s[i].pair < s[j].pair // 3. the order in the index, just to keep the order stable
		}
		return s[i].sumPrefix > s[j].sumPrefix
	}

	return s[i].nMasks > s[j].nMasks
}
func (s PairResults) Swap(i, j int) { s[i], s[j] = s[j], s[i] }

type DensePairResults struct {
	pairs []uint32
	stats *DensePairStats
}

func (s DensePairResults) Len() int { return len(s.pairs) }
func (s DensePairResults) Less(i, j int) bool {
	idx1, idx2 := s.pairs[i], s.pairs[j]
	matches1, matches2 := s.stats.matches[idx1], s.stats.matches[idx2]
	if matches1 == matches2 {
		sum1, sum2 := s.stats.sumPrefixes[idx1], s.stats.sumPrefixes[idx2]
		if sum1 == sum2 {
			return idx1 < idx2
		}
		return sum1 > sum2
	}
	return matches1 > matches2
}
func (s DensePairResults) Swap(i, j int) { s.pairs[i], s.pairs[j] = s.pairs[j], s.pairs[i] }
