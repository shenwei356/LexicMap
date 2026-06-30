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
	"strings"
	"sync"
	"time"

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
  1. Genomes stored in multiple chunks are not evaluated as a whole.

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
		masks := make(map[uint64]struct{}, len(lh.Masks))
		if nMasks == 0 {
			for _, mask := range lh.Masks {
				masks[mask] = struct{}{}
			}
		} else {
			maskPrefix = int(math.Log2(float64(nMasks)) / 2)
			m := make(map[uint64]struct{}, maskPrefix)
			for _, mask := range lh.Masks {
				prefix := mask >> (uint64(lh.K-maskPrefix) << 1)
				if _, ok := m[prefix]; !ok {
					masks[mask] = struct{}{}

					m[prefix] = struct{}{}
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
			bar = pbs.AddBar(int64(len(masks)),
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

		// Active pairs tracking for probabilistic pruning
		activePairs := make(map[uint64]int, 10240) // pair -> number of matches

		// Global statistics across all masks - accumulate prefix length sums
		globalCounts := make(map[uint64]uint32, 10240) // gid1<<32|gid2 -> sum of prefix lengths

		// Calculate threshold for minimum prefix length
		// threshold = 1 << ((k - minPrefix) * 2)
		k := int(info.K)
		kMinus32 := k - 32 // precompute to avoid repeated calculation
		threshold := uint64(1) << ((k - minPrefix) * 2)
		minPrefixU8 := uint8(minPrefix) // convert to uint8 for comparison

		totalMasks := len(masks) // len(lh.Masks)
		requiredMatches := int(minMaskFraction * float64(totalMasks))

		if outputLog {
			log.Infof("  minimum prefix length between k-mers captured by a mask: %d", minPrefix)
			log.Infof("  total masks: %d, required matches: %d (%.1f%%)", totalMasks, requiredMatches, minMaskFraction*100)
		}

		// -------------------------------------------------------------------------
		// collect counting results

		type Result struct {
			Counts    *map[uint64]uint8
			StartTime time.Time
		}

		ch := make(chan *Result, opt.NumCPUs)
		done := make(chan int)
		go func() {
			var processedMasks, remaining int
			remaining = totalMasks
			var pair uint64
			var ok, shouldAddNewPair bool
			var matches int
			var prefixLen uint8

			for result := range ch {
				maskCounts := result.Counts

				processedMasks++
				remaining--

				if maskCounts == nil { // no k-mers
					if showProgressBar {
						chDuration <- time.Duration(float64(time.Since(result.StartTime)) / fcpus)
					}
					continue
				}

				if len(*maskCounts) == 0 { // no data
					poolMaskCounts.Put(maskCounts)

					if showProgressBar {
						chDuration <- time.Duration(float64(time.Since(result.StartTime)) / fcpus)
					}
					continue
				}

				if probThreshold == 0 { //  no pruning
					// Simply accumulate all pairs
					for pair, prefixLen = range *maskCounts {
						activePairs[pair]++
						globalCounts[pair] += uint32(prefixLen)
					}

				} else {

					// Check if new pairs can still reach the threshold
					if 1+remaining >= requiredMatches {
						// Pre-compute probability check for new pairs (count=1)
						shouldAddNewPair = shouldKeepPair(processedMasks, 1, minMaskFraction, totalMasks, probThreshold)
					}

					// Update match counts for pairs that matched in this mask
					for pair, prefixLen = range *maskCounts {
						matches, ok = activePairs[pair]
						if !ok {
							// New pair: check if it passes probability check
							if shouldAddNewPair {
								activePairs[pair] = 1
								globalCounts[pair] += uint32(prefixLen)
							}
						} else {
							// Existing pair: increment count
							activePairs[pair] = matches + 1
							globalCounts[pair] += uint32(prefixLen)
						}
					}

					// Probabilistic pruning: check all active pairs to remove impossible ones early
					if processedMasks < totalMasks && processedMasks&7 == 0 {
						for pair, matches = range activePairs {
							if matches > 1 && !shouldKeepPair(processedMasks, matches, minMaskFraction, totalMasks, probThreshold) {
								delete(activePairs, pair)
								delete(globalCounts, pair)
							}
						}
					}
				}

				clear(*maskCounts)
				poolMaskCounts.Put(maskCounts)

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
		tokens := make(chan int, opt.NumCPUs)

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

				r := bufio.NewReaderSize(fh, 4096)

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
				_, _, indexes, _, _, _, err := kv.ReadKVIndex(filepath.Clean(fileSeeds) + kv.KVIndexFileExt)
				if err != nil {
					checkError(fmt.Errorf("failed to read kv-data index file: %s", err))
				}

				// -------------------------------
				// data of all masks

				buf := make([]byte, 64)
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
				var i int

				var mask uint64

				for iMask := 0; iMask < nMasks; iMask++ {
					mask = lh.Masks[iFirstMask+iMask]
					if _, ok := masks[mask]; !ok { // not wanted
						continue
					}

					timeStart := time.Now()

					if len(indexes[iMask]) == 0 { // no k-mers
						ch <- &Result{
							Counts:    nil,
							StartTime: timeStart,
						}
						continue
					}

					// genome id list
					genomes := poolGenomes.Get().(*[]uint32)

					// Sliding window for all-to-all comparison
					window := poolKmerWindow.Get().(*[]*KmerRecord)

					// Per-mask tracking: which genomes appear in this mask
					// local counts for this mask (max prefix per pair)
					maskCounts := poolMaskCounts.Get().(*map[uint64]uint8)

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

						// ------------------ values for kmer1 -------------------

						*genomes = (*genomes)[:0] // reuse slice
						for j = 0; j < lenVal1; j++ {
							nReaded, err = io.ReadFull(r, buf[:bytesPos])
							if nReaded < bytesPos {
								checkError(kv.ErrBrokenFile)
							}

							v = fUint64(buf[:bytesPos])
							if v&MASK_REVERSE == 1 {
								continue // skip reverse complement
							}
							// Extract genome ID (batchID + refID)
							batchIDAndRefID = (v >> BITS_NONE_IDX) & 4294967295
							*genomes = append(*genomes, uint32(batchIDAndRefID))
						}

						// Process kmer1 with sliding window
						if len(*genomes) > 0 {
							processKmerWithWindow(kmer1, genomes, window, maskCounts, threshold, kMinus32, minPrefixU8) // , blacklist)
						}

						if lastPair && !hasKmer2 {
							break
						}

						// ------------------ values for kmer2 -------------------

						*genomes = (*genomes)[:0] // reuse slice
						for j = 0; j < lenVal2; j++ {
							nReaded, err = io.ReadFull(r, buf[:bytesPos])
							if nReaded < bytesPos {
								checkError(kv.ErrBrokenFile)
							}

							v = fUint64(buf[:bytesPos])
							if v&MASK_REVERSE == 1 {
								continue // skip reverse complement
							}
							batchIDAndRefID = (v >> BITS_NONE_IDX) & 4294967295
							*genomes = append(*genomes, uint32(batchIDAndRefID))
						}

						// Process kmer2 with sliding window
						if len(*genomes) > 0 {
							processKmerWithWindow(kmer2, genomes, window, maskCounts, threshold, kMinus32, minPrefixU8) // , blacklist)
						}

						if lastPair {
							break
						}
					}

					// recycle objects and send result
					poolGenomes.Put(genomes)

					for i = range *window {
						(*window)[i].genomes = (*window)[i].genomes[:0]
						poolKmerRecord.Put((*window)[i])
					}
					*window = (*window)[:0]
					poolKmerWindow.Put(window)

					ch <- &Result{
						Counts:    maskCounts,
						StartTime: timeStart,
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

		id2name, err := readGenomeMapIdx2Name(filepath.Join(dbDir, FileGenomeIndex))
		if err != nil {
			checkError(fmt.Errorf("failed to read %s: %s", filepath.Join(dbDir, FileGenomeIndex), err))
		}

		results := make([]PairResult, 0, len(globalCounts))
		var matchedMasks int
		for pair, sumPrefix := range globalCounts {
			// Only output pairs that meet the required threshold
			matchedMasks = activePairs[pair]
			if matchedMasks >= requiredMatches {
				results = append(results, PairResult{
					pair:      pair,
					nMasks:    matchedMasks,
					sumPrefix: sumPrefix,
				})
			}
		}
		if outputLog {
			log.Info()
			log.Infof("total genome pairs: %d", len(results))
		}

		// Sort by nMasks (then sumPrefix) in descending order
		sorts.Quicksort(PairResults(results))

		// Write header
		outfh.WriteString("genome1\tgenome2\tminPrefix\tfracMasks\tnMasks\tsumPrefix\tavgPrefix\n")

		// Write sorted results
		var gid1, gid2 uint64
		for _, result := range results {
			gid1 = result.pair >> 32
			gid2 = result.pair & 0xFFFFFFFF

			fracMasks := float64(result.nMasks) / float64(totalMasks)
			fmt.Fprintf(outfh, "%s\t%s\t%d\t%.4f\t%d\t%d\t%.2f\n",
				id2name[gid1], id2name[gid2], minPrefix, fracMasks, result.nMasks,
				result.sumPrefix, float64(result.sumPrefix)/float64(result.nMasks))
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
		formatFlagUsage(`Number of LexicHash masks to use. It should be 0 (for all masks in the index) or power of 4 (needs to be >= 64, e.g., 64, 256, 1024, 4096, 16384).`))

	pairCmd.Flags().StringP("out-file", "o", "-",
		formatFlagUsage(`Out file, supports and recommends a ".gz" suffix ("-" for stdout).`))

	pairCmd.SetUsageTemplate(usageTemplate("-d <index path> [-o out.tsv.gz]"))

	pairCmd.Flags().IntP("min-prefix", "p", 21,
		formatFlagUsage(`Minimum prefix length between k-mers captured by a mask.`))

	pairCmd.Flags().Float64P("min-mask-fraction", "f", 0.25,
		formatFlagUsage(`Minimum fraction of masks that must match for a genome pair to be reported.`))

	pairCmd.Flags().Float64P("prob-threshold", "s", 0.001,
		formatFlagUsage(`Probabilistic threshold for early termination heuristic (lower = more aggressive pruning， 0 = disable pruning).`))

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

var poolKmerWindow = &sync.Pool{New: func() interface{} {
	tmp := make([]*KmerRecord, 0, WindowInitialSize)
	return &tmp
}}

var poolMaskCounts = &sync.Pool{New: func() interface{} {
	tmp := make(map[uint64]uint8, 4096)
	return &tmp
}}

// processKmerWithWindow processes a k-mer against the sliding window
func processKmerWithWindow(currentCode uint64, currentGenomes *[]uint32, window *[]*KmerRecord, counts *map[uint64]uint8, threshold uint64, kMinus32 int, minPrefix uint8) { // , blacklist *sync.Map) {
	// only move elements when waste exceeds this threshold
	// const moveThreshold = 8

	// Clean up window: remove k-mers that are too far away
	windowStart := 0
	for windowStart < len(*window) && currentCode-(*window)[windowStart].code >= threshold {
		// Return KmerRecord to pool
		(*window)[windowStart].genomes = (*window)[windowStart].genomes[:0]
		poolKmerRecord.Put((*window)[windowStart])
		windowStart++
	}

	// Move valid elements to the front only when waste exceeds threshold
	// if windowStart >= moveThreshold {
	// 	if windowStart < len(*window) {
	// 		copy(*window, (*window)[windowStart:])
	// 	}
	// 	*window = (*window)[:len(*window)-windowStart]
	// } else
	if windowStart > 0 {
		// Just trim the slice without copying
		*window = (*window)[windowStart:]
	}

	// Compare with all k-mers in the window
	var key uint64
	var g1, g2 uint32
	// var ok bool
	var prefixLen uint8
	for i := range *window {
		// Calculate exact prefix length using XOR and leading zeros
		// prefixLen = (bits.LeadingZeros64(kmer1^kmer2) >> 1) + kMinus32
		prefixLen = uint8((bits.LeadingZeros64(currentCode^(*window)[i].code) >> 1) + kMinus32)

		// Skip if prefix length is less than minimum
		if prefixLen < minPrefix {
			continue
		}

		// Cartesian product of genome IDs
		for _, g1 = range (*window)[i].genomes {
			for _, g2 = range *currentGenomes {
				if g1 == g2 {
					continue // skip self-comparison
				}

				// Ensure gid1 < gid2 for consistent key
				if g1 < g2 {
					key = uint64(g1)<<32 | uint64(g2)
				} else {
					key = uint64(g2)<<32 | uint64(g1)
				}

				// // Skip impossible pairs
				// if _, ok = blacklist.Load(key); ok {
				// 	continue
				// }

				// Keep maximum prefix length within this mask
				if prefixLen > (*counts)[key] {
					(*counts)[key] = prefixLen
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
		for i = 0; i < n; i++ {
			for j = i + 1; j < n; j++ {
				g1, g2 = (*currentGenomes)[i], (*currentGenomes)[j]
				if g1 == g2 {
					continue // skip self-comparison
				}

				if g1 < g2 {
					key = uint64(g1)<<32 | uint64(g2)
				} else {
					key = uint64(g2)<<32 | uint64(g1)
				}

				// // skip impossible pairs
				// if _, ok = blacklist.Load(key); ok {
				// 	continue
				// }

				if prefixLen > (*counts)[key] { // it's definitely the longest
					(*counts)[key] = prefixLen
				}
			}
		}
	}

	// Get a KmerRecord from pool
	record := poolKmerRecord.Get().(*KmerRecord)
	record.code = currentCode
	record.genomes = append(record.genomes, (*currentGenomes)...)

	if cap(*window) == len(*window) { // need to resize
		if len(*window) == 0 {
			*window = make([]*KmerRecord, 0, WindowInitialSize)
		} else {
			tmp := make([]*KmerRecord, len(*window), WindowInitialSize)
			copy(tmp, *window)
			*window = tmp
		}
	}

	// Add current k-mer to window
	*window = append(*window, record)
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
