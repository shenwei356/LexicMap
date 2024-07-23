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
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"sync"
	"time"

	"github.com/dustin/go-humanize"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/seedposition"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/lexichash"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
	"github.com/vbauerster/mpb/v8"
	"github.com/vbauerster/mpb/v8/decor"
	"gonum.org/v1/gonum/stat"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
)

var seedPosCmd = &cobra.Command{
	Use:   "seed-pos",
	Short: "Extract and plot seed positions via reference name(s)",
	Long: `Extract and plot seed positions via reference name(s)

Attention:
  0. This command requires the index to be created with the flag --save-seed-pos in lexicmap index.
  1. Seed/K-mer positions (column pos) are 1-based.
     For reference genomes with multiple sequences, the sequences were
     concatenated to a single sequence with intervals of N's.
     So values of column pos_gnm and pos_seq might be different.
     The positions can be used to extract subsequence with 'lexicmap utils subseq'.
  2. All degenerate bases in reference genomes were converted to the lexicographic first bases.
     E.g., N was converted to A. Therefore, consecutive A's in output might be N's in the genomes.
	
Extra columns:
  Using -v/--verbose will output more columns:
     len_aaa,  length of consecutive A's.
     seq,      sequence between the previous and current seed.

Figures:
  Using -O/--plot-dir will write plots into given directory:
    - Histograms of seed distances.
	- Histograms of numbers of seeds in sliding windows.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		seq.ValidateSeq = false

		timeStart := time.Now()

		// ------------------------------

		dbDir := getFlagString(cmd, "index")
		if dbDir == "" {
			checkError(fmt.Errorf("flag -d/--index needed"))
		}

		allGenomes := getFlagBool(cmd, "all-refs")

		refnames := getFlagStringSlice(cmd, "ref-name")
		if !allGenomes && len(refnames) == 0 {
			checkError(fmt.Errorf("flag -n/--ref-name needed, or use -a/--all-refs for all ref genomes"))
		}

		outFile := getFlagString(cmd, "out-file")
		moreColumns := getFlagBool(cmd, "verbose")
		minDist := getFlagNonNegativeInt(cmd, "min-dist")
		maxOpenFiles := getFlagPositiveInt(cmd, "max-open-files")
		step := uint32(getFlagPositiveInt(cmd, "slid-step"))
		window := uint32(getFlagPositiveInt(cmd, "slid-window"))

		// ------------------------------

		plotDir := getFlagString(cmd, "plot-dir")
		force := getFlagBool(cmd, "force")

		outputPlotDir := plotDir != ""
		if outputPlotDir {
			makeOutDir(plotDir, force, "plot-dir")
		}

		bins := getFlagPositiveInt(cmd, "bins")

		colorIndex := getFlagPositiveInt(cmd, "color-index")
		if colorIndex > 7 {
			checkError(fmt.Errorf("unsupported color index"))
		}

		// percentiles := getFlagBool(cmd, "percentiles")
		width := vg.Length(getFlagPositiveFloat64(cmd, "width"))
		height := vg.Length(getFlagPositiveFloat64(cmd, "height"))
		plotExt := getFlagString(cmd, "plot-ext")
		if plotExt == "" {
			checkError(fmt.Errorf("the value of --plot-ext should not be empty"))
		}
		// plotTitle := getFlagBool(cmd, "plot-title")

		// ---------------------------------------------------------------

		// info file for the number of genome batches
		fileInfo := filepath.Join(dbDir, FileInfo)
		info, err := readIndexInfo(fileInfo)
		if err != nil {
			checkError(fmt.Errorf("failed to read info file: %s", err))
		}

		fileSeedLoc := filepath.Join(dbDir, DirGenomes, batchDir(0), FileSeedPositions)
		ok, err := pathutil.Exists(fileSeedLoc)
		if err != nil {
			checkError(fmt.Errorf("check index file structure: %s", err))
		}
		if !ok {
			log.Warningf("no seed position file detected in %s, which was not built with --save-seed-pos.", dbDir)
			return
		}

		// genomes.map file for mapping index to genome id
		m, err := readGenomeMapName2Idx(filepath.Join(dbDir, FileGenomeIndex))
		if err != nil {
			checkError(fmt.Errorf("failed to read genomes index mapping file: %s", err))
		}

		if allGenomes {
			refnames = refnames[:0]
			for ref := range m {
				refnames = append(refnames, ref)
			}
		}

		// pool of genome reader
		var hasGenomeRdrs bool
		var poolGenomeRdrs []chan *genome.Reader

		{
			openFiles := info.Chunks + 2
			if maxOpenFiles < openFiles {
				checkError(fmt.Errorf("invalid max open files: %d, should be >= %d", maxOpenFiles, openFiles))
			}

			// we can create genome reader pools
			// n := (maxOpenFiles - info.Chunks) / info.GenomeBatches
			// if n < 2 {
			// } else {
			// 	n >>= 1
			// 	if n > opt.NumCPUs {
			// 		n = opt.NumCPUs
			// 	}
			n := 1
			if opt.Verbose || opt.Log2File {
				log.Infof("creating genome reader pools, each batch with %d readers...", n)
			}
			poolGenomeRdrs = make([]chan *genome.Reader, info.GenomeBatches)
			for i := 0; i < info.GenomeBatches; i++ {
				poolGenomeRdrs[i] = make(chan *genome.Reader, n)
			}

			// parallelize it
			var wg sync.WaitGroup
			tokens := make(chan int, opt.NumCPUs)
			for i := 0; i < info.GenomeBatches; i++ {
				for j := 0; j < n; j++ {
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

			hasGenomeRdrs = true
			// }
		}

		// ---------------------------------------------------------------

		showProgressBar := len(refnames) > 1 && opt.Verbose

		// process bar
		var pbs *mpb.Progress
		var bar *mpb.Bar
		var chDuration chan time.Duration
		var doneDuration chan int
		if showProgressBar {
			pbs = mpb.New(mpb.WithWidth(40), mpb.WithOutput(os.Stderr))
			bar = pbs.AddBar(int64(len(refnames)),
				mpb.PrependDecorators(
					decor.Name("processed files: ", decor.WC{W: len("processed files: "), C: decor.DindentRight}),
					decor.Name("", decor.WCSyncSpaceR),
					decor.CountersNoUnit("%d / %d", decor.WCSyncWidth),
				),
				mpb.AppendDecorators(
					decor.Name("ETA: ", decor.WC{W: len("ETA: ")}),
					decor.EwmaETA(decor.ET_STYLE_GO, 5),
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

		fmt.Fprintf(outfh, "ref\tseqid\tpos_gnm\tpos_seq\tstrand\tdistance")
		if moreColumns {
			fmt.Fprintf(outfh, "\tlen_aaa\tseq")
		}
		fmt.Fprintln(outfh)

		type Ref2Locs struct {
			Ref         string
			GenomeBatch int
			GenomeIdx   int
			Locs        *[]uint32
			StartTime   time.Time

			Genome *genome.Genome
		}

		poolRef2Locs := &sync.Pool{New: func() interface{} {
			tmp := make([]uint32, 0, 40<<20)
			return &Ref2Locs{
				Locs: &tmp,
			}
		}}

		// seed position readers
		readers := make([]*seedposition.Reader, 0, info.GenomeBatches) // save for closing them in the end
		readerPools := make([]*sync.Pool, info.GenomeBatches)          // one for each genome batch
		for batch := 0; batch < info.GenomeBatches; batch++ {
			_batch := batch
			readerPools[batch] = &sync.Pool{New: func() interface{} {
				fileSeedLoc := filepath.Join(dbDir, DirGenomes, batchDir(_batch), FileSeedPositions)
				rdr, err := seedposition.NewReader(fileSeedLoc)
				if err != nil {
					checkError(fmt.Errorf("failed to read seed position data file: %s", err))
				}
				readers = append(readers, rdr)
				return rdr
			}}
		}

		var wg sync.WaitGroup
		ch := make(chan *Ref2Locs, opt.NumCPUs)
		tokens := make(chan int, opt.NumCPUs)
		done := make(chan int)

		// 2. receive and output
		var n int
		var nPlot int
		go func() {
			var pos2str, pos, pre uint32
			var dist int
			var refname string
			var v plotter.Values  // for seed distance
			var v2 plotter.Values // for number of seed
			var filePlot string
			var p *plot.Plot
			threadsFloat := float64(opt.NumCPUs)

			if outputPlotDir {
				v = make(plotter.Values, 0, 40<<20)
				v2 = make(plotter.Values, 0, 40<<10)
			}

			var tSeq *genome.Genome
			var rdr *genome.Reader
			var genomeIdx int

			seqRegions := poolSkipRegions.Get().(*[][2]int)
			var interval int
			var _n, _n2 int // len of concatenated seqs
			var slen int
			var ri, rs int
			var sseqid []byte

			var seqRegion [2]int
			var ws, we, end uint32 // start position of current window
			var nSeeds int
			var locs *[]uint32
			var nPos, pi, ps int
			k32p1 := uint32(info.K) - 1
			var first bool

			// var kp1 int = int(info.K) - 1

			for ref2locs := range ch {
				if len(*ref2locs.Locs) == 0 {
					if showProgressBar {
						chDuration <- time.Duration(float64(time.Since(ref2locs.StartTime)) / threadsFloat)
					}
					genome.RecycleGenome(ref2locs.Genome)
					poolRef2Locs.Put(ref2locs)
				}

				n++
				pre = 0
				refname = ref2locs.Ref
				genomeIdx = ref2locs.GenomeIdx

				if hasGenomeRdrs {
					rdr = <-poolGenomeRdrs[ref2locs.GenomeBatch]
				} else {
					fileGenome := filepath.Join(dbDir, DirGenomes, batchDir(ref2locs.GenomeBatch), FileGenomes)
					rdr, err = genome.NewReader(fileGenome)
					if err != nil {
						checkError(fmt.Errorf("failed to read genome data file: %s", err))
					}
				}

				*seqRegions = (*seqRegions)[:0]
				if ref2locs.Genome.NumSeqs > 1 {
					interval = (ref2locs.Genome.Len - ref2locs.Genome.GenomeSize) / (ref2locs.Genome.NumSeqs - 1)
				}
				_n = 0
				for _, slen = range ref2locs.Genome.SeqSizes {
					*seqRegions = append(*seqRegions, [2]int{_n, _n + slen - 1})
					_n += slen + interval
				}

				v = v[:0]
				pre = 0

				ri = 0                    // index of seq region
				rs = (*seqRegions)[ri][0] // end of that region
				sseqid = *ref2locs.Genome.SeqIDs[ri]
				_n2 = ref2locs.Genome.SeqSizes[ri]
				for _, pos2str = range *ref2locs.Locs {
					pos = pos2str >> 2

					if pos2str&1 > 0 && int(pos) > _n2 { // this is the first pos after an interval region
						ri++
						rs = (*seqRegions)[ri][0] // end of that region
						sseqid = *ref2locs.Genome.SeqIDs[ri]
						_n2 += interval + ref2locs.Genome.SeqSizes[ri]

						pre = uint32(rs)
					}
					dist = int(pos - pre)

					if dist < minDist {
						pre = pos
						continue
					}

					if outputPlotDir {
						v = append(v, float64(dist))
					}

					fmt.Fprintf(outfh, "%s\t%s\t%d\t%d\t%c\t%d", refname, sseqid, pos+1, int(pos)-rs+1, lexichash.Strands[pos2str>>1&1], dist)

					if moreColumns {
						tSeq, err = rdr.SubSeq(genomeIdx, int(pre), int(pos)-1)
						if err != nil {
							checkError(fmt.Errorf("failed to read subsequence: %s", err))
						}

						fmt.Fprintf(outfh, "\t%d\t%s", lengthAAs(tSeq.Seq), tSeq.Seq)

						genome.RecycleGenome(tSeq)
					}

					fmt.Fprintln(outfh)

					pre = pos
				}

				// sliding window
				v2 = v2[:0]
				locs = ref2locs.Locs
				nPos = len(*locs)
				ps = 0
				for _, seqRegion = range *seqRegions { // each sequence
					end = uint32(seqRegion[1]) - window
					// fmt.Printf("seq: %d-%d\n", seqRegion[0], end)

					nSeeds = 0
					first = true
					for ws = uint32(seqRegion[0]); ws <= end; ws += step { // each window
						we = ws + window - k32p1
						// fmt.Printf("  window: %d-%d\n", ws, we)

						for pi = ps; pi < nPos; pi++ {
							pos = (*locs)[pi] >> 2
							// fmt.Printf("    pi:%d, pos:%d\n", pi, pos)

							if first && pos >= ws+step {
								first = false
								ps = pi // start index of next window
								// fmt.Printf("      next ps:%d\n", ps)
							}

							if pos < we {
								nSeeds++
								// fmt.Printf("      in range. nSeeds:%d\n", nSeeds)
							} else { // out of current window
								v2 = append(v2, float64(nSeeds))
								// fmt.Printf("      out range. nSeeds:%d\n", nSeeds)
								nSeeds = 0
								first = true
								break
							}
						}
					}
				}

				if !outputPlotDir || len(v) == 0 { // no distance > -D
					// ---------------------------------------------------------

					if showProgressBar {
						chDuration <- time.Duration(float64(time.Since(ref2locs.StartTime)) / threadsFloat)
					}
					genome.RecycleGenome(ref2locs.Genome)
					poolRef2Locs.Put(ref2locs)

					// return or close genome reader
					if hasGenomeRdrs {
						poolGenomeRdrs[ref2locs.GenomeBatch] <- rdr
					} else {
						err = rdr.Close()
						if err != nil {
							checkError(fmt.Errorf("failed to close genome data file: %s", err))
						}
					}
					continue
				}

				nPlot++

				title := fmt.Sprintf("%s\n%s bp, %d contig(s), %s seeds",
					refname, humanize.Comma(int64(ref2locs.Genome.GenomeSize)),
					len(ref2locs.Genome.SeqSizes), humanize.Comma(int64(len(*ref2locs.Locs))))
				// -------------------------------------------------------------

				p = plot.New()

				h, err := plotter.NewHist(v, bins)
				if err != nil {
					checkError(err)
				}

				// h.Normalize(1)
				h.FillColor = plotutil.Color(0)
				p.Add(h)

				p.Title.Text = title
				p.Title.TextStyle.Font.Size = 16
				sort.Float64s(v)
				p.X.Label.Text = fmt.Sprintf("%s\n99th pctl=%.0f, 99.9th pctl=%.0f, median=%.0f, max=%.0f\n",
					"Seed distance (bp)", getPercentile(0.99, v), getPercentile(0.999, v), getPercentile(0.5, v), v[len(v)-1])

				p.Y.Label.Text = "Frequency"
				p.X.Label.TextStyle.Font.Size = 14
				p.Y.Label.TextStyle.Font.Size = 14
				p.X.Width = 1.5
				p.Y.Width = 1.5
				p.X.Tick.Width = 1.5
				p.Y.Tick.Width = 1.5
				p.X.Tick.Label.Font.Size = 12
				p.Y.Tick.Label.Font.Size = 12

				// Save image
				filePlot = filepath.Join(plotDir, refname+plotExt)
				checkError(p.Save(width*vg.Inch, height*vg.Inch, filePlot))

				// -------------------------------------------------------------

				p = plot.New()

				h, err = plotter.NewHist(v2, bins)
				if err != nil {
					checkError(err)
				}

				// h.Normalize(1)
				h.FillColor = plotutil.Color(0)
				p.Add(h)

				p.Title.Text = title
				p.Title.TextStyle.Font.Size = 16
				sort.Float64s(v2)
				// p.X.Label.Text = fmt.Sprintf("%s\n99th pctl=%.0f, 99.9th pctl=%.0f, median=%.0f, max=%.0f\n",
				// 	fmt.Sprintf("Number of seeds in %d-bp sliding windows", window), getPercentile(0.99, v2), getPercentile(0.999, v2), getPercentile(0.5, v), v2[len(v2)-1])
				m, M := minAndMax(v2)
				p.X.Label.Text = fmt.Sprintf("Number of seeds in %d-bp sliding windows of step %d-bp\nmin=%.0f, median=%.0f, max=%.0f",
					window, step, m, getPercentile(0.5, v2), M)

				p.Y.Label.Text = "Frequency"
				p.X.Label.TextStyle.Font.Size = 14
				p.Y.Label.TextStyle.Font.Size = 14
				p.X.Width = 1.5
				p.Y.Width = 1.5
				p.X.Tick.Width = 1.5
				p.Y.Tick.Width = 1.5
				p.X.Tick.Label.Font.Size = 12
				p.Y.Tick.Label.Font.Size = 12

				// Save image
				filePlot = filepath.Join(plotDir, refname+".seed_number"+plotExt)
				checkError(p.Save(width*vg.Inch, height*vg.Inch, filePlot))

				// ---------------------------------------------------------

				if showProgressBar {
					chDuration <- time.Duration(float64(time.Since(ref2locs.StartTime)) / threadsFloat)
				}
				genome.RecycleGenome(ref2locs.Genome)
				poolRef2Locs.Put(ref2locs)

				// return or close genome reader
				if hasGenomeRdrs {
					poolGenomeRdrs[ref2locs.GenomeBatch] <- rdr
				} else {
					err = rdr.Close()
					if err != nil {
						checkError(fmt.Errorf("failed to close genome data file: %s", err))
					}
				}
			}

			poolSkipRegions.Put(seqRegions)
			done <- 1
		}()

		// 1. extract data
		for _, refname := range refnames {
			tokens <- 1
			wg.Add(1)

			go func(refname string) {
				defer func() {
					wg.Done()
					<-tokens
				}()

				ref2locs := poolRef2Locs.Get().(*Ref2Locs)
				ref2locs.Ref = refname
				ref2locs.StartTime = time.Now()

				var batchIDAndRefID uint64
				var ok bool
				if batchIDAndRefID, ok = m[refname]; !ok {
					log.Warningf("reference name not found: %s", refname)
					ch <- ref2locs
					return
				}

				genomeBatch := int(batchIDAndRefID >> 17)
				genomeIdx := int(batchIDAndRefID & 131071)

				// seed positions ------------------------------------------------------------

				ref2locs.GenomeBatch = genomeBatch
				ref2locs.GenomeIdx = genomeIdx

				rdr := readerPools[genomeBatch].Get().(*seedposition.Reader)
				err = rdr.SeedPositions(genomeIdx, ref2locs.Locs)
				if err != nil {
					checkError(fmt.Errorf("failed to read seed position for %s: %s", refname, err))
				}

				readerPools[genomeBatch].Put(rdr)

				// sequence length ----------------------------------------------------------

				var gRdr *genome.Reader
				if hasGenomeRdrs {
					gRdr = <-poolGenomeRdrs[ref2locs.GenomeBatch]
				} else {
					fileGenome := filepath.Join(dbDir, DirGenomes, batchDir(ref2locs.GenomeBatch), FileGenomes)
					gRdr, err = genome.NewReader(fileGenome)
					if err != nil {
						checkError(fmt.Errorf("failed to read genome data file: %s", err))
					}
				}
				g, err := gRdr.GenomeInfo(genomeIdx)
				if err != nil {
					checkError(fmt.Errorf("failed to read genome infomation: %s", err))
				}
				ref2locs.Genome = g

				if hasGenomeRdrs {
					poolGenomeRdrs[ref2locs.GenomeBatch] <- gRdr
				} else {
					err = gRdr.Close()
					if err != nil {
						checkError(fmt.Errorf("failed to close genome data file: %s", err))
					}
				}

				ch <- ref2locs

			}(refname)
		}

		wg.Wait()
		close(ch)
		<-done

		// close all readers
		for _, rdr := range readers {
			checkError(rdr.Close())
		}

		// genome reader
		if hasGenomeRdrs {
			var _err error
			var wg sync.WaitGroup
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
			checkError(_err)
		}

		if showProgressBar {
			close(chDuration)
			<-doneDuration
			pbs.Wait()
		}

		if opt.Verbose {
			log.Infof("seed positions of %d genomes(s) saved to %s", n, outFile)
			if outputPlotDir {
				log.Infof("histograms of %d genomes(s) saved to %s", nPlot, plotDir)
			}
			log.Info()
			log.Infof("elapsed time: %s", time.Since(timeStart))
			log.Info()
		}
	},
}

func init() {
	utilsCmd.AddCommand(seedPosCmd)

	seedPosCmd.Flags().StringP("index", "d", "",
		formatFlagUsage(`Index directory created by "lexicmap index".`))

	seedPosCmd.Flags().StringSliceP("ref-name", "n", []string{},
		formatFlagUsage(`Reference name(s).`))

	seedPosCmd.Flags().BoolP("all-refs", "a", false,
		formatFlagUsage(`Output for all reference genomes. This would take a long time for an index with a lot of genomes.`))

	seedPosCmd.Flags().StringP("out-file", "o", "-",
		formatFlagUsage(`Out file, supports and recommends a ".gz" suffix ("-" for stdout).`))

	seedPosCmd.Flags().BoolP("verbose", "v", false,
		formatFlagUsage(`Show more columns including position of the previous seed and sequence between the two seeds. `+
			`Warning: it's slow to extract the sequences, recommend set -D 1000 or higher values to filter results `))

	seedPosCmd.Flags().IntP("min-dist", "D", 0,
		formatFlagUsage(`Only output records with seed distance >= this value.`))

	seedPosCmd.Flags().IntP("max-open-files", "", 512,
		formatFlagUsage(`Maximum opened files, used for extracting sequences.`))

	seedPosCmd.Flags().StringP("plot-dir", "O", "",
		formatFlagUsage(`Output directory for 1) histograms of seed distances, 2) histograms of numbers of seeds in sliding windows.`))

	seedPosCmd.Flags().IntP("slid-step", "s", 200,
		formatFlagUsage(`The step size of sliding windows for counting the number of seeds`))

	seedPosCmd.Flags().IntP("slid-window", "w", 500,
		formatFlagUsage(`The window size of sliding windows for counting the number of seeds`))

	seedPosCmd.Flags().BoolP("force", "", false,
		formatFlagUsage(`Overwrite existing output directory.`))

	// for histogram
	seedPosCmd.Flags().IntP("bins", "b", 100,
		formatFlagUsage(`Number of bins in histograms.`))
	seedPosCmd.Flags().IntP("color-index", "", 1,
		formatFlagUsage(`Color index (1-7).`))
	// seedPosCmd.Flags().BoolP("percentiles", "p", false,
	// 	formatFlagUsage(`Calculate percentiles`))
	seedPosCmd.Flags().Float64P("width", "", 6,
		formatFlagUsage(`Histogram width (unit: inch).`))
	seedPosCmd.Flags().Float64P("height", "", 4,
		formatFlagUsage(`Histogram height (unit: inch).`))
	seedPosCmd.Flags().StringP("plot-ext", "", ".png",
		formatFlagUsage(`Histogram plot file extention.`))
	// seedPosCmd.Flags().BoolP("plot-title", "t", false,
	// 	formatFlagUsage(`Plot genome ID and size as the title.`))

	seedPosCmd.SetUsageTemplate(usageTemplate(""))
}

func getPercentile(percentile float64, vals []float64) (p float64) {
	p = stat.Quantile(percentile, stat.Empirical, vals, nil)
	return
}

func lengthAAs(s []byte) int {
	var p, b byte
	var n int
	for _, b = range s {
		if b == 'A' && b == p {
			n++
		}

		p = b
	}
	return n
}

func minAndMax(vals []float64) (m, M float64) {
	m = math.MaxFloat64
	M = 0
	for _, v := range vals {
		if v > M {
			M = v
		}
		if v < m {
			m = v
		}
	}
	return m, M
}
