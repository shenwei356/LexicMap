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
	"sort"
	"strings"
	"sync"
	"time"

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
	Short: "Extract seed positions via reference names",
	Long: `Extract seed positions via reference names

Attentions:
  1. Seed/K-mer positions (column pos) are 1-based.
     For reference genomes with multiple sequences, the sequences were
     concatenated to a single sequence with intervals of (k-1) N's.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		seq.ValidateSeq = false

		// ------------------------------

		dbDir := getFlagString(cmd, "index")
		if dbDir == "" {
			checkError(fmt.Errorf("flag -d/--index needed"))
		}

		allGenomes := getFlagBool(cmd, "all-refs")

		refnames := getFlagStringSlice(cmd, "ref-name")
		if !allGenomes && len(refnames) == 0 {
			checkError(fmt.Errorf("flag -n/--ref-name needed"))
		}

		outFile := getFlagString(cmd, "out-file")

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
		plotTitle := getFlagBool(cmd, "plot-title")

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

		fmt.Fprintf(outfh, "ref\tpos\tstrand\tdistance\n")

		type Ref2Locs struct {
			Ref       string
			Locs      *[]uint32
			StartTime time.Time
		}

		poolRef2Locs := &sync.Pool{New: func() interface{} {
			tmp := make([]uint32, 0, 40<<20)
			return &Ref2Locs{
				Locs: &tmp,
			}
		}}

		// readers
		var readers []*seedposition.Reader                    // save for closing them in the end
		readerPools := make([]*sync.Pool, info.GenomeBatches) // one for each genome batch
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
		go func() {
			var pos2str, pos, pre uint32
			var refname string
			v := make(plotter.Values, 0, 40<<20)
			var filePlot string
			var p *plot.Plot
			threadsFloat := float64(opt.NumCPUs)

			for ref2locs := range ch {
				if len(*ref2locs.Locs) == 0 {
					if showProgressBar {
						chDuration <- time.Duration(float64(time.Since(ref2locs.StartTime)) / threadsFloat)
					}
					poolRef2Locs.Put(ref2locs)
				}

				n++
				pre = 0
				refname = ref2locs.Ref
				for _, pos2str = range *ref2locs.Locs {
					pos = pos2str >> 1
					fmt.Fprintf(outfh, "%s\t%d\t%c\t%d\n", refname, pos+1, lexichash.Strands[pos2str&1], pos-pre)
					pre = pos
				}

				if !outputPlotDir {
					if showProgressBar {
						chDuration <- time.Duration(float64(time.Since(ref2locs.StartTime)) / threadsFloat)
					}
					poolRef2Locs.Put(ref2locs)
					continue
				}

				// ---------------------------------------------------------
				// plot histogram

				v = v[:0]
				pre = 0
				for _, pos2str = range *ref2locs.Locs {
					pos = pos2str >> 1
					v = append(v, float64(pos-pre))
					pre = pos
				}

				p = plot.New()

				h, err := plotter.NewHist(v, bins)
				if err != nil {
					checkError(err)
				}

				// h.Normalize(1)
				h.FillColor = plotutil.Color(0)
				p.Add(h)

				if plotTitle {
					p.Title.Text = refname
				} else {
					p.Title.Text = ""
				}
				p.Title.TextStyle.Font.Size = 16
				//if percentiles {
				sort.Float64s(v)
				p.X.Label.Text = fmt.Sprintf("%s\n99th pctl=%.0f, 99.9th pctl=%.0f, median=%.0f, max=%.0f\n",
					"Seed distance (bp)", getPercentile(0.99, v), getPercentile(0.999, v), getPercentile(0.5, v), v[len(v)-1])
				// } else {
				// 	p.X.Label.Text = "Seed distance (bp)"
				// }
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

				// ---------------------------------------------------------

				if showProgressBar {
					chDuration <- time.Since(ref2locs.StartTime)
				}
				poolRef2Locs.Put(ref2locs)
			}
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

				rdr := readerPools[genomeBatch].Get().(*seedposition.Reader)
				err = rdr.SeedPositions(genomeIdx, ref2locs.Locs)
				if err != nil {
					checkError(fmt.Errorf("failed to read seed position for %s: %s", refname, err))
				}

				readerPools[genomeBatch].Put(rdr)

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

		if showProgressBar {
			close(chDuration)
			<-doneDuration
			pbs.Wait()
		}

		if opt.Verbose {
			log.Infof("seed positions of %d genomes(s) saved to %s", n, outFile)
			log.Infof("histograms of %d genomes(s) saved to %s", n, plotDir)
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

	seedPosCmd.Flags().StringP("plot-dir", "O", "",
		formatFlagUsage(`Output directory for histograms of seed distances.`))

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
	seedPosCmd.Flags().BoolP("plot-title", "t", false,
		formatFlagUsage(`Plot genome ID as the title.`))

	seedPosCmd.SetUsageTemplate(usageTemplate(""))
}

func getPercentile(percentile float64, vals []float64) (p float64) {
	p = stat.Quantile(percentile, stat.Empirical, vals, nil)
	return
}
