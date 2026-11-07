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
	"fmt"
	"math"
	"strconv"
	"strings"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
	"github.com/spf13/cobra"
	"github.com/spf13/pflag"
	"gonum.org/v1/plot/vg/vgimg"
)

const (
	// Resource counts allocate workers, queues, or tables directly. These CLI
	// limits reject accidental enormous allocations without limiting byte budgets.
	maxCLIWorkers      = 1024
	maxCLITableEntries = 1 << 20
	// Scanner lengths must also fit int on 32-bit platforms.
	maxScannerBufferSize = math.MaxInt32
	// Bound both sides and total raster work, using the plotting backend's DPI.
	maxPlotSidePixels = 1 << 16
	maxPlotPixels     = 1 << 26
)

// checkCommandOptions rejects invalid numbers before logs, outputs, or runtime
// worker pools are opened, including numeric flags unused by the selected mode.
func checkCommandOptions(cmd *cobra.Command, args []string) error {
	if err := checkNumericFlags(cmd); err != nil {
		return err
	}
	return checkCommandFiles(cmd, args)
}

// checkNumericFlags applies common representation and resource limits. Local
// command validation still checks relationships and index-dependent constraints.
func checkNumericFlags(cmd *cobra.Command) error {
	var invalid error
	cmd.Flags().VisitAll(func(flag *pflag.Flag) {
		if invalid != nil {
			return
		}
		name, text := flag.Name, flag.Value.String()
		switch flag.Value.Type() {
		case "float64", "float32":
			value, err := strconv.ParseFloat(text, 64)
			if err != nil || math.IsNaN(value) || math.IsInf(value, 0) {
				invalid = fmt.Errorf("--%s must be finite, got %s", name, text)
				return
			}
			upper := math.MaxFloat64
			switch name {
			case "min-af", "min-ani", "align-min-match-pident", "min-qcov-per-genome", "min-qcov-per-hsp":
				upper = 100
			case "min-mask-fraction", "prob-threshold":
				upper = 1
			case "width", "height":
				upper = float64(maxPlotSidePixels) / vgimg.DefaultDPI
				if value <= 0 {
					invalid = fmt.Errorf("--%s must be positive", name)
					return
				}
			}
			if value < 0 || value > upper {
				invalid = fmt.Errorf("--%s must be in [0, %g], got %s", name, upper, text)
			}
		case "int":
			value, err := strconv.ParseInt(text, 10, 64)
			if err != nil {
				invalid = fmt.Errorf("invalid --%s: %w", name, err)
				return
			}
			lower, upper := int64(0), int64(math.MaxInt32)
			switch name {
			case "seed", "rand-seed":
				lower, upper = 1, math.MaxInt // seeds do not size allocations
			case "min-seq-len":
				lower = -1
			case "threads", "max-query-conc", "seed-data-threads", "max-seed-matching-conc":
				upper = maxCLIWorkers
			case "partitions":
				upper = kv.MaxIndexPartitions
			case "masks", "max-open-files", "bins", "windows":
				upper = maxCLITableEntries
			case "kmer", "min-prefix", "seed-min-prefix", "seed-min-single-prefix", "prefix":
				upper = 32
			case "max-genome":
				upper = MAX_GENOME_SIZE
			case "batch-size":
				upper = 1 << BITS_GENOME_IDX
			case "max-dense-genomes":
				upper = maxSupportedDenseGenomes
			case "max-genome-size", "max-subject-genome-size":
				upper = min(upper, int64(math.MaxInt)/1_000_000) // Mb to bases
			case "upstream", "downstream", "align-ext-len":
				upper = min(upper, int64(math.MaxInt)/2) // endpoint additions
			}
			if value < lower || value > upper {
				invalid = fmt.Errorf("--%s must be in [%d, %d], got %s", name, lower, upper, text)
			}
		case "string":
			switch name {
			case "buffer-size", "max-seed-memory", "max-align-result-memory", "max-genome-cache-memory", "seed-index2-threshold", "threshold":
				value, err := ParseByteSize(text)
				if err != nil {
					invalid = fmt.Errorf("invalid --%s: %w", name, err)
					return
				}
				if name == "buffer-size" && (value < 1 || value > maxScannerBufferSize) {
					invalid = fmt.Errorf("--buffer-size must be in [1, %d] bytes", maxScannerBufferSize)
				}
			case "region":
				if text == "" {
					return
				}
				parts := strings.Split(text, ":")
				if len(parts) != 2 {
					invalid = fmt.Errorf("--region must be begin:end")
					return
				}
				for _, part := range parts {
					value, err := strconv.ParseInt(part, 10, 32)
					if err != nil || value < 1 || value > int64(math.MaxInt)/2 {
						invalid = fmt.Errorf("--region positions must be in [1, %d]", min(int64(math.MaxInt)/2, int64(math.MaxInt32)))
						return
					}
				}
			}
		}
	})
	if invalid != nil {
		return invalid
	}
	if cmd.Flags().Lookup("width") != nil && cmd.Flags().Lookup("height") != nil {
		width, _ := cmd.Flags().GetFloat64("width")
		height, _ := cmd.Flags().GetFloat64("height")
		if width*height*vgimg.DefaultDPI*vgimg.DefaultDPI > maxPlotPixels {
			return fmt.Errorf("--width and --height exceed %d pixels at %d DPI", maxPlotPixels, vgimg.DefaultDPI)
		}
	}
	return nil
}
