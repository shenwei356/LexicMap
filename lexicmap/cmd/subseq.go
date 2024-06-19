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
	"path/filepath"
	"regexp"
	"strconv"
	"strings"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
	"github.com/shenwei356/bio/seq"
	"github.com/spf13/cobra"
)

var subseqCmd = &cobra.Command{
	Use:   "subseq",
	Short: "Extract subsequence via reference name, sequence ID, position and strand",
	Long: `Exextract subsequence via reference name, sequence ID, position and strand

Attention:
  1. The option -s/--seq-id is optional.
     1) If given, the positions are these in the original sequence.
     2) If not given, the positions are these in the concatenated sequence.
  2. All degenerate bases in reference genomes were converted to the lexicographic first bases.
     E.g., N was converted to A. Therefore, consecutive A's in output might be N's in the genomes.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		seq.ValidateSeq = false

		// ------------------------------

		dbDir := getFlagString(cmd, "index")
		if dbDir == "" {
			checkError(fmt.Errorf("flag -d/--index needed"))
		}

		refname := getFlagString(cmd, "ref-name")
		if refname == "" {
			checkError(fmt.Errorf("flag -n/--ref-name needed"))
		}

		seqid := getFlagString(cmd, "seq-id")
		var concatenatedPositions bool
		if seqid == "" {
			concatenatedPositions = true
		}

		var reRegion = regexp.MustCompile(`\-?\d+:\-?\d+`)

		region := getFlagString(cmd, "region")
		if region == "" {
			checkError(fmt.Errorf("flag -r/--region needed"))
		}
		revcom := getFlagBool(cmd, "revcom")

		lineWidth := getFlagNonNegativeInt(cmd, "line-width")

		if !reRegion.MatchString(region) {
			checkError(fmt.Errorf(`invalid region: %s. type "lexicmap utils subseq -h" for more examples`, region))
		}
		var start, end int
		var err error

		r := strings.Split(region, ":")
		start, err = strconv.Atoi(r[0])
		checkError(err)
		end, err = strconv.Atoi(r[1])
		checkError(err)
		if start <= 0 || end <= 0 {
			checkError(fmt.Errorf("both begin and end position should not be <= 0"))
		}
		if start > end {
			checkError(fmt.Errorf("begin position should be < end position"))
		}

		outFile := getFlagString(cmd, "out-file")

		// ---------------------------------------------------------------

		// genomes.map file for mapping index to genome id
		m, err := readGenomeMapName2Idx(filepath.Join(dbDir, FileGenomeIndex))
		if err != nil {
			checkError(fmt.Errorf("failed to read genomes index mapping file: %s", err))
		}

		var batchIDAndRefID uint64
		var ok bool
		if batchIDAndRefID, ok = m[refname]; !ok {
			checkError(fmt.Errorf("reference name not found: %s", refname))
		}

		genomeBatch := int(batchIDAndRefID >> 17)
		genomeIdx := int(batchIDAndRefID & 131071)

		fileGenome := filepath.Join(dbDir, DirGenomes, batchDir(genomeBatch), FileGenomes)
		rdr, err := genome.NewReader(fileGenome)
		if err != nil {
			checkError(fmt.Errorf("failed to read genome data file: %s", err))
		}

		var tSeq *genome.Genome
		if concatenatedPositions {
			tSeq, err = rdr.SubSeq(genomeIdx, start-1, end-1)
		} else {
			tSeq, end, err = rdr.SubSeq2(genomeIdx, []byte(seqid), start-1, end-1)
			end++ // returned end is 0-based.
		}
		if err != nil {
			checkError(fmt.Errorf("failed to read subsequence: %s", err))
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

		s, err := seq.NewSeq(seq.DNAredundant, tSeq.Seq)
		checkError(err)
		if revcom {
			s.RevComInplace()
		}

		if concatenatedPositions {
			fmt.Fprintf(outfh, ">%s:%d-%d\n", refname, start, end)
		} else {
			fmt.Fprintf(outfh, ">%s:%d-%d\n", seqid, start, end)
		}
		outfh.Write(s.FormatSeq(lineWidth))
		outfh.WriteByte('\n')

		genome.RecycleGenome(tSeq)
		checkError(rdr.Close())
	},
}

func init() {
	utilsCmd.AddCommand(subseqCmd)

	subseqCmd.Flags().StringP("index", "d", "",
		formatFlagUsage(`Index directory created by "lexicmap index".`))

	subseqCmd.Flags().StringP("ref-name", "n", "",
		formatFlagUsage(`Reference name.`))

	subseqCmd.Flags().StringP("seq-id", "s", "",
		formatFlagUsage(`Sequence ID. If the value is empty, the positions in the region are treated as that in the concatenated sequence.`))

	subseqCmd.Flags().StringP("out-file", "o", "-",
		formatFlagUsage(`Out file, supports the ".gz" suffix ("-" for stdout).`))

	subseqCmd.Flags().StringP("region", "r", "",
		formatFlagUsage(`Region of the subsequence (1-based).`))

	subseqCmd.Flags().BoolP("revcom", "R", false,
		formatFlagUsage("Extract subsequence on the negative strand."))

	subseqCmd.Flags().IntP("line-width", "w", 60,
		formatFlagUsage("Line width of sequence (0 for no wrap)."))

	subseqCmd.SetUsageTemplate(usageTemplate(""))
}
