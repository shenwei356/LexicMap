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
	"bytes"
	"fmt"
	"path/filepath"
	"strings"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
	"github.com/shenwei356/bio/seq"
	"github.com/spf13/cobra"
)

var gseqCmd = &cobra.Command{
	Use:   "genome-seqs",
	Short: "Extract all sequences of a given genome",
	Long: `Extract all sequences of a given genome

Attention:
  1. All degenerate bases in reference genomes were converted to the lexicographic first bases.
     E.g., N was converted to A. Therefore, consecutive A's in output might be N's in the genomes.
  2. Large genomes fragmented into multiple chunks during indexing (total size > 15 Mb by default,
     configurable with -g/--max-genome in 'lexicmap index'), such as many fungal genomes,
     may have their sequence order rearranged relative to the original input files.

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
		lineWidth := getFlagNonNegativeInt(cmd, "line-width")

		refname := getFlagString(cmd, "ref-name")
		if refname == "" {
			checkError(fmt.Errorf("flag -n/--ref-name needed"))
		}

		// ---------------------------------------------------------------

		// info file
		fileInfo := filepath.Join(dbDir, FileInfo)
		info, err := readIndexInfo(fileInfo)
		if err != nil {
			checkError(fmt.Errorf("failed to read info file: %s", err))
		}
		if info.MainVersion != MainVersion {
			checkError(fmt.Errorf("index main versions do not match: %d (index) != %d (tool). please re-create the index", info.MainVersion, MainVersion))
		}

		// genomes.map file for mapping index to genome id
		m, err := readGenomeMapName2Idx(filepath.Join(dbDir, FileGenomeIndex))
		if err != nil {
			checkError(fmt.Errorf("failed to read genomes index mapping file: %s", err))
		}

		// ---------------------------------------------------------------

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

		// ---------------------------------------------------------------

		var batchIDAndRefIDs *[]uint64

		var ok bool
		if batchIDAndRefIDs, ok = m[refname]; !ok {
			checkError(fmt.Errorf("reference name not found: %s", refname))
		}

		var fileGenome string
		var g *genome.Genome
		var genomeBatch, genomeIdx int
		var rdr *genome.Reader
		var i int
		var s *[]byte
		var text []byte
		var buffer *bytes.Buffer

		for _, batchIDAndRefID := range *batchIDAndRefIDs {
			genomeBatch = int(batchIDAndRefID >> BITS_GENOME_IDX)
			genomeIdx = int(batchIDAndRefID & MASK_GENOME_IDX)

			fileGenome = filepath.Join(dbDir, DirGenomes, batchDir(genomeBatch), FileGenomes)
			rdr, err = genome.NewReader(fileGenome)
			if err != nil {
				checkError(fmt.Errorf("failed to read genome data file: %s", err))
			}

			g, err = rdr.Seqs(genomeIdx)
			checkError(err)

			for i, s = range g.Seqs {
				outfh.Write(_mark_fasta)
				outfh.Write(*g.SeqIDs[i])
				outfh.Write(_mark_newline)

				text, buffer = wrapByteSlice(*s, lineWidth, buffer)
				outfh.Write(text)
				outfh.Write(_mark_newline)
			}

			genome.RecycleGenome(g)
		}

		checkError(rdr.Close())
	},
}

func init() {
	utilsCmd.AddCommand(gseqCmd)

	gseqCmd.Flags().StringP("index", "d", "",
		formatFlagUsage(`Index directory created by "lexicmap index".`))

	gseqCmd.Flags().StringP("ref-name", "n", "",
		formatFlagUsage(`Reference name.`))

	gseqCmd.Flags().StringP("out-file", "o", "-",
		formatFlagUsage(`Out file, supports the ".gz" suffix ("-" for stdout).`))

	gseqCmd.Flags().IntP("line-width", "w", 60,
		formatFlagUsage("Line width of sequence (0 for no wrap)."))

	gseqCmd.SetUsageTemplate(usageTemplate(""))
}
