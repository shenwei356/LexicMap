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
	"bufio"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"

	"github.com/shenwei356/bio/seq"
	"github.com/spf13/cobra"
)

var genomesCmd = &cobra.Command{
	Use:   "genomes",
	Short: "View genome IDs in the index",
	Long: `View genome IDs in the index

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
		// read genome chunks data if existed
		genomeChunks, err := readGenomeChunksMapBig2Small(filepath.Join(dbDir, FileGenomeChunks))
		if err != nil {
			checkError(fmt.Errorf("failed to read genome chunk file: %s", err))
		}
		var hasGenomeChunks bool
		if len(genomeChunks) > 0 {
			hasGenomeChunks = true
		}

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
		var batchIDAndRefID uint64
		var ok bool

		outfh.WriteString("ref\tchunked\n")
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

			if hasGenomeChunks {
				if _, ok = genomeChunks[batchIDAndRefID]; ok {
					fmt.Fprintf(outfh, "%s\t%s\n", id, "yes")
				} else {
					fmt.Fprintf(outfh, "%s\t\n", id)
				}
			} else {
				fmt.Fprintf(outfh, "%s\t\n", id)
			}

		}
	},
}

func init() {
	utilsCmd.AddCommand(genomesCmd)

	genomesCmd.Flags().StringP("index", "d", "",
		formatFlagUsage(`Index directory created by "lexicmap index".`))

	genomesCmd.Flags().StringP("out-file", "o", "-",
		formatFlagUsage(`Out file, supports the ".gz" suffix ("-" for stdout).`))

	genomesCmd.SetUsageTemplate(usageTemplate(""))
}
