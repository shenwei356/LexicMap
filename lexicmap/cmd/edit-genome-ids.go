// Copyright © 2023-2025 Wei Shen <shenwei356@gmail.com>
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
	"fmt"
	"io"
	"os"
	"path/filepath"
	"regexp"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/util/pathutil"
	"github.com/spf13/cobra"
)

var editGenomeIDsCmd = &cobra.Command{
	Use:   "edit-genome-ids",
	Short: "Edit genome IDs in the index via a regular expression",
	Long: `Edit genome IDs in the index via a regular expression

Use cases:
  In the 'lexicmap index' command, users might forget to use the flag
  -N/--ref-name-regexp to extract the genome ID from the sequence file.
  A genome file from NCBI looks like:

    GCF_009818595.1_ASM981859v1_genomic.fna.gz

  In this case, the genome ID would be GCF_009818595.1_ASM981859v1_genomic,
  which is too long. So we can use this command to extract the assembly
  accession via:

    lexicmap utils edit-genome-ids -d t.lmi/ -p '^(\w{3}_\d{9}\.\d+).*' -r '$1'

Tips:
  - A backup file (genomes.map.bin.bak) will be created on the first run.

`,
	Run: func(cmd *cobra.Command, args []string) {
		// opt := getOptions(cmd)
		seq.ValidateSeq = false

		pattern := getFlagString(cmd, "pattern")
		replacement := []byte(getFlagString(cmd, "replacement"))

		// ------------------------------

		dbDir := getFlagString(cmd, "index")
		if dbDir == "" {
			checkError(fmt.Errorf("flag -d/--index needed"))
		}

		if pattern == "" {
			checkError(fmt.Errorf("flags -p (--pattern) needed"))
		}

		patternRegexp, err := regexp.Compile(pattern)
		checkError(err)

		// ---------------------------------------------------------------

		// new genomes.map file
		fileGenomeIndex := filepath.Join(dbDir, FileGenomeIndex+".new")
		fhGI, err := os.Create(fileGenomeIndex)
		if err != nil {
			checkError(fmt.Errorf("%s", err))
		}
		bw := bufio.NewWriter(fhGI)

		// ----------------------------------

		// old genomes.map file
		fileGenomeIndex0 := filepath.Join(dbDir, FileGenomeIndex)
		fh, err := os.Open(fileGenomeIndex0)
		if err != nil {
			checkError(fmt.Errorf("failed to read genome index mapping file: %s", err))
		}
		defer fh.Close()

		r := bufio.NewReader(fh)

		buf := make([]byte, 8)
		buf2 := make([]byte, 8)
		var n, lenID int
		// var batchIDAndRefID uint64
		var id2 []byte

		N := 0
		_n := 0

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

			// batchIDAndRefID = be.Uint64(buf)

			// replace id with the regular expression
			id2 = patternRegexp.ReplaceAll(id, replacement)
			if len(id2) == 0 {
				checkError(fmt.Errorf("genome ID is rename to empty: %s", id))
			}

			N++
			if !bytes.Equal(id, id2) {
				_n++
			}

			// output
			be.PutUint16(buf2[:2], uint16(len(id2)))
			bw.Write(buf2[:2])
			bw.Write(id2)
			bw.Write(buf)

		}

		bw.Flush()
		checkError(fhGI.Close())

		backFile := fileGenomeIndex0 + ".bak"

		var hasBackup bool
		hasBackup, err = pathutil.Exists(backFile)
		if err != nil {
			checkError(fmt.Errorf("failed to check backup file %s: %s", backFile, err))
		}
		if hasBackup {
			log.Infof("found the backup of genome index mapping file: %s", backFile)
		} else {
			err = os.Rename(fileGenomeIndex0, backFile)
			if err != nil {
				checkError(fmt.Errorf("failed to create backup file %s: %s", backFile, err))
			}
			log.Infof("created a backup of genome index mapping file: %s", backFile)
		}

		if _n == 0 {
			err = os.RemoveAll(fileGenomeIndex)
			if err != nil {
				checkError(fmt.Errorf("failed to delete tmp file %s: %s", fileGenomeIndex, err))
			}
		} else {
			err = os.Rename(fileGenomeIndex, fileGenomeIndex0)
			if err != nil {
				checkError(fmt.Errorf("failed to update %s: %s", fileGenomeIndex0, err))
			}
		}

		log.Infof("%d of %d genome IDs are changed", _n, N)
	},
}

func init() {
	utilsCmd.AddCommand(editGenomeIDsCmd)

	editGenomeIDsCmd.Flags().StringP("index", "d", "",
		formatFlagUsage(`Index directory created by "lexicmap index".`))

	editGenomeIDsCmd.Flags().StringP("pattern", "p", "",
		formatFlagUsage(`Search regular expression".`))

	editGenomeIDsCmd.Flags().StringP("replacement", "r", "",
		formatFlagUsage("Replacement. Supporting capture variables. "+
			" e.g. $1 represents the text of the first submatch. "+
			"ATTENTION: for *nix OS, use SINGLE quote NOT double quotes or "+
			`use the \ escape character.`))

	editGenomeIDsCmd.SetUsageTemplate(usageTemplate(""))
}
