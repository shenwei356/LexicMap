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
	"bytes"
	"fmt"
	"io"
	"regexp"
	"strconv"
	"strings"

	"github.com/pkg/errors"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/spf13/cobra"
)

var subseqCmd = &cobra.Command{
	Use:   "subseq",
	Short: "extract subsequence via position and strand",
	Long: `extract subsequence via position and strand

Function:
  For reference genomes with multiple sequences, the sequences were
  concatenated to a single sequence with intervals of (k-1) N's.
  This command is for extracting subsequences via the pseudo positions
  outputed by 'utils kmers/kmer-locations'.

Attention:
  1. You need to set the same -k/--kmer value used in 'lexicmap index'.
  2. You need to set the same -B/--seq-name-filter values used in
     'lexicmap index' to filter out unwanted sequences like plasmids.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		seq.ValidateSeq = false

		// ------------------------------

		var reRegion = regexp.MustCompile(`\-?\d+:\-?\d+`)

		region := getFlagString(cmd, "region")
		if region == "" {
			checkError(fmt.Errorf("flag -r/--region needed"))
		}
		revcom := getFlagBool(cmd, "revcom")

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
		if start == 0 || end == 0 {
			checkError(fmt.Errorf("both start and end should not be 0"))
		}
		if start < 0 && end > 0 {
			checkError(fmt.Errorf("when start < 0, end should not > 0"))
		}

		// ----------

		k := getFlagPositiveInt(cmd, "kmer")
		if k < 4 || k > 32 {
			checkError(fmt.Errorf("the value of flag -k/--kmer should be in range of [4, 32]"))
		}
		nnn := bytes.Repeat([]byte{'N'}, k-1)

		outFile := getFlagString(cmd, "out-file")

		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)
		if len(files) > 1 {
			checkError(fmt.Errorf("no more than one file should be given"))
		}

		reSeqNameStrs := getFlagStringSlice(cmd, "seq-name-filter")
		reSeqNames := make([]*regexp.Regexp, 0, len(reSeqNameStrs))
		for _, kw := range reSeqNameStrs {
			if !reIgnoreCase.MatchString(kw) {
				kw = reIgnoreCaseStr + kw
			}
			re, err := regexp.Compile(kw)
			if err != nil {
				checkError(errors.Wrapf(err, "failed to parse regular expression for matching sequence header: %s", kw))
			}
			reSeqNames = append(reSeqNames, re)
		}
		filterNames := len(reSeqNames) > 0

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

		// --------------------------------

		var record *fastx.Record
		var fastxReader *fastx.Reader

		file := files[0]
		fastxReader, err = fastx.NewReader(nil, file, "")
		checkError(err)

		var ignoreSeq bool
		var re *regexp.Regexp
		var s *seq.Seq
		var i int = 0
		for {
			record, err = fastxReader.Read()
			if err != nil {
				if err == io.EOF {
					break
				}
				checkError(err)
				break
			}

			// filter out sequences shorter than k
			if len(record.Seq.Seq) < k {
				continue
			}

			// filter out sequences with names in the blast list
			if filterNames {
				ignoreSeq = false
				for _, re = range reSeqNames {
					if re.Match(record.Name) {
						ignoreSeq = true
						break
					}
				}
				if ignoreSeq {
					continue
				}
			}
			if i == 0 {
				s = record.Seq.Clone()
			} else {
				s.Seq = append(s.Seq, nnn...)
				s.Seq = append(s.Seq, record.Seq.Seq...)
			}

			i++
		}

		if s == nil {
			log.Warningf("skipping %s: no valid sequences", file)
			return
		}

		s.SubSeqInplace(start, end)
		if revcom {
			s.RevComInplace()
		}
		outfh.Write(s.Seq)
		outfh.WriteByte('\n')
	},
}

func init() {
	utilsCmd.AddCommand(subseqCmd)

	subseqCmd.Flags().StringP("out-file", "o", "-",
		formatFlagUsage(`Out file, supports and recommends a ".gz" suffix ("-" for stdout).`))

	subseqCmd.Flags().IntP("kmer", "k", 31,
		formatFlagUsage(`Maximum k-mer size. K needs to be <= 32.`))

	subseqCmd.Flags().StringSliceP("seq-name-filter", "B", []string{},
		formatFlagUsage(`List of regular expressions for filtering out sequences by header/name, case ignored.`))

	subseqCmd.Flags().StringP("region", "r", "",
		formatFlagUsage(`region of the subsequence`))

	subseqCmd.Flags().BoolP("revcom", "R", false,
		formatFlagUsage("extract subsequence on the negative strand"))

	subseqCmd.SetUsageTemplate(usageTemplate(""))
}
