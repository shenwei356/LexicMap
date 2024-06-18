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
	"strconv"
	"strings"

	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
)

var toBlastCmd = &cobra.Command{
	Use:   "2blast",
	Short: "Convert the default search output to blast-style format",
	Long: `Convert the default search output to blast-style format

Input:
   - Output of 'lexicmap search' with the flag -a/--all.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)

		outFile := getFlagString(cmd, "out-file")

		bufferSizeS := getFlagString(cmd, "buffer-size")
		if bufferSizeS == "" {
			checkError(fmt.Errorf("value of buffer size. supported unit: K, M, G"))
		}

		bufferSize, err := ParseByteSize(bufferSizeS)
		if err != nil {
			checkError(fmt.Errorf("invalid value of buffer size. supported unit: K, M, G"))
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

		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)

		buf := make([]byte, bufferSize)
		var fh *xopen.Reader
		var line string
		var scanner *bufio.Scanner

		ncols := 21
		items := make([]string, ncols)

		var query, qlen, hits, sgenome, sseqid, qcovGnm, hsp, qcovHSP, alenHSP, pident, gaps, qstart, qend, sstart, send, sstr, slen string
		var cigar, qseq, sseq, align string

		var headerLine bool
		var iGenome, iSeq int
		var preQuery, preGenome, preSeq string
		var rows, i, j, end int
		var _qstart, _qend, _sstart, _send int
		var _qstart2, _qend2, _sstart2, _send2 int
		var posW int
		var fA, fQ, fT string
		var _strand string
		var q, t string
		var rc bool

		for _, file := range files {
			fh, err = xopen.Ropen(file)
			checkError(err)

			headerLine = true

			scanner = bufio.NewScanner(fh)
			scanner.Buffer(buf, int(bufferSize))
			for scanner.Scan() {
				line = strings.TrimRight(scanner.Text(), "\r\n")
				if line == "" {
					continue
				}
				if headerLine {
					headerLine = false
					continue
				}

				stringSplitNByByte(line, '\t', ncols, &items)
				if len(items) < ncols {
					checkError(fmt.Errorf("the input has only %d columns, did you forgot to add -a/--all for 'lexicmap search'?", len(items)))
				}

				query = items[0]
				qlen = items[1]
				hits = items[2]
				sgenome = items[3]
				sseqid = items[4]
				qcovGnm = items[5]
				hsp = items[6]
				qcovHSP = items[7]
				alenHSP = items[8]
				pident = items[9]
				gaps = items[10]
				qstart = items[11]
				qend = items[12]
				sstart = items[13]
				send = items[14]
				sstr = items[15]
				slen = items[16]
				cigar = items[17]
				qseq = items[18]
				sseq = items[19]
				align = items[20]

				_qstart, _ = strconv.Atoi(qstart)
				_qend, _ = strconv.Atoi(qend)
				_sstart, _ = strconv.Atoi(sstart)
				_send, _ = strconv.Atoi(send)

				_ = _qend

				if sstr == "+" {
					_strand = "Plus"
					rc = false
				} else {
					_strand = "Minus"
					rc = true
				}

				posW = max(len(qend), len(send))
				fQ = fmt.Sprintf("Query  %%-%dd  %%s  %%d\n", posW)
				fA = fmt.Sprintf("       %%%ds  %%s\n", posW)
				fT = fmt.Sprintf("Sbjct  %%-%dd  %%s  %%d\n", posW)

				_ = cigar

				if preQuery != query {
					iGenome = 0
					fmt.Fprintf(outfh, "Query = %s\nLength = %s\n\n", query, qlen)
				}
				if preGenome != sgenome {
					iGenome++
					fmt.Fprintf(outfh, "[Subject genome #%d/%s] = %s\nQuery coverage per genome = %s%%\n\n", iGenome, hits, sgenome, qcovGnm)
				}
				if preSeq != sseqid {
					iSeq = 1
					fmt.Fprintf(outfh, ">%s\nLength = %s\n\n", sseqid, slen)
				}

				fmt.Fprintf(outfh, " HSP #%s\n", hsp)
				fmt.Fprintf(outfh, " Query coverage per seq = %s%%, Aligned length = %s, Identities = %s%%, Gaps = %s\n",
					qcovHSP, alenHSP, pident, gaps)
				fmt.Fprintf(outfh, " Query range = %s-%s, Subject range = %s-%s, Strand = Plus/%s\n\n",
					qstart, qend, sstart, send, _strand)

				rows = (len(qseq) + 59) / 60

				_qstart2 = _qstart
				if rc {
					_sstart2 = _send
				} else {
					_sstart2 = _sstart
				}
				for i = 0; i < rows; i++ {
					j = i * 60
					if i < rows-1 {
						end = j + 60
					} else {
						end = len(qseq)
					}
					q, t = qseq[j:end], sseq[j:end]

					_qend2 = _qstart2 + len(q) - ngaps(q) - 1
					if rc {
						_send2 = _sstart2 - (len(t) - ngaps(t)) + 1
					} else {
						_send2 = _sstart2 + len(t) - ngaps(t) - 1
					}

					fmt.Fprintf(outfh, fQ, _qstart2, q, _qend2)
					fmt.Fprintf(outfh, fA, " ", align[j:end])
					fmt.Fprintf(outfh, fT, _sstart2, t, _send2)
					outfh.WriteByte('\n')

					_qstart2 = _qend2 + 1
					if rc {
						_sstart2 = _send2 - 1
					} else {
						_sstart2 = _send2 + 1
					}
				}
				outfh.WriteByte('\n')

				iSeq++
				preQuery = query
				preGenome = sgenome
				preSeq = sseqid

			}

			checkError(scanner.Err())
			checkError(fh.Close())
		}
	},
}

func init() {
	utilsCmd.AddCommand(toBlastCmd)

	toBlastCmd.Flags().StringP("out-file", "o", "-",
		formatFlagUsage(`Out file, supports and recommends a ".gz" suffix ("-" for stdout).`))

	toBlastCmd.Flags().StringP("buffer-size", "b", "20M",
		formatFlagUsage(`Size of buffer, supported unit: K, M, G. You need increase the value when "bufio.Scanner: token too long" error reported`))

	toBlastCmd.SetUsageTemplate(usageTemplate(""))
}

func ngaps(s string) int {
	if len(s) == 0 {
		return 0
	}
	n := 0
	for _, b := range s {
		if b == '-' {
			n++
		}
	}
	return n
}