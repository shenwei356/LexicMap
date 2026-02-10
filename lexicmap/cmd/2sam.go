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
	"fmt"
	"math"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
)

var toSamCmd = &cobra.Command{
	Use:   "2sam",
	Short: "Convert the default search output to SAM format",
	Long: `Convert the default search output to SAM format

Input:
   - Output of 'lexicmap search' with the flag -a/--all.
   - Do not support STDIN.

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

		concatSgenomeAndSseqid := getFlagBool(cmd, "concat-sgenome-sseqid")
		separater := getFlagString(cmd, "separater")

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

		ncols := 24
		items := make([]string, ncols)

		var query, qlen, sgenome, sseqid, qcovHSP, alenHSP, pident, gaps, qstart, qend, sstart, send, sstr, slen, bitscore string
		var cigar, qseq string

		var headerLine bool
		var preQuery string
		var _sstart, _send int

		// ---------------------------------------------------------------
		// 1st round: extract subject sequence id and length
		timeStart := time.Now()
		if opt.Verbose {
			log.Info("round 1/2: extracting subject sequence IDs and lengths ...")
		}

		var key string
		var ok bool
		refs := make([]string, 0, 2048)
		m := make(map[string]struct{}, 1024)
		for _, file := range files {
			if file == "-" {
				checkError(fmt.Errorf("stdin not supported"))
			}
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
					checkError(fmt.Errorf("the input has only %d columns (<%d), did you forget to add -a/--all for 'lexicmap search'?", ncols, len(items)))
				}

				slen = items[17] // let the compiler to reduce boundary checking
				sgenome = items[3]
				sseqid = items[4]

				// reference sequence name -> sequence length
				if concatSgenomeAndSseqid {
					key = sgenome + separater + sseqid
				} else {
					key = sseqid
				}
				if _, ok = m[key]; !ok {
					refs = append(refs, key)
					refs = append(refs, slen)

					m[key] = struct{}{}
				}
			}

			checkError(scanner.Err())
			checkError(fh.Close())
		}

		fmt.Fprintf(outfh, "@HD\tVN:1.6\tSO:unsorted\tGO:query\n")
		e := len(refs) - 2
		for i := 0; i <= e; i += 2 {
			fmt.Fprintf(outfh, "@SQ\tSN:%s\tLN:%s\n", refs[i], refs[i+1])
		}
		fmt.Fprintf(outfh, "@PG\tID:lexicmap\tPN:lexicmap\tVN:%s\n", VERSION)

		clear(m)
		refs = refs[:0]
		if opt.Verbose {
			log.Infof("elapsed time: %s", time.Since(timeStart))
		}

		// ---------------------------------------------------------------
		// 2nd round: convert to sam
		timeStart = time.Now()
		if opt.Verbose {
			log.Info("round 2/2: converting to SAM format ...")
		}

		// from blastn_values_2_3 in ncbi-blast-2.15.0+-src/c++/src/algo/blast/core/blast_stat.c
		lambda := 0.625
		lnK := math.Log(0.41)
		// bitScore := (lambda*float64(_score) - lnK) / math.Ln2

		var _bitscore, _alenHSP, algnScore, _gaps, _qlen, _qstart, _qend int
		var _qcovHSP, _pident float64
		var flag uint32
		aligns := make([]*SearchResultOfASequence2, 0, 1024)
		var a, b *SearchResultOfASequence2
		var i, maxScore, maxI int
		var mapq float64
		var clip5, clip3 string
		for _, file := range files {
			fh, err = xopen.Ropen(file)
			checkError(err)

			headerLine = true
			preQuery = "shenwei356"
			aligns = aligns[:0]
			clear(m)

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
					checkError(fmt.Errorf("the input has only %d columns (<%d), did you forget to add -a/--all for 'lexicmap search'?", ncols, len(items)))
				}

				cigar = items[20] // let the compiler to reduce boundary checking
				query = items[0]
				qlen = items[1]
				sgenome = items[3]
				sseqid = items[4]
				qcovHSP = items[8]
				alenHSP = items[9]
				pident = items[10]
				gaps = items[11]
				qstart = items[12]
				qend = items[13]
				sstart = items[14]
				send = items[15]
				sstr = items[16]
				slen = items[17]
				bitscore = items[19]
				qseq = items[21]

				_qlen, _ = strconv.Atoi(qlen)
				_qcovHSP, _ = strconv.ParseFloat(qcovHSP, 64)
				_alenHSP, _ = strconv.Atoi(alenHSP)
				_pident, _ = strconv.ParseFloat(pident, 64)
				_gaps, _ = strconv.Atoi(gaps)
				_qstart, _ = strconv.Atoi(qstart)
				_qend, _ = strconv.Atoi(qend)
				_sstart, _ = strconv.Atoi(sstart)
				_send, _ = strconv.Atoi(send)
				_bitscore, _ = strconv.Atoi(bitscore)
				algnScore = int((float64(_bitscore)*math.Ln2 + lnK) / lambda)

				flag = 0
				if sstr == "-" {
					flag |= 0x10 // SEQ being reverse complemented
				}

				r := poolSearchResultOfASequence2.Get().(*SearchResultOfASequence2)
				r.sgenome = sgenome
				r.score = algnScore
				r.qcovHSP = _qcovHSP
				r.pident = _pident
				r.gaps = float64(_gaps)
				r.alenHSP = float64(_alenHSP)

				r.FLAG = flag
				if concatSgenomeAndSseqid {
					r.RNAME = sgenome + separater + sseqid
				} else {
					r.RNAME = sseqid
				}
				r.POS = sstart
				r.MAPQ = 0 // compute later
				if _qstart > 1 {
					clip5 = fmt.Sprintf("%dS", _qstart-1)
				}
				if _qend < _qlen {
					clip3 = fmt.Sprintf("%dS", _qlen-_qend)
				}
				if clip5 != "" || clip3 != "" {
					// cigar = clip5 + cigar + clip3
				}
				r.CIGAR = cigar // todo: append S
				r.RNEXT = "*"
				r.PNEXT = "0"
				r.TLEN = _send - _sstart + 1
				r.SEQ = strings.ReplaceAll(qseq, "-", "")
				r.QUAL = "*"
				r.NM = int(float64(_alenHSP) * (1 - _pident/100))
				r.AS = algnScore

				if preQuery != query && len(aligns) > 0 {
					if len(aligns) == 1 {
						aligns[0].MAPQ = 60
					} else {
						a = aligns[0]
						b = aligns[1]
						maxScore, maxI = b.score, 1
						b.FLAG |= 0x100 // secondary alignment
						b.SEQ = "*"
						for i, b = range aligns[2:] {
							b.FLAG |= 0x100 // secondary alignment
							b.SEQ = "*"
							if b.score > maxScore {
								maxScore, maxI = b.score, i
							}
						}
						b = aligns[maxI]

						mapq = 40 * float64(a.score-b.score) / float64(a.score)
						mapq *= a.qcovHSP / 100                           // cov_factor
						mapq *= (a.pident / 100) * (1 - a.gaps/a.alenHSP) // qual_factor

						a.MAPQ = uint32(min(60, max(0, int(mapq))))
					}

					for _, a = range aligns {
						fmt.Fprintf(outfh, "%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%d\t%s\t%s\tNM:i:%d\tAS:i:%d\n",
							preQuery, a.FLAG, a.RNAME, a.POS, a.MAPQ, a.CIGAR, a.RNEXT, a.PNEXT,
							a.TLEN, a.SEQ, a.QUAL, a.NM, a.AS)

						poolSearchResultOfASequence2.Put(a)
					}

					aligns = aligns[:0]
					clear(m)
				}

				aligns = append(aligns, r)
				preQuery = query
			}
			if len(aligns) > 0 {
				if len(aligns) == 1 {
					aligns[0].MAPQ = 60
				} else {
					a = aligns[0]
					b = aligns[1]
					maxScore, maxI = b.score, 1
					b.FLAG |= 0x100 // secondary alignment
					b.SEQ = "*"
					for i, b = range aligns[2:] {
						b.FLAG |= 0x100 // secondary alignment
						b.SEQ = "*"
						if b.score > maxScore {
							maxScore, maxI = b.score, i
						}
					}
					b = aligns[maxI]

					mapq = 40 * float64(a.score-b.score) / float64(a.score)
					mapq *= a.qcovHSP / 100                           // cov_factor
					mapq *= (a.pident / 100) * (1 - a.gaps/a.alenHSP) // qual_factor

					a.MAPQ = uint32(min(60, max(0, int(mapq))))
				}

				for _, a = range aligns {
					fmt.Fprintf(outfh, "%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%d\t%s\t%s\tNM:i:%d\tAS:i:%d\n",
						preQuery, a.FLAG, a.RNAME, a.POS, a.MAPQ, a.CIGAR, a.RNEXT, a.PNEXT,
						a.TLEN, a.SEQ, a.QUAL, a.NM, a.AS)

					poolSearchResultOfASequence2.Put(a)
				}
			}

			checkError(scanner.Err())
			checkError(fh.Close())
		}
	},
}

type SearchResultOfASequence2 struct {
	sgenome string
	score   int
	qcovHSP float64
	pident  float64
	gaps    float64
	alenHSP float64

	FLAG  uint32
	RNAME string
	POS   string
	MAPQ  uint32
	CIGAR string
	RNEXT string
	PNEXT string
	TLEN  int
	SEQ   string
	QUAL  string
	NM    int
	AS    int
}

var poolSearchResultOfASequence2 = &sync.Pool{New: func() interface{} {
	return &SearchResultOfASequence2{}
}}

func init() {
	utilsCmd.AddCommand(toSamCmd)

	toSamCmd.Flags().StringP("out-file", "o", "-",
		formatFlagUsage(`Out file, supports and recommends a ".gz" suffix ("-" for stdout).`))

	toSamCmd.Flags().StringP("buffer-size", "b", "20M",
		formatFlagUsage(`Size of buffer, supported unit: K, M, G. You need increase the value when "bufio.Scanner: token too long" error reported`))

	toSamCmd.Flags().BoolP("concat-sgenome-sseqid", "c", false,
		formatFlagUsage(`Concatenate sgenome and sseqid to make sure the reference sequence names are distinct.`))
	toSamCmd.Flags().StringP("separater", "s", "~",
		formatFlagUsage(`Separater between sgenome and sseqid`))

	toSamCmd.SetUsageTemplate(usageTemplate(""))
}
