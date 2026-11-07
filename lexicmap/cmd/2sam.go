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
	"bufio"
	"fmt"
	"io"
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
   - Output file of 'lexicmap search' with the flag -a/--all.
     Results with --show-avg-qual are also supported.
   - Each reference name must have one consistent length across all input files.
     Use -c/--concat-sgenome-sseqid to distinguish sequence IDs shared by different genomes.
   - Names must meet SAM character restrictions. Query names are limited to 254
     printable ASCII characters and cannot contain @.
   - Records for each query must be contiguous, including across input files,
     and must have the same query length. A group may span adjacent files.
   - Do not support STDIN.

Output:
   - Clipped regions in SEQ are represented as N's,
     as the input contains only the aligned portion of the sequences.
   - Different from SAM files produced by Minimap2,
     'X' (mismatch) in CIGAR is not converted to 'M' (match).
   - NM (edit distance) and AS (alignment score) fields are produced.
     AS uses match +2, mismatch -3, and gap -(5 + 2 * length), excluding clipping.
     MAPQ uses the same raw alignment scores.
     With alternative alignments, a nonpositive primary score gives MAPQ=0.
   - Custom integer length tags are emitted in this order:
     ZA:i: alignment length (alenHSP, without the optional quality suffix).
     ZQ:i: query alignment span (qend - qstart + 1).
     ZR:i: reference alignment span (send - sstart + 1).
     Both spans exclude clipping; ZQ includes insertions, ZR includes deletions.
   - Reverse-strand alignments use reverse-complemented SEQ and reversed
     CIGAR operation order, including soft clipping.
   - Records are unpaired, with RNEXT=*, PNEXT=0, and TLEN=0.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)

		outFile := getFlagString(cmd, "out-file")
		defer setupCommandLog(opt, outFile)()
		outputLog := opt.Verbose || opt.Log2File

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
		outfh, gw, w, err := outStream(outFile, strings.HasSuffix(outFile, ".gz") || opt.ForceGzipOutput, opt.CompressionLevel)
		checkError(err)
		defer func() {
			checkError(closeOutStream(outfh, gw, w))
		}()

		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)
		for _, file := range files {
			if file == "-" {
				checkError(fmt.Errorf("stdin not supported, as the input needs to be read twice"))
			}
		}

		// Grow for long records instead of eagerly allocating an unusually large limit.
		buf := make([]byte, min(bufferSize, 20<<20))
		var fh *xopen.Reader
		var line string
		var scanner *bufio.Scanner

		ncols := 24
		items := make([]string, ncols)

		var query, sgenome, sseqid, qcovHSP, alenHSP, pident, gaps, sstart, send, sstr, slen string
		var cigar, qseq string
		var qlen, qstart, qend string

		var headerLine bool
		var lineNo int
		var preQuery string

		// ---------------------------------------------------------------
		// 1st round: extract subject sequence id and length
		timeStart := time.Now()
		if outputLog {
			log.Info("round 1/2: extracting subject sequence IDs and lengths ...")
		}

		var key string
		refs := make([]string, 0, 2048)
		m := make(map[string]int, 1024)                // reference name -> validated sequence length
		seenQueries := make(map[string]struct{}, 1024) // completed/current query groups
		var previousQueryLength int
		for _, file := range files {
			fh, err = xopen.Ropen(file)
			checkError(err)

			headerLine = true
			lineNo = 0

			scanner = bufio.NewScanner(fh)
			scanner.Buffer(buf, int(bufferSize))
			for scanner.Scan() {
				lineNo++
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
					checkError(fmt.Errorf("the input has only %d columns (<%d), did you forget to add -a/--all for 'lexicmap search'?", len(items), ncols))
				}
				query = items[0]
				if !validSAMQueryName(query) {
					checkError(fmt.Errorf("%s: line %d: invalid SAM query name %q (expected 1-254 printable ASCII characters excluding @)", file, lineNo, query))
				}
				queryLength, _, _, coordErr := parseSearchResultCoordinates("query", items[1], items[12], items[13])
				if coordErr != nil {
					checkError(fmt.Errorf("%s: line %d: %w", file, lineNo, coordErr))
				}
				if query != preQuery {
					if _, exists := seenQueries[query]; exists {
						checkError(fmt.Errorf("%s: line %d: noncontiguous records for SAM query %q; group all records for each query before conversion", file, lineNo, query))
					}
					// Do not retain the entire alignment row for each query ID.
					seenQueries[strings.Clone(query)] = struct{}{}
					preQuery, previousQueryLength = query, queryLength
				} else if queryLength != previousQueryLength {
					checkError(fmt.Errorf("%s: line %d: inconsistent lengths for SAM query %q: %d and %d", file, lineNo, query, previousQueryLength, queryLength))
				}

				slen = items[17] // let the compiler to reduce boundary checking
				sgenome = items[3]
				sseqid = items[4]
				referenceLength, _, _, coordErr := parseSearchResultCoordinates("subject", slen, items[14], items[15])
				if coordErr != nil {
					checkError(fmt.Errorf("%s: line %d: %w", file, lineNo, coordErr))
				}

				// reference sequence name -> sequence length
				if concatSgenomeAndSseqid {
					key = sgenome + separater + sseqid
				} else {
					key = sseqid
				}
				if previous, ok := m[key]; ok {
					if previous != referenceLength {
						hint := ""
						if !concatSgenomeAndSseqid {
							hint = "; use -c/--concat-sgenome-sseqid for sequence IDs shared by different genomes"
						}
						checkError(fmt.Errorf("%s: line %d: inconsistent lengths for SAM reference %q: %d and %d%s", file, lineNo, key, previous, referenceLength, hint))
					}
				} else {
					if !validSAMReferenceName(key) {
						checkError(fmt.Errorf("%s: line %d: invalid SAM reference name %q", file, lineNo, key))
					}
					// These short fields outlive the row, whose aligned sequences
					// can be large. Detach both, including slen in concatenated mode.
					key = strings.Clone(key)
					refs = append(refs, key, strings.Clone(slen))
					m[key] = referenceLength
				}
			}

			checkError(scanner.Err())
			checkError(fh.Close())
		}
		seenQueries = nil // validation is complete; no query-name table is needed in round 2

		fmt.Fprintf(outfh, "@HD\tVN:1.6\tSO:unsorted\tGO:query\n")
		e := len(refs) - 2
		for i := 0; i <= e; i += 2 {
			fmt.Fprintf(outfh, "@SQ\tSN:%s\tLN:%s\n", refs[i], refs[i+1])
		}
		fmt.Fprintf(outfh, "@PG\tID:lexicmap\tPN:lexicmap\tVN:%s\n", VERSION)

		clear(m)
		clear(refs)
		refs = refs[:0]
		if outputLog {
			log.Infof("  elapsed time: %s", time.Since(timeStart))
		}

		// ---------------------------------------------------------------
		// 2nd round: convert to sam
		timeStart = time.Now()
		if outputLog {
			log.Info("round 2/2: converting to SAM format ...")
		}

		var _alenHSP, algnScore, editDistance, _gaps int
		var _qcovHSP, _pident float64
		var flag uint32
		aligns := make([]*SearchResultOfASequence2, 0, 1024)
		var _qlen int
		var clip5, clip3 string
		preQuery = ""
		for _, file := range files {
			fh, err = xopen.Ropen(file)
			checkError(err)

			headerLine = true
			lineNo = 0

			scanner = bufio.NewScanner(fh)
			scanner.Buffer(buf, int(bufferSize))
			for scanner.Scan() {
				lineNo++
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
					checkError(fmt.Errorf("the input has only %d columns (<%d), did you forget to add -a/--all for 'lexicmap search'?", len(items), ncols))
				}

				qseq = items[21] // let the compiler to reduce boundary checking
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
				cigar = items[20]

				// MAPQ needs finite percentages and a valid gap count; parse
				// failures must not silently lower or inflate confidence.
				_qcovHSP, err = strconv.ParseFloat(qcovHSP, 64)
				if err != nil || math.IsNaN(_qcovHSP) || math.IsInf(_qcovHSP, 0) || _qcovHSP < 0 || _qcovHSP > 100 {
					checkError(fmt.Errorf("%s: line %d: invalid qcovHSP %q (expected a finite percentage in [0, 100])", file, lineNo, qcovHSP))
				}
				// --show-avg-qual decorates the length as length:quality.
				length, _, _ := strings.Cut(alenHSP, ":")
				_alenHSP, err = strconv.Atoi(length)
				if err != nil || _alenHSP <= 0 {
					checkError(fmt.Errorf("%s: line %d: invalid alenHSP %q", file, lineNo, alenHSP))
				}
				_pident, err = strconv.ParseFloat(pident, 64)
				if err != nil || math.IsNaN(_pident) || math.IsInf(_pident, 0) || _pident < 0 || _pident > 100 {
					checkError(fmt.Errorf("%s: line %d: invalid pident %q (expected a finite percentage in [0, 100])", file, lineNo, pident))
				}
				_gaps, err = strconv.Atoi(gaps)
				if err != nil || _gaps < 0 || _gaps > _alenHSP {
					checkError(fmt.Errorf("%s: line %d: invalid gaps %q (expected an integer in [0, %d])", file, lineNo, gaps, _alenHSP))
				}
				var _qstart, _qend int
				_qlen, _qstart, _qend, err = parseSearchResultCoordinates("query", qlen, qstart, qend)
				if err != nil {
					checkError(fmt.Errorf("%s: line %d: %w", file, lineNo, err))
				}
				_, start, end, coordErr := parseSearchResultCoordinates("subject", slen, sstart, send)
				if coordErr != nil {
					checkError(fmt.Errorf("%s: line %d: %w", file, lineNo, coordErr))
				}
				if sstr != "+" && sstr != "-" {
					checkError(fmt.Errorf("%s: line %d: invalid strand %q (expected + or -)", file, lineNo, sstr))
				}
				querySpan, referenceSpan := _qend-_qstart+1, end-start+1
				editDistance, algnScore, err = alignmentStatsFromCIGAR(cigar, querySpan, referenceSpan, _alenHSP)
				if err != nil {
					checkError(fmt.Errorf("%s: line %d: %w", file, lineNo, err))
				}
				if n := len(qseq) - strings.Count(qseq, "-"); n != querySpan {
					checkError(fmt.Errorf("%s: line %d: query sequence length %d does not match query span %d", file, lineNo, n, querySpan))
				}

				flag = 0
				if sstr == "-" {
					flag |= 0x10 // SEQ being reverse complemented
				}

				r := poolSearchResultOfASequence2.Get().(*SearchResultOfASequence2)

				// data for later use
				r.score = algnScore
				r.qcovHSP = _qcovHSP
				r.pident = _pident
				r.gaps = float64(_gaps)
				r.alenHSP = _alenHSP

				// sam data
				r.FLAG = flag
				if concatSgenomeAndSseqid {
					r.RNAME = sgenome + separater + sseqid
				} else {
					r.RNAME = strings.Clone(sseqid)
				}
				// These fields survive until the query is complete. Detach them
				// from the row's large aligned sequences, including secondary hits.
				r.POS = strings.Clone(sstart)
				r.MAPQ = 0 // compute later
				clip5, clip3 = "", ""
				if _qstart > 1 {
					clip5 = fmt.Sprintf("%dS", _qstart-1)
				}
				if _qend < _qlen {
					clip3 = fmt.Sprintf("%dS", _qlen-_qend)
				}
				if clip5 != "" || clip3 != "" {
					cigar = clip5 + cigar + clip3
				}
				r.CIGAR = strings.Clone(cigar)
				r.RNEXT = "*"
				r.PNEXT = "0"
				r.TLEN = 0 // SAM single-segment templates have no paired template length.
				// Keep the former TLEN value in a custom tag with explicit span semantics.
				r.ReferenceSpan = referenceSpan
				r.QuerySpan = querySpan
				if len(aligns) == 0 || preQuery != query {
					r.SEQ = strings.Repeat("N", _qstart-1) + strings.ReplaceAll(qseq, "-", "") + strings.Repeat("N", _qlen-_qend)
				} else {
					r.SEQ = "*"
				}
				// Search TSV follows the query strand; SAM follows the reference strand.
				// Reverse whole operations (including clipping), never their digits.
				if flag&0x10 != 0 {
					r.CIGAR = reverseSAMCIGAR(r.CIGAR)
					if r.SEQ != "*" {
						r.SEQ = string(RC([]byte(r.SEQ)))
					}
				}
				r.QUAL = "*"
				r.NM = editDistance
				r.AS = algnScore

				if preQuery != query && len(aligns) > 0 {
					checkError(writeSAMQuery(outfh, preQuery, aligns))
					clear(aligns)
					aligns = aligns[:0]
				}

				aligns = append(aligns, r)
				preQuery = query
			}
			checkError(scanner.Err())
			checkError(fh.Close())
		}
		// An input file boundary does not end a query that continues in the next file.
		if len(aligns) > 0 {
			checkError(writeSAMQuery(outfh, preQuery, aligns))
			clear(aligns)
			aligns = aligns[:0]
		}

		if outputLog {
			log.Infof("  elapsed time: %s", time.Since(timeStart))
		}
	},
}

// validSAMQueryName checks the SAM v1.6 QNAME length and ASCII restrictions.
func validSAMQueryName(name string) bool {
	if len(name) == 0 || len(name) > 254 {
		return false
	}
	for i := 0; i < len(name); i++ {
		if name[i] < '!' || name[i] > '~' || name[i] == '@' {
			return false
		}
	}
	return true
}

// validSAMReferenceName checks SN/RNAME delimiters and reserved leading symbols.
// See https://samtools.github.io/hts-specs/SAMv1.pdf, section 1.2.1.
func validSAMReferenceName(name string) bool {
	if len(name) == 0 || name[0] == '*' || name[0] == '=' {
		return false
	}
	for i := 0; i < len(name); i++ {
		if name[i] < '!' || name[i] > '~' {
			return false
		}
		switch name[i] {
		case '\\', ',', '"', '\'', '`', '(', ')', '[', ']', '{', '}', '<', '>':
			return false
		}
	}
	return true
}

// reverseSAMCIGAR changes from query-forward to reference-forward operation order.
func reverseSAMCIGAR(cigar string) string {
	var reversed strings.Builder
	reversed.Grow(len(cigar))
	for end := len(cigar); end > 0; {
		start := end - 1 // operation letter
		for start > 0 && cigar[start-1] >= '0' && cigar[start-1] <= '9' {
			start--
		}
		reversed.WriteString(cigar[start:end])
		end = start
	}
	return reversed.String()
}

// writeSAMQuery applies the same MAPQ and secondary flags at query boundaries
// and EOF. Keep the existing primary alignment, and find its strongest competitor
// directly rather than translating indexes from a subslice.
func writeSAMQuery(w io.Writer, query string, aligns []*SearchResultOfASequence2) error {
	if len(aligns) == 0 {
		return nil
	}
	a := aligns[0]
	if len(aligns) == 1 {
		a.MAPQ = 60
	} else {
		maxScore := aligns[1].score
		for _, other := range aligns[1:] {
			other.FLAG |= 0x100
			maxScore = max(maxScore, other.score)
		}
		// A nonpositive denominator can invert the score difference or divide
		// by zero. Confidence also requires beating every alternative.
		a.MAPQ = 0
		if a.score > 0 && a.score > maxScore {
			mapq := 40 * float64(a.score-maxScore) / float64(a.score)
			mapq *= a.qcovHSP / 100
			mapq *= (a.pident / 100) * (1 - a.gaps/float64(a.alenHSP))
			a.MAPQ = uint32(min(60, max(0, int(mapq))))
		}
	}
	for _, a := range aligns {
		if _, err := fmt.Fprintf(w, "%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%d\t%s\t%s\tNM:i:%d\tAS:i:%d\tZA:i:%d\tZQ:i:%d\tZR:i:%d\n",
			query, a.FLAG, a.RNAME, a.POS, a.MAPQ, a.CIGAR, a.RNEXT, a.PNEXT,
			a.TLEN, a.SEQ, a.QUAL, a.NM, a.AS, a.alenHSP, a.QuerySpan, a.ReferenceSpan); err != nil {
			return err
		}
		recycleSearchResultOfASequence2(a)
	}
	return nil
}

// SearchResultOfASequence2 holds one SAM record and the inputs to its MAPQ.
type SearchResultOfASequence2 struct {
	score   int
	qcovHSP float64
	pident  float64
	gaps    float64
	alenHSP int // Alignment length used by MAPQ and the ZA tag.

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

	QuerySpan     int // ZQ: inclusive query span, excluding clipping.
	ReferenceSpan int // ZR: inclusive reference span, independent of strand and query clipping.
}

var poolSearchResultOfASequence2 = &sync.Pool{New: func() interface{} {
	return &SearchResultOfASequence2{}
}}

func recycleSearchResultOfASequence2(r *SearchResultOfASequence2) {
	*r = SearchResultOfASequence2{}
	poolSearchResultOfASequence2.Put(r)
}

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

// alignmentStatsFromCIGAR uses search's +2/-3 and gap-open 5 + extension 2*k
// scoring. Search emits M/=, X, I, and D; clipping is added by the converter.
// Rounded bitscores cannot recover the raw score, especially odd scores.
func alignmentStatsFromCIGAR(cigar string, querySpan, referenceSpan, alignmentLen int) (distance, score int, err error) {
	// A single-base gap costs seven, the largest absolute cost per column.
	if cigar == "" || alignmentLen <= 0 || alignmentLen > int(^uint(0)>>1)/7 {
		return 0, 0, fmt.Errorf("invalid CIGAR %q or alignment length %d", cigar, alignmentLen)
	}
	var query, reference, aligned int
	for i := 0; i < len(cigar); {
		start := i
		for i < len(cigar) && cigar[i] >= '0' && cigar[i] <= '9' {
			i++
		}
		if start == i || i == len(cigar) {
			return 0, 0, fmt.Errorf("invalid CIGAR %q", cigar)
		}
		n, parseErr := strconv.Atoi(cigar[start:i])
		if parseErr != nil || n <= 0 || n > alignmentLen-aligned {
			return 0, 0, fmt.Errorf("invalid CIGAR %q for alignment length %d", cigar, alignmentLen)
		}
		switch cigar[i] {
		case 'M', '=':
			query += n
			reference += n
			score += 2 * n
		case 'X':
			query += n
			reference += n
			distance += n
			score -= 3 * n
		case 'I', 'D':
			if cigar[i] == 'I' {
				query += n
			} else {
				reference += n
			}
			distance += n
			score -= 5 + 2*n
		default:
			return 0, 0, fmt.Errorf("invalid CIGAR %q: unsupported search operation %q", cigar, cigar[i])
		}
		aligned += n
		i++
	}
	if query != querySpan || reference != referenceSpan || aligned != alignmentLen {
		return 0, 0, fmt.Errorf("CIGAR spans (%d, %d, %d) do not match query, subject, and alignment lengths (%d, %d, %d)", query, reference, aligned, querySpan, referenceSpan, alignmentLen)
	}
	return distance, score, nil
}
