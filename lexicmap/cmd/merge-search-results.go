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
	"cmp"
	"container/heap"
	"errors"
	"fmt"
	"io"
	"math"
	"slices"
	"strconv"
	"strings"
	"sync"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
)

var mergeCmd = &cobra.Command{
	Use:   "merge-search-results",
	Short: "Merge a query's search results from multiple indexes",
	Long: `Merge a query's search results from multiple indexes

Attention:
  1. These search results should come from the same ONE query.
     If not, please specify one query with the flag -q/--query.
     Query selection also applies when only one input file is given.
     Merging without -q rejects mixed query IDs, including within one subject group.
     All records for the selected query must have the same positive query length.
     Hit counts must be positive integers, identities finite values in [0,100],
     and bit scores integers (negative scores are supported).
  2. We assume that genome IDs are distinct across all indexes.
  3. One or more input files are accepted, via positional parameters
     and/or a file list via the flag -X/--infile-list.
  4. Both the default 20- and 24-column formats are supported,
     and the first retained result selects the output format.
     In 24-column output, missing alignment fields are left empty.
     In 20-column output, alignment fields are omitted.
  5. Flag -n/--top-n-genomes ranks genomes by their highest alignment
     bitscore * pident, retaining all ties at the cutoff score.
     Unlike 'lexicmap search -n', this does not use chaining scores.
     Ranking uses the rounded TSV values, which can change the original search order.
     Only the top N genomes and cutoff ties are buffered for final output.
  6. Reader buffers start at up to 64 KiB and grow for long alignment rows up to
     -b/--buffer-size (default 20M). Increase -b if a row exceeds this limit.
  7. Empty inputs are skipped with a warning naming each file (unless --quiet).
     If all inputs are empty, the command succeeds without writing results.
     A completed search writes a header even when no matches are found.
     Header-only inputs are accepted.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		seq.ValidateSeq = false

		outFile := getFlagString(cmd, "out-file")
		defer setupCommandLog(opt, outFile)()
		outputLog := opt.Verbose || opt.Log2File

		// ------------------------------

		query := getFlagString(cmd, "query")
		topN := getFlagNonNegativeInt(cmd, "top-n-genomes")

		bufferSizeS := getFlagString(cmd, "buffer-size")
		if bufferSizeS == "" {
			checkError(fmt.Errorf("value of buffer size. supported unit: K, M, G"))
		}
		bufferSize, err := ParseByteSize(bufferSizeS)
		if err != nil {
			checkError(fmt.Errorf("invalid value of buffer size. supported unit: K, M, G"))
		}

		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)

		if outputLog {
			log.Infof("merging search results from %d files...", len(files))
			if query == "" {
				log.Infof("  query (--query): not specified (input must contain one query)")
			} else {
				log.Infof("  query (--query): %s", query)
			}
			if topN == 0 {
				log.Infof("  top-n-genomes: 0 (all genome hits)")
			} else {
				log.Infof("  top-n-genomes: %d (alignment bitscore * pident, including cutoff ties)", topN)
			}
			if outFile == "-" {
				log.Infof("  output file: - (stdout)")
			} else {
				log.Infof("  output file: %s", outFile)
			}
		}

		// output file handler
		outfh, gw, w, err := outStream(outFile, strings.HasSuffix(outFile, ".gz") || opt.ForceGzipOutput, opt.CompressionLevel)
		checkError(err)
		defer func() {
			checkError(closeOutStream(outfh, gw, w))
		}()

		var retained []*SearchResultOfAGenome
		selected := topSearchResults{SearchResultsHeap: SearchResultsHeap{entries: &retained}, limit: topN}

		// readers

		readers := make(map[int]*SearchResultReader, len(files))
		for i, file := range files {
			reader, err := NewSearchResultReader(file, query, bufferSize)
			// Only clean empty streams can be skipped; read errors remain fatal.
			if errors.Is(err, xopen.ErrNoContent) {
				log.Warningf("skipping %s", err)
				continue
			}
			checkError(err)
			readers[i] = reader
		}
		// Preserve a header-only single input after it passes through the reader.
		// Next synchronizes access to its header by draining the result channel.
		var singleInput *SearchResultReader
		if len(files) == 1 && topN == 0 && query == "" {
			singleInput = readers[0]
		}

		entries := make([]*SearchResultOfAGenome, 0, len(files))
		results := SearchResultsHeap{entries: &entries}

		fillBuffer := func() error {
			var rGnm *SearchResultOfAGenome
			for i, reader := range readers {
				reader = readers[i]

				rGnm = reader.Next()
				if rGnm == nil {
					delete(readers, i)
					continue
				}

				rGnm.idx = i // the index of reader
				heap.Push(results, rGnm)
			}

			return nil
		}

		// read

		var reader *SearchResultReader
		var rGnm *SearchResultOfAGenome
		var n int
		var idx int
		checkColumns := true
		checkQuery := query == ""
		var qlen string
		var queryLength uint64 // numerical length shared by all selected-query records
		var moreColumns bool
		var idx0 int
		countTotalHits := true
		var hits int
		var order int // traversal order preserves ties after the selected heap is sorted

		for {
			if len(*(results.entries)) == 0 {
				checkError(fillBuffer())
				if countTotalHits {
					countTotalHits = false

					for _, r := range *results.entries {
						hits += r.Hits
					}
				}
			}
			if len(*(results.entries)) == 0 {
				break
			}

			rGnm = heap.Pop(results).(*SearchResultOfAGenome)

			// -------------------------------------------------

			// output
			if checkColumns { // for once
				checkColumns = false

				if query == "" {
					query = rGnm.Query
					if outputLog {
						log.Infof("  query inferred from input: %s", query)
					}
				}
				qlen = rGnm.Qlen
				idx0 = rGnm.idx

				moreColumns = rGnm.Records[0].Extra != ""

				if topN == 0 {
					writeMergedSearchHeader(outfh, moreColumns)
				}
			}

			if checkQuery && query != rGnm.Query {
				checkError(fmt.Errorf("inconsistent queries: '%s' in file '%s' and '%s' in file '%s. Please specify one query with flag -q/--query",
					query, files[idx0], rGnm.Query, files[rGnm.idx],
				))
			}
			// Each reader validates its own rows; check consistency across files too.
			length, _ := strconv.ParseUint(rGnm.Qlen, 10, 64)
			if queryLength == 0 {
				queryLength = length
			} else if length != queryLength {
				checkError(fmt.Errorf("inconsistent query lengths in %s and %s: %s != %s", files[idx0], files[rGnm.idx], qlen, rGnm.Qlen))
			}

			idx = rGnm.idx
			if topN == 0 {
				n++
				writeMergedGenomeResult(outfh, rGnm, query, qlen, hits, moreColumns)
				RecycleSearchResultOfAGenome(rGnm)
			} else {
				// Rounded TSV identities can invert the full-precision search order.
				// Inspect every genome, retaining only the best N and cutoff ties.
				rGnm.order = order
				order++
				selected.Add(rGnm)
			}

			// -------------------------------------------------

			reader = readers[idx]
			if reader != nil {
				rGnm = reader.Next()
				if rGnm == nil {
					delete(readers, idx)
					continue
				}

				rGnm.idx = idx
				heap.Push(results, rGnm)
			}
		}

		if checkColumns && singleInput != nil && singleInput.header != "" {
			fmt.Fprintln(outfh, singleInput.header)
		}
		if topN > 0 && len(retained) > 0 {
			retained = selected.Results()
			n = len(retained)
			moreColumns = retained[0].Records[0].Extra != ""
			writeMergedSearchHeader(outfh, moreColumns)
		}
		for i, r := range retained {
			writeMergedGenomeResult(outfh, r, query, qlen, n, moreColumns)
			RecycleSearchResultOfAGenome(r)
			retained[i] = nil
		}

		if outputLog {
			log.Infof("%d genome hits merged from %d files for query: %s", n, len(files), query)
		}
	},
}

func init() {
	utilsCmd.AddCommand(mergeCmd)

	mergeCmd.Flags().StringP("out-file", "o", "-",
		formatFlagUsage(`Out file, supports the ".gz" suffix ("-" for stdout).`))

	mergeCmd.Flags().StringP("query", "q", "",
		formatFlagUsage(`Query ID to merge`))

	mergeCmd.Flags().IntP("top-n-genomes", "n", 0,
		formatFlagUsage(`Keep the top N genome matches by the highest alignment bitscore * pident per genome, including all matches tied at the cutoff score (0 for all). This filters alignment results, not chaining scores.`))

	mergeCmd.Flags().StringP("buffer-size", "b", "20M",
		formatFlagUsage(`Size of buffer, supported unit: K, M, G. You need increase the value when "bufio.Scanner: token too long" error reported`))

	mergeCmd.SetUsageTemplate(usageTemplate(""))
}

// writeMergedSearchHeader waits for Top-N selection to determine the output format.
func writeMergedSearchHeader(out io.Writer, moreColumns bool) {
	fmt.Fprint(out, "query\tqlen\thits\tsgenome\tsseqid\tqcovGnm\tcls\thsp\tqcovHSP\talenHSP\tpident\tgaps\tqstart\tqend\tsstart\tsend\tsstr\tslen\tevalue\tbitscore")
	if moreColumns {
		fmt.Fprint(out, "\tcigar\tqseq\tsseq\talign")
	}
	fmt.Fprintln(out)
}

func writeMergedGenomeResult(out io.Writer, r *SearchResultOfAGenome, query, qlen string, hits int, moreColumns bool) {
	for _, record := range r.Records {
		fmt.Fprintf(out, "%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
			query, qlen, hits, r.Sgenome, record.Sseqid, r.QcovGnm,
			record.Cls, record.Hsp, record.QcovHSP, record.AlenHSP, record.Pident,
			record.Gaps, record.Qstart, record.Qend, record.Sstart, record.Send,
			record.Sstr, record.Slen, record.Evalue, record.Bitscore)
		if moreColumns {
			if record.Extra == "" {
				// Preserve four optional columns when merging a 20-column row.
				fmt.Fprint(out, "\t\t\t\t")
			} else {
				fmt.Fprintf(out, "\t%s", record.Extra)
			}
		}
		fmt.Fprintln(out)
	}
}

type SearchResultOfAGenome struct {
	idx   int // the index of reader, only for merge
	order int // traversal order for stable Top-N ties

	Hits    int
	Query   string
	Qlen    string
	Sgenome string
	QcovGnm string

	Score float64 // BitScore * PIdent; three-decimal TSV ties are preserved

	Records []*SearchResultOfASequence
}

var poolSearchResultOfAGenome = &sync.Pool{New: func() interface{} {
	r := &SearchResultOfAGenome{
		Records: make([]*SearchResultOfASequence, 0, 8),
	}
	return r
}}

func RecycleSearchResultOfAGenome(r *SearchResultOfAGenome) {
	r.order = 0
	r.Hits = 0
	r.Query = ""
	r.Qlen = ""
	r.Sgenome = ""
	r.QcovGnm = ""
	r.Score = 0
	for _, record := range r.Records {
		// Parsed strings share the input line, which may include sequences.
		*record = SearchResultOfASequence{}
		poolSearchResultOfASequence.Put(record)
	}
	clear(r.Records)
	r.Records = r.Records[:0]
	poolSearchResultOfAGenome.Put(r)
}

type SearchResultOfASequence struct {
	Sseqid   string
	Cls      string // should be smaller than int
	Hsp      string // should be smaller than int
	QcovHSP  string
	AlenHSP  string
	Pident   string
	Gaps     string
	Qstart   string
	Qend     string
	Sstart   string
	Send     string
	Sstr     string
	Slen     string
	Evalue   string
	Bitscore string

	Extra string
}

var poolSearchResultOfASequence = &sync.Pool{New: func() interface{} {
	return &SearchResultOfASequence{}
}}

type SearchResultsHeap struct {
	entries *[]*SearchResultOfAGenome
}

func (h SearchResultsHeap) Len() int { return len(*(h.entries)) }

func (h SearchResultsHeap) Less(i, j int) bool {
	return (*(h.entries))[i].Score > (*(h.entries))[j].Score
}

func (h SearchResultsHeap) Swap(i, j int) {
	(*(h.entries))[i], (*(h.entries))[j] = (*(h.entries))[j], (*(h.entries))[i]
}

func (h SearchResultsHeap) Push(x interface{}) {
	*(h.entries) = append(*(h.entries), x.(*SearchResultOfAGenome))
}

func (h SearchResultsHeap) Pop() interface{} {
	n := len(*(h.entries))
	x := (*(h.entries))[n-1]
	(*(h.entries))[n-1] = nil
	*(h.entries) = (*(h.entries))[:n-1]
	return x
}

// topSearchResults owns a min-heap of at most limit genomes, plus all extra
// genomes tied at its minimum score. A rising cutoff releases obsolete ties.
type topSearchResults struct {
	SearchResultsHeap
	limit int
	ties  []*SearchResultOfAGenome
}

func (h topSearchResults) Less(i, j int) bool {
	a, b := (*h.entries)[i], (*h.entries)[j]
	return a.Score < b.Score || a.Score == b.Score && a.order > b.order
}

func (h *topSearchResults) Add(r *SearchResultOfAGenome) {
	if h.Len() < h.limit {
		heap.Push(h, r)
		return
	}
	cutoff := (*h.entries)[0].Score
	if r.Score < cutoff {
		RecycleSearchResultOfAGenome(r)
		return
	}
	if r.Score == cutoff {
		h.ties = append(h.ties, r)
		return
	}
	old := heap.Pop(h).(*SearchResultOfAGenome)
	heap.Push(h, r)
	if (*h.entries)[0].Score == cutoff {
		h.ties = append(h.ties, old)
		return
	}
	RecycleSearchResultOfAGenome(old)
	for _, tied := range h.ties {
		RecycleSearchResultOfAGenome(tied)
	}
	clear(h.ties)
	h.ties = h.ties[:0]
}

// Results transfers retained records in score order, preserving traversal ties.
func (h *topSearchResults) Results() []*SearchResultOfAGenome {
	results := append(*h.entries, h.ties...)
	clear(h.ties)
	h.ties = nil
	*h.entries = results
	slices.SortFunc(results, func(a, b *SearchResultOfAGenome) int {
		if order := cmp.Compare(b.Score, a.Score); order != 0 {
			return order
		}
		return cmp.Compare(a.order, b.order)
	})
	return results
}

type SearchResultReader struct {
	query  string
	header string // published before the first result or the channel close

	fh *xopen.Reader

	ch chan *SearchResultOfAGenome
}

// openSearchResultFile reports empty input with its source path.
func openSearchResultFile(file string) (*xopen.Reader, error) {
	fh, err := xopen.Ropen(file)
	if errors.Is(err, xopen.ErrNoContent) {
		return nil, fmt.Errorf("empty search result file %q: %w", file, err)
	}
	if err != nil {
		return nil, fmt.Errorf("failed to open search result file %q: %w", file, err)
	}
	return fh, nil
}

func NewSearchResultReader(file string, query string, bufferSize int64) (*SearchResultReader, error) {
	if bufferSize < 1 || bufferSize > maxScannerBufferSize {
		return nil, fmt.Errorf("buffer size must be in [1, %d] bytes", maxScannerBufferSize)
	}
	fh, err := openSearchResultFile(file)
	if err != nil {
		return nil, err
	}

	r := &SearchResultReader{
		query: query,
		fh:    fh,
		ch:    make(chan *SearchResultOfAGenome),
	}

	// ------------------------------------------

	scanner := bufio.NewScanner(fh)
	// Start small per input; Scanner grows for long -a rows up to bufferSize.
	buf := make([]byte, min(bufferSize, 64<<10))
	scanner.Buffer(buf, int(bufferSize))
	ncols := 20

	go func() {
		defer fh.Close()
		var line string
		var lineNumber uint64 // Physical line number, including headers and blank lines.
		headerLine := true
		var query, qlen, hits, sgenome, sseqid, qcovGnm, cls, hsp, qcovHSP, alenHSP string
		var pident, gaps, qstart, qend, sstart, send, sstr, slen, evalue, bitscore string
		var extra string
		var _pident float64
		var _bitscore int
		var _hits int

		var sgenomePre, firstQuery string
		var firstQueryLength uint64

		items := make([]string, ncols)

		checkQuery := r.query != ""

		rGnm := poolSearchResultOfAGenome.Get().(*SearchResultOfAGenome)
		var score float64

		for scanner.Scan() {
			lineNumber++
			line = strings.TrimRight(scanner.Text(), "\r\n")
			if line == "" {
				continue
			}
			if headerLine {
				headerLine = false
				r.header = line
				continue
			}
			stringSplitNByByte(line, '\t', ncols, &items)
			if len(items) < ncols {
				checkError(fmt.Errorf("the input has only %d columns (<%d), please use output from 'lexicmap search'", len(items), ncols))
			}

			query = items[0]

			if checkQuery && query != r.query {
				continue
			}
			// Validate before subject grouping can discard a row's query ID.
			// An explicit -q still filters a multi-query file as documented.
			if !checkQuery {
				if firstQuery == "" {
					firstQuery = query
				} else if query != firstQuery {
					checkError(fmt.Errorf("multiple queries in %s: %s and %s; use -q/--query", file, firstQuery, query))
				}
			}

			qlen = items[1]
			length, parseErr := strconv.ParseUint(qlen, 10, 64)
			if parseErr != nil || length == 0 {
				checkError(fmt.Errorf("invalid query length %q in %s", qlen, file))
			}
			// Validate before subject grouping can hide a row's query length.
			if firstQueryLength == 0 {
				firstQueryLength = length
			} else if length != firstQueryLength {
				checkError(fmt.Errorf("inconsistent query lengths for %q in %s: %d != %d", query, file, firstQueryLength, length))
			}
			hits = items[2]
			sgenome = items[3]
			sseqid = items[4]
			qcovGnm = items[5]
			cls = items[6]
			hsp = items[7]
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
			evalue = items[18]

			// Optional alignment fields belong only to the current row.
			extra = ""
			if i := strings.IndexByte(items[19], '\t'); i > 0 { // the search result has >=20 columns
				bitscore = items[19][:i]
				extra = items[19][i+1:]
			} else {
				bitscore = items[19]
			}

			_bitscore, err = strconv.Atoi(bitscore)
			if err != nil {
				checkError(fmt.Errorf("%s:%d: invalid bitscore %q: expected an integer", file, lineNumber, bitscore))
			}

			_hits, err = strconv.Atoi(hits)
			if err != nil || _hits <= 0 {
				checkError(fmt.Errorf("%s:%d: invalid hits %q: expected a positive integer", file, lineNumber, hits))
			}

			_pident, err = strconv.ParseFloat(pident, 64)
			// NaN breaks heap comparisons; other invalid identities distort ranking.
			if err != nil || math.IsNaN(_pident) || math.IsInf(_pident, 0) || _pident < 0 || _pident > 100 {
				checkError(fmt.Errorf("%s:%d: invalid pident %q: expected a finite value in [0,100]", file, lineNumber, pident))
			}

			// -------

			if sgenome != sgenomePre { // new one
				if rGnm.Sgenome != "" {
					r.ch <- rGnm

					rGnm = poolSearchResultOfAGenome.Get().(*SearchResultOfAGenome)
				}

				rGnm.Hits = _hits
				rGnm.Query = query
				rGnm.Qlen = qlen
				rGnm.Sgenome = sgenome
				rGnm.QcovGnm = qcovGnm

				sgenomePre = sgenome
			}

			rSeq := poolSearchResultOfASequence.Get().(*SearchResultOfASequence)

			rSeq.Sseqid = sseqid
			rSeq.Cls = cls
			rSeq.Hsp = hsp
			rSeq.QcovHSP = qcovHSP
			rSeq.AlenHSP = alenHSP
			rSeq.Pident = pident
			rSeq.Gaps = gaps
			rSeq.Qstart = qstart
			rSeq.Qend = qend
			rSeq.Sstart = sstart
			rSeq.Send = send
			rSeq.Sstr = sstr
			rSeq.Slen = slen
			rSeq.Evalue = evalue
			rSeq.Bitscore = bitscore
			rSeq.Extra = extra

			rGnm.Records = append(rGnm.Records, rSeq)

			score = searchResultScore(_bitscore, _pident)
			// Seed the maximum from a real HSP, including a negative score.
			if len(rGnm.Records) == 1 || score > rGnm.Score {
				rGnm.Score = score
			}
		}
		if err = scanner.Err(); err != nil {
			checkError(fmt.Errorf("failed to scan file %s: %s", file, err))
		}

		if rGnm.Sgenome != "" {
			r.ch <- rGnm // do not forget the last one
		} else {
			RecycleSearchResultOfAGenome(rGnm)
		}

		close(r.ch)
	}()

	return r, nil
}

func (r *SearchResultReader) Next() *SearchResultOfAGenome {
	return <-r.ch
}

// searchResultScore preserves exact ties for search's three-decimal identities
// by multiplying integer thousandths. Keep additional input precision when present.
func searchResultScore(bitscore int, pident float64) float64 {
	if identity := math.Round(pident * 1000); identity/1000 == pident {
		return float64(bitscore) * identity / 1000
	}
	return float64(bitscore) * pident
}
