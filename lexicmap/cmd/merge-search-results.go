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
	"container/heap"
	"fmt"
	"io"
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
     If not, please specify one query with the flag -q/--query
  2. We assume that genome IDs are distinct across all indexes.
  3. One or more input files are accepted, via positional parameters
     and/or a file list via the flag -X/--infile-list.
  4. Both the default 20- and 24-column formats are supported,
     and formats better be consistent across all input files.
     If not, the output format would be the one with a valid record.

`,
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		seq.ValidateSeq = false

		// ------------------------------

		outFile := getFlagString(cmd, "out-file")

		query := getFlagString(cmd, "query")

		bufferSizeS := getFlagString(cmd, "buffer-size")
		if bufferSizeS == "" {
			checkError(fmt.Errorf("value of buffer size. supported unit: K, M, G"))
		}
		bufferSize, err := ParseByteSize(bufferSizeS)
		if err != nil {
			checkError(fmt.Errorf("invalid value of buffer size. supported unit: K, M, G"))
		}

		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)

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

		if len(files) == 1 {
			if opt.Verbose {
				log.Infof("only one input file '%s' is given, just copy data to '%s'", files[0], outFile)
			}
			fh, err := xopen.Ropen(files[0])
			checkError(err)

			_, err = io.Copy(outfh, fh)
			if err != nil {
				checkError(fmt.Errorf("failed to copy data from %s to %s: %s", files[0], outFile, err))
			}
			return
		}

		// readers

		readers := make(map[int]*SearchResultReader, len(files))
		for i, file := range files {
			reader, err := NewSearchResultReader(file, query, bufferSize)
			checkError(err)
			readers[i] = reader
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
		var rSeq *SearchResultOfASequence
		var n int
		var idx int
		checkColumns := true
		checkQuery := query == ""
		var qlen string
		var moreColumns bool
		var idx0 int
		countTotalHits := true
		var hits int

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
				}
				qlen = rGnm.Qlen
				idx0 = rGnm.idx

				moreColumns = rGnm.Records[0].Extra != ""

				fmt.Fprintf(outfh, "query\tqlen\thits\tsgenome\tsseqid\tqcovGnm\tcls\thsp\tqcovHSP\talenHSP\tpident\tgaps\tqstart\tqend\tsstart\tsend\tsstr\tslen\tevalue\tbitscore")
				if moreColumns {
					fmt.Fprintf(outfh, "\tcigar\tqseq\tsseq\talign")
				}
				fmt.Fprintln(outfh)
			}

			if checkQuery && query != rGnm.Query {
				checkError(fmt.Errorf("inconsistent queries: '%s' in file '%s' and '%s' in file '%s. Please specify one query with flag -q/--query",
					query, files[idx0], rGnm.Query, files[rGnm.idx],
				))
			}

			n++

			for _, rSeq = range rGnm.Records {
				fmt.Fprintf(outfh, "%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
					query,
					qlen,
					hits,
					rGnm.Sgenome,
					rSeq.Sseqid,
					rGnm.QcovGnm,
					rSeq.Cls,
					rSeq.Hsp,
					rSeq.QcovHSP,
					rSeq.AlenHSP,
					rSeq.Pident,
					rSeq.Gaps,
					rSeq.Qstart,
					rSeq.Qend,
					rSeq.Sstart,
					rSeq.Send,
					rSeq.Sstr,
					rSeq.Slen,
					rSeq.Evalue,
					rSeq.Bitscore,
				)
				if moreColumns {
					fmt.Fprintf(outfh, "\t%s", rSeq.Extra)
				}
				fmt.Fprintln(outfh)
			}

			RecycleSearchResultOfAGenome(rGnm)

			// -------------------------------------------------

			idx = rGnm.idx
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

		if opt.Verbose {
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

	mergeCmd.Flags().StringP("buffer-size", "b", "20M",
		formatFlagUsage(`Size of buffer, supported unit: K, M, G. You need increase the value when "bufio.Scanner: token too long" error reported`))

	mergeCmd.SetUsageTemplate(usageTemplate(""))
}

type SearchResultOfAGenome struct {
	idx int // the index of reader, only for merge

	Hits    int
	Query   string
	Qlen    string
	Sgenome string
	QcovGnm string

	Score float64 // float64(BitScore) * PIdent

	Records []*SearchResultOfASequence
}

var poolSearchResultOfAGenome = &sync.Pool{New: func() interface{} {
	r := &SearchResultOfAGenome{
		Records: make([]*SearchResultOfASequence, 0, 8),
	}
	return r
}}

func RecycleSearchResultOfAGenome(r *SearchResultOfAGenome) {
	r.Hits = 0
	r.Sgenome = ""
	r.QcovGnm = ""
	r.Score = 0
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
	*(h.entries) = (*(h.entries))[:n-1]
	return x
}

type SearchResultReader struct {
	query string

	fh *xopen.Reader

	ch chan *SearchResultOfAGenome
}

func NewSearchResultReader(file string, query string, bufferSize int64) (*SearchResultReader, error) {
	fh, err := xopen.Ropen(file)
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
	buf := make([]byte, bufferSize)
	scanner.Buffer(buf, int(bufferSize))
	ncols := 20

	go func() {
		var line string
		headerLine := true
		var query, qlen, hits, sgenome, sseqid, qcovGnm, cls, hsp, qcovHSP, alenHSP string
		var pident, gaps, qstart, qend, sstart, send, sstr, slen, evalue, bitscore string
		var extra string
		var _pident float64
		var _bitscore int
		var _hits int

		var sgenomePre string

		items := make([]string, ncols)

		checkQuery := r.query != ""

		rGnm := poolSearchResultOfAGenome.Get().(*SearchResultOfAGenome)
		var score float64

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
				checkError(fmt.Errorf("the input has only %d columns (<%d), please use output from 'lexicmap search'", ncols, len(items)))
			}

			query = items[0]

			if checkQuery && query != r.query {
				continue
			}

			qlen = items[1]
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

			if i := strings.IndexByte(items[19], '\t'); i > 0 { // the search result has >=20 columns
				bitscore = items[19][:i]
				extra = items[19][i+1:]
			} else {
				bitscore = items[19]
			}

			_bitscore, err = strconv.Atoi(bitscore)
			if err != nil {
				checkError(fmt.Errorf("failed to parse bitsocre: %s", pident))
			}

			_hits, err = strconv.Atoi(hits)
			if err != nil {
				checkError(fmt.Errorf("failed to parse hits: %s", hits))
			}

			_pident, err = strconv.ParseFloat(pident, 64)
			if err != nil {
				checkError(fmt.Errorf("failed to parse pident: %s", pident))
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

			score = float64(_bitscore) * _pident // upate score
			if score > rGnm.Score {
				rGnm.Score = score
			}
		}
		if err = scanner.Err(); err != nil {
			checkError(fmt.Errorf("failed to scan file %s: %s", file, err))
		}

		if rGnm.Sgenome != "" {
			r.ch <- rGnm // do not forget the last one
		}

		close(r.ch)
	}()

	return r, nil
}

func (r *SearchResultReader) Next() *SearchResultOfAGenome {
	return <-r.ch
}
