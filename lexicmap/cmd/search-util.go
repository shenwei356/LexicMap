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
	"strconv"
	"strings"
	"sync"

	"github.com/shenwei356/xopen"
)

// Strands could be used to output strand for a reverse complement flag
var Strands = [2]byte{'+', '-'}

// Query is an object for each query sequence, it also contains the query result.
type Query struct {
	seqID  []byte
	seq    []byte
	result *[]*SearchResult
}

// Reset reset the data for next round of using
func (q *Query) Reset() {
	q.seqID = q.seqID[:0]
	q.seq = q.seq[:0]
	q.result = nil
}

var poolQuery = &sync.Pool{New: func() interface{} {
	return &Query{
		seqID: make([]byte, 0, 128),     // the id should be not too long
		seq:   make([]byte, 0, 100<<10), // initialize with 100K
	}
}}

func parseTaxids(taxdumpDir string, genome2taxidFile string, taxidsStr []string, taxidFile string) ([]uint32, []uint32) {
	var err error

	var taxids, negativeTaxids []uint32

	var m, negativeM map[uint32]interface{}
	var ok bool
	var v uint32

	if len(taxidsStr) > 0 {
		if !(taxdumpDir != "" && genome2taxidFile != "") {
			checkError(fmt.Errorf("flags -T/--taxdump and -G/--genome2taxid are need if -t/--taxids is given"))
		}
		m = make(map[uint32]interface{}, len(taxidsStr))
		negativeM = make(map[uint32]interface{}, len(taxidsStr))
		taxids = make([]uint32, 0, len(taxidsStr))
		negativeTaxids = make([]uint32, 0, len(taxidsStr))

		var val int64
		for _, tmp := range taxidsStr {
			val, err = strconv.ParseInt(tmp, 10, 32)
			if err != nil {
				checkError(fmt.Errorf("invalid TaxId: %s", tmp))
			}

			if val > 0 {
				v = uint32(val)
				if _, ok = m[v]; !ok {
					taxids = append(taxids, v)
					m[v] = struct{}{}
				}
			} else if val < 0 {
				v = uint32(-val)
				if _, ok = negativeM[v]; !ok {
					negativeTaxids = append(negativeTaxids, v)
					negativeM[v] = struct{}{}
				}
			}
		}
	}
	if taxidFile != "" {
		if m == nil {
			m = make(map[uint32]interface{}, len(taxidsStr))
		}

		fh, err := xopen.Ropen(taxidFile)
		if err != nil {
			checkError(fmt.Errorf("failed to read taxid file: %s", taxidFile))
		}

		scanner := bufio.NewScanner(fh)
		var line string
		var val int64
		for scanner.Scan() {
			line = strings.TrimSpace(strings.TrimRight(scanner.Text(), "\r\n"))
			if line == "" {
				continue
			}

			val, err = strconv.ParseInt(line, 10, 32)
			if err != nil {
				checkError(fmt.Errorf("invalid TaxId: %s", line))
			}

			if val > 0 {
				v = uint32(val)
				if _, ok = m[v]; !ok {
					taxids = append(taxids, v)
					m[v] = struct{}{}
				}
			} else if val < 0 {
				v = uint32(-val)
				if _, ok = negativeM[v]; !ok {
					negativeTaxids = append(negativeTaxids, v)
					negativeM[v] = struct{}{}
				}
			}
		}
		if err = scanner.Err(); err != nil {
			checkError(fmt.Errorf("failed to read taxid file: %s", taxidFile))
		}
	}
	// } else if taxdumpDir != "" {
	// 	checkError(fmt.Errorf("the flag -T/--taxdump is given, but -t/--taxids is not"))
	// } else if genome2taxidFile != "" {
	// 	checkError(fmt.Errorf("the flag -G/--genome2taxid is given, but -t/--taxids is not"))
	// }

	return taxids, negativeTaxids
}
