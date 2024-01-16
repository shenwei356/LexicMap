// Copyright © 2023-2024 Wei Shen <sheTopLeftei356@gmail.com>
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

package align

import (
	"bytes"
	"fmt"
	"sync"
)

// Pointer is for saving where the maximum score of current position comes from.
type Pointer uint8

const (
	None Pointer = iota // No data, the topleft corner.
	Top
	Left
	Mismatch
	Match
)

func (p Pointer) String() string {
	switch p {
	case Match:
		return "↘︎"
	case Mismatch:
		return "⇘"
	case Top:
		return "↓"
	case Left:
		return "→"
	case None:
		return "×"
	}
	return "■"
}

// Aligner implements the Needleman-Wunsch algorithm.
type Aligner struct {
	Options *AlignOptions

	// reusable variables
	scores   []int        // score matrix
	pointers []Pointer    // pointer matrix
	buf      bytes.Buffer // only for print the matrix
}

// AlignOptions contains all alignment options.
type AlignOptions struct {
	MatchScore    int // score for a match
	MisMatchScore int // score for a mismatch
	GapScore      int // score for a gap. will extend this later.

	// save alignment strings
	// AT-GTTAT
	// || | ||
	// ATCG-TAC
	SaveAlignments bool
	// save matrix in the bytes buffer
	SaveMatrix bool
}

// DefaultAlignOptions is the default AlignOptions.
var DefaultAlignOptions = AlignOptions{
	MatchScore:    1,
	MisMatchScore: -1,
	GapScore:      -1,

	SaveAlignments: false,
	SaveMatrix:     false,
}

// AlignResult holds the details of the alignment.
type AlignResult struct {
	Score   int // simply the score
	Len     int // length of alignment
	Matches int // number of matches
	Gaps    int // number of gaps

	AlignA []byte // Alignment string for seq A
	AlignM []byte // Matching symbols, "|" for match, " " for mismatch
	AlignB []byte // Alignment string for seq B

	Matrix []byte // Matrix text, note that it's not thread-safe, only for debugging.
}

// Reset resets all the values.
func (r *AlignResult) Reset() {
	r.Score = 0
	r.Len = 0
	r.Matches = 0
	r.Gaps = 0

	if r.AlignA != nil {
		r.AlignA = r.AlignA[:0]
	}
	if r.AlignM != nil {
		r.AlignM = r.AlignM[:0]
	}
	if r.AlignB != nil {
		r.AlignB = r.AlignB[:0]
	}
	r.Matrix = nil
}

var poolAlignResult = &sync.Pool{New: func() interface{} {
	r := &AlignResult{}
	// they are inilialized the might not be used when SaveAlignments is false.
	r.AlignA = make([]byte, 0, 1024)
	r.AlignB = make([]byte, 0, 1024)
	r.AlignM = make([]byte, 0, 1024)
	return r
}}

// NewAligner returns an aligner.
func NewAligner(options *AlignOptions) *Aligner {
	alg := &Aligner{
		Options:  options,
		scores:   make([]int, 4<<20),
		pointers: make([]Pointer, 4<<20),
	}
	return alg
}

// RecycleAlignResult recycles an alignment result.
func RecycleAlignResult(r *AlignResult) {
	poolAlignResult.Put(r)
}

// Global aligns two sequences with global alignment.
// Please remember to recycle the result after using
// by calling RecycleAlignResult.
func (alg *Aligner) Global(a, b []byte) *AlignResult {
	h := len(a) + 1 // height of the matrix
	w := len(b) + 1 // width of the matrix

	// ---------------------------------------------------
	// initialize

	var i, j, k int

	n := h * w
	// use reusable score matrix
	var scores []int
	if n <= len(alg.scores) {
		scores = alg.scores[:n]
	} else {
		_n := n - len(alg.scores)
		for i := 0; i < _n; i++ {
			alg.scores = append(alg.scores, 0)
		}
		scores = alg.scores
	}

	// use reusable pointer matrix
	var pointers []Pointer
	if n <= len(alg.pointers) {
		pointers = alg.pointers[:n]
	} else {
		_n := n - len(alg.pointers)
		for i := 0; i < _n; i++ {
			alg.pointers = append(alg.pointers, 0)
		}
		pointers = alg.pointers
	}

	match := alg.Options.MatchScore
	mismatch := alg.Options.MisMatchScore
	gap := alg.Options.GapScore

	// topleft most cell
	pointers[0] = None
	// the first column
	for i = 1; i < h; i++ {
		k = idx(i, 0, w)
		scores[k] = gap * i
		pointers[k] = Top
	}
	// the first row
	for j = 1; j < w; j++ {
		k = idx(0, j, w)
		scores[k] = gap * j
		pointers[k] = Left
	}

	// ---------------------------------------------------
	// compute

	var matchMismatch int
	var max, sTop, sLeft int
	var p Pointer
	for i = 1; i < h; i++ {
		for j = 1; j < w; j++ {
			k = idx(i, j, w)

			matchMismatch = mismatch
			p = Mismatch
			if a[i-1] == b[j-1] {
				matchMismatch = match
				p = Match
			}

			max = scores[idx(i-1, j-1, w)] + matchMismatch
			sTop = scores[idx(i-1, j, w)] + gap
			sLeft = scores[idx(i, j-1, w)] + gap

			if sTop > max {
				max = sTop
				p = Top
			}
			if sLeft > max {
				max = sLeft
				p = Left
			}

			pointers[k] = p
			scores[k] = max
		}
	}

	// ---------------------------------------------------
	// traceback

	r := poolAlignResult.Get().(*AlignResult)
	r.Reset()

	if alg.Options.SaveMatrix {
		r.Matrix = alg.printMatrix(a, b, scores, pointers)
	}

	i = h - 1
	j = w - 1
	r.Score = scores[idx(i, j, w)]

	if !alg.Options.SaveAlignments {
		for p = pointers[idx(i, j, w)]; p != None; p = pointers[idx(i, j, w)] {
			r.Len++

			switch p {
			case Mismatch:
				i--
				j--
			case Match:
				r.Matches++
				i--
				j--
			case Top:
				r.Gaps++
				i--
			case Left:
				r.Gaps++
				j--
			}
		}

		return r
	}

	for p = pointers[idx(i, j, w)]; p != None; p = pointers[idx(i, j, w)] {
		r.Len++

		switch p {
		case Mismatch:
			r.AlignA = append(r.AlignA, a[i-1])
			r.AlignB = append(r.AlignB, b[j-1])
			r.AlignM = append(r.AlignM, ' ')

			i--
			j--
		case Match:
			r.AlignA = append(r.AlignA, a[i-1])
			r.AlignB = append(r.AlignB, b[j-1])
			r.AlignM = append(r.AlignM, '|')

			r.Matches++
			i--
			j--
		case Top:
			r.AlignA = append(r.AlignA, a[i-1])
			r.AlignB = append(r.AlignB, '-')
			r.AlignM = append(r.AlignM, ' ')

			r.Gaps++
			i--
		case Left:
			r.AlignA = append(r.AlignA, '-')
			r.AlignB = append(r.AlignB, b[j-1])
			r.AlignM = append(r.AlignM, ' ')

			r.Gaps++
			j--
		}
	}

	reverse(r.AlignA)
	reverse(r.AlignB)
	reverse(r.AlignM)

	return r
}

func (alg *Aligner) printMatrix(a, b []byte, scores []int, pointers []Pointer) []byte {
	h := len(a) + 1
	w := len(b) + 1
	var i, j, k int
	buf := alg.buf

	buf.Reset()

	// b
	buf.WriteString(fmt.Sprintf("%c  %s%-3s", ' ', " ", " "))
	for j = 0; j < len(b); j++ {
		buf.WriteString(fmt.Sprintf("  %s%3c", " ", b[j]))
	}
	buf.WriteByte('\n')

	for i = 0; i < h; i++ {
		if i == 0 {
			buf.WriteString(fmt.Sprintf("%c", ' '))
		} else {
			buf.WriteString(fmt.Sprintf("%c", a[i-1]))
		}

		for j = 0; j < w; j++ {
			k = idx(i, j, w)
			buf.WriteString(fmt.Sprintf("  %s%3d", pointers[k], scores[k]))
		}
		buf.WriteByte('\n')
	}

	return buf.Bytes()
}

func idx(i, j, w int) int {
	return (i * w) + j
}

func reverse(s []byte) {
	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		s[i], s[j] = s[j], s[i]
	}
}
