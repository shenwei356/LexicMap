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

package twobit

import (
	"bytes"
	"fmt"
	"os"
	"testing"
)

func TestGenomeWritingAndSeqExtraction(t *testing.T) {
	_seq := []byte("ACTAGACGACGTACGCGTACGTAGTACGATGCTCGA")
	var s, s2 []byte
	var b2 *[]byte
	var err error
	for n := 1; n < len(_seq); n++ {
		s = _seq[:n]
		b2 = Seq2TwoBit(s)
		s2, err = TwoBit2Seq(*b2, n)
		if err != nil {
			t.Error(err)
			return
		}
		if !bytes.Equal(s, s2) {
			t.Errorf("expected: %s, results: %s\n", s, s2)
			return
		}
		RecycleTwoBit(b2)
	}
}

func TestReadAndWrite(t *testing.T) {
	file := "t.2bit"

	// ----------------------- write --------------

	w, err := NewWriter(file, 1)
	if err != nil {
		t.Error(err)
		return
	}

	_seqs := [][]byte{
		[]byte("A"),
		[]byte("C"),
		[]byte("CA"),
		[]byte("CAT"),
		[]byte("CATG"),
		[]byte("CATGC"),
		[]byte("CATGCC"),
		[]byte("CATGCCA"),
		[]byte("CATGCCAC"),
		[]byte("CATGCCACG"),
		[]byte("ACCCTCGAGCGACTAG"),
		[]byte("ACTAGACGACGTACGCGTACGTAGTACGATGCTCGA"),
		[]byte("ACGCAGTCGTCATCATGCGTGTCGCATGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACATGCTGCATGCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGCTGTGATGCGTCTCAGTAGATGAT"),
	}

	for i, s := range _seqs {
		id := []byte(fmt.Sprintf("seq_%d", i+1))

		g := PoolGenome.Get().(*Genome)
		g.Reset()
		g.ID = append(g.ID, id...)
		g.Seq = append(g.Seq, s...)
		g.GenomeSize = len(s)
		g.Len = len(s)
		g.NumSeqs = 1
		g.SeqSizes = append(g.SeqSizes, len(s))

		err = w.Write(g)
		if err != nil {
			t.Error(err)
			return
		}

		RecycleGenome(g)
	}

	err = w.Close()
	if err != nil {
		t.Error(err)
		return
	}

	// ----------------------- read --------------

	r, err := NewReader(file)
	if err != nil {
		t.Error(err)
		return
	}

	var start, end int
	var s1 []byte
	var s2 *Genome
	for i, s := range _seqs {
		// subseq
		for start = 0; start < len(s); start++ {
			for end = start; end < len(s); end++ {
				s2, err = r.SubSeq(i, start, end)
				if err != nil {
					t.Error(err)
					return
				}
				s1 = s[start : end+1]
				if !bytes.Equal(s1, s2.Seq) {
					t.Errorf("idx: %d:%d-%d, expected: %s, results: %s",
						i, start, end, s1, s2.Seq)
					return
				}
				RecycleGenome(s2)
			}
		}

		// whole seq
		s2, err = r.Seq(i)
		if err != nil {
			t.Error(err)
			return
		}
		if !bytes.Equal(s, s2.Seq) {
			t.Errorf("idx: %d not matched", i)
		}
		RecycleGenome(s2)
	}

	r.Close()

	// clean up

	err = os.RemoveAll(file)
	if err != nil {
		t.Error(err)
		return
	}

	err = os.RemoveAll(file + GenomeIndexFileExt)
	if err != nil {
		t.Error(err)
		return
	}
}
