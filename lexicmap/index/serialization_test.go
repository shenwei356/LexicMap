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

package index

import (
	"io"
	"os"
	"runtime"
	"strconv"
	"testing"
	"time"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/lexichash"
)

func TestSerialization(t *testing.T) {
	k := 21
	nMasks := 1000
	var seed int64 = 1

	idx, err := NewIndexWithSeed(k, nMasks, seed, 0)
	if err != nil {
		t.Error(err)
		return
	}

	sTime := time.Now()
	t.Logf("starting to build the index ...")

	input, done := idx.BatchInsert()

	seq.ValidateSeq = false
	var record *fastx.Record
	fastxReader, err := fastx.NewReader(nil, "test_data/hairpin.fasta", "")
	if err != nil {
		t.Error(err)
		return
	}
	var nSeqs int
	for {
		record, err = fastxReader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			t.Error(err)
			return
		}

		nSeqs++

		_seq := make([]byte, len(record.Seq.Seq))
		copy(_seq, record.Seq.Seq)
		input <- &RefSeq{
			ID:  []byte(string(record.ID)),
			Seq: _seq,
		}
	}
	close(input) // wait BatchInsert
	<-done       // wait BatchInsert

	t.Logf("finished to build the index in %s from %d sequences with %d masks",
		time.Since(sTime), nSeqs, nMasks)

	// ---------------------------------------------------

	outDir := "testindex"
	threads := runtime.NumCPU()

	sTime = time.Now()

	err = idx.WriteToPath(outDir, true, threads)
	if err != nil {
		t.Logf("%s", err)
	}

	t.Logf("finished to dump the index in %s", time.Since(sTime))

	// ---------------------------------------------------

	sTime = time.Now()

	idx2, err := NewFromPath(outDir, threads)
	if err != nil {
		t.Logf("%s", err)
	}

	t.Logf("finished to read the index in %s", time.Since(sTime))

	if idx.K() != idx2.K() {
		t.Errorf("Ks unmatched: %d vs %d", idx.K(), idx2.K())
		return
	}

	if len(idx.lh.Masks) != len(idx2.lh.Masks) {
		t.Errorf("number of masks unmatched: %d vs %d", len(idx.lh.Masks), len(idx2.lh.Masks))
		return
	}

	var m2 uint64
	for i, m := range idx.lh.Masks {
		if m != idx2.lh.Masks[i] {
			t.Errorf("unmatched mask: %d vs %d", m, m2)
			return
		}
	}

	t.Logf("checking loaded index passed")

	// ---------------------------------------------------

	queries, err := fastx.GetSeqs("test_data/hairpin.query.fasta", nil, 8, 100, "")
	if err != nil {
		t.Error(err)
		return
	}
	if len(queries) == 0 {
		t.Error(err)
		return
	}

	minLen := 13

	decoder := lexichash.MustDecoder()
	idx.SetSearchingOptions(&SearchOptions{
		MinPrefix: uint8(minLen),
		TopN:      0,
	})
	for _, s := range queries {
		sr, err := idx2.Search(s.Seq.Seq)
		if err != nil {
			t.Log(err)
			return
		}
		if sr == nil {
			continue
		}
		t.Log()
		t.Logf("query: %s, targets: %d\n", s.ID, len(*sr))

		for i, r := range *sr {
			t.Logf("%4s %s\n", "#"+strconv.Itoa(i+1), idx.IDs[r.IdIdx])
			for _, v := range *r.Subs {
				t.Logf("     (%3d,%3d) vs (%3d,%3d) %3d %s\n",
					v.QBegin+1, v.QBegin+v.Len,
					v.TBegin+1, v.TBegin+v.Len,
					v.Len, decoder(v.Code, uint8(v.Len)))
			}
		}
		idx.RecycleSearchResults(sr)
	}

	// ---------------------------------------------------

	idx.ExtractKmerLocations()

	err = idx.WriteKmerLocations()
	if err != nil {
		t.Log(err)
		return
	}

	old := idx.KmerLocations

	err = idx.ReadKmerLocations()
	if err != nil {
		t.Log(err)
		return
	}

	if len(old) != len(idx.KmerLocations) {
		t.Log("number of mask location records error")
		return
	}

	var vs2 []uint64
	var v, v2 uint64
	var j int
	for i, vs := range old {
		vs2 = idx.KmerLocations[i]
		if len(vs) != len(vs2) {
			t.Logf("number of mask locations error for ref id: %d", i)
			return
		}
		for j, v = range vs {
			v2 = vs2[j]
			if v != v2 {
				t.Logf("mask location records error for id: %d", i)
				return
			}
		}
	}

	// ---------------------------------------------------

	if os.RemoveAll(outDir) != nil {
		t.Errorf("failed to remove the directory: %s", outDir)
		return
	}
}
