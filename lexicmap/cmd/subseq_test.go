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
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
)

func TestSubseqGenomeReadersWithoutPool(t *testing.T) {
	dbDir := t.TempDir()
	for batch := 0; batch < 2; batch++ {
		writeSubseqTestGenome(t, dbDir, batch)
	}

	readers, err := newSubseqGenomeReaders(dbDir, 2, 1, 2)
	if err != nil {
		t.Fatal(err)
	}
	defer func() {
		if err := readers.close(); err != nil {
			t.Errorf("close readers: %v", err)
		}
	}()
	if readers.pooled() {
		t.Fatal("expected on-demand readers when the file budget cannot create one reader per batch")
	}

	errs := make(chan error, 2)
	for batch := 0; batch < 2; batch++ {
		go func(batch int) {
			rdr, err := readers.acquire(batch)
			if err != nil {
				errs <- err
				return
			}
			g, err := rdr.SubSeq(0, 0, 3)
			if err == nil {
				genome.RecycleGenome(g)
			}
			if releaseErr := readers.release(batch, rdr); err == nil {
				err = releaseErr
			}
			errs <- err
		}(batch)
	}

	for range 2 {
		select {
		case err := <-errs:
			if err != nil {
				t.Fatal(err)
			}
		case <-time.After(2 * time.Second):
			t.Fatal("timed out waiting for an on-demand genome reader")
		}
	}
}

func writeSubseqTestGenome(t *testing.T, dbDir string, batch int) {
	t.Helper()
	dir := filepath.Join(dbDir, DirGenomes, batchDir(batch))
	if err := os.MkdirAll(dir, 0o755); err != nil {
		t.Fatal(err)
	}
	wtr, err := genome.NewWriter(filepath.Join(dir, FileGenomes), uint32(batch))
	if err != nil {
		t.Fatal(err)
	}

	g := genome.PoolGenome.Get().(*genome.Genome)
	g.Reset()
	g.ID = append(g.ID, fmt.Sprintf("genome-%d", batch)...)
	g.Seq = append(g.Seq, "ACGTACGT"...)
	g.GenomeSize = len(g.Seq)
	g.Len = len(g.Seq)
	g.NumSeqs = 1
	g.SeqSizes = append(g.SeqSizes, len(g.Seq))
	seqID := []byte("seq")
	g.SeqIDs = append(g.SeqIDs, &seqID)
	if err := wtr.Write(g); err != nil {
		genome.RecycleGenome(g)
		_ = wtr.Close()
		t.Fatal(err)
	}
	genome.RecycleGenome(g)
	if err := wtr.Close(); err != nil {
		t.Fatal(err)
	}
}

func TestSubseqConcatenatedCoordinates(t *testing.T) {
	db := writeDetailsFixture(t, 1)
	for _, tc := range []struct {
		region, start, end, sequence string
		reverse                      bool
		extra                        []string
	}{
		{"1:4", "1", "4", "ACGT", false, nil},
		{"1:20", "1", "8", "ACGTACGT", false, nil},
		{"2:5", "2", "5", "TACG", true, nil},
		{"3:4", "1", "7", "ACGTACG", false, []string{"-U", "2", "-D", "3"}},
		{"3:4", "1", "6", "GTACGT", true, []string{"-U", "2", "-D", "3"}},
	} {
		args := []string{"utils", "subseq", "--quiet", "-d", db, "-n", "g0", "-r", tc.region}
		strand := "+"
		if tc.reverse {
			args = append(args, "-R")
			strand = "-"
		}
		output := runSeedIndexCommand(t, append(args, tc.extra...)...)
		want := fmt.Sprintf(">g0:%s-%s:%s\n%s\n", tc.start, tc.end, strand, tc.sequence)
		if string(output) != want {
			t.Fatalf("got %s, want %s", output, want)
		}
	}
}

func TestSubseqDecoratedSequenceIDs(t *testing.T) {
	dir := t.TempDir()
	// Exact IDs that resemble prefixes must win over the stripped fallback.
	a, b := indexTestDNA(1000), indexTestDNA(2000)[1000:]
	ref := writeTestInput(t, filepath.Join(dir, "ref.fa"), ">a:b\n"+a+"\n>c1/1:s1/2:a:b\n"+b+"\n")
	db := filepath.Join(dir, "db")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-c", "2", "--partitions", "4", "-O", db, ref)
	query := writeTestInput(t, filepath.Join(dir, "query.fa"), ">q\n"+a+"\n")
	result := runSeedIndexCommand(t, "search", "--quiet", "-j", "1", "-d", db, "-n", "0", "-N", "0", "--show-sseq-idx", query)
	lines := strings.Split(strings.TrimSpace(string(result)), "\n")
	if len(lines) != 2 {
		t.Fatalf("expected one exact hit: %s", result)
	}
	file := writeTestInput(t, filepath.Join(dir, "result.tsv"), string(result))
	// The decorated spelling is also a real ID: exact lookup is intentional.
	got := string(runSeedIndexCommand(t, "utils", "subseq", "--quiet", "-j", "1", "-d", db, "-f", file, "-w", "0"))
	if !strings.HasPrefix(got, ">c1/1:s1/2:a:b:1-1000:+ ") || !strings.HasSuffix(got, "\n"+b+"\n") {
		t.Fatalf("exact prefix-like ID lost: %s", got)
	}
	// A nested display prefix around that real ID must fall back once.
	fields := strings.Split(lines[1], "\t")
	fields[4] = "c1/1:s2/2:" + fields[4]
	writeTestInput(t, file, lines[0]+"\n"+strings.Join(fields, "\t")+"\n")
	got = string(runSeedIndexCommand(t, "utils", "subseq", "--quiet", "-j", "1", "-d", db, "-f", file, "-w", "0"))
	if !strings.HasSuffix(got, "\n"+b+"\n") {
		t.Fatalf("decorated ID not resolved: %s", got)
	}
	// An unambiguous displayed ID must extract the same bases as its raw ID.
	fields[4] = "c1/1:s1/2:a:b"
	fields[3] = "ref"
	// Use a distinct prefix count to avoid the exact spelling above.
	fields[4] = "c1/1:s1/3:a:b"
	writeTestInput(t, file, lines[0]+"\n"+strings.Join(fields, "\t")+"\n")
	got = string(runSeedIndexCommand(t, "utils", "subseq", "--quiet", "-j", "1", "-d", db, "-f", file, "-w", "0"))
	if !strings.HasPrefix(got, ">a:b:1-1000:+ ") || !strings.HasSuffix(got, "\n"+a+"\n") {
		t.Fatalf("plain ID fallback failed: %s", got)
	}
}
