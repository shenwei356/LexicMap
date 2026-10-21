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
	"encoding/binary"
	"fmt"
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"
	"time"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
)

func TestResolveSearchGenomeID(t *testing.T) {
	a, b, literal := []uint64{1, 2}, []uint64{3}, []uint64{4}
	refs := map[string]*[]uint64{"a:b": &a, "b": &b, "species:a:b": &literal}
	for _, tc := range []struct {
		name string
		want *[]uint64
	}{
		{"a:b", &a}, {"species:a:b", &literal}, {"species:b", &b},
		{"name:with:colons:b", &b},
	} {
		got, err := resolveSearchGenomeID(tc.name, refs)
		if err != nil || !reflect.DeepEqual(got, tc.want) {
			t.Fatalf("%q: got %v, %v; want %v", tc.name, got, err, tc.want)
		}
	}
	for _, tc := range []struct{ name, message string }{
		{"other:a:b", "ambiguous"}, {"absent", "not found"},
	} {
		if _, err := resolveSearchGenomeID(tc.name, refs); err == nil || !strings.Contains(err.Error(), tc.message) {
			t.Fatalf("%q: expected %s, got %v", tc.name, tc.message, err)
		}
	}
}

func TestSubseqTaxonomicGenomeNames(t *testing.T) {
	dir := t.TempDir()
	dna := indexTestDNA(1000)
	ref := writeTestInput(t, filepath.Join(dir, "a:b.fa"), ">seq\n"+dna+"\n")
	db := filepath.Join(dir, "db")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-c", "1", "--partitions", "4", "-O", db, ref)
	tax := filepath.Join(dir, "taxdump")
	if err := os.Mkdir(tax, 0700); err != nil {
		t.Fatal(err)
	}
	writeTestInput(t, filepath.Join(tax, "nodes.dmp"), "1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tspecies\t|\n")
	writeTestInput(t, filepath.Join(tax, "names.dmp"), "1\t|\troot\t|\t\t|\tscientific name\t|\n2\t|\tSynthetic:species\t|\t\t|\tscientific name\t|\n")
	writeTestInput(t, filepath.Join(tax, "taxid.map"), "a:b\t2\n")
	for _, flag := range []string{"--show-genome-name", "--show-species-name"} {
		result := runSeedIndexCommand(t, "search", "--quiet", "-j", "1", "-n", "0", "-N", "0", "-d", db, "-T", tax, flag, ref)
		if !strings.Contains(string(result), "\tSynthetic:species:a:b\t") {
			t.Fatalf("missing taxonomic prefix: %s", result)
		}
		output, stderr, err := runAuditCommand(t, string(result), "utils", "subseq", "--quiet", "-j", "1", "-d", db, "-f", "-", "-w", "0")
		if err != nil || !strings.HasSuffix(string(output), "\n"+dna+"\n") || strings.Count(string(output), "\n") != 2 {
			t.Fatalf("%s: error=%v stdout=%s stderr=%s", flag, err, output, stderr)
		}
	}
}

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

// TestSubseqInvalidStart rejects contig and concatenated requests without FASTA.
func TestSubseqInvalidStart(t *testing.T) {
	db := writeDetailsFixture(t, 1)
	for _, tc := range []struct {
		seqid  string
		length int
	}{
		{"", 8}, {"contig1", 3}, {"contig2", 5},
	} {
		for _, start := range []int{tc.length + 1, tc.length + 8} {
			args := []string{"utils", "subseq", "--quiet", "-d", db, "-n", "g0", "-s", tc.seqid, "-r", fmt.Sprintf("%d:%d", start, start+20)}
			out, stderr, err := runAuditCommand(t, "", args...)
			want := fmt.Sprintf("start (%d) out of range: [0, %d]", start-1, tc.length-1)
			if err == nil || len(out) != 0 || !strings.Contains(string(stderr), want) || strings.Contains(string(stderr), "panic:") {
				t.Fatalf("seqid=%q start=%d: error=%v stdout=%s stderr=%s", tc.seqid, start, err, out, stderr)
			}
		}
	}
	// A valid start at the final base still permits an oversized end.
	for _, tc := range []struct{ seqid, region, want string }{
		{"contig1", "3:20", ">contig1:3-3:+\nG\n"},
		{"contig2", "5:20", ">contig2:5-5:+\nT\n"},
	} {
		out := runSeedIndexCommand(t, "utils", "subseq", "--quiet", "-d", db, "-n", "g0", "-s", tc.seqid, "-r", tc.region)
		if string(out) != tc.want {
			t.Fatalf("got %s, want %s", out, tc.want)
		}
	}
}

// TestSubseqSearchResultInvalidStart checks fail-fast and skip-record behavior.
func TestSubseqSearchResultInvalidStart(t *testing.T) {
	db := writeDetailsFixture(t, 1)
	row := func(start, end int) string {
		return fmt.Sprintf("q%d\t3\t1\tg0\tcontig1\t100\t1\t1\t100\t3\t100\t0\t1\t3\t%d\t%d\t+\t3\t0\t10\n", start, start, end)
	}
	args := []string{"utils", "subseq", "--quiet", "-j", "2", "-d", db, "-f", "-", "-H"}
	out, stderr, err := runAuditCommand(t, row(4, 9), args...)
	if err == nil || len(out) != 0 || !strings.Contains(string(stderr), "failed to extract subsequence") || strings.Contains(string(stderr), "panic:") {
		t.Fatalf("error=%v stdout=%s stderr=%s", err, out, stderr)
	}
	// Invalid records surround valid reads, including one whose end is clipped.
	out, stderr, err = runAuditCommand(t, row(4, 9)+row(1, 3)+row(6, 10)+row(3, 9), append(args, "--ignore-err")...)
	if err != nil {
		t.Fatalf("--ignore-err failed: %v\n%s", err, stderr)
	}
	lines := strings.Split(strings.TrimSpace(string(out)), "\n")
	if len(lines) != 4 || !strings.HasPrefix(lines[0], ">contig1:1-3:+ query=q1 ") || lines[1] != "ACG" ||
		!strings.HasPrefix(lines[2], ">contig1:3-3:+ query=q3 ") || lines[3] != "G" {
		t.Fatalf("--ignore-err did not retain only valid records in order: %s", out)
	}
}

// TestSubseqExactIDExtractionErrors preserves exact-ID precedence across chunks.
func TestSubseqExactIDExtractionErrors(t *testing.T) {
	db := writeDetailsFixture(t, 2)
	const exact = "c1/1:s1/2:a:b"
	var mapping []byte
	// The earlier chunk has only the longer fallback ID; the exact ID is later.
	for batch, contig := range []struct{ id, sequence string }{{"a:b", "ACGTACGT"}, {exact, "TGA"}} {
		writer, err := genome.NewWriter(filepath.Join(db, DirGenomes, batchDir(batch), FileGenomes), uint32(batch))
		if err != nil {
			t.Fatal(err)
		}
		id := []byte(contig.id)
		g := &genome.Genome{ID: []byte("g0"), Seq: []byte(contig.sequence), GenomeSize: len(contig.sequence),
			Len: len(contig.sequence), NumSeqs: 1, SeqSizes: []int{len(contig.sequence)}, SeqIDs: []*[]byte{&id}}
		if err := writer.Write(g); err != nil {
			t.Fatal(err)
		}
		if err := writer.Close(); err != nil {
			t.Fatal(err)
		}
		mapping = binary.BigEndian.AppendUint16(mapping, 2)
		mapping = append(mapping, "g0"...)
		mapping = binary.BigEndian.AppendUint64(mapping, uint64(batch)<<BITS_GENOME_IDX)
	}
	if err := os.WriteFile(filepath.Join(db, FileGenomeIndex), mapping, 0600); err != nil {
		t.Fatal(err)
	}
	// Direct requests must pass an earlier miss and preserve a later range error.
	args := []string{"utils", "subseq", "--quiet", "-d", db, "-n", "g0", "-s", exact}
	if got := string(runSeedIndexCommand(t, append(args, "-r", "1:2")...)); got != ">"+exact+":1-2:+\nTG\n" {
		t.Fatalf("exact ID in later chunk lost: %s", got)
	}
	out, stderr, err := runAuditCommand(t, "", append(args, "-r", "4:5")...)
	if err == nil || len(out) != 0 || !strings.Contains(string(stderr), "start (3) out of range: [0, 2]") || strings.Contains(string(stderr), "panic:") {
		t.Fatalf("direct exact-ID range error lost: error=%v stdout=%s stderr=%s", err, out, stderr)
	}
	row := func(query, seqid string, start, end int) string {
		return fmt.Sprintf("%s\t3\t1\tg0\t%s\t100\t1\t1\t100\t3\t100\t0\t1\t3\t%d\t%d\t+\t3\t0\t10\n", query, seqid, start, end)
	}
	input := row("before", exact, 1, 2) + row("invalid", exact, 4, 5) + row("after", exact, 3, 20) + row("fallback", "c1/1:s1/3:a:b", 4, 5)
	args = []string{"utils", "subseq", "--quiet", "-j", "2", "-d", db, "-f", "-", "-H"}
	out, stderr, err = runAuditCommand(t, input, args...)
	if err == nil || !strings.Contains(string(stderr), "failed to extract subsequence: "+exact+":4-5:+") || strings.Contains(string(stderr), "panic:") {
		t.Fatalf("search-result exact-ID range error lost: error=%v stdout=%s stderr=%s", err, out, stderr)
	}
	out, stderr, err = runAuditCommand(t, input, append(args, "--ignore-err")...)
	if err != nil {
		t.Fatalf("--ignore-err failed: %v\n%s", err, stderr)
	}
	lines := strings.Split(strings.TrimSpace(string(out)), "\n")
	if len(lines) != 6 || !strings.HasPrefix(lines[0], ">"+exact+":1-2:+ query=before ") || lines[1] != "TG" ||
		!strings.HasPrefix(lines[2], ">"+exact+":3-3:+ query=after ") || lines[3] != "A" ||
		!strings.HasPrefix(lines[4], ">a:b:4-5:+ query=fallback ") || lines[5] != "TA" {
		t.Fatalf("--ignore-err lost exact-ID precedence, fallback or record order: %s", out)
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
