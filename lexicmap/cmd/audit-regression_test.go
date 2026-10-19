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
	"bytes"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// runAuditCommand isolates CLI exits and global command state, including failures.
func runAuditCommand(t *testing.T, stdin string, args ...string) ([]byte, []byte, error) {
	t.Helper()
	command := exec.Command(os.Args[0], append([]string{"-test.run=^TestSeedIndexCommand$", "--"}, args...)...)
	command.Env = append(os.Environ(), "LEXICMAP_SEED_INDEX_TEST_HELPER=1")
	command.Stdin = strings.NewReader(stdin)
	var out, stderr bytes.Buffer
	command.Stdout, command.Stderr = &out, &stderr
	err := command.Run()
	return out.Bytes(), stderr.Bytes(), err
}

func TestEditGenomeIDsNoChange(t *testing.T) {
	for _, backup := range []bool{false, true} {
		for _, pattern := range []string{"missing", "^g1$"} {
			t.Run(fmt.Sprintf("backup=%v/pattern=%s", backup, pattern), func(t *testing.T) {
				dir := t.TempDir()
				file := filepath.Join(dir, FileGenomeIndex)
				var original bytes.Buffer
				binary.Write(&original, binary.BigEndian, uint16(2))
				original.WriteString("g1")
				binary.Write(&original, binary.BigEndian, uint64(7))
				if err := os.WriteFile(file, original.Bytes(), 0600); err != nil {
					t.Fatal(err)
				}
				if backup {
					if err := os.WriteFile(file+".bak", []byte("older backup"), 0600); err != nil {
						t.Fatal(err)
					}
				}
				cache := filepath.Join(dir, FileGenomeDetails)
				if err := os.WriteFile(cache, []byte("cached details"), 0600); err != nil {
					t.Fatal(err)
				}
				runSeedIndexCommand(t, "utils", "edit-genome-ids", "--quiet", "-d", dir, "-p", pattern, "-r", "g1")
				got, err := os.ReadFile(file)
				if err != nil || !bytes.Equal(got, original.Bytes()) {
					t.Fatalf("original changed: %q, %v", got, err)
				}
				got, err = os.ReadFile(file + ".bak")
				if backup {
					if err != nil || string(got) != "older backup" {
						t.Fatalf("backup changed: %q, %v", got, err)
					}
				} else if !os.IsNotExist(err) {
					t.Fatalf("unexpected backup: %v", err)
				}
				if got, err := os.ReadFile(cache); err != nil || string(got) != "cached details" {
					t.Fatalf("no-op edit changed cached details: %q, %v", got, err)
				}
			})
		}
	}
}

func TestEditGenomeIDsInvalidatesDetails(t *testing.T) {
	dir := t.TempDir()
	input := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(input, []byte(">seq\n"+strings.Repeat("ACGT", 100)), 0600); err != nil {
		t.Fatal(err)
	}
	index := filepath.Join(dir, "db.lmi")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-k", "16", "--partitions", "4", "-O", index, input)
	oldID := "ref"
	for _, newID := range []string{"renamed", "renamed-again"} {
		before := runSeedIndexCommand(t, "utils", "genome-details", "--quiet", "-j", "1", "-i", "-e", "-d", index)
		if !bytes.Contains(before, []byte("\n"+oldID+"\t")) {
			t.Fatalf("unexpected initial details: %s", before)
		}
		runSeedIndexCommand(t, "utils", "edit-genome-ids", "--quiet", "-d", index, "-p", "^"+oldID+"$", "-r", newID)
		if _, err := os.Stat(filepath.Join(index, FileGenomeDetails)); !os.IsNotExist(err) {
			t.Fatalf("stale details retained after rename: %v", err)
		}
		after := runSeedIndexCommand(t, "utils", "genome-details", "--quiet", "-j", "1", "-i", "-e", "-d", index)
		if !bytes.Equal(after, bytes.ReplaceAll(before, []byte("\n"+oldID+"\t"), []byte("\n"+newID+"\t"))) {
			t.Fatalf("rebuilt details disagree with renamed mapping: %s", after)
		}
		oldID = newID
	}
}

func TestEditGenomeIDsInvalidationFailurePreservesMapping(t *testing.T) {
	dir := t.TempDir()
	file := filepath.Join(dir, FileGenomeIndex)
	var original bytes.Buffer
	binary.Write(&original, binary.BigEndian, uint16(2))
	original.WriteString("g1")
	binary.Write(&original, binary.BigEndian, uint64(0))
	if err := os.WriteFile(file, original.Bytes(), 0600); err != nil {
		t.Fatal(err)
	}
	// A nonempty directory makes cache removal fail without relying on permissions.
	cache := filepath.Join(dir, FileGenomeDetails)
	if err := os.Mkdir(cache, 0700); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(cache, "keep"), []byte("keep"), 0600); err != nil {
		t.Fatal(err)
	}
	_, stderr, err := runAuditCommand(t, "", "utils", "edit-genome-ids", "--quiet", "-d", dir, "-p", "g1", "-r", "g2")
	if err == nil || !bytes.Contains(stderr, []byte("failed to invalidate genome details")) {
		t.Fatalf("invalidation failure ignored: %v, %s", err, stderr)
	}
	got, err := os.ReadFile(file)
	if err != nil || !bytes.Equal(got, original.Bytes()) {
		t.Fatalf("mapping changed despite invalidation failure: %q, %v", got, err)
	}
	if _, err := os.Stat(file + ".bak"); !os.IsNotExist(err) {
		t.Fatalf("backup changed despite invalidation failure: %v", err)
	}
}

// failingOutput checks both flush errors and errors reported only on close.
type failingOutput struct {
	writeErr, closeErr error
	closed             bool
}

func (w *failingOutput) Write(p []byte) (int, error) {
	if w.writeErr != nil {
		return 0, w.writeErr
	}
	return len(p), nil
}
func (w *failingOutput) Close() error { w.closed = true; return w.closeErr }

func TestCloseOutStreamErrors(t *testing.T) {
	flushErr, closeErr := errors.New("flush failed"), errors.New("compression close failed")
	w := &failingOutput{writeErr: flushErr, closeErr: closeErr}
	buffered := bufio.NewWriter(w)
	buffered.WriteString("buffered result")
	file, err := os.CreateTemp(t.TempDir(), "output")
	if err != nil {
		t.Fatal(err)
	}
	err = closeOutStream(buffered, w, file)
	if !errors.Is(err, flushErr) || !errors.Is(err, closeErr) || !w.closed {
		t.Fatalf("lost error/cleanup: %v, closed=%v", err, w.closed)
	}
	if _, err = file.WriteString("closed"); !errors.Is(err, os.ErrClosed) {
		t.Fatalf("file not closed: %v", err)
	}
	// Some compressors discover a disk error only while writing their footer.
	w = &failingOutput{closeErr: closeErr}
	if err = closeOutStream(bufio.NewWriter(w), w, nil); !errors.Is(err, closeErr) {
		t.Fatalf("lost footer error: %v", err)
	}
}

func TestIndexIdentifierLimits(t *testing.T) {
	maxBatches, maxRecords := 1<<BITS_BATCH_IDX, 1<<BITS_GENOME_IDX
	for _, tt := range []struct {
		files, size, want int
		bad               bool
	}{
		{0, 1, 0, false}, {maxBatches, 1, maxBatches, false}, {maxBatches + 1, 1, 0, true},
		{maxBatches * 2, 2, maxBatches, false}, {maxBatches*2 + 1, 2, 0, true}, {1, 0, 0, true},
	} {
		got, err := indexBatchCount(tt.files, tt.size)
		if (err != nil) != tt.bad || (!tt.bad && got != tt.want) {
			t.Fatalf("batch count %+v: %d, %v", tt, got, err)
		}
	}
	for _, tt := range []struct {
		batch  int
		record uint64
		bad    bool
	}{
		{0, 0, false}, {maxBatches - 1, uint64(maxRecords - 1), false}, {maxBatches, 0, true}, {0, uint64(maxRecords), true}, {-1, 0, true},
	} {
		got, err := packGenomeIndex(tt.batch, tt.record)
		if (err != nil) != tt.bad {
			t.Fatalf("pack %+v: %d, %v", tt, got, err)
		}
		if !tt.bad {
			seed := got << BITS_NONE_IDX
			decoded := seed >> BITS_NONE_IDX
			if decoded>>BITS_GENOME_IDX != uint64(tt.batch) || decoded&MASK_GENOME_IDX != tt.record {
				t.Fatalf("identifier lost bits: %+v", tt)
			}
		}
	}
	// Exercise the real build entry point: oversized batches must fail before
	// masks, file I/O or fixed build buffers are allocated.
	err := BuildIndex(filepath.Join(t.TempDir(), "index"), make([]string, maxBatches+1), &IndexBuildingOptions{GenomeBatchSize: 1})
	if err == nil || !strings.Contains(err.Error(), "batches") {
		t.Fatalf("missing early batch validation: %v", err)
	}
}

const auditSearchHeader = "query\tqlen\thits\tsgenome\tsseqid\tqcovGnm\tcls\thsp\tqcovHSP\talenHSP\tpident\tgaps\tqstart\tqend\tsstart\tsend\tsstr\tslen\tevalue\tbitscore"

func TestMergeSingleFileQuerySelection(t *testing.T) {
	dir := t.TempDir()
	file := filepath.Join(dir, "input.tsv")
	row := func(q string) string {
		return q + "\t100\t1\tg\ts\t100\t1\t1\t100\t100\t100\t0\t1\t100\t1\t100\t+\t100\t0\t100\n"
	}
	if err := os.WriteFile(file, []byte(auditSearchHeader+"\n"+row("q1")+row("q2")), 0600); err != nil {
		t.Fatal(err)
	}
	for _, query := range []string{"q1", "q2", "missing"} {
		out := runSeedIndexCommand(t, "utils", "merge-search-results", "--quiet", "-q", query, file)
		lines := strings.Split(strings.TrimSpace(string(out)), "\n")
		if query == "missing" {
			if len(out) != 0 {
				t.Fatalf("missing query returned data: %s", out)
			}
			continue
		}
		if len(lines) != 2 || strings.Split(lines[1], "\t")[0] != query {
			t.Fatalf("wrong selected query: %s", out)
		}
	}
}

func TestSAMReferenceOrientation(t *testing.T) {
	file := filepath.Join(t.TempDir(), "input.tsv")
	var input strings.Builder
	input.WriteString(auditSearchHeader + "\tcigar\tqseq\tsseq\talign\n")
	for _, q := range []struct{ id, strand string }{{"forward", "+"}, {"reverse", "-"}, {"reverse", "-"}} {
		fmt.Fprintf(&input, "%s\t12\t1\tg\ts\t58.3\t1\t1\t58.3\t8\t75\t2\t3\t9\t10\t16\t%s\t100\t0\t100\t2M1D3M1I1M\tAC-GTACG\tACTGTA-G\t|| ||| |\n", q.id, q.strand)
	}
	if err := os.WriteFile(file, []byte(input.String()), 0600); err != nil {
		t.Fatal(err)
	}
	out := runSeedIndexCommand(t, "utils", "2sam", "--quiet", file)
	var records [][]string
	for _, line := range strings.Split(strings.TrimSpace(string(out)), "\n") {
		if !strings.HasPrefix(line, "@") {
			records = append(records, strings.Split(line, "\t"))
		}
	}
	if len(records) != 3 {
		t.Fatalf("unexpected SAM: %s", out)
	}
	for i, want := range []struct{ flag, cigar, seq string }{
		{"0", "2S2M1D3M1I1M3S", "NNACGTACGNNN"},
		{"16", "3S1M1I3M1D2M2S", "NNNCGTACGTNN"},
		{"272", "3S1M1I3M1D2M2S", "*"},
	} {
		r := records[i]
		if r[1] != want.flag || r[3] != "10" || r[5] != want.cigar || r[9] != want.seq || r[11] != "NM:i:2" {
			t.Fatalf("SAM record %d: %v", i, r)
		}
	}
}

func TestSAMMappingQualityAtQueryBoundaryAndEOF(t *testing.T) {
	file := filepath.Join(t.TempDir(), "input.tsv")
	var input strings.Builder
	input.WriteString(auditSearchHeader + "\tcigar\tqseq\tsseq\talign\n")
	for _, query := range []string{"first", "last"} {
		for i, score := range []int{100, 10, 90} {
			// Two HSPs from the best genome, then the second-best genome.
			genome := "g1"
			if i == 2 {
				genome = "g2"
			}
			fmt.Fprintf(&input, "%s\t100\t2\t%s\ts\t100\t1\t%d\t100\t100\t100\t0\t1\t100\t1\t100\t+\t100\t0\t%d\t100M\t%s\t%s\t%s\n", query, genome, i+1, score, strings.Repeat("A", 100), strings.Repeat("A", 100), strings.Repeat("|", 100))
		}
	}
	if err := os.WriteFile(file, []byte(input.String()), 0600); err != nil {
		t.Fatal(err)
	}
	out := runSeedIndexCommand(t, "utils", "2sam", "--quiet", file)
	primary := 0
	for _, line := range strings.Split(strings.TrimSpace(string(out)), "\n") {
		if strings.HasPrefix(line, "@") {
			continue
		}
		r := strings.Split(line, "\t")
		if r[1] == "0" {
			primary++
			if r[4] != "4" {
				t.Fatalf("wrong MAPQ: %s", line)
			}
		} else if r[1] != "256" {
			t.Fatalf("wrong secondary flag: %s", line)
		}
	}
	if primary != 2 {
		t.Fatalf("primary alignments: %d", primary)
	}
}

// Check delayed errors without relying on device files on non-Linux platforms.
var _ io.WriteCloser = (*failingOutput)(nil)
