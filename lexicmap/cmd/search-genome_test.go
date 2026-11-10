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
	"bytes"
	"fmt"
	"math/rand"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
)

// Screening uses seeds, so ANI fragment sizes must not discard short queries.
func TestGenomeScreenShortQueries(t *testing.T) {
	dir := t.TempDir()
	rng := rand.New(rand.NewSource(608))
	sequence := make([]byte, 600)
	for i := range sequence {
		sequence[i] = "ACGT"[rng.Intn(4)]
	}
	var refs []string
	for i := range 2 {
		refs = append(refs, writeTestInput(t, filepath.Join(dir, fmt.Sprintf("g%d.fa", i)), ">s\n"+string(sequence)+"\n"))
	}
	db := filepath.Join(dir, "db")
	runSeedIndexCommand(t, append([]string{"index", "--quiet", "-j", "1", "-m", "64", "-k", "31", "-c", "1", "--partitions", "4", "-O", db}, refs...)...)
	var forward []byte
	for _, reverse := range []bool{false, true} {
		q := append([]byte(nil), sequence...)
		if reverse {
			RC(q)
		}
		query := writeTestInput(t, filepath.Join(dir, "query.fa"), ">q\n"+string(q)+"\n")
		var expected []byte
		for _, fragment := range []int{100, 600, 601, 0, 2000} {
			t.Run(fmt.Sprintf("reverse=%t/fragment=%d", reverse, fragment), func(t *testing.T) {
				args := []string{"genome", "search", "--quiet", "-j", "2", "-d", db, "-S", "--extra", "-m", "0", "-p", "11", "-n", "1"}
				if fragment != 0 { // zero denotes the default --frag-size=1020
					args = append(args, "--frag-size", strconv.Itoa(fragment))
				}
				got := runSeedIndexCommand(t, append(args, query)...)
				lines := strings.Split(strings.TrimSpace(string(got)), "\n")
				if len(lines) != 3 { // both identical references tie at Top-1
					t.Fatalf("got %d screening hits, want two", len(lines)-1)
				}
				if expected == nil {
					expected = got
				} else if !bytes.Equal(got, expected) {
					t.Fatal("ANI fragment size changed screening fields or order")
				}
			})
		}
		if !reverse {
			forward = expected
		} else if !bytes.Equal(expected, forward) {
			t.Fatal("reverse complement changed exact-match screening results")
		}
		// The length guard remains in force when ANI is requested.
		got := runSeedIndexCommand(t, "genome", "search", "--quiet", "-j", "1", "-d", db, "-m", "0", "-p", "11", query)
		if bytes.Count(got, []byte("\n")) != 1 {
			t.Fatal("short-query ANI behavior changed")
		}
		// Even when their combined length exceeds k, short contigs cannot seed.
		for _, contigs := range []string{
			">a\n" + string(q[:30]) + "\n",
			">a\n" + string(q[:30]) + "\n>b\n" + string(q[30:60]) + "\n",
		} {
			short := writeTestInput(t, filepath.Join(dir, "short.fa"), contigs)
			got = runSeedIndexCommand(t, "genome", "search", "--quiet", "-j", "1", "-d", db, "-S", "-m", "0", "-p", "11", short)
			if bytes.Count(got, []byte("\n")) != 1 {
				t.Fatal("contigs shorter than k produced screening hits")
			}
		}
	}
}

// Sparse shared seeds must be filtered before ANI work, while -M 0 preserves
// all cutoff ties.
func TestGenomeSearchMinMaskFrac(t *testing.T) {
	dir := t.TempDir()
	rng := rand.New(rand.NewSource(9012))
	querySeq, refSeq := make([]byte, 16000), make([]byte, 16000)
	for i := range querySeq {
		querySeq[i], refSeq[i] = "ACGT"[rng.Intn(4)], "ACGT"[rng.Intn(4)]
	}
	copy(refSeq[8000:8128], querySeq[5000:5128])
	query := writeTestInput(t, filepath.Join(dir, "query.fa"), ">q\n"+string(querySeq)+"\n")
	var refs []string
	for i := range 2 {
		refs = append(refs, writeTestInput(t, filepath.Join(dir, fmt.Sprintf("ref%d.fa", i)), ">s\n"+string(refSeq)+"\n"))
	}
	db := filepath.Join(dir, "db")
	runSeedIndexCommand(t, append([]string{"index", "--quiet", "-j", "2", "-m", "4096", "-c", "2", "--partitions", "4", "-O", db}, refs...)...)
	args := []string{"genome", "search", "--quiet", "-j", "2", "-d", db, "-p", "21", "-n", "1", query}
	baseline := runSeedIndexCommand(t, append(args, "-S", "-M", "0")...)
	lines := strings.Split(strings.TrimSpace(string(baseline)), "\n")
	if len(lines) != 3 {
		t.Fatalf("expected two sparse cutoff ties, got %s", baseline)
	}
	hits, err := strconv.Atoi(strings.Split(lines[1], "\t")[4])
	if err != nil || hits < 1 || hits >= 41 {
		t.Fatalf("fixture needs 1-40 matched masks, got %d (%v)", hits, err)
	}
	for _, mode := range [][]string{{"-S"}, {}, {"-O"}} {
		for _, cutoff := range [][]string{nil, {"--min-mask-frac", "0.01"}} {
			t.Run(fmt.Sprint(mode, cutoff), func(t *testing.T) {
				command := append(append(append([]string(nil), args...), mode...), cutoff...)
				got := runSeedIndexCommand(t, command...)
				if bytes.Count(got, []byte("\n")) != 1 {
					t.Fatalf("sparse candidates passed the mask fraction cutoff: %s", got)
				}
			})
		}
	}
	// Compare the unrounded fraction: a threshold exactly on the boundary keeps
	// both ties, even when its four-decimal representation is rounded upward.
	fraction := strconv.FormatFloat(float64(hits)/4096, 'g', -1, 64)
	got := runSeedIndexCommand(t, append(args, "-S", "-M", fraction)...)
	if !bytes.Equal(got, baseline) {
		t.Fatalf("boundary cutoff changed screening results:\n%s\n%s", got, baseline)
	}
}

func TestGenomeSearchRejectsInvalidMinMaskFrac(t *testing.T) {
	for _, fraction := range []string{"-0.01", "1.01", "NaN", "+Inf", "-Inf"} {
		t.Run(fraction, func(t *testing.T) {
			out := filepath.Join(t.TempDir(), "out.tsv")
			command := exec.Command(os.Args[0], "-test.run=^TestSeedIndexCommand$", "--",
				"genome", "search", "--quiet", "-d", "unused-index", "-M", fraction, "-o", out)
			command.Env = append(os.Environ(), "LEXICMAP_SEED_INDEX_TEST_HELPER=1")
			data, err := command.CombinedOutput()
			if err == nil || !bytes.Contains(data, []byte("min-mask-frac")) {
				t.Fatalf("invalid fraction %s: error=%v, output=%s", fraction, err, data)
			}
			if _, err := os.Stat(out); !os.IsNotExist(err) {
				t.Fatalf("invalid fraction created output: %v", err)
			}
		})
	}
}
