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
	"math/rand"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
)

// TestGenomeANIOutputPrecisionAndTerminalAF checks output precision and the
// whole-genome AF denominator when a terminal fragment cannot align.
func TestGenomeANIOutputPrecisionAndTerminalAF(t *testing.T) {
	dir := t.TempDir()
	rng := rand.New(rand.NewSource(1))
	sequence := make([]byte, 4980) // four full fragments plus a 900-bp tail
	for i := range sequence {
		sequence[i] = "ACGT"[rng.Intn(4)]
	}
	query := filepath.Join(dir, "query.fa")
	reference := filepath.Join(dir, "reference.fa")
	fastq := filepath.Join(dir, "query.fq")
	if err := os.WriteFile(query, append([]byte(">seq\n"), sequence...), 0600); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(fastq, []byte("@seq\n"+string(sequence)+"\n+\n"+strings.Repeat("I", len(sequence))+"\n"), 0600); err != nil {
		t.Fatal(err)
	}
	// A substitution gives a non-integer ANI so the test also covers rounding.
	sequence[1040] = "ACGT"[(strings.IndexByte("ACGT", sequence[1040])+1)%4]
	// The query's 900-bp tail has no counterpart with usable seeds in the
	// reference, so AF must be close to 4080/4980 instead of 4080/4080.
	for i := 4080; i < len(sequence); i++ {
		sequence[i] = 'N'
	}
	if err := os.WriteFile(reference, append([]byte(">seq\n"), sequence...), 0600); err != nil {
		t.Fatal(err)
	}
	db := filepath.Join(dir, "index.lmi")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-c", "2", "-m", "64", "-k", "16", "--partitions", "4", "-O", db, reference)
	for _, command := range []struct {
		name string
		args []string
	}{
		{"search", []string{"genome", "search", "-d", db, "-m", "0", "-p", "12", query}},
		{"search OrthoANI", []string{"genome", "search", "-O", "-d", db, "-m", "0", "-p", "12", query}},
		{"search adjusted ANI", []string{"genome", "search", "--adjust-ani-by-quality", "-d", db, "-m", "0", "-p", "12", fastq}},
		{"compare", []string{"genome", "compare", query, reference}},
		{"compare OrthoANI", []string{"genome", "compare", "-O", query, reference}},
	} {
		t.Run(command.name, func(t *testing.T) {
			output := runSeedIndexCommand(t, append(command.args, "--quiet", "-j", "1")...)
			lines := strings.Split(strings.TrimSuffix(string(output), "\n"), "\n")
			if len(lines) < 2 {
				t.Fatalf("expected genome matches, got %s", output)
			}
			header := strings.Split(lines[0], "\t")
			for _, line := range lines[1:] {
				fields := strings.Split(line, "\t")
				if len(fields) != len(header) {
					t.Fatalf("unexpected column count: %s", line)
				}
				for i, name := range header {
					decimals := 0
					if strings.Contains(name, "ANI") {
						decimals = 4
					} else if strings.Contains(name, "AF") {
						decimals = 3
					} else {
						continue
					}
					parts := strings.Split(fields[i], ".")
					if len(parts) != 2 || len(parts[1]) != decimals {
						t.Fatalf("%s: got %q, want %d decimal places", name, fields[i], decimals)
					}
					if strings.Contains(name, "AF") {
						value, err := strconv.ParseFloat(fields[i], 64)
						if err != nil || value < 75 || value > 90 {
							t.Fatalf("%s: got %q, want about 82%% using the full genome size", name, fields[i])
						}
					}
				}
			}
		})
	}
}

// TestGenomeANISeedlessTailSelfAF checks that a sub-k terminal piece is
// aligned with its preceding fragment in both search and compare modes.
func TestGenomeANISeedlessTailSelfAF(t *testing.T) {
	dir := t.TempDir()
	rng := rand.New(rand.NewSource(2))
	sequence := make([]byte, 4*1020+7)
	for i := range sequence {
		sequence[i] = "ACGT"[rng.Intn(4)]
	}
	query := filepath.Join(dir, "query.fa")
	reference := filepath.Join(dir, "reference.fa")
	data := append([]byte(">seq\n"), sequence...)
	for _, file := range []string{query, reference} {
		if err := os.WriteFile(file, data, 0600); err != nil {
			t.Fatal(err)
		}
	}
	db := filepath.Join(dir, "index.lmi")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-c", "2", "-m", "64", "-k", "16", "--partitions", "4", "-O", db, reference)
	for _, command := range []struct {
		name string
		args []string
	}{
		{"search", []string{"genome", "search", "-d", db, "-m", "0", "-p", "12", query}},
		{"search OrthoANI", []string{"genome", "search", "-O", "-d", db, "-m", "0", "-p", "12", query}},
		{"compare", []string{"genome", "compare", query, reference}},
		{"compare OrthoANI", []string{"genome", "compare", "-O", query, reference}},
	} {
		t.Run(command.name, func(t *testing.T) {
			output := runSeedIndexCommand(t, append(command.args, "--quiet", "-j", "1")...)
			lines := strings.Split(strings.TrimSpace(string(output)), "\n")
			if len(lines) != 2 {
				t.Fatalf("expected one genome match, got %s", output)
			}
			header, fields := strings.Split(lines[0], "\t"), strings.Split(lines[1], "\t")
			for i, name := range header {
				if strings.Contains(name, "AF") && fields[i] != "100.000" {
					t.Fatalf("%s: got %s, want 100.000 with a 7-bp tail", name, fields[i])
				}
			}
		})
	}
}

// TestGenomeANIExactFragmentEndCoverage catches a self-matching fragment
// whose heuristic extension previously stopped before its exact right end.
func TestGenomeANIExactFragmentEndCoverage(t *testing.T) {
	genome, err := ReadGenomeFromFile("../../demo/refs/GCF_000006945.2.fa.gz", nil, false)
	if err != nil {
		t.Fatal(err)
	}
	defer RecycleGQuery(genome)
	const start = 1114 * 1020
	seq := *genome.seqs[0]
	dir := t.TempDir()
	query := filepath.Join(dir, "query.fa")
	subject := filepath.Join(dir, "subject.fa")
	for file, fragment := range map[string][]byte{
		query:   seq[start : start+1020],
		subject: seq[start-100 : start+1120],
	} {
		if err := os.WriteFile(file, append([]byte(">seq\n"), fragment...), 0600); err != nil {
			t.Fatal(err)
		}
	}
	output := runSeedIndexCommand(t, "genome", "compare", query, subject, "--quiet", "-j", "1")
	lines := strings.Split(strings.TrimSpace(string(output)), "\n")
	if len(lines) != 2 {
		t.Fatalf("expected one genome match, got %s", output)
	}
	fields := strings.Split(lines[1], "\t")
	if fields[5] != "100.000" {
		t.Fatalf("AF1: got %s, want 100.000 for an exact query fragment", fields[5])
	}
}

// TestGenomeANIExactRepetitiveSelf checks that sampled exact hits survive
// repetitive-region chaining in both genome search and genome compare.
func TestGenomeANIExactRepetitiveSelf(t *testing.T) {
	reference := "../../demo/refs/GCF_002949675.1.fa.gz"
	db := filepath.Join(t.TempDir(), "index.lmi")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-c", "2", "-m", "64", "-k", "16", "--partitions", "4", "-O", db, reference)

	for _, tc := range []struct {
		name string
		args []string
		want []string
	}{
		{"search", []string{"genome", "search", "-d", db, "-m", "0", "-p", "12", reference, "--quiet", "-j", "1"}, []string{"ANI", "qAF", "sAF"}},
		{"compare", []string{"genome", "compare", reference, reference, "--quiet", "-j", "1"}, []string{"tANI", "ANI1", "ANI2", "AF1", "AF2"}},
	} {
		t.Run(tc.name, func(t *testing.T) {
			lines := strings.Split(strings.TrimSpace(string(runSeedIndexCommand(t, tc.args...))), "\n")
			if len(lines) != 2 {
				t.Fatalf("expected one self match, got %q", lines)
			}
			header, fields := strings.Split(lines[0], "\t"), strings.Split(lines[1], "\t")
			for _, name := range tc.want {
				want := "100.000"
				if strings.Contains(name, "ANI") {
					want = "100.0000"
				}
				i := -1
				for j, h := range header {
					if h == name {
						i = j
						break
					}
				}
				if i < 0 || i >= len(fields) || fields[i] != want {
					t.Fatalf("%s: got %q, want %s", name, lines[1], want)
				}
			}
		})
	}
}

// TestGenomeANIMutatedRepetitiveFragment checks that a renamed, imperfect
// query uses its actual alignment instead of an exact self-hit shortcut.
func TestGenomeANIMutatedRepetitiveFragment(t *testing.T) {
	reference := "../../demo/refs/GCF_002949675.1.fa.gz"
	genome, err := ReadGenomeFromFile(reference, nil, false)
	if err != nil {
		t.Fatal(err)
	}
	defer RecycleGQuery(genome)
	const start = 1297 * 1020
	fragment := append([]byte(nil), (*genome.seqs[0])[start:start+1020]...)
	for i := 50; i < len(fragment); i += 100 {
		fragment[i] = "ACGT"[(strings.IndexByte("ACGT", fragment[i])+1)%4]
	}
	query := filepath.Join(t.TempDir(), "unknown.fa")
	if err := os.WriteFile(query, append([]byte(">unknown\n"), fragment...), 0600); err != nil {
		t.Fatal(err)
	}
	output := runSeedIndexCommand(t, "genome", "compare", query, reference, "--quiet", "-j", "1")
	lines := strings.Split(strings.TrimSpace(string(output)), "\n")
	if len(lines) != 2 {
		t.Fatalf("expected one match, got %s", output)
	}
	header, fields := strings.Split(lines[0], "\t"), strings.Split(lines[1], "\t")
	for i, name := range header {
		if name == "AF1" && fields[i] != "100.000" {
			t.Fatalf("AF1: got %s, want 100.000", fields[i])
		}
		if name == "ANI1" {
			ani, err := strconv.ParseFloat(fields[i], 64)
			if err != nil || ani < 98 || ani >= 100 {
				t.Fatalf("ANI1: got %s, want an imperfect alignment", fields[i])
			}
		}
	}
}

// TestGenomeANIShortSampledTailRC checks that a terminal piece long enough
// for k-mers but too short to yield sampled seeds stays covered after RC.
func TestGenomeANIShortSampledTailRC(t *testing.T) {
	rng := rand.New(rand.NewSource(3))
	sequence := make([]byte, 2*1020+41)
	for i := range sequence {
		sequence[i] = "ACGT"[rng.Intn(4)]
	}
	dir := t.TempDir()
	reference := filepath.Join(dir, "reference.fa")
	query := filepath.Join(dir, "unknown-rc.fa")
	if err := os.WriteFile(reference, append([]byte(">reference\n"), sequence...), 0600); err != nil {
		t.Fatal(err)
	}
	RC(sequence)
	if err := os.WriteFile(query, append([]byte(">unknown\n"), sequence...), 0600); err != nil {
		t.Fatal(err)
	}
	db := filepath.Join(dir, "index.lmi")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-c", "2", "-m", "64", "-k", "16", "--partitions", "4", "-O", db, reference)
	for _, tc := range []struct {
		name string
		args []string
	}{
		{"search", []string{"genome", "search", "-d", db, "-m", "0", "-p", "12", query, "--quiet", "-j", "1"}},
		{"compare", []string{"genome", "compare", query, reference, "--quiet", "-j", "1"}},
	} {
		t.Run(tc.name, func(t *testing.T) {
			output := runSeedIndexCommand(t, tc.args...)
			lines := strings.Split(strings.TrimSpace(string(output)), "\n")
			if len(lines) != 2 {
				t.Fatalf("expected one match, got %s", output)
			}
			header, fields := strings.Split(lines[0], "\t"), strings.Split(lines[1], "\t")
			for i, name := range header {
				if strings.Contains(name, "ANI") && fields[i] != "100.0000" || strings.Contains(name, "AF") && fields[i] != "100.000" {
					t.Fatalf("%s: got %s, want 100%%", name, fields[i])
				}
			}
		})
	}
}
