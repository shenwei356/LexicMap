//go:build linux

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
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
)

func TestSearchResourceFileBudget(t *testing.T) {
	db, query := buildMultiBatchSearchFixture(t)
	for _, command := range [][]string{{"search", "-P", "12", "-n", "0", "-N", "0"}, {"genome", "search", "-m", "0", "-n", "10"}} {
		args := append(append([]string{}, command...), "--quiet", "-j", "16", "-J", "32", "-d", db, "-p", "12", query)
		want := runSeedIndexCommand(t, args...)
		if bytes.Count(want, []byte("\n")) < 2 {
			t.Fatalf("fixture produced no hits: %s", want)
		}
		for _, budget := range []int{32, 15, 11} {
			limited := append(append([]string{}, args...), "--max-open-files", strconv.Itoa(budget))
			if budget == 11 {
				limited = append(limited, "-w")
			}
			got, stderr, err := runResourceLimitedCommand(t, "open-files="+strconv.Itoa(budget), limited...)
			if err != nil {
				t.Fatalf("budget %d: %v, %s", budget, err, stderr)
			}
			if !bytes.Equal(got, want) {
				t.Fatalf("budget %d changed output:\ngot: %s\nwant: %s", budget, got, want)
			}
		}
	}
}

// TestSearchReaderDescriptorAccounting checks real handles, including the
// in-memory path and the minimum one-reader pool for every genome batch.
func TestSearchReaderDescriptorAccounting(t *testing.T) {
	db, _ := buildMultiBatchSearchFixture(t)
	count := func() int {
		t.Helper()
		entries, err := os.ReadDir("/proc/self/fd")
		if err != nil {
			t.Fatal(err)
		}
		n := 0
		for _, entry := range entries {
			target, err := os.Readlink(filepath.Join("/proc/self/fd", entry.Name()))
			if err == nil && strings.HasPrefix(target, db+string(os.PathSeparator)) {
				n++
			}
		}
		return n
	}
	for _, budget := range []int{32, 15, 11} {
		opt := DefaultIndexSearchingOptions
		opt.NumCPUs, opt.MaxOpenFiles, opt.MaxSeedSearchingConcurrency = 16, budget, 16
		opt.MinPrefix, opt.MinSinglePrefix = 12, 12
		opt.InMemorySearch = budget == 11
		idx, err := NewIndexSearcher(db, &opt)
		if err != nil {
			t.Fatal(err)
		}
		if !idx.hasGenomeRdrs {
			t.Fatalf("budget %d: missing persistent genome readers", budget)
		}
		check := func() {
			t.Helper()
			if got := count(); got != len(idx.openFileTokens) || got > budget-reservedOpenFiles {
				t.Fatalf("budget %d: %d descriptors, %d tokens", budget, got, len(idx.openFileTokens))
			}
		}
		check()
		for batch := 0; batch < 2; batch++ {
			rdr, err := idx.acquireGenomeReader(batch)
			if err != nil {
				t.Fatal(err)
			}
			check()
			if err := idx.releaseGenomeReader(batch, rdr); err != nil {
				t.Fatal(err)
			}
			check()
		}
		if err := idx.Close(); err != nil {
			t.Fatal(err)
		}
		if got := count(); got != 0 {
			t.Fatalf("Close leaked %d descriptors", got)
		}
	}
}

func TestSearchRecountFileBudget(t *testing.T) {
	db, _ := buildMultiBatchSearchFixture(t)
	file := filepath.Join(db, FileInfo)
	info, err := readIndexInfo(file)
	if err != nil {
		t.Fatal(err)
	}
	want := info.InputBases
	info.InputBases = 0
	if err := writeIndexInfo(file, info); err != nil {
		t.Fatal(err)
	}
	budget := reservedOpenFiles + info.Chunks + info.GenomeBatches
	logFile := filepath.Join(t.TempDir(), "search.log")
	_, stderr, err := runResourceLimitedCommand(t, "open-files="+strconv.Itoa(budget), "search", "--quiet", "-j", "64", "-d", db,
		"--max-open-files", strconv.Itoa(budget), "--log", logFile, "-o", filepath.Join(t.TempDir(), "out.tsv"), "/dev/null")
	if err != nil {
		log, _ := os.ReadFile(logFile)
		t.Fatalf("recount exceeded the budget: %v, %s\n%s", err, stderr, log)
	}
	info, err = readIndexInfo(file)
	if err != nil || info.InputBases != want {
		t.Fatalf("incorrect recount: %+v, %v; want %d", info, err, want)
	}
}
