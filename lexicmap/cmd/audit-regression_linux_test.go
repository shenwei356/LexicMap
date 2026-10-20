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
	"context"
	"math/rand"
	"os"
	"os/exec"
	"os/signal"
	"path/filepath"
	"runtime/debug"
	"strconv"
	"strings"
	"syscall"
	"testing"
	"time"
)

// Do not rely on os.File finalizers to release a batch before the next opens.
func TestUpdateInputBasesClosesReaders(t *testing.T) {
	previous := debug.SetGCPercent(-1)
	defer debug.SetGCPercent(previous)
	dir := t.TempDir()
	info := writeRecountFixture(t, dir, 48)
	for _, broken := range []bool{false, true} {
		if broken {
			if err := os.Truncate(filepath.Join(dir, DirGenomes, batchDir(20), FileGenomes), 0); err != nil {
				t.Fatal(err)
			}
		}
		before, err := os.ReadDir("/proc/self/fd")
		if err != nil {
			t.Fatal(err)
		}
		_, err = updateInputBases(info, dir, 3)
		if (err != nil) != broken {
			t.Fatalf("broken=%v: unexpected recount error: %v", broken, err)
		}
		after, err := os.ReadDir("/proc/self/fd")
		if err != nil || len(after) != len(before) {
			t.Fatalf("broken=%v: open files grew from %d to %d: %v", broken, len(before), len(after), err)
		}
	}
}

func TestCommandOutputWriteFailures(t *testing.T) {
	dir := t.TempDir()
	input := filepath.Join(dir, "input.tsv")
	if err := os.WriteFile(input, []byte("buffered result\n"), 0600); err != nil {
		t.Fatal(err)
	}
	for _, compressed := range []bool{false, true} {
		output := "/dev/full"
		if compressed {
			output = filepath.Join(dir, "full.gz")
			if err := os.Symlink("/dev/full", output); err != nil {
				t.Fatal(err)
			}
		}
		_, stderr, err := runAuditCommand(t, "", "utils", "merge-search-results", "--quiet", input, "-o", output)
		if err == nil || !bytes.Contains(stderr, []byte("no space left")) {
			t.Fatalf("lost output error: %v, %s", err, stderr)
		}
	}
	// Verify the real search output path, including the header-only case.
	rng := rand.New(rand.NewSource(91))
	sequence := make([]byte, 600)
	for i := range sequence {
		sequence[i] = "ACGT"[rng.Intn(4)]
	}
	ref := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(ref, append([]byte(">seq\n"), sequence...), 0600); err != nil {
		t.Fatal(err)
	}
	db := filepath.Join(dir, "db")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-c", "2", "-m", "64", "-k", "16", "--partitions", "4", "-O", db, ref)
	for _, query := range []string{">match\n" + string(sequence), ">short\nACGT\n"} {
		_, stderr, err := runAuditCommand(t, query, "search", "--quiet", "-j", "1", "-d", db, "-p", "12", "-P", "12", "-n", "0", "-N", "0", "-o", "/dev/full")
		if err == nil || !bytes.Contains(stderr, []byte("no space left")) {
			t.Fatalf("lost search output error: %v, %s", err, stderr)
		}
	}
}

func TestRedirectedStdinOutputCollision(t *testing.T) {
	file := filepath.Join(t.TempDir(), "input.tsv")
	data := []byte("original input\n")
	if err := os.WriteFile(file, data, 0600); err != nil {
		t.Fatal(err)
	}
	input, err := os.Open(file)
	if err != nil {
		t.Fatal(err)
	}
	defer input.Close()
	command := exec.Command(os.Args[0], "-test.run=^TestSeedIndexCommand$", "--", "utils", "merge-search-results", "--quiet", "-o", file)
	command.Env = append(os.Environ(), "LEXICMAP_SEED_INDEX_TEST_HELPER=1")
	command.Stdin = input
	out, err := command.CombinedOutput()
	if err == nil || !strings.Contains(string(out), "conflicts with stdin") {
		t.Fatalf("stdin alias accepted: %v, %s", err, out)
	}
	got, err := os.ReadFile(file)
	if err != nil || !bytes.Equal(got, data) {
		t.Fatalf("stdin overwritten: %q, %v", got, err)
	}
}

// TestResourceLimitedCommand applies limits only in a child, after Go startup.
func TestResourceLimitedCommand(t *testing.T) {
	kind := os.Getenv("LEXICMAP_RESOURCE_LIMIT")
	if kind == "" {
		return
	}
	resource, limit := syscall.RLIMIT_NOFILE, uint64(64)
	if kind == "file-size" {
		resource, limit = syscall.RLIMIT_FSIZE, 8
		signal.Ignore(syscall.SIGXFSZ)
	}
	var old syscall.Rlimit
	if err := syscall.Getrlimit(resource, &old); err != nil {
		t.Fatal(err)
	}
	old.Cur = limit
	if err := syscall.Setrlimit(resource, &old); err != nil {
		t.Fatal(err)
	}
	for i, arg := range os.Args {
		if arg == "--" {
			RootCmd.SetArgs(os.Args[i+1:])
			if err := RootCmd.Execute(); err != nil {
				os.Exit(1)
			}
			os.Exit(0)
		}
	}
	t.Fatal("missing command arguments")
}

// runResourceLimitedCommand catches both descriptor exhaustion and budget deadlocks.
func runResourceLimitedCommand(t *testing.T, kind string, args ...string) ([]byte, []byte, error) {
	t.Helper()
	ctx, cancel := context.WithTimeout(context.Background(), 30*time.Second)
	defer cancel()
	command := exec.CommandContext(ctx, os.Args[0], append([]string{"-test.run=^TestResourceLimitedCommand$", "--"}, args...)...)
	command.Env = append(os.Environ(), "LEXICMAP_RESOURCE_LIMIT="+kind)
	var output, stderr bytes.Buffer
	command.Stdout, command.Stderr = &output, &stderr
	err := command.Run()
	if ctx.Err() != nil {
		t.Fatalf("command timed out: %v", args)
	}
	return output.Bytes(), stderr.Bytes(), err
}

func TestGenomeDetailsResourceLimits(t *testing.T) {
	db := writeDetailsFixture(t, 64)
	args := []string{"utils", "genome-details", "--quiet", "-j", "1", "-i", "-e", "-d", db}
	_, stderr, err := runResourceLimitedCommand(t, "file-size", args...)
	if err == nil || !bytes.Contains(stderr, []byte("file too large")) {
		t.Fatalf("write failure not propagated: %v, %s", err, stderr)
	}
	if _, err := os.Stat(filepath.Join(db, FileGenomeDetails)); !os.IsNotExist(err) {
		t.Fatalf("failed write published cache: %v", err)
	}
	if files, err := filepath.Glob(filepath.Join(db, "."+FileGenomeDetails+"-*")); err != nil || len(files) != 0 {
		t.Fatalf("failed write left temporary cache: %v, %v", files, err)
	}
	got, stderr, err := runResourceLimitedCommand(t, "open-files", args...)
	if err != nil {
		t.Fatalf("bounded-reader extraction failed: %v, %s", err, stderr)
	}
	if bytes.Count(got, []byte("\n")) != 65 {
		t.Fatalf("retry lost genomes: %s", got)
	}
	want := runSeedIndexCommand(t, args...)
	if !bytes.Equal(got, want) {
		t.Fatal("cached details differ from extraction")
	}
}

func TestSearchResourceFileBudget(t *testing.T) {
	db, query := buildMultiBatchSearchFixture(t)
	for _, command := range [][]string{{"search", "-P", "12", "-n", "0", "-N", "0"}, {"genome", "search", "-m", "0", "-n", "10"}} {
		args := append(append([]string{}, command...), "--quiet", "-j", "16", "-J", "32", "-d", db, "-p", "12", query)
		want := runSeedIndexCommand(t, args...)
		if bytes.Count(want, []byte("\n")) < 2 {
			t.Fatalf("fixture produced no hits: %s", want)
		}
		for _, budget := range []int{32, 13, 9} {
			limited := append(append([]string{}, args...), "--max-open-files", strconv.Itoa(budget))
			if budget == 9 {
				limited = append(limited, "-w")
			}
			got, stderr, err := runResourceLimitedCommand(t, "open-files", limited...)
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
// in-memory path and on-demand readers when not every batch can have a pool.
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
	for _, budget := range []int{32, 13, 9} {
		opt := DefaultIndexSearchingOptions
		opt.NumCPUs, opt.MaxOpenFiles, opt.MaxSeedSearchingConcurrency = 16, budget, 16
		opt.MinPrefix, opt.MinSinglePrefix = 12, 12
		opt.InMemorySearch = budget == 9
		idx, err := NewIndexSearcher(db, &opt)
		if err != nil {
			t.Fatal(err)
		}
		check := func() {
			t.Helper()
			if got := count(); got != len(idx.openFileTokens) || got > budget-searchReservedFiles {
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
