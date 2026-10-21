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
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"slices"
	"strings"
	"testing"
	"time"
)

// Slow output keeps completed queries' spill files alive while later queries run.
func TestAlignmentPayloadResourceFileBudget(t *testing.T) {
	dir := t.TempDir()
	sequence := indexTestDNA(6000)
	ref := writeTestInput(t, filepath.Join(dir, "g.fa"), ">g\n"+sequence+"\n")
	db := filepath.Join(dir, "db")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-c", "1", "--partitions", "4", "-O", db, ref)
	var data strings.Builder
	for i := 0; i < 100; i++ {
		fmt.Fprintf(&data, ">q%d\n%s\n", i, sequence)
	}
	query := writeTestInput(t, filepath.Join(dir, "queries.fa"), data.String())
	args := []string{"search", "--quiet", "-j", "32", "-J", "32", "-d", db, "-a", "-n", "0", "-N", "0", query}
	want := strings.Split(string(runSeedIndexCommand(t, append(args, "--max-align-result-memory", "0")...)), "\n")
	if len(want) != 102 {
		t.Fatalf("expected 100 matches, got %d lines", len(want))
	}
	slices.Sort(want)
	for _, tc := range []struct {
		seedMemory, fileBudget string
		inMemory               bool
	}{
		{"0", "24", false}, {"1M", "24", false},
		{"0", "11", false}, {"1M", "14", false},
		{"0", "10", true}, {"1M", "13", true},
	} {
		ctx, cancel := context.WithTimeout(context.Background(), 30*time.Second)
		limited := append(append([]string{}, args...), "--max-open-files", tc.fileBudget, "--max-align-result-memory", "1", "--max-seed-memory", tc.seedMemory)
		if tc.inMemory {
			limited = append(limited, "-w")
		}
		command := exec.CommandContext(ctx, os.Args[0], append([]string{"-test.run=^TestResourceLimitedCommand$", "--"}, limited...)...)
		temp := t.TempDir()
		command.Env = append(os.Environ(), "LEXICMAP_RESOURCE_LIMIT=open-files=24", "TMPDIR="+temp)
		var stderr bytes.Buffer
		command.Stderr = &stderr
		output, err := command.StdoutPipe()
		if err != nil {
			cancel()
			t.Fatal(err)
		}
		if err := command.Start(); err != nil {
			cancel()
			t.Fatal(err)
		}
		// Backpressure must neither exhaust descriptors nor deadlock the shared readers.
		time.Sleep(600 * time.Millisecond)
		out, readErr := io.ReadAll(output)
		err = command.Wait()
		timedOut := ctx.Err() != nil
		cancel()
		got := strings.Split(string(out), "\n")
		slices.Sort(got)
		if err != nil || readErr != nil || timedOut || !slices.Equal(got, want) {
			t.Fatalf("%+v: wait=%v read=%v timeout=%v stderr=%s output differs=%v", tc, err, readErr, timedOut, stderr.String(), !slices.Equal(got, want))
		}
		entries, err := os.ReadDir(temp)
		if err != nil || len(entries) != 0 {
			t.Fatalf("temporary files leaked: %v, %v", entries, err)
		}
	}
}
