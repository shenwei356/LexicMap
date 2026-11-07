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
	"os"
	"os/exec"
	"path/filepath"
	"syscall"
	"testing"
	"time"
)

// TestSAMRejectsStreamsBeforeOutput never supplies a FIFO producer: validation must not open it.
func TestSAMRejectsStreamsBeforeOutput(t *testing.T) {
	for _, kind := range []string{"fifo", "symlink", "list", "stdin", "stdin-alias"} {
		t.Run(kind, func(t *testing.T) {
			dir := t.TempDir()
			fifo := filepath.Join(dir, "input.fifo")
			if err := syscall.Mkfifo(fifo, 0600); err != nil {
				t.Fatal(err)
			}
			input := fifo
			if kind == "symlink" {
				input = filepath.Join(dir, "alias")
				if err := os.Symlink(fifo, input); err != nil {
					t.Fatal(err)
				}
			}
			if kind == "stdin" {
				input = "-"
			}
			if kind == "stdin-alias" {
				input = "/dev/stdin"
			}
			output, logFile := filepath.Join(dir, "out.sam"), filepath.Join(dir, "run.log")
			sentinel := []byte("preserve existing content\n")
			for _, file := range []string{output, logFile} {
				if err := os.WriteFile(file, sentinel, 0600); err != nil {
					t.Fatal(err)
				}
			}
			args := []string{"-test.run=^TestSeedIndexCommand$", "--", "utils", "2sam", "--quiet", "-o", output, "--log", logFile}
			if kind == "list" {
				list := filepath.Join(dir, "inputs.txt")
				if err := os.WriteFile(list, []byte(fifo+"\n"), 0600); err != nil {
					t.Fatal(err)
				}
				args = append(args, "-X", list)
			} else {
				args = append(args, input)
			}
			ctx, cancel := context.WithTimeout(context.Background(), 3*time.Second)
			defer cancel()
			command := exec.CommandContext(ctx, os.Args[0], args...)
			command.Env = append(os.Environ(), "LEXICMAP_SEED_INDEX_TEST_HELPER=1")
			reader, writer, err := os.Pipe()
			if err != nil {
				t.Fatal(err)
			}
			defer reader.Close()
			defer writer.Close()
			command.Stdin = reader
			stderr, err := command.CombinedOutput()
			if ctx.Err() != nil {
				t.Fatal("waited for a stream instead of rejecting it")
			}
			if err == nil || !bytes.Contains(stderr, []byte("read twice")) {
				t.Fatalf("stream accepted: %v %s", err, stderr)
			}
			for _, file := range []string{output, logFile} {
				got, err := os.ReadFile(file)
				if err != nil || !bytes.Equal(got, sentinel) {
					t.Fatalf("output modified before validation: %s: %q, %v", file, got, err)
				}
			}
		})
	}
}

// TestSAMRegularFileSymlink keeps replayable symlinks usable after rejecting streams.
func TestSAMRegularFileSymlink(t *testing.T) {
	dir := t.TempDir()
	file := writeTestInput(t, filepath.Join(dir, "input.tsv"), auditSearchHeader+"\tcigar\tqseq\tsseq\talign\n"+
		"q\t4\t1\tg\ts\t100\t1\t1\t100\t4\t100\t0\t1\t4\t1\t4\t+\t4\t0\t8\t4M\tACGT\tACGT\t||||\n")
	link := filepath.Join(dir, "alias")
	if err := os.Symlink(file, link); err != nil {
		t.Fatal(err)
	}
	want := runSeedIndexCommand(t, "utils", "2sam", "--quiet", file)
	got := runSeedIndexCommand(t, "utils", "2sam", "--quiet", link)
	if !bytes.Equal(got, want) {
		t.Fatal("regular-file alias changed SAM output")
	}
}
