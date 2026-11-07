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
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"slices"
	"strings"
	"syscall"
	"testing"
	"time"
)

// checkStreamingOutputFailure keeps a pipe open after a complete record so a
// successful test requires error propagation before input EOF.
func checkStreamingOutputFailure(t *testing.T, input string, args ...string) {
	t.Helper()
	for _, mode := range []string{"plain", "gzip", "gzip_after_header"} {
		t.Run(mode, func(t *testing.T) {
			checkStreamingOutputFailureMode(t, input, mode, args...)
		})
	}
}

// The late-close case lets gzip write its header successfully before breaking
// the pipe. Flushing only the outer buffer cannot detect this failure.
func checkStreamingOutputFailureMode(t *testing.T, input, mode string, args ...string) {
	t.Helper()
	for _, output := range []string{"/dev/full", "-"} {
		if mode == "gzip_after_header" && output != "-" {
			continue
		}
		for _, kind := range []string{"stdin", "fifo"} {
			t.Run(output+"/"+kind, func(t *testing.T) {
				ctx, cancel := context.WithTimeout(context.Background(), 10*time.Second)
				defer cancel()
				arguments := slices.Clone(args)
				if mode != "plain" {
					arguments = append(arguments, "--force-gzip-output")
				}
				var inputWriter io.WriteCloser
				if kind == "fifo" {
					file := filepath.Join(t.TempDir(), "input.pipe")
					if err := syscall.Mkfifo(file, 0600); err != nil {
						t.Fatal(err)
					}
					// Keep both ends open in the parent until the child reports failure.
					writer, err := os.OpenFile(file, os.O_RDWR, 0600)
					if err != nil {
						t.Fatal(err)
					}
					inputWriter = writer
					for i, arg := range arguments {
						if arg == "-" {
							arguments[i] = file
						}
					}
				}
				argv := append([]string{"-test.run=^TestSeedIndexCommand$", "--"}, arguments...)
				command := exec.CommandContext(ctx, os.Args[0], append(argv, "-o", output)...)
				command.Env = append(os.Environ(), "LEXICMAP_SEED_INDEX_TEST_HELPER=1")
				if inputWriter == nil {
					writer, err := command.StdinPipe()
					if err != nil {
						t.Fatal(err)
					}
					inputWriter = writer
				}
				defer inputWriter.Close()
				var stderr bytes.Buffer
				command.Stderr = &stderr
				var outputReader *os.File
				defer func() {
					if outputReader != nil {
						outputReader.Close()
					}
				}()
				if output == "-" {
					reader, writer, err := os.Pipe()
					if err != nil {
						t.Fatal(err)
					}
					if mode == "gzip_after_header" {
						outputReader = reader
						if err := reader.SetReadDeadline(time.Now().Add(5 * time.Second)); err != nil {
							t.Fatal(err)
						}
					} else {
						reader.Close()
					}
					defer writer.Close()
					command.Stdout = writer
				}
				if err := command.Start(); err != nil {
					t.Fatal(err)
				}
				if _, err := io.WriteString(inputWriter, input); err != nil {
					t.Fatal(err)
				}
				if outputReader != nil {
					var header [10]byte
					if _, err := io.ReadFull(outputReader, header[:]); err != nil {
						t.Fatalf("failed to read gzip header: %v", err)
					}
					if header[0] != 0x1f || header[1] != 0x8b {
						t.Fatal("missing gzip header")
					}
					outputReader.Close()
					outputReader = nil
					// The first flush may already have noticed the closed pipe.
					_, _ = io.WriteString(inputWriter, input)
				}
				err := command.Wait()
				if ctx.Err() != nil {
					t.Fatal("output failure waited for input EOF")
				}
				want := "no space left on device"
				if output == "-" {
					want = "broken pipe"
				}
				if err == nil || !strings.Contains(stderr.String(), want) || strings.Contains(stderr.String(), "panic:") {
					t.Fatalf("expected %s, got %v: %s", want, err, &stderr)
				}
			})
		}
	}
}

func TestBLASTStopsOnOutputFailure(t *testing.T) {
	input := auditSearchHeader + "\tcigar\tqseq\tsseq\talign\n" +
		"q\t4\t1\tg\ts\t100\t1\t1\t100\t4\t100\t0\t1\t4\t1\t4\t+\t4\t0\t8\t4M\tACGT\tACGT\t||||\n"
	checkStreamingOutputFailure(t, input, "utils", "2blast", "--quiet", "-")
}
