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
	"compress/gzip"
	"context"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
	"time"
)

// Output failures must clean the current payload and drain concurrent results.
func TestSearchOutputFailureCleansSpills(t *testing.T) {
	dir := t.TempDir()
	sequence := indexTestDNA(40000)
	ref := writeTestInput(t, filepath.Join(dir, "ref.fa"), ">s\n"+sequence+"\n")
	db := filepath.Join(dir, "db")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-k", "21", "-c", "1", "--partitions", "16", "-O", db, ref)
	base := []string{"search", "--quiet", "-j", "4", "-J", "8", "-d", db, "-a", "-n", "0", "-N", "0", "--max-align-result-memory", "1", "--max-seed-memory", "64K", "--max-open-files", "24"}
	want := runSeedIndexCommand(t, append(append([]string{}, base...), ref)...)
	if len(want) <= maxSearchOutputRowBytes || bytes.Count(want, []byte{'\n'}) < 2 {
		t.Fatal("fixture must produce an alignment row larger than the output buffer")
	}
	short := writeTestInput(t, filepath.Join(dir, "short.fa"), ">short\n"+sequence[:400]+"\n")
	shortOutput := runSeedIndexCommand(t, append(append([]string{}, base...), short)...)
	if len(shortOutput) >= IOBufferSize || bytes.Count(shortOutput, []byte{'\n'}) < 2 {
		t.Fatal("short fixture must fail during flush rather than row writing")
	}
	for _, tc := range []struct {
		name, query string
		queued      bool
		streaming   bool
		compressed  bool
	}{
		{"write", ref, false, false, false},
		{"flush", short, false, false, false},
		{"broken_pipe", ref, true, false, false},
		{"streaming_flush", short, false, true, false},
		{"streaming_broken_pipe", ref, true, true, false},
		{"streaming_gzip", ref, true, true, true},
	} {
		t.Run(tc.name, func(t *testing.T) {
			temp := t.TempDir()
			query := tc.query
			output, diagnostic := "/dev/full", "no space left"
			var pipe io.ReadCloser
			if tc.queued {
				var input strings.Builder
				for i := 0; i < 32; i++ {
					fmt.Fprintf(&input, ">q%d\n%s\n", i, sequence)
				}
				query = writeTestInput(t, filepath.Join(t.TempDir(), "queries.fa"), input.String())
				output, diagnostic = "-", "broken pipe"
			}
			ctx, cancel := context.WithTimeout(context.Background(), 30*time.Second)
			defer cancel()
			if tc.streaming {
				query = "-"
			}
			args := append(append([]string{}, base...), query, "-o", output)
			command := exec.CommandContext(ctx, os.Args[0], append([]string{"-test.run=^TestSeedIndexCommand$", "--"}, args...)...)
			command.Env = append(os.Environ(), "LEXICMAP_SEED_INDEX_TEST_HELPER=1", "TMPDIR="+temp)
			var stderr bytes.Buffer
			command.Stderr = &stderr
			var stdin io.WriteCloser
			if tc.streaming {
				var err error
				stdin, err = command.StdinPipe()
				if err != nil {
					t.Fatal(err)
				}
				defer stdin.Close()
			}
			if tc.queued {
				var err error
				pipe, err = command.StdoutPipe()
				if err != nil {
					t.Fatal(err)
				}
				defer pipe.Close()
			}
			if err := command.Start(); err != nil {
				t.Fatal(err)
			}
			if tc.streaming {
				inputSequence := sequence
				if !tc.queued {
					inputSequence = sequence[:400]
				}
				// The next header releases q2, then leaves the input reader
				// waiting for data. Keep stdin open until the command exits.
				var input io.Writer = stdin
				var compressor *gzip.Writer
				if tc.compressed {
					compressor = gzip.NewWriter(stdin)
					input = compressor
				}
				if _, err := fmt.Fprintf(input, ">q1\n%s\n>q2\n%s\n>q3\n", inputSequence, inputSequence); err != nil {
					t.Fatal(err)
				}
				if compressor != nil {
					// Release the decoder's 1 MiB read-ahead block while leaving
					// q3 and the gzip stream unfinished, awaiting upstream data.
					if _, err := io.WriteString(compressor, strings.Repeat("N", 1<<20)); err != nil {
						t.Fatal(err)
					}
					if err := compressor.Flush(); err != nil {
						t.Fatal(err)
					}
				}
			}
			if tc.queued {
				// Withhold reads until the blocked writer and another query own
				// spill files. Closing stdout's reader then triggers SIGPIPE.
				for {
					files, err := filepath.Glob(filepath.Join(temp, "lexicmap-alignment-*.bin"))
					if err != nil {
						t.Fatal(err)
					}
					if len(files) >= 2 {
						break
					}
					if ctx.Err() != nil {
						_ = command.Wait()
						t.Fatalf("concurrent payloads were not created: %s", stderr.String())
					}
					time.Sleep(10 * time.Millisecond)
				}
				if err := pipe.Close(); err != nil {
					t.Fatal(err)
				}
			}
			err := command.Wait()
			if err == nil || ctx.Err() != nil {
				t.Fatalf("output failure: %v, timeout: %v, stderr: %s", err, ctx.Err(), stderr.String())
			}
			if !strings.Contains(stderr.String(), diagnostic) {
				t.Errorf("missing output error %q: %v, stderr: %s", diagnostic, err, stderr.String())
			}
			entries, err := os.ReadDir(temp)
			if err != nil || len(entries) != 0 {
				t.Fatalf("temporary files leaked after output failure: %v, %v", entries, err)
			}
		})
	}
}
