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
	"os"
	"os/exec"
	"os/signal"
	"path/filepath"
	"strings"
	"syscall"
	"testing"
)

// A subprocess confines the file-size limit to the command under test.
func TestEditGenomeIDsWriteFailureHelper(t *testing.T) {
	dir := os.Getenv("LEXICMAP_EDIT_WRITE_FAILURE_DIR")
	if dir == "" {
		return
	}
	signal.Ignore(syscall.SIGXFSZ)
	if err := syscall.Setrlimit(syscall.RLIMIT_FSIZE, &syscall.Rlimit{Cur: 128, Max: 128}); err != nil {
		t.Fatal(err)
	}
	RootCmd.SetArgs([]string{"utils", "edit-genome-ids", "--quiet", "-d", dir, "-p", "^g$", "-r", strings.Repeat("r", 256)})
	if err := RootCmd.Execute(); err != nil {
		t.Fatal(err)
	}
}

// Cover a buffered Flush failure and a Write failure after the buffer fills.
func TestEditGenomeIDsCleansWriteFailures(t *testing.T) {
	for _, records := range []int{1, 2000} {
		t.Run(fmt.Sprint(records), func(t *testing.T) {
			dir := t.TempDir()
			record := []byte{0, 1, 'g', 0, 0, 0, 0, 0, 0, 0, 0}
			original := bytes.Repeat(record, records)
			mapping := filepath.Join(dir, FileGenomeIndex)
			if err := os.WriteFile(mapping, original, 0600); err != nil {
				t.Fatal(err)
			}
			command := exec.Command(os.Args[0], "-test.run=^TestEditGenomeIDsWriteFailureHelper$")
			command.Env = append(os.Environ(), "LEXICMAP_EDIT_WRITE_FAILURE_DIR="+dir)
			output, err := command.CombinedOutput()
			if err == nil || !bytes.Contains(output, []byte("file too large")) {
				t.Fatalf("expected a write error: %v\n%s", err, output)
			}
			got, err := os.ReadFile(mapping)
			if err != nil || !bytes.Equal(got, original) {
				t.Fatalf("original mapping changed: %v", err)
			}
			entries, err := os.ReadDir(dir)
			if err != nil || len(entries) != 1 || entries[0].Name() != FileGenomeIndex {
				t.Fatalf("staging file or backup remains: %v, %v", entries, err)
			}
		})
	}
}
