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
	"context"
	"encoding/binary"
	"errors"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"syscall"
	"testing"
	"time"
)

// TestGenomeDetailsRejectsOversizedCount bounds the child so the old reader cannot exhaust host memory.
func TestGenomeDetailsRejectsOversizedCount(t *testing.T) {
	if file := os.Getenv("LEXICMAP_DETAILS_COUNT_TEST"); file != "" {
		// Leave one GiB above the current mappings, including race-detector
		// reservations, while preventing a disk count from requesting 16 GiB.
		statm, err := os.ReadFile("/proc/self/statm")
		if err != nil {
			t.Fatal(err)
		}
		pages, err := strconv.ParseUint(strings.Fields(string(statm))[0], 10, 64)
		if err != nil {
			t.Fatal(err)
		}
		var limit syscall.Rlimit
		if err := syscall.Getrlimit(syscall.RLIMIT_AS, &limit); err != nil {
			t.Fatal(err)
		}
		limit.Cur = min(limit.Cur, pages*uint64(os.Getpagesize())+(1<<30))
		if err := syscall.Setrlimit(syscall.RLIMIT_AS, &limit); err != nil {
			t.Fatal(err)
		}
		for _, extra := range []bool{false, true} {
			if err := readGenomeDetails(file, bufio.NewWriter(io.Discard), extra, 1); !errors.Is(err, ErrBrokenFile) {
				t.Fatalf("expected broken cache: %v", err)
			}
		}
		return
	}
	// One genome ID and one chunk header; the declared length table is absent.
	data := make([]byte, 31)
	binary.BigEndian.PutUint16(data[8:10], 1)
	binary.BigEndian.PutUint32(data[10:14], 1)
	data[14] = 'g'
	binary.BigEndian.PutUint32(data[23:27], 1)
	binary.BigEndian.PutUint32(data[27:31], ^uint32(0))
	file := filepath.Join(t.TempDir(), "details.bin")
	if err := os.WriteFile(file, data, 0600); err != nil {
		t.Fatal(err)
	}
	ctx, cancel := context.WithTimeout(context.Background(), 10*time.Second)
	defer cancel()
	command := exec.CommandContext(ctx, os.Args[0], "-test.run=^TestGenomeDetailsRejectsOversizedCount$")
	command.Env = append(os.Environ(), "LEXICMAP_DETAILS_COUNT_TEST="+file)
	if output, err := command.CombinedOutput(); err != nil {
		t.Fatalf("invalid count was not rejected safely: %v\n%s", err, output)
	}
}
