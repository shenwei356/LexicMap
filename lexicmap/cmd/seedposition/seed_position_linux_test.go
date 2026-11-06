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

package seedposition

import (
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"runtime/debug"
	"syscall"
	"testing"
)

// Count only this failure device's descriptors, excluding unrelated runtime files.
func fullDeviceDescriptors(t *testing.T) int {
	t.Helper()
	entries, err := os.ReadDir("/proc/self/fd")
	if err != nil {
		t.Fatal(err)
	}
	n := 0
	for _, entry := range entries {
		target, _ := os.Readlink(filepath.Join("/proc/self/fd", entry.Name()))
		if target == "/dev/full" {
			n++
		}
	}
	return n
}

// Header and final-flush failures must release both data and index handles.
func TestWriterDiskFullCleanup(t *testing.T) {
	gc := debug.SetGCPercent(-1) // finalizers must not hide leaked descriptors
	defer debug.SetGCPercent(gc)
	originalBufferSize := BufferSize
	defer func() { BufferSize = originalBufferSize }()
	for _, phase := range []string{"data-header", "data-flush", "index"} {
		for _, size := range []int{1, 8, 16, 65536} {
			if phase == "data-header" && size >= 16 || phase == "data-flush" && size < 16 {
				continue
			}
			t.Run(fmt.Sprintf("%s/buffer=%d", phase, size), func(t *testing.T) {
				BufferSize = size
				file := filepath.Join(t.TempDir(), "data.bin")
				target := file
				if phase == "index" {
					target += PositionsIndexFileExt
				}
				if err := os.Symlink("/dev/full", target); err != nil {
					t.Fatal(err)
				}
				before := fullDeviceDescriptors(t)
				w, err := NewWriter(file, 0)
				if phase != "data-header" {
					if err != nil {
						t.Fatal(err)
					}
					defer w.fh.Close()
					err = w.Close()
				}
				if !errors.Is(err, syscall.ENOSPC) {
					t.Errorf("expected disk-full failure, got %v", err)
				}
				if after := fullDeviceDescriptors(t); after != before {
					t.Errorf("disk-full handles grew from %d to %d", before, after)
				}
			})
		}
	}
}
