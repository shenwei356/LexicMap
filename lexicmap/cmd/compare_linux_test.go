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
	"testing"
)

// A near-full batch pool must leave room for the process and command output.
func TestCompareResourceFileBudget(t *testing.T) {
	db := writeDetailsFixture(t, 22)
	args := []string{"genome", "compare", "--quiet", "-j", "4", "-d", db, "g0", "g1"}
	want := runSeedIndexCommand(t, args...)
	for _, budget := range []string{"24", "9"} {
		got, stderr, err := runResourceLimitedCommand(t, "open-files=24", append(args, "--max-open-files", budget)...)
		if err != nil || !bytes.Equal(got, want) {
			t.Fatalf("budget %s: error=%v stdout=%s stderr=%s", budget, err, got, stderr)
		}
	}
}
