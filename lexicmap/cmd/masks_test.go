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
	"testing"
)

func TestMasksSelectionBounds(t *testing.T) {
	dir := t.TempDir()
	input := writeTestInput(t, filepath.Join(dir, "g.fa"), ">g\n"+indexTestDNA(500)+"\n")
	db := filepath.Join(dir, "db")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-c", "2", "--partitions", "4", "-O", db, input)
	all := runSeedIndexCommand(t, "utils", "masks", "--quiet", "-d", db)
	last := runSeedIndexCommand(t, "utils", "masks", "--quiet", "-d", db, "-m", "64")
	if bytes.Count(all, []byte("\n")) != 64 || !bytes.HasSuffix(all, last) || !bytes.HasPrefix(last, []byte("64\t")) {
		t.Fatalf("incorrect valid selection: %s", last)
	}
	output := writeTestInput(t, filepath.Join(dir, "out.tsv"), "keep existing output\n")
	for _, mask := range []string{"0", "65"} {
		_, stderr, err := runAuditCommand(t, "", "utils", "masks", "--quiet", "-d", db, "-m", mask, "-o", output)
		if err == nil || bytes.Contains(stderr, []byte("panic:")) {
			t.Fatalf("expected normal parameter error for %s: %v\n%s", mask, err, stderr)
		}
		got, err := os.ReadFile(output)
		if err != nil || string(got) != "keep existing output\n" {
			t.Fatalf("invalid selection changed output: %q, %v", got, err)
		}
	}
}
