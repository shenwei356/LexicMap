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
	"testing"
)

func TestSubseqClosesUnmatchedChunkReaders(t *testing.T) {
	db, sequence := writeChunkExtractionFixture(t)
	t.Setenv("GOGC", "off")
	output, stderr, err := runResourceLimitedCommand(t, "open-files", "utils", "subseq", "--quiet", "-d", db, "-n", "many", "-s", "s95", "-r", "1:10", "-w", "0")
	if err != nil {
		t.Fatalf("last chunk extraction failed: %v %s", err, stderr)
	}
	if string(output) != ">s95:1-10:+\n"+sequence[:10]+"\n" {
		t.Fatalf("wrong sequence: %s", output)
	}
}
