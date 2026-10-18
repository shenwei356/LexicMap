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
	"path/filepath"
	"testing"

	"github.com/shenwei356/lexichash"
)

func TestLoadIndexMasksUsesCustomMaskCount(t *testing.T) {
	lh, err := lexichash.NewWithSeed(31, 64, 1, 0)
	if err != nil {
		t.Fatal(err)
	}

	var data bytes.Buffer
	for _, mask := range lh.Masks {
		if _, err := fmt.Fprintf(&data, "%s\n", lexichash.MustDecode(mask, uint8(lh.K))); err != nil {
			t.Fatal(err)
		}
	}

	file := filepath.Join(t.TempDir(), "masks.txt")
	if err := os.WriteFile(file, data.Bytes(), 0o644); err != nil {
		t.Fatal(err)
	}

	opt := &IndexBuildingOptions{MaskFile: file, K: 21, Masks: 20_000, Chunks: 16}
	got, err := loadIndexMasks(opt)
	if err != nil {
		t.Fatal(err)
	}
	if got.K != 31 || opt.K != 31 {
		t.Fatalf("k-mer size: loaded=%d option=%d, want 31", got.K, opt.K)
	}
	if len(got.Masks) != 64 || opt.Masks != 64 {
		t.Fatalf("mask count: loaded=%d option=%d, want 64", len(got.Masks), opt.Masks)
	}
}
