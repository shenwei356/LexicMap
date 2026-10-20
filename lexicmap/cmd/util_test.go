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
	"os"
	"path/filepath"
	"testing"
)

func TestMakeOutDirRejectsCurrentDirectory(t *testing.T) {
	cwd := t.TempDir()
	previous, err := os.Getwd()
	if err != nil {
		t.Fatal(err)
	}
	if err := os.Chdir(cwd); err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() {
		if err := os.Chdir(previous); err != nil {
			t.Error(err)
		}
	})

	sentinel := writeTestInput(t, filepath.Join(cwd, "sentinel"), "preserve")
	alias := filepath.Join(t.TempDir(), "cwd")
	if err := os.Symlink(cwd, alias); err != nil {
		t.Fatal(err)
	}
	for _, path := range []string{".", "./", cwd, "child" + string(filepath.Separator) + "..", alias} {
		for _, force := range []bool{false, true} {
			forceName := "without-force"
			if force {
				forceName = "with-force"
			}
			t.Run(filepath.Base(path)+"/"+forceName, func(t *testing.T) {
				if err := makeOutDir(path, force, "out-dir", false); err == nil {
					t.Fatalf("accepted current directory %q", path)
				}
				if got, err := os.ReadFile(sentinel); err != nil || string(got) != "preserve" {
					t.Fatalf("sentinel changed: %q, %v", got, err)
				}
			})
		}
	}
}
