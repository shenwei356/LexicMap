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
	"strings"
	"testing"
)

func TestEditGenomeIDsPreservesTemporaryAlias(t *testing.T) {
	for _, hard := range []bool{false, true} {
		t.Run(map[bool]string{false: "symlink", true: "hardlink"}[hard], func(t *testing.T) {
			db := writeSeedPosFixture(t, 1)
			mapping := filepath.Join(db, FileGenomeIndex)
			if err := os.Chmod(mapping, 0640); err != nil {
				t.Fatal(err)
			}
			original, err := os.ReadFile(mapping)
			if err != nil {
				t.Fatal(err)
			}
			link := os.Symlink
			if hard {
				link = os.Link
			}
			if err := link(mapping, mapping+".new"); err != nil {
				t.Fatal(err)
			}
			runSeedIndexCommand(t, "utils", "edit-genome-ids", "--quiet", "-d", db, "-p", "^g0$", "-r", "renamed")
			if got, err := os.ReadFile(mapping + ".bak"); err != nil || !bytes.Equal(got, original) {
				t.Fatalf("original mapping lost: %v", err)
			}
			if info, err := os.Stat(mapping); err != nil || info.Mode().Perm() != 0640 {
				t.Fatalf("mapping permissions changed: %v, %v", info, err)
			}
			ids, err := readGenomeMapName2Idx(mapping)
			if err != nil || ids["renamed"] == nil || len(ids) != 1 {
				t.Fatalf("rename failed: %v, %v", ids, err)
			}
			if _, err := os.Lstat(mapping + ".new"); err != nil {
				t.Fatalf("unrelated entry removed: %v", err)
			}
		})
	}
}

// Failures after buffered writes must remove staging files and preserve the mapping.
func TestEditGenomeIDsCleansFailedRewrite(t *testing.T) {
	for _, failure := range []string{"empty replacement", "truncated mapping", "missing mapping", "cache removal"} {
		t.Run(failure, func(t *testing.T) {
			dir := t.TempDir()
			file := filepath.Join(dir, FileGenomeIndex)
			var data bytes.Buffer
			var length [2]byte
			var code [8]byte
			for i := range 20000 {
				id := fmt.Sprintf("g%05d", i)
				be.PutUint16(length[:], uint16(len(id)))
				be.PutUint64(code[:], uint64(i))
				data.Write(length[:])
				data.WriteString(id)
				data.Write(code[:])
			}
			if failure == "truncated mapping" {
				data.WriteByte(0)
			}
			if failure != "missing mapping" {
				if err := os.WriteFile(file, data.Bytes(), 0600); err != nil {
					t.Fatal(err)
				}
			}
			cache := filepath.Join(dir, FileGenomeDetails)
			if failure == "cache removal" {
				if err := os.Mkdir(cache, 0700); err != nil {
					t.Fatal(err)
				}
				writeTestInput(t, filepath.Join(cache, "keep"), "keep")
			} else {
				writeTestInput(t, cache, "keep cached details")
			}
			replacement := "renamed"
			if failure == "empty replacement" {
				replacement = ""
			}
			_, stderr, err := runAuditCommand(t, "", "utils", "edit-genome-ids", "--quiet", "-d", dir, "-p", "^g19999$", "-r", replacement)
			if err == nil {
				t.Fatalf("accepted %s: %s", failure, stderr)
			}
			wantError := map[string]string{
				"empty replacement": "rename to empty", "truncated mapping": "failed to read genome index mapping file",
				"missing mapping": "failed to read genome index mapping file", "cache removal": "failed to invalidate genome details",
			}[failure]
			if !strings.Contains(string(stderr), wantError) {
				t.Fatalf("wrong error: %s", stderr)
			}
			if failure != "missing mapping" {
				got, err := os.ReadFile(file)
				if err != nil || !bytes.Equal(got, data.Bytes()) {
					t.Fatalf("original mapping changed: %v", err)
				}
			}
			entries, err := filepath.Glob(filepath.Join(dir, "."+FileGenomeIndex+"-*"))
			if err != nil || len(entries) != 0 {
				t.Errorf("temporary files remain: %v, %v", entries, err)
			}
			if _, err := os.Stat(file + ".bak"); !os.IsNotExist(err) {
				t.Errorf("backup created on failure: %v", err)
			}
			if failure != "cache removal" {
				got, err := os.ReadFile(cache)
				if err != nil || string(got) != "keep cached details" {
					t.Errorf("cache changed on failure: %v", err)
				}
			}
		})
	}
}
