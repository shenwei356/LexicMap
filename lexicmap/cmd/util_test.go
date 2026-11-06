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
	"compress/bzip2"
	"compress/gzip"
	"encoding/base64"
	"io"
	"os"
	"path/filepath"
	"reflect"
	"regexp"
	"strings"
	"testing"
)

func TestGenomeNamesFromCompressedFiles(t *testing.T) {
	// A bzip2-compressed synthetic 600-base FASTA record, with an ID-only header.
	compressed, err := base64.StdEncoding.DecodeString(
		"QlpoOTFBWSZTWRstSnwAAAVPgHAQAAEogAQAAgAoADABmACGptqmgeoQNA0Bp4qUmkH5vw/Rd/ftuv58cKCVAWcFpagNGt9n3g333vPjLEtJo5RFTS50ZDPVkXhi10M5lyHvXeMHDu9arNXluPFZsXaQHZpyA8wt" +
			"jRrP5+/L3wLQRhHAijoM4pkKxDYsKiVzzJnSYjnQBmt9+/PnvhRznOrBBIV70WV9jDUCEw9gbRJ00IT1mysgay56F6hKNS82aCcU8k0LkAbIuBmboYoNNE60hTrgxcVbBEjnhg5yOEUFACaTEJkk+LuSKcKEgNlqU+A=")
	if err != nil {
		t.Fatal(err)
	}
	plain, err := io.ReadAll(bzip2.NewReader(bytes.NewReader(compressed)))
	if err != nil {
		t.Fatal(err)
	}
	var gz bytes.Buffer
	zw := gzip.NewWriter(&gz)
	if _, err := zw.Write(plain); err != nil {
		t.Fatal(err)
	}
	if err := zw.Close(); err != nil {
		t.Fatal(err)
	}
	for _, format := range []struct {
		ext  string
		data []byte
	}{{".fa", plain}, {".fa.gz", gz.Bytes()}, {".fa.bz2", compressed}, {".FA.BZ2", compressed}} {
		for _, pattern := range []struct{ regexp, want string }{
			{"", "sample"}, {`missing-(.+)`, "sample"}, {`(sam)ple`, "sam"},
		} {
			t.Run(format.ext+"/"+pattern.regexp, func(t *testing.T) {
				dir := t.TempDir()
				file := writeTestInput(t, filepath.Join(dir, "sample"+format.ext), string(format.data))
				re, err := compileRefNameRegexp(pattern.regexp)
				if err != nil {
					t.Fatal(err)
				}
				// Both genome commands and index construction use the same fallback.
				for _, read := range []func() (*GQuery, error){
					func() (*GQuery, error) { return ReadGenomeFromFile(file, re, false) },
					func() (*GQuery, error) { return NewGenomeReader(31, re).Read(file, false, false) },
				} {
					q, err := read()
					if err != nil || q == nil {
						t.Fatalf("read genome: %v", err)
					}
					if string(q.id) != pattern.want || q.genomeSize != 600 {
						t.Errorf("genome ID/size = %q/%d, want %q/600", q.id, q.genomeSize, pattern.want)
					}
					RecycleGQuery(q)
				}
				db := filepath.Join(dir, "db")
				runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-c", "1", "--partitions", "4",
					"--ref-name-regexp", pattern.regexp, "-O", db, file)
				got := runSeedIndexCommand(t, "utils", "genomes", "--quiet", "-d", db)
				if string(got) != "ref\tchunked\n"+pattern.want+"\t\n" {
					t.Errorf("indexed genome names = %q", got)
				}
			})
		}
	}
}

func TestRefNameRegexpCaptureValidation(t *testing.T) {
	for _, pattern := range []string{`(?:.*)`, `(?i:.*)`, `\(.+\)`, `[(].*[)]`, `(`} {
		if _, err := compileRefNameRegexp(pattern); err == nil {
			t.Errorf("accepted a regexp without a usable capture: %q", pattern)
		}
	}
	for _, pattern := range []string{`(ref)`, `(?i)(ref)`, `(?P<id>ref)`, `(?:prefix)?(ref)`} {
		re, err := compileRefNameRegexp(pattern)
		if err != nil || !re.MatchString("REF") {
			t.Fatalf("valid case-insensitive capture %q: %v", pattern, err)
		}
	}
	if re, err := compileRefNameRegexp(""); re != nil || err != nil {
		t.Fatalf("disabled extraction: %v, %v", re, err)
	}
	// Library entry points must reject the same regexp before opening inputs.
	re := regexp.MustCompile(`(?:.*)`)
	if _, err := ReadGenomeFromFile("missing.fa", re, false); err == nil || !strings.Contains(err.Error(), "capturing group") {
		t.Fatalf("file reader validation: %v", err)
	}
	if _, err := NewGenomeReader(31, re).Read("missing.fa", false, false); err == nil || !strings.Contains(err.Error(), "capturing group") {
		t.Fatalf("genome reader validation: %v", err)
	}
	if err := CheckIndexBuildingOptions(&IndexBuildingOptions{ReRefName: re}); err == nil || !strings.Contains(err.Error(), "capturing group") {
		t.Fatalf("index option validation: %v", err)
	}
	if err := BuildIndex(filepath.Join(t.TempDir(), "db"), nil, &IndexBuildingOptions{ReRefName: re}); err == nil || !strings.Contains(err.Error(), "capturing group") {
		t.Fatalf("index builder validation: %v", err)
	}
}

func TestCommandsRejectRefNameRegexpWithoutCapture(t *testing.T) {
	dir := t.TempDir()
	ref := writeTestInput(t, filepath.Join(dir, "ref.fa"), ">seq\n"+indexTestDNA(1000)+"\n")
	db := filepath.Join(dir, "db")
	for _, args := range [][]string{
		{"index", "-O", db, ref},
		{"genome", "search", "-d", db, ref},
		{"genome", "compare", ref, ref},
	} {
		args = append(args, "--quiet", "--ref-name-regexp", "(?:.*)")
		out, stderr, err := runAuditCommand(t, "", args...)
		if err == nil || !strings.Contains(string(stderr), "capturing group") || strings.Contains(string(stderr), "panic:") || len(out) != 0 {
			t.Fatalf("%v: error=%v stdout=%s stderr=%s", args, err, out, stderr)
		}
	}
	if _, err := os.Stat(db); !os.IsNotExist(err) {
		t.Fatalf("invalid regexp created an index directory: %v", err)
	}
}

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

// Reused split buffers must recover their full length after skipped short rows.
func TestSplitNReusesShortRows(t *testing.T) {
	for _, split := range []func(string, *[]string){
		func(s string, a *[]string) { stringSplitN(s, "\t", 2, a) },
		func(s string, a *[]string) { stringSplitNByByte(s, '\t', 2, a) },
	} {
		fields := make([]string, 2)
		for _, row := range []string{"short", "a\tb", "short", "c\td\te"} {
			split(row, &fields)
			if want := strings.SplitN(row, "\t", 2); !reflect.DeepEqual(fields, want) {
				t.Fatalf("row %q: got %v, want %v", row, fields, want)
			}
		}
	}
	fields := make([][]byte, 2)
	for _, row := range []string{"short", "a\tb", "short", "c\td\te"} {
		bytesSplitN([]byte(row), []byte("\t"), 2, &fields)
		want := strings.SplitN(row, "\t", 2)
		if len(fields) != len(want) {
			t.Fatalf("row %q: %q", row, fields)
		}
		for i := range fields {
			if string(fields[i]) != want[i] {
				t.Fatalf("row %q: %q", row, fields)
			}
		}
	}
}

func TestMappingReadersSkipShortRows(t *testing.T) {
	file := writeTestInput(t, filepath.Join(t.TempDir(), "map.tsv"), "short\nA\t12\n\nshort\nB\t34\n")
	got, err := readKVs(file, true)
	if err != nil || !reflect.DeepEqual(got, map[string]string{"a": "12", "b": "34"}) {
		t.Fatalf("string map: %v, %v", got, err)
	}
	numbers, err := readKVsUint32(file, true, true)
	if err != nil || !reflect.DeepEqual(numbers, map[string]uint32{"a": 12, "b": 34}) {
		t.Fatalf("numeric map: %v, %v", numbers, err)
	}
}
