// Copyright © 2023-2024 Wei Shen <shenwei356@gmail.com>
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
	"bytes"
	"fmt"
	"math"
	"math/rand"
	"os"
	"path/filepath"
	"regexp"
	"runtime"
	"strconv"
	"strings"

	"github.com/iafan/cwalk"
	"github.com/pkg/errors"
	"github.com/shenwei356/util/pathutil"
	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
	"github.com/twotwotwo/sorts"
)

var mapInitSize = 1 << 20 // 1M

// the minimum k value
var minK = 5

// Options contains the global flags
type Options struct {
	NumCPUs int
	Verbose bool

	LogFile  string
	Log2File bool

	Compress         bool
	CompressionLevel int
}

func getOptions(cmd *cobra.Command) *Options {
	threads := getFlagNonNegativeInt(cmd, "threads")
	if threads == 0 {
		threads = runtime.NumCPU()
	}

	sorts.MaxProcs = threads
	runtime.GOMAXPROCS(threads)

	logfile := getFlagString(cmd, "log")
	return &Options{
		NumCPUs: threads,
		Verbose: !getFlagBool(cmd, "quiet"),

		LogFile:  logfile,
		Log2File: logfile != "",

		Compress:         true,
		CompressionLevel: -1,
	}
}

func checkFileSuffix(suffix string, files ...string) {
	for _, file := range files {
		if isStdin(file) {
			continue
		}

		if suffix != "" && !strings.HasSuffix(file, suffix) {
			checkError(fmt.Errorf("input should be stdin or %s files: %s", suffix, file))
		}
	}
}

func makeOutDir(outDir string, force bool, logname string, verbose bool) {
	pwd, _ := os.Getwd()
	if outDir != "./" && outDir != "." && pwd != filepath.Clean(outDir) {
		existed, err := pathutil.DirExists(outDir)
		checkError(errors.Wrap(err, outDir))
		if existed {
			empty, err := pathutil.IsEmpty(outDir)
			checkError(errors.Wrap(err, outDir))
			if !empty {
				if force {
					if verbose {
						log.Infof("removing old output directory: %s", outDir)
					}
					checkError(os.RemoveAll(outDir))
				} else {
					checkError(fmt.Errorf("%s not empty: %s, use --force to overwrite", logname, outDir))
				}
			} else {
				checkError(os.RemoveAll(outDir))
			}
		}
		checkError(os.MkdirAll(outDir, 0777))
	} else {
		log.Errorf("%s should not be current directory", logname)
	}
}

func getFileListFromDir(path string, pattern *regexp.Regexp, threads int) ([]string, error) {
	files := make([]string, 0, 512)
	ch := make(chan string, threads)
	done := make(chan int)
	go func() {
		for file := range ch {
			files = append(files, file)
		}
		done <- 1
	}()

	cwalk.NumWorkers = threads
	err := cwalk.WalkWithSymlinks(path, func(_path string, info os.FileInfo, err error) error {
		if err != nil {
			return err
		}
		if !info.IsDir() && pattern.MatchString(info.Name()) {
			ch <- filepath.Join(path, _path)
		}
		return nil
	})
	close(ch)
	<-done
	if err != nil {
		return nil, err
	}

	return files, err
}

var defaultExts = []string{".gz", ".xz", ".zst", ".bz"}

func filepathTrimExtension(file string, suffixes []string) (string, string, string) {
	if suffixes == nil {
		suffixes = defaultExts
	}

	var e, e1, e2 string
	f := strings.ToLower(file)
	for _, s := range suffixes {
		// e = strings.ToLower(s)
		e = s
		if strings.HasSuffix(f, e) {
			e2 = e
			file = file[0 : len(file)-len(e)]
			break
		}
	}

	e1 = filepath.Ext(file)
	name := file[0 : len(file)-len(e1)]

	return name, e1, e2
}

func stringSplitN(s string, sep string, n int, a *[]string) {
	if a == nil {
		tmp := make([]string, n)
		a = &tmp
	}

	n--
	i := 0
	for i < n {
		m := strings.Index(s, sep)
		if m < 0 {
			break
		}
		(*a)[i] = s[:m]
		s = s[m+len(sep):]
		i++
	}
	(*a)[i] = s

	(*a) = (*a)[:i+1]
}

func bytesSplitN(s []byte, sep []byte, n int, a *[][]byte) {
	if a == nil {
		tmp := make([][]byte, n)
		a = &tmp
	}

	n--
	i := 0
	for i < n {
		m := bytes.Index(s, sep)
		if m < 0 {
			break
		}
		(*a)[i] = s[:m]
		s = s[m+len(sep):]
		i++
	}
	(*a)[i] = s

	(*a) = (*a)[:i+1]
}

func stringSplitNByByte(s string, sep byte, n int, a *[]string) {
	if a == nil {
		tmp := make([]string, n)
		a = &tmp
	}

	n--
	i := 0
	for i < n {
		m := strings.IndexByte(s, sep)
		if m < 0 {
			break
		}
		(*a)[i] = s[:m]
		s = s[m+1:]
		i++
	}
	(*a)[i] = s

	(*a) = (*a)[:i+1]
}

type Uint64Slice []uint64

func (s Uint64Slice) Len() int           { return len(s) }
func (s Uint64Slice) Less(i, j int) bool { return s[i] < s[j] }
func (s Uint64Slice) Swap(i, j int)      { s[i], s[j] = s[j], s[i] }

func (s *Uint64Slice) Push(x interface{}) {
	*s = append(*s, x.(uint64))
}

func (s *Uint64Slice) Pop() interface{} {
	old := *s
	n := len(old)
	x := old[n-1]
	*s = old[0 : n-1]
	return x
}

// Note: set should not have duplicates
func Combinations2(set []uint64) [][2]uint64 {
	if len(set) < 2 {
		return nil
	}
	var i, j int
	var n = len(set)
	var np1 = n - 1
	comb := make([][2]uint64, 0, n*(n-1)/2)
	for i = 0; i < np1; i++ {
		for j = i + 1; j < n; j++ {
			comb = append(comb, [2]uint64{set[i], set[j]})
		}
	}
	return comb
}

func IntSlice2StringSlice(vals []int) []string {
	s := make([]string, len(vals))
	for i, v := range vals {
		s[i] = strconv.Itoa(v)
	}
	return s
}

func MeanStdev(values []float64) (float64, float64) {
	n := len(values)

	if n == 0 {
		return 0, 0
	}

	if n == 1 {
		return values[0], 0
	}

	var sum float64
	for _, v := range values {
		sum += v
	}

	mean := sum / float64(n)

	var variance float64
	for _, v := range values {
		variance += (v - mean) * (v - mean)
	}

	return mean, math.Sqrt(variance / float64(n))
}

func wrapByteSlice(s []byte, width int, buffer *bytes.Buffer) ([]byte, *bytes.Buffer) {
	if width < 1 {
		return s, buffer
	}
	l := len(s)
	if l == 0 {
		return s, buffer
	}

	var lines int
	if l%width == 0 {
		lines = l/width - 1
	} else {
		lines = int(l / width)
	}

	if buffer == nil {
		buffer = bytes.NewBuffer(make([]byte, 0, l+lines))
	} else {
		buffer.Reset()
	}

	var start, end int
	for i := 0; i <= lines; i++ {
		start = i * width
		end = (i + 1) * width
		if end > l {
			end = l
		}

		buffer.Write(s[start:end])
		if i < lines {
			buffer.Write(_mark_newline)
		}
	}
	return buffer.Bytes(), buffer
}

var _mark_fasta = []byte{'>'}
var _mark_fastq = []byte{'@'}
var _mark_plus_newline = []byte{'+', '\n'}
var _mark_newline = []byte{'\n'}

func randInts(start int, end int, n int) []int {
	if start > end {
		start, end = end, start
	}
	ns := make([]int, n)
	w := end - start
	for i := 0; i < n; i++ {
		ns[i] = rand.Intn(w) + start
	}
	return ns
}

func readKVs(file string, ignoreCase bool) (map[string]string, error) {
	fh, err := xopen.Ropen(file)
	if err != nil {
		return nil, err
	}

	m := make(map[string]string, 1024)

	items := make([]string, 2)
	scanner := bufio.NewScanner(fh)
	var line string
	for scanner.Scan() {
		line = strings.TrimRight(scanner.Text(), "\r\n")
		if line == "" {
			continue
		}

		stringSplitNByByte(line, '\t', 2, &items)
		if len(items) < 2 {
			continue
		}

		if ignoreCase {
			m[strings.ToLower(items[0])] = items[1]
		} else {
			m[items[0]] = items[1]
		}
	}
	if err = scanner.Err(); err != nil {
		return nil, err
	}

	return m, fh.Close()
}
