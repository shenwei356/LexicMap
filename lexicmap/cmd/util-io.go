// Copyright © 2023-2025 Wei Shen <shenwei356@gmail.com>
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
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"

	gzip "github.com/klauspost/pgzip"
)

// BufferSize is size of buffer
var BufferSize = 65536 // os.Getpagesize()

func outStream(file string, gzipped bool, level int) (*bufio.Writer, io.WriteCloser, *os.File, error) {
	var w *os.File
	if file == "-" {
		w = os.Stdout
	} else {
		dir := filepath.Dir(file)
		fi, err := os.Stat(dir)
		if err == nil && !fi.IsDir() {
			return nil, nil, nil, fmt.Errorf("can not write file into a non-directory path: %s", dir)
		}
		if os.IsNotExist(err) {
			os.MkdirAll(dir, 0755)
		}

		w, err = os.Create(file)
		if err != nil {
			return nil, nil, nil, fmt.Errorf("fail to write %s: %s", file, err)
		}
	}

	if gzipped {
		// gw := gzip.NewWriter(w)
		gw, err := gzip.NewWriterLevel(w, level)
		if err != nil {
			return nil, nil, nil, fmt.Errorf("fail to write %s: %s", file, err)
		}
		return bufio.NewWriterSize(gw, BufferSize), gw, w, nil
	}
	return bufio.NewWriterSize(w, BufferSize), nil, w, nil
}

func inStream(file string) (*bufio.Reader, *os.File, bool, error) {
	var err error
	var r *os.File
	var gzipped bool
	if file == "-" {
		if !detectStdin() {
			return nil, nil, gzipped, errors.New("stdin not detected")
		}
		r = os.Stdin
	} else {
		r, err = os.Open(file)
		if err != nil {
			return nil, nil, gzipped, fmt.Errorf("fail to read %s: %s", file, err)
		}
	}

	br := bufio.NewReaderSize(r, BufferSize)

	if gzipped, err = isGzip(br); err != nil {
		return nil, nil, gzipped, fmt.Errorf("fail to check is file (%s) gzipped: %s", file, err)
	} else if gzipped {
		// gr, err := gzip.NewReader(br)
		gr, err := gzip.NewReaderN(br, 65536, 8)
		if err != nil {
			return nil, r, gzipped, fmt.Errorf("fail to create gzip reader for %s: %s", file, err)
		}
		br = bufio.NewReaderSize(gr, BufferSize)
	}
	return br, r, gzipped, nil
}

func isGzip(b *bufio.Reader) (bool, error) {
	return checkBytes(b, []byte{0x1f, 0x8b})
}

func checkBytes(b *bufio.Reader, buf []byte) (bool, error) {
	m, err := b.Peek(len(buf))
	if err != nil {
		return false, fmt.Errorf("no content")
	}
	for i := range buf {
		if m[i] != buf[i] {
			return false, nil
		}
	}
	return true, nil
}

func detectStdin() bool {
	// http://stackoverflow.com/a/26567513
	stat, err := os.Stdin.Stat()
	if err != nil {
		return false
	}
	return (stat.Mode() & os.ModeCharDevice) == 0
}
