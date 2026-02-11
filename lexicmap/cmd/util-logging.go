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
	"fmt"
	"io"
	"os"
	"runtime"

	"github.com/mattn/go-colorable"
	"github.com/shenwei356/go-logging"
)

var log *logging.Logger

var logFormat = logging.MustStringFormatter(
	`%{time:15:04:05.000} %{color}[%{level:.4s}]%{color:reset} %{message}`,
)

var backendFormatter logging.Backend

func init() {
	var stderr io.Writer = os.Stderr
	if runtime.GOOS == "windows" {
		stderr = colorable.NewColorableStderr()
	}
	backend := logging.NewLogBackend(stderr, "", 0)
	backendFormatter = logging.NewBackendFormatter(backend, logFormat)

	logging.SetBackend(backendFormatter)

	log = logging.MustGetLogger("lexicmap")
}

func addLog(file string, verbose bool) *os.File {
	w, err := os.Create(file)
	if err != nil {
		checkError(fmt.Errorf("failed to write log file %s: %s", file, err))
	}

	var logFormat2 = logging.MustStringFormatter(
		`%{time:15:04:05.000} [%{level:.4s}] %{message}`,
	)
	backend := logging.NewLogBackend(w, "", 0)
	backendFormatter2 := logging.NewBackendFormatter(backend, logFormat2)

	if !verbose {
		logging.SetBackend(backendFormatter2)
	} else {
		logging.SetBackend(backendFormatter, backendFormatter2)
	}

	log = logging.MustGetLogger("lexicmap")

	return w
}
