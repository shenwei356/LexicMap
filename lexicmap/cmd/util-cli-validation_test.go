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
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/spf13/cobra"
	"github.com/spf13/pflag"
)

// Exercise every registered numeric flag through the actual CLI preflight.
// Invalid options must fail before opening an existing log or result file.
func TestNumericFlagsAcrossCommands(t *testing.T) {
	var walk func(*cobra.Command)
	walk = func(cmd *cobra.Command) {
		if cmd.Run != nil || cmd.RunE != nil {
			cmd.Flags().VisitAll(func(flag *pflag.Flag) {
				var values []string
				switch flag.Value.Type() {
				case "float64", "float32":
					values = []string{"NaN", "+Inf", "-Inf", "1e309"}
				case "int":
					if flag.Name != "seed" && flag.Name != "rand-seed" {
						values = []string{"9223372036854775807"}
					}
				case "string":
					switch flag.Name {
					case "buffer-size", "max-seed-memory", "max-align-result-memory", "max-genome-cache-memory", "seed-index2-threshold", "threshold":
						values = []string{"NaN", "Inf", "8E18K", "-1G"}
					}
				}
				for _, value := range values {
					args := append(strings.Fields(cmd.CommandPath())[1:], "--"+flag.Name+"="+value)
					t.Run(strings.Join(args, "/"), func(t *testing.T) {
						dir := t.TempDir()
						logFile := writeTestInput(t, filepath.Join(dir, "run.log"), "preserve\n")
						args = append(args, "--quiet", "--log", logFile)
						files := []string{logFile}
						if cmd.Flags().Lookup("out-file") != nil {
							output := writeTestInput(t, filepath.Join(dir, "output"), "preserve\n")
							args = append(args, "--out-file", output)
							files = append(files, output)
						}
						out, stderr, err := runAuditCommand(t, "", args...)
						message := string(out) + string(stderr)
						if err == nil || !strings.Contains(message, "--"+flag.Name) || strings.Contains(message, "panic:") {
							t.Errorf("invalid number not rejected: %v, %s", err, message)
						}
						for _, file := range files {
							if got, err := os.ReadFile(file); err != nil || string(got) != "preserve\n" {
								t.Errorf("file changed before validation: %s: %q, %v", file, got, err)
							}
						}
					})
				}
			})
		}
		for _, child := range cmd.Commands() {
			walk(child)
		}
	}
	walk(RootCmd)
}

func TestNumericFlagResourceLimits(t *testing.T) {
	for _, args := range [][]string{
		{"version", "--threads", "1025"},
		{"search", "--max-query-conc", "1025"},
		{"index", "--seed-data-threads", "1025"},
		{"index", "--masks", "1048577"},
		{"index", "--partitions", "1048577"},
		{"search", "--max-open-files", "1048577"},
		{"genome", "search", "--windows", "1048577"},
		{"utils", "seed-pos", "--bins", "1048577"},
		{"utils", "seed-pos", "--slid-step", "4294967296"},
		{"utils", "seed-pos", "--width", "1e100"},
		{"utils", "seed-pos", "--width", "100", "--height", "100"},
		{"utils", "2sam", "--buffer-size", "2G"},
		{"utils", "2blast", "--buffer-size", "0"},
		{"genome", "compare", "--min-ani", "100.1"},
		{"genome", "search", "--min-af", "100.1"},
		{"genome", "pair", "--prob-threshold", "1.01"},
		{"utils", "subseq", "--region", "1:9223372036854775807"},
	} {
		t.Run(strings.Join(args, "/"), func(t *testing.T) {
			logFile := writeTestInput(t, filepath.Join(t.TempDir(), "run.log"), "preserve\n")
			out, stderr, err := runAuditCommand(t, "", append(args, "--quiet", "--log", logFile)...)
			if err == nil || strings.Contains(string(out)+string(stderr), "panic:") {
				t.Errorf("invalid resource limit accepted: %v, %s, %s", err, out, stderr)
			}
			if got, err := os.ReadFile(logFile); err != nil || string(got) != "preserve\n" {
				t.Errorf("log changed: %q, %v", got, err)
			}
		})
	}
}

// Large budgets and random seeds are valid; they do not allocate that many bytes
// or workers. Zero retains its documented auto/disabled meaning.
func TestNumericFlagsRetainLargeBudgets(t *testing.T) {
	cmd := &cobra.Command{}
	cmd.Flags().Int("threads", 0, "")
	cmd.Flags().Int("max-query-conc", 0, "")
	cmd.Flags().Int("rand-seed", math.MaxInt, "")
	cmd.Flags().Float64("max-evalue", math.MaxFloat64, "")
	cmd.Flags().Float64("min-af", 100, "")
	cmd.Flags().Float64("min-mask-fraction", 1, "")
	cmd.Flags().String("max-genome-cache-memory", strconv.FormatInt(math.MaxInt64, 10), "")
	cmd.Flags().String("max-seed-memory", "0", "")
	cmd.Flags().String("buffer-size", "2147483647", "")
	if err := checkNumericFlags(cmd); err != nil {
		t.Fatal(err)
	}
	// A positive Mb value must be rejected before multiplication can overflow int.
	cmd.Flags().Int("max-genome-size", math.MaxInt/1_000_000+1, "")
	if err := checkNumericFlags(cmd); err == nil {
		t.Fatal("overflowing Mb-to-base conversion accepted")
	}
}
