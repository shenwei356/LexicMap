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
	"context"
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"strings"

	"github.com/spf13/cobra"
)

// commandFileListKey scopes a preflight-read file list to one CLI invocation.
// In particular, -X - must consume stdin only once.
type commandFileListKey struct{}

// commandFileList retains the parsed paths until the command consumes them.
type commandFileList struct {
	path  string
	files []string
}

// checkCommandFiles runs before logging or output creation. Paths named as
// genome IDs are not files; explicit auxiliary inputs and file lists are.
func checkCommandFiles(cmd *cobra.Command, args []string) error {
	ctx := cmd.Context()
	if ctx == nil {
		ctx = context.Background()
	}
	// Commands may be executed repeatedly in-process; do not reuse an old list.
	cmd.SetContext(context.WithValue(ctx, commandFileListKey{}, commandFileList{}))
	flag := func(name string) string {
		if cmd.Flags().Lookup(name) == nil {
			return ""
		}
		value, _ := cmd.Flags().GetString(name)
		return value
	}
	type outputPath struct {
		path, absolute string
		info           os.FileInfo
		directory      bool
	}
	var outputs []outputPath
	for _, name := range []string{"out-file", "log", "big-genomes", "out-dir"} {
		path := flag(name)
		if path == "" || (path == "-" && name == "out-file") {
			continue
		}
		absolute, err := filepath.Abs(path)
		if err != nil {
			return err
		}
		info, err := os.Stat(path)
		if err != nil && !os.IsNotExist(err) {
			return err
		}
		if real, err := filepath.EvalSymlinks(absolute); err == nil {
			absolute = real
		}
		outputs = append(outputs, outputPath{path, absolute, info, name == "out-dir"})
	}
	if len(outputs) == 0 {
		return nil
	}
	for i, a := range outputs {
		for _, b := range outputs[:i] {
			if a.directory || b.directory {
				continue
			}
			if a.absolute == b.absolute || (a.info != nil && b.info != nil && os.SameFile(a.info, b.info)) {
				return fmt.Errorf("output paths should not be the same: %s and %s", a.path, b.path)
			}
		}
	}

	check := func(input string) error {
		if input == "" {
			return nil
		}
		if input == "-" {
			info, err := os.Stdin.Stat()
			if err != nil {
				return err
			}
			for _, output := range outputs {
				if output.info != nil && os.SameFile(output.info, info) {
					return fmt.Errorf("output path %s conflicts with stdin", output.path)
				}
			}
			return nil
		}
		absolute, err := filepath.Abs(input)
		if err != nil {
			return err
		}
		var info os.FileInfo
		statDone := false
		for _, output := range outputs {
			same := absolute == output.absolute
			if !same && output.directory {
				real := absolute
				if resolved, err := filepath.EvalSymlinks(real); err == nil {
					real = resolved
				}
				rel, err := filepath.Rel(output.absolute, real)
				same = err == nil && rel != ".." && !strings.HasPrefix(rel, ".."+string(filepath.Separator))
			}
			if !same && output.info != nil {
				if !statDone {
					info, err = os.Stat(input)
					if err != nil && !os.IsNotExist(err) {
						return err
					}
					statDone = true
				}
				same = info != nil && os.SameFile(output.info, info)
			}
			if same {
				return fmt.Errorf("output path %s conflicts with input %s", output.path, input)
			}
		}
		return nil
	}

	for _, name := range []string{"infile-list", "mask-file", "taxid-file", "genome2taxid", "kv-file-seq", "kv-file-genome", "search-result"} {
		if err := check(flag(name)); err != nil {
			return err
		}
	}
	if cmd.Flags().Lookup("pair-file") != nil {
		pairs, _ := cmd.Flags().GetStringSlice("pair-file")
		for _, path := range pairs {
			if err := check(path); err != nil {
				return err
			}
		}
	}

	// Only these commands interpret positional arguments and -X entries as files.
	usesList, positionalFiles := false, false
	switch cmd.Name() {
	case "index", "search", "2sam", "2blast", "merge-search-results":
		usesList, positionalFiles = true, true
	case "compare":
		usesList, positionalFiles = true, flag("index") == ""
	}
	if positionalFiles {
		if len(args) == 0 && flag("infile-list") == "" && flag("in-dir") == "" {
			if err := check("-"); err != nil {
				return err
			}
		}
		for _, path := range args {
			if err := check(path); err != nil {
				return err
			}
		}
	}
	if list := flag("infile-list"); usesList && list != "" {
		files, err := getFileListFromFile(list, false)
		if err != nil {
			return err
		}
		ctx := cmd.Context()
		if ctx == nil {
			ctx = context.Background()
		}
		cmd.SetContext(context.WithValue(ctx, commandFileListKey{}, commandFileList{list, files}))
		if positionalFiles {
			for _, path := range files {
				if err := check(path); err != nil {
					return err
				}
			}
		}
	}

	// Protect index metadata, including files created lazily by utilities.
	index := flag("index")
	if index != "" {
		for _, name := range []string{FileInfo, FileMasks, FileGenomeIndex, FileGenomeChunks, FileGenomeDetails, FileGenomeIndex + ".bak", FileGenomeIndex + ".new"} {
			if err := check(filepath.Join(index, name)); err != nil {
				return err
			}
		}
	}
	if cmd.Flags().Lookup("taxdump") != nil {
		dir := flag("taxdump")
		if dir == "" && index != "" {
			dir = filepath.Join(index, "taxdump")
		}
		if dir != "" {
			for _, name := range []string{"nodes.dmp", "names.dmp"} {
				if err := check(filepath.Join(dir, name)); err != nil {
					return err
				}
			}
			if flag("genome2taxid") == "" {
				if err := check(filepath.Join(dir, "taxid.map")); err != nil {
					return err
				}
			}
		}
	}
	// Missing outputs cannot alias existing directory inputs. Avoid traversing
	// large indexes on the usual first write; existing outputs need identity checks
	// to catch hard links whose names lie outside the input directory.
	existingOutput := false
	for _, output := range outputs {
		existingOutput = existingOutput || output.info != nil
	}
	if !existingOutput {
		return nil
	}
	walk := func(dir string, pattern *regexp.Regexp) error {
		if dir == "" {
			return nil
		}
		root, err := filepath.EvalSymlinks(dir)
		if os.IsNotExist(err) {
			return nil
		} // normal command validation reports missing inputs
		if err != nil {
			return err
		}
		// Follow linked batch directories just as input discovery does.
		files, err := getFileListFromDir(root, pattern, 1)
		if err != nil {
			return err
		}
		for _, file := range files {
			if err := check(file); err != nil {
				return err
			}
		}
		return nil
	}
	if index != "" {
		allFiles := regexp.MustCompile(".*")
		for _, name := range []string{DirSeeds, DirGenomes} {
			if err := walk(filepath.Join(index, name), allFiles); err != nil {
				return err
			}
		}
	}
	if dir := flag("in-dir"); cmd.Name() == "index" && dir != "" {
		pattern, err := regexp.Compile("(?i)" + flag("file-regexp"))
		if err != nil {
			return err
		}
		if err = walk(dir, pattern); err != nil {
			return err
		}
	}
	return nil
}
