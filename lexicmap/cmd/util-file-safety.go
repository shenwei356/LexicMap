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

// resolveOutputPath resolves existing parent symlinks even before a log or its
// containing output directory exists.
func resolveOutputPath(path string) (string, error) {
	absolute, err := filepath.Abs(path)
	if err != nil {
		return "", err
	}
	parent, tail := absolute, ""
	for {
		real, err := filepath.EvalSymlinks(parent)
		if err == nil {
			return filepath.Join(real, tail), nil
		}
		if !os.IsNotExist(err) || filepath.Dir(parent) == parent {
			return "", err
		}
		tail = filepath.Join(filepath.Base(parent), tail)
		parent = filepath.Dir(parent)
	}
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
	var replacedDirs []string
	addOutput := func(path string, directory bool) error {
		if directory {
			if err := checkOutputDirectory(path, path); err != nil {
				return err
			}
			replacedDirs = append(replacedDirs, path)
		}
		absolute, err := resolveOutputPath(path)
		if err != nil {
			return err
		}
		info, err := os.Stat(path)
		if err != nil && !os.IsNotExist(err) {
			return err
		}
		outputs = append(outputs, outputPath{path, absolute, info, directory})
		return nil
	}
	for _, name := range []string{"out-file", "log", "big-genomes", "out-dir", "plot-dir"} {
		path := flag(name)
		if path == "" || (path == "-" && name == "out-file") {
			continue
		}
		if err := addOutput(path, name == "out-dir" || name == "plot-dir"); err != nil {
			return err
		}
	}
	// Multi-batch builds replace this directory before reading their inputs.
	// Protect it even when the final output directory does not exist yet.
	if dir := flag("out-dir"); cmd.Name() == "index" && dir != "" {
		if err := addOutput(filepath.Clean(dir)+ExtTmpDir, true); err != nil {
			return err
		}
	}
	// Remerge replaces --index and removes its temporary batches after publishing.
	// Those batches are also inputs, so validate these directories separately.
	if dir := flag("index"); cmd.Name() == "remerge" && dir != "" {
		for _, path := range []string{dir, filepath.Clean(dir) + ExtTmpDir} {
			if err := checkOutputDirectory(path, path); err != nil {
				return err
			}
			replacedDirs = append(replacedDirs, path)
		}
	}
	// Logs open before directory cleanup. Reject containment before either the
	// log or an existing directory can be modified, including symlink aliases.
	// Preserve the more specific input-alias errors when both checks apply.
	var logDirectoryErr error
	if path := flag("log"); path != "" {
		logPath, err := resolveOutputPath(path)
		if err != nil {
			return err
		}
		for _, dir := range replacedDirs {
			real, err := resolveOutputPath(dir)
			if err != nil {
				return err
			}
			rel, err := filepath.Rel(real, logPath)
			if err == nil && rel != ".." && !strings.HasPrefix(rel, ".."+string(filepath.Separator)) {
				logDirectoryErr = fmt.Errorf("log path %s must be outside output directory %s", path, dir)
				break
			}
		}
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

	// Count actual stdin consumers, not genome IDs named "-". Run even when
	// every output is stdout, and reject explicit reuse before reading -X -.
	stdinInputs := 0
	checkInput := func(path string) error {
		if isStdin(path) {
			stdinInputs++
			if stdinInputs > 1 {
				return fmt.Errorf("stdin (-) may only be used once across input files and auxiliary inputs")
			}
		}
		return check(path)
	}
	for _, name := range []string{"infile-list", "mask-file", "taxid-file", "genome2taxid", "kv-file-seq", "kv-file-genome", "search-result"} {
		if err := checkInput(flag(name)); err != nil {
			return err
		}
	}
	if cmd.Flags().Lookup("pair-file") != nil {
		pairs, _ := cmd.Flags().GetStringSlice("pair-file")
		for _, path := range pairs {
			if err := checkInput(path); err != nil {
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
		for _, path := range args {
			if err := checkInput(path); err != nil {
				return err
			}
		}
	}
	var listedFiles []string
	if list := flag("infile-list"); usesList && list != "" {
		files, err := getFileListFromFile(list, false)
		if err != nil {
			return err
		}
		listedFiles = files
		ctx := cmd.Context()
		if ctx == nil {
			ctx = context.Background()
		}
		cmd.SetContext(context.WithValue(ctx, commandFileListKey{}, commandFileList{list, files}))
		if positionalFiles {
			for _, path := range files {
				if err := checkInput(path); err != nil {
					return err
				}
			}
		}
	}
	// An empty file list falls back to implicit stdin, just like no list.
	if positionalFiles && len(args) == 0 && len(listedFiles) == 0 && flag("in-dir") == "" {
		if err := checkInput("-"); err != nil {
			return err
		}
	}

	// Protect index metadata, including files created lazily by utilities.
	index := flag("index")
	var inputIndexes []string
	if index != "" {
		inputIndexes = append(inputIndexes, index)
	}
	if cmd.Name() == "remerge" && index != "" {
		// Recovery reads the temporary batches, before the final index exists.
		// Discover these inputs before setupCommandLog can truncate any of them.
		tmpDir := filepath.Clean(index) + ExtTmpDir
		entries, err := os.ReadDir(tmpDir)
		if err != nil && !os.IsNotExist(err) {
			return err
		}
		batchPattern := regexp.MustCompile(`^batch_\d+$`)
		for _, entry := range entries {
			if entry.IsDir() && batchPattern.MatchString(entry.Name()) {
				inputIndexes = append(inputIndexes, filepath.Join(tmpDir, entry.Name()))
			}
		}
	}
	for _, dir := range inputIndexes {
		// Directory outputs such as --plot-dir must not remove an input index.
		if err := check(dir); err != nil {
			return err
		}
		for _, name := range []string{FileInfo, FileMasks, FileGenomeIndex, FileGenomeChunks, FileGenomeDetails, FileGenomeIndex + ".bak", FileGenomeIndex + ".new"} {
			if err := check(filepath.Join(dir, name)); err != nil {
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
		return logDirectoryErr
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
	for _, dir := range inputIndexes {
		allFiles := regexp.MustCompile(".*")
		for _, name := range []string{DirSeeds, DirGenomes} {
			if err := walk(filepath.Join(dir, name), allFiles); err != nil {
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
	return logDirectoryErr
}
