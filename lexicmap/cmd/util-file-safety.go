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
	"runtime"
	"strings"

	"github.com/spf13/cobra"
)

// commandFileListKey scopes a preflight-read file list to one CLI invocation.
// In particular, -X - must consume stdin only once.
type commandFileListKey struct{}

// commandDirectoryListKey freezes directory discovery before outputs are created.
// Otherwise a new log matching --file-regexp would become an input genome.
type commandDirectoryListKey struct{}

// commandFileList retains the parsed paths until the command consumes them.
type commandFileList struct {
	path  string
	files []string
}

// resolveOutputPath resolves existing parent symlinks even before a log or its
// containing output directory exists.
func resolveOutputPath(path string) (string, error) {
	paths, err := filePathVariants(path)
	if err != nil {
		return "", err
	}
	return paths[len(paths)-1], nil
}

// cleanDirectoryPath resolves .. after following parent symlinks. Retain the
// final entry so --force still removes a leaf symlink rather than its target.
func cleanDirectoryPath(path string) (string, error) {
	separator := func(c rune) bool { return c < 128 && os.IsPathSeparator(uint8(c)) }
	hasParent := false
	for _, part := range strings.FieldsFunc(path, separator) {
		hasParent = hasParent || part == ".."
	}
	if !hasParent {
		return filepath.Clean(path), nil
	}
	path = strings.TrimRightFunc(path, separator)
	last := strings.LastIndexFunc(path, separator)
	leaf := path[last+1:]
	if leaf == "." || leaf == ".." {
		return resolveOutputPath(path)
	}
	parent, err := resolveOutputPath(path[:last+1])
	if err != nil {
		return "", err
	}
	return filepath.Join(parent, leaf), nil
}

// directoryEntryPath resolves the parent but preserves the entry that RemoveAll
// removes. A final symlink and its target need separate containment checks.
func directoryEntryPath(path string) (string, error) {
	clean, err := cleanDirectoryPath(path)
	if err != nil {
		return "", err
	}
	parent, err := resolveOutputPath(filepath.Dir(clean))
	if err != nil {
		return "", err
	}
	return filepath.Join(parent, filepath.Base(clean)), nil
}

// filePathVariants preserves every symlink entry traversed on the way to a file.
// Directory cleanup can remove an intermediate entry even if the target is outside.
func filePathVariants(path string) ([]string, error) {
	absolute, err := filepath.Abs(path)
	if err != nil {
		return nil, err
	}
	paths := []string{absolute}
	// Keep .. components until preceding links have been followed.
	if filepath.IsAbs(path) {
		absolute = path
	} else if filepath.VolumeName(path) == "" {
		cwd, err := os.Getwd()
		if err != nil {
			return nil, err
		}
		absolute = cwd + string(filepath.Separator) + path
	}
	volume := filepath.VolumeName(absolute)
	parent := volume + string(filepath.Separator)
	components := func(p string) []string {
		return strings.FieldsFunc(p, func(c rune) bool { return c < 128 && os.IsPathSeparator(uint8(c)) })
	}
	pending := components(absolute[len(volume):])
	links := 0
	for len(pending) > 0 {
		part := pending[0]
		pending = pending[1:]
		if part == "." {
			continue
		}
		if part == ".." {
			parent = filepath.Dir(parent)
			continue
		}
		named := filepath.Join(parent, part)
		info, err := os.Lstat(named)
		if os.IsNotExist(err) {
			// Auxiliary outputs may not exist yet; their existing parents matter.
			// Continue through .. and later symlinks after a missing component.
			parent = named
			continue
		}
		if err != nil {
			return nil, err
		}
		if info.Mode()&os.ModeSymlink == 0 {
			if len(pending) > 0 && !info.IsDir() {
				return nil, fmt.Errorf("path component %s is not a directory", named)
			}
			parent = named
			continue
		}
		paths = append(paths, named)
		links++
		if links > 255 {
			return nil, fmt.Errorf("too many symbolic links in %s", path)
		}
		target, err := os.Readlink(named)
		if err != nil {
			return nil, err
		}
		// Resolve one link at a time; EvalSymlinks on the parent would hide
		// entries in chains such as alias/db/inputs -> outside -> elsewhere.
		v := filepath.VolumeName(target)
		if filepath.IsAbs(target) {
			parent = v + string(filepath.Separator)
		} else if len(target) > len(v) && os.IsPathSeparator(target[len(v)]) {
			parent = filepath.VolumeName(parent) + string(filepath.Separator)
		}
		pending = append(components(target[len(v):]), pending...)
	}
	paths = append(paths, parent)
	return paths, nil
}

// checkFileOutsideDirectories protects auxiliary files from directory cleanup.
func checkFileOutsideDirectories(path, label string, dirs ...string) error {
	if path == "" {
		return nil
	}
	paths, err := filePathVariants(path)
	if err != nil {
		return err
	}
	for _, dir := range dirs {
		absoluteDir, err := directoryEntryPath(dir)
		if err != nil {
			return err
		}
		realDir, err := resolveOutputPath(dir)
		if err != nil {
			return err
		}
		for _, variant := range paths {
			for _, directory := range []string{absoluteDir, realDir} {
				rel, err := filepath.Rel(directory, variant)
				if err == nil && rel != ".." && !strings.HasPrefix(rel, ".."+string(filepath.Separator)) {
					return fmt.Errorf("%s path %s must be outside output directory %s", label, path, dir)
				}
			}
		}
	}
	return nil
}

// checkCommandFiles runs before logging or output creation. Paths named as
// genome IDs are not files; explicit auxiliary inputs and file lists are.
func checkCommandFiles(cmd *cobra.Command, args []string) error {
	// All later Clean/Join calls must use the same directory as preflight. Do
	// not change Flag.Changed: commands use it to distinguish defaults from input.
	var originalInputDirs []string
	for _, name := range []string{"out-dir", "plot-dir", "index", "in-dir", "taxdump"} {
		f := cmd.Flags().Lookup(name)
		if f == nil || f.Value.String() == "" {
			continue
		}
		if name == "index" || name == "in-dir" || name == "taxdump" {
			originalInputDirs = append(originalInputDirs, f.Value.String())
		}
		path, err := cleanDirectoryPath(f.Value.String())
		if err != nil {
			return err
		}
		if err := f.Value.Set(path); err != nil {
			return err
		}
	}
	ctx := cmd.Context()
	if ctx == nil {
		ctx = context.Background()
	}
	// Commands may be executed repeatedly in-process; do not reuse an old list.
	ctx = context.WithValue(ctx, commandFileListKey{}, commandFileList{})
	cmd.SetContext(context.WithValue(ctx, commandDirectoryListKey{}, commandFileList{}))
	flag := func(name string) string {
		if cmd.Flags().Lookup(name) == nil {
			return ""
		}
		value, _ := cmd.Flags().GetString(name)
		return value
	}
	type outputPath struct {
		path, absolute, entry string
		info                  os.FileInfo
		directory             bool
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
		entry := absolute
		if directory {
			entry, err = directoryEntryPath(path)
			if err != nil {
				return err
			}
		}
		outputs = append(outputs, outputPath{path, absolute, entry, info, directory})
		return nil
	}
	for _, name := range []string{"out-file", "file", "log", "big-genomes", "out-dir", "plot-dir"} {
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
	// Output and temporary directories must remain independent throughout
	// cleanup and publication, even when either is reached through a symlink.
	for i, dir := range replacedDirs {
		for _, other := range replacedDirs[:i] {
			if err := checkFileOutsideDirectories(dir, "output directory", other); err != nil {
				return err
			}
			if err := checkFileOutsideDirectories(other, "output directory", dir); err != nil {
				return err
			}
		}
	}
	// Logs and skipped-genome reports open before directory cleanup. Reject
	// containment before any file or directory can be modified.
	// Preserve the more specific input-alias errors when both checks apply.
	var outputDirectoryErr error
	for _, name := range []string{"log", "big-genomes"} {
		if err := checkFileOutsideDirectories(flag(name), name, replacedDirs...); err != nil {
			outputDirectoryErr = err
			break
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
				if output.directory && output.info != nil && info.Mode().IsRegular() {
					// A redirected file remains readable after unlink, so comparing
					// only the output directory's inode would silently lose the input.
					// Walk does not follow links: cleanup removes those entries only.
					err := filepath.Walk(output.path, func(path string, entry os.FileInfo, err error) error {
						if err != nil {
							return err
						}
						if entry.Mode().IsRegular() && os.SameFile(info, entry) {
							return fmt.Errorf("output directory %s conflicts with stdin file %s", output.path, path)
						}
						return nil
					})
					if err != nil {
						return err
					}
				}
			}
			return nil
		}
		paths, err := filePathVariants(input)
		if err != nil {
			return err
		}
		var info os.FileInfo
		statDone := false
		for _, output := range outputs {
			same := false
			for _, path := range paths {
				if path == output.absolute {
					same = true
					break
				}
			}
			if !same && output.directory {
				for _, path := range paths {
					for _, dir := range []string{output.entry, output.absolute} {
						rel, err := filepath.Rel(dir, path)
						if err == nil && rel != ".." && !strings.HasPrefix(rel, ".."+string(filepath.Separator)) {
							same = true
							break
						}
					}
					if same {
						break
					}
				}
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

	// Normalization must not hide an input symlink removed during cleanup.
	for _, dir := range originalInputDirs {
		if err := check(dir); err != nil {
			return err
		}
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
	// Trust ordinary input paths with -S; keep auxiliary and stdin checks.
	// Directory discovery is explicit via --in-dir.
	skipFileCheck, _ := cmd.Flags().GetBool("skip-file-check")
	checkFileInput := func(path string) error {
		if skipFileCheck && !isStdin(path) {
			return nil
		}
		if err := checkInput(path); err != nil {
			return err
		}
		if cmd.Name() == "index" && !isStdin(path) {
			if info, err := os.Stat(path); err == nil && info.IsDir() {
				return fmt.Errorf("input %s is a directory; use -I/--in-dir to read sequence files from a directory", path)
			}
		}
		return nil
	}
	if positionalFiles {
		for _, path := range args {
			if err := checkFileInput(path); err != nil {
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
				if err := checkFileInput(path); err != nil {
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
	var remergeOutputs []string // directories replaced before recovery inputs are consumed
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
		remergeOutputs = append(remergeOutputs, index)
		batchPattern := regexp.MustCompile(`^batch_\d+$`)
		roundPattern := regexp.MustCompile(`^r\d+_b\d+$`)
		for _, entry := range entries {
			if roundPattern.MatchString(entry.Name()) {
				remergeOutputs = append(remergeOutputs, filepath.Join(tmpDir, entry.Name()))
			}
			if entry.IsDir() && batchPattern.MatchString(entry.Name()) {
				inputIndexes = append(inputIndexes, filepath.Join(tmpDir, entry.Name()))
			}
		}
	}
	checkIndexInput := func(path, sourceDir string) error {
		if sourceDir != index {
			if err := checkFileOutsideDirectories(path, "input index", remergeOutputs...); err != nil {
				return err
			}
		}
		return check(path)
	}
	for _, dir := range inputIndexes {
		// Directory outputs such as --plot-dir must not remove an input index.
		if err := checkIndexInput(dir, dir); err != nil {
			return err
		}
		for _, name := range []string{FileInfo, FileMasks, FileGenomeIndex, FileGenomeChunks, FileGenomeDetails, FileGenomeIndex + ".bak", FileGenomeIndex + ".new"} {
			if err := checkIndexInput(filepath.Join(dir, name), dir); err != nil {
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
	// Scan genome inputs before logging or directory cleanup and reuse this exact
	// list in index. New logs/reports must never be discovered as input genomes.
	if dir := flag("in-dir"); cmd.Name() == "index" && dir != "" {
		pattern, err := regexp.Compile("(?i)" + flag("file-regexp"))
		if err != nil {
			return err
		}
		threads, _ := cmd.Flags().GetInt("threads")
		if threads == 0 {
			threads = runtime.NumCPU()
		}
		files, err := getFileListFromDir(dir, pattern, threads)
		if err != nil {
			return err
		}
		if !skipFileCheck {
			for _, path := range files {
				if err := check(path); err != nil {
					return err
				}
			}
		}
		cmd.SetContext(context.WithValue(cmd.Context(), commandDirectoryListKey{}, commandFileList{dir, files}))
	}
	// Missing outputs cannot alias existing directory inputs. Avoid traversing
	// large indexes on the usual first write; existing outputs need identity checks
	// to catch hard links whose names lie outside the input directory.
	existingOutput := false
	for _, output := range outputs {
		existingOutput = existingOutput || output.info != nil
	}
	if !existingOutput && len(remergeOutputs) == 0 {
		return outputDirectoryErr
	}
	walk := func(dir string, pattern *regexp.Regexp, sourceDir string) error {
		if dir == "" {
			return nil
		}
		if err := checkIndexInput(dir, sourceDir); err != nil {
			return err
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
			if err := checkIndexInput(file, sourceDir); err != nil {
				return err
			}
		}
		return nil
	}
	for _, dir := range inputIndexes {
		allFiles := regexp.MustCompile(".*")
		for _, name := range []string{DirSeeds, DirGenomes} {
			if err := walk(filepath.Join(dir, name), allFiles, dir); err != nil {
				return err
			}
		}
	}

	return outputDirectoryErr
}
