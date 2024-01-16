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

package index

import (
	"bufio"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"io/fs"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"sync"

	"github.com/shenwei356/LexicMap/lexicmap/index/twobit"
	"github.com/shenwei356/LexicMap/lexicmap/tree"
	"github.com/shenwei356/lexichash"
	"github.com/shenwei356/util/pathutil"
	"github.com/shenwei356/xopen"
)

// var Magic = [8]byte{'l', 'e', 'x', 'i', 'c', 'i', 'd', 'x'}

// MainVersion is use for checking compatibility
var MainVersion uint8 = 0

// MinorVersion is less important
var MinorVersion uint8 = 1

// ErrInvalidFileFormat means invalid file format.
var ErrInvalidFileFormat = errors.New("lexichash index: invalid binary format")

// ErrBrokenFile means the file is not complete.
var ErrBrokenFile = errors.New("lexichash index: broken file")

// ErrBrokenGenomeInfoFile means the genome info file is not complete.
var ErrBrokenGenomeInfoFile = errors.New("lexichash index: broken genome info file")

// ErrVersionMismatch means version mismatch between files and program.
var ErrVersionMismatch = errors.New("lexichash index: version mismatch")

// ErrDirNotEmpty means the output directory is not empty.
var ErrDirNotEmpty = errors.New("lexichash index: output directory not empty")

// ErrDirNotEmpty means the output directory is not empty.
var ErrPWDAsOutDir = errors.New("lexichash index: current directory cant't be the output dir")

// ErrInvalidIndexDir means the path is not a valid index directory.
var ErrInvalidIndexDir = errors.New("lexichash index: invalid index directory")

// ErrTreeFileMissing means some tree files are missing.
var ErrTreeFileMissing = errors.New("lexichash index: some tree files missing")

// IDListFile defines the name of the ID list file.
// Users can edit the file to show different names,
// but please do not change the order of IDs.
const IDListFile = "IDs.txt"

// MaskFile is the name of the mask file.
const MaskFile = "masks.bin"

// InfoFile contains some summary infomation.
const InfoFile = "info.txt"

// TreeDir is the name of director of trees.
const TreeDir = "trees"

// TreeFileExt is the file extension of the tree file.
// It would be ".gz", ".xz", ".zst" or ".bz2",
// but they are not recommended when saving a lot of trees,
// as it would assume a lot of RAM.
const TreeFileExt = ".bin"

// GenomeInfoFile stores some basic information of the indexed sequences.
const GenomeInfoFile = "genome_info.bin"

// TwoBitFile stores the 2bit-packed reference sequences
const TwoBitFile = "seqs.2bit"

// WriteToPath writes an index to a directory.
//
// Files:
//
//	Mask file, binary
//	ID list file, plain text
//	Trees directory, binary. Files numbers: 1-5000
//	...
//
// Note that, if you plan to save reference sequences, then call SetOutputPath()
// first. Next, you need to set 'overwrite' false when calling this method.
func (idx *Index) WriteToPath(outDir string, overwrite bool, threads int) error {
	pwd, _ := os.Getwd()
	if outDir != "./" && outDir != "." && pwd != filepath.Clean(outDir) {
		existed, err := pathutil.DirExists(outDir)
		if err != nil {
			return err
		}
		if existed {
			empty, err := pathutil.IsEmpty(outDir)
			if err != nil {
				return err
			}

			if !empty {
				if overwrite {
					err = os.RemoveAll(outDir)
					if err != nil {
						return err
					}
				} else {
					files, err := os.ReadDir(outDir)
					if err != nil {
						return err
					}

					// check if the directory is created by SetOutputPath,
					// which only two files: seqs.2bit and seqs.2bit.idx.
					flag := false
					for _, file := range files {
						if file.Name() == TwoBitFile {
							flag = true
						}
					}
					if !(flag && len(files) < 5) {
						return ErrDirNotEmpty
					}
				}
			} else {
				err = os.RemoveAll(outDir)
				if err != nil {
					return err
				}
			}
		}
		err = os.MkdirAll(outDir, 0777)
		if err != nil {
			return err
		}
	} else {
		return ErrPWDAsOutDir
	}

	// ID list file
	err := idx.writeIDlist(filepath.Join(outDir, IDListFile))
	if err != nil {
		return err
	}

	// Genome info file
	err = idx.writeGenomeInfo(filepath.Join(outDir, GenomeInfoFile))
	if err != nil {
		return err
	}

	// Mask file
	_, err = idx.lh.WriteToFile(filepath.Join(outDir, MaskFile))
	if err != nil {
		return err
	}

	if threads <= 0 {
		threads = runtime.NumCPU()
	}

	// Trees
	var wg sync.WaitGroup
	tokens := make(chan int, threads)

	var idStr, subDir, file string
	for i, t := range idx.Trees {
		idStr = fmt.Sprintf("%04d", i)
		subDir = idStr[len(idStr)-2:]
		file = filepath.Join(outDir, TreeDir, subDir, idStr+TreeFileExt)

		wg.Add(1)
		tokens <- 1

		go func(t *tree.Tree, file string) {
			defer func() {
				wg.Done()
				<-tokens
			}()
			t.WriteToFile(file)
			_, _err := t.WriteToFile(file)
			if _err != nil {
				err = _err
			}
		}(t, file)
	}
	wg.Wait()

	// Info file
	err = idx.writeInfo(filepath.Join(outDir, InfoFile))
	if err != nil {
		return err
	}

	if idx.saveTwoBit && idx.path != outDir {
		return fmt.Errorf("please set the same path to that used in SetOutputPath()")
	} else {
		idx.path = outDir
	}

	return err
}

// NewFromPath reads an index from a directory.
func NewFromPath(outDir string, threads int) (*Index, error) {
	// ------------- checking directory structure -----------

	ok, err := pathutil.DirExists(outDir)
	if err != nil {
		return nil, err
	}
	if !ok {
		return nil, fmt.Errorf("index path not found: %s", outDir)
	}

	// Mask file
	fileMask := filepath.Join(outDir, MaskFile)
	ok, err = pathutil.Exists(fileMask)
	if err != nil {
		return nil, err
	}
	if !ok {
		return nil, fmt.Errorf("mask file not found: %s", fileMask)
	}

	// ID list file
	fileIDList := filepath.Join(outDir, IDListFile)
	ok, err = pathutil.Exists(fileIDList)
	if err != nil {
		return nil, err
	}
	if !ok {
		return nil, fmt.Errorf("ID list file not found: %s", fileIDList)
	}

	// Genome info file
	fileGenomeiInfo := filepath.Join(outDir, GenomeInfoFile)
	ok, err = pathutil.Exists(fileIDList)
	if err != nil {
		return nil, err
	}
	if !ok {
		return nil, fmt.Errorf("genome info file not found: %s", fileIDList)
	}

	// Trees
	dirTrees := filepath.Join(outDir, TreeDir)
	ok, err = pathutil.DirExists(dirTrees)
	if err != nil {
		return nil, err
	}
	if !ok {
		return nil, fmt.Errorf("trees path not found: %s", dirTrees)
	}

	// ------------- parsing -----------

	idx := &Index{path: outDir}

	// Mask file
	idx.lh, err = lexichash.NewFromFile(fileMask)
	if err != nil {
		return nil, err
	}
	idx.k = uint8(idx.lh.K)

	// ID list file
	err = idx.readIDlist(fileIDList)
	if err != nil {
		return nil, err
	}

	// Genome info file
	err = idx.readGenomeInfo(fileGenomeiInfo)
	if err != nil {
		return nil, err
	}

	// Trees
	nMasks := len(idx.lh.Masks)
	idx.Trees = make([]*tree.Tree, nMasks)

	if threads <= 0 {
		threads = runtime.NumCPU()
	}

	treePaths := make([]string, 0, nMasks)
	fs.WalkDir(os.DirFS(dirTrees), ".", func(p string, d fs.DirEntry, err error) error {
		if filepath.Ext(p) == TreeFileExt {
			treePaths = append(treePaths, filepath.Join(dirTrees, p))
		}
		return nil
	})
	if len(treePaths) != nMasks {
		return nil, ErrTreeFileMissing
	}

	var wg sync.WaitGroup
	tokens := make(chan int, threads)
	for _, file := range treePaths {
		wg.Add(1)
		tokens <- 1
		go func(file string) {
			defer func() {
				wg.Done()
				<-tokens
			}()

			// idx of tree
			base := filepath.Base(file)
			i, _err := strconv.Atoi(base[0 : len(base)-len(TreeFileExt)])
			if _err != nil {
				err = _err
				return
			}

			t, _err := tree.NewFromFile(file)
			if _err != nil {
				err = _err
				return
			}

			idx.Trees[i] = t
		}(file)
	}
	wg.Wait()

	if err != nil {
		return nil, err
	}

	for i, t := range idx.Trees {
		if t == nil {
			return nil, fmt.Errorf("tree missing: %d", i)
		}
	}

	// searching
	idx.SetSearchingOptions(&DefaultSearchOptions)
	// idx.SetAlignOptions(&align.DefaultAlignOptions)
	idx.SetSeqCompareOptions(&DefaultSeqComparatorOptions)

	// 2bit file

	fileTwoBit := filepath.Join(outDir, TwoBitFile)
	ok, err = pathutil.Exists(fileTwoBit)
	if err != nil {
		return nil, err
	}
	if ok {
		idx.twobitReaders = make(chan *twobit.Reader, threads)
		for i := 0; i < threads; i++ {
			rdr, err := twobit.NewReader(fileTwoBit)
			if err != nil {
				return nil, err
			}
			idx.twobitReaders <- rdr
		}
		idx.saveTwoBit = true
	}

	return idx, nil
}

// Close closes all.
func (idx *Index) Close() error {
	if idx.saveTwoBit {
		close(idx.twobitReaders)
		var err error
		for rdr := range idx.twobitReaders {
			err = rdr.Close()
			if err != nil {
				return err
			}
		}
	}
	return nil
}

func (idx *Index) writeInfo(file string) error {
	outfh, err := xopen.Wopen(file)
	if err != nil {
		return err
	}
	defer outfh.Close()

	fmt.Fprintf(outfh, "index format version: v%d.%d\n", MainVersion, MinorVersion)
	fmt.Fprintf(outfh, "k-mer size: %d\n", idx.k)
	fmt.Fprintf(outfh, "number of masks: %d\n", len(idx.lh.Masks))
	fmt.Fprintf(outfh, "seed for generating masks: %d\n", idx.lh.Seed)
	fmt.Fprintf(outfh, "number of indexed sequences: %d\n", len(idx.IDs))

	return nil
}

// --------------------------------------------------------------

func (idx *Index) writeIDlist(file string) error {
	outfh, err := xopen.Wopen(file)
	if err != nil {
		return err
	}
	defer outfh.Close()

	for _, id := range idx.IDs {
		outfh.Write(id)
		outfh.WriteByte('\n')
	}

	return nil
}

func (idx *Index) readIDlist(file string) error {
	fh, err := xopen.Ropen(file)
	if err != nil {
		return err
	}
	defer fh.Close()

	idx.IDs, err = ReadIDlistFromFile(file)
	return err
}

// ReadIDlistFromFile read ID list from a file
func ReadIDlistFromFile(file string) ([][]byte, error) {
	fh, err := xopen.Ropen(file)
	if err != nil {
		return nil, err
	}
	defer fh.Close()

	ids := make([][]byte, 0, 1024)

	scanner := bufio.NewScanner(fh)
	for scanner.Scan() {
		ids = append(ids, []byte(scanner.Text()))
	}
	return ids, scanner.Err()
}

// --------------------------------------------------------------

var be = binary.BigEndian

func (idx *Index) writeGenomeInfo(file string) error {
	outfh, err := xopen.Wopen(file)
	if err != nil {
		return err
	}
	defer outfh.Close()

	err = binary.Write(outfh, be, uint64(len(idx.RefSeqInfos)))
	if err != nil {
		return err
	}
	if len(idx.RefSeqInfos) == 0 {
		return fmt.Errorf("no genome info to write")
	}

	var i int
	var size int
	for _, info := range idx.RefSeqInfos {
		err = binary.Write(outfh, be, [3]uint64{uint64(info.GenomeSize), uint64(info.Len), uint64(info.NumSeqs)})
		if err != nil {
			return err
		}

		data := make([]byte, 8*info.NumSeqs)
		i = 0
		for _, size = range info.SeqSizes {
			be.PutUint64(data[i:i+8], uint64(size))
			i += 8
		}
		_, err = outfh.Write(data)
		if err != nil {
			return err
		}
	}

	return nil
}

func (idx *Index) readGenomeInfo(file string) error {
	fh, err := xopen.Ropen(file)
	if err != nil {
		return err
	}
	defer fh.Close()

	buf := make([]byte, 8)

	// the number of records
	n, err := io.ReadFull(fh, buf)
	if err != nil {
		return err
	}
	if n < 8 {
		return fmt.Errorf("broken genome info file")
	}
	N := int(be.Uint64(buf))
	if N == 0 {
		return fmt.Errorf("no genome info to read")
	}

	idx.RefSeqInfos = make([]RefSeqInfo, N)

	var i, j int
	for i = 0; i < N; i++ {
		info := RefSeqInfo{}
		// genome size
		n, err := io.ReadFull(fh, buf)
		if err != nil {
			return err
		}
		if n < 8 {
			return ErrBrokenGenomeInfoFile
		}
		info.GenomeSize = int(be.Uint64(buf))

		// len of concatenated seqs
		n, err = io.ReadFull(fh, buf)
		if err != nil {
			return err
		}
		if n < 8 {
			return ErrBrokenGenomeInfoFile
		}
		info.Len = int(be.Uint64(buf))

		// number of sequences
		n, err = io.ReadFull(fh, buf)
		if err != nil {
			return err
		}
		if n < 8 {
			return ErrBrokenGenomeInfoFile
		}
		info.NumSeqs = int(be.Uint64(buf))

		info.SeqSizes = make([]int, info.NumSeqs)
		for j = 0; j < info.NumSeqs; j++ {
			n, err = io.ReadFull(fh, buf)
			if err != nil {
				return err
			}
			if n < 8 {
				return ErrBrokenGenomeInfoFile
			}
			info.SeqSizes[j] = int(be.Uint64(buf))
		}
		idx.RefSeqInfos[i] = info
	}
	return nil
}

// --------------------------------------------------------------
