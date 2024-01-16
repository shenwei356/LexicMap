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
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OFTestSerializationTestSerialization ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

package index

import (
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"sync"

	"github.com/shenwei356/LexicMap/lexicmap/index/twobit"
	"github.com/shenwei356/LexicMap/lexicmap/tree"
	"github.com/shenwei356/lexichash"
	"github.com/shenwei356/util/pathutil"
)

// Strands could be used to output strand for a reverse complement flag
var Strands = [2]byte{'+', '-'}

// ErrKConcurrentInsert occurs when calling Insert during calling BatchInsert.
var ErrKConcurrentInsert = errors.New("index: concurrent insertion")

// Index creates LexicHash index for mutitple reference sequences
// and supports searching with a query sequence.
type Index struct {
	lh          *lexichash.LexicHash
	k           uint8
	batchInsert bool

	// each record of the k-mer value is an uint64
	//  ref idx: 26 bits
	//  pos:     36 bits (0-based position)
	//  strand:   2 bits
	Trees []*tree.Tree

	IDs         [][]byte     // IDs of the reference genomes
	RefSeqInfos []RefSeqInfo // Reference sequence basic information
	i           uint32       // curent index, for inserting a new ref seq

	// ------------- optional -------------

	path string // path of the index directory

	KmerLocations [][]uint64 // mask locations of each reference genomes

	// for writing and reading 2bit-packed sequences
	saveTwoBit    bool
	twobitWriter  *twobit.Writer
	twobitReaders chan *twobit.Reader // reader pool

	// for searching
	chainingOptions *ChainingOptions
	searchOptions   *SearchOptions
	poolChainers    *sync.Pool

	// for sequence comparing
	seqCompareOption  *SeqComparatorOptions
	poolSeqComparator *sync.Pool

	// for sequence alignment
	// alignOptions *align.AlignOptions
	// poolAligner  *sync.Pool
}

// NewIndex ceates a new Index.
// nMasks better be >= 1024 and better be power of 4,
// i.e., 4, 16, 64, 256, 1024, 4096 ...
// p is the length of mask k-mer prefixes which need to be checked for low-complexity.
// p == 0 for no checking.
func NewIndex(k int, nMasks int, p int) (*Index, error) {
	return NewIndexWithSeed(k, nMasks, 1, p)
}

// NewIndexWithSeed ceates a new Index with given seed.
// nMasks better be >= 1024 and better be power of 4,
// i.e., 4, 16, 64, 256, 1024, 4096 ...
// p is the length of mask k-mer prefixes which need to be checked for low-complexity.
// p == 0 for no checking.
func NewIndexWithSeed(k int, nMasks int, seed int64, p int) (*Index, error) {
	lh, err := lexichash.NewWithSeed(k, nMasks, seed, p)
	if err != nil {
		return nil, err
	}

	// create a tree for each mask
	trees := make([]*tree.Tree, len(lh.Masks))
	for i := range trees {
		trees[i] = tree.New(uint8(k))
	}

	idx := &Index{
		lh:    lh,
		k:     uint8(k),
		Trees: trees,
		IDs:   make([][]byte, 0, 128),
		i:     0,

		RefSeqInfos: make([]RefSeqInfo, 0, 128),
	}

	idx.SetSearchingOptions(&DefaultSearchOptions)
	// lidx.SetAlignOptions(&align.DefaultAlignOptions)
	idx.SetSeqCompareOptions(&DefaultSeqComparatorOptions)

	return idx, nil
}

// SetOutputPath sets the output directory of index.
// If you want to save the reference sequences in 2bit-packed binary format,
// Please calls this method right after creating a new index with NewIndex
// or NewIndexWithSeed.
// And later please call WriteToPath() to
// save other data into the same output directory.
// Why is it designed to this? Because sequences, which might be tens of GB,
// need to be written to file during index building. While index data
// are optionally saved.
func (idx *Index) SetOutputPath(outDir string, overwrite bool) error {
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
					return ErrDirNotEmpty
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

	idx.path = outDir

	twobitFile := filepath.Join(outDir, TwoBitFile)
	var err error
	idx.twobitWriter, err = twobit.NewWriter(twobitFile)
	if err != nil {
		return fmt.Errorf("fail to write 2bit seq file: %s", twobitFile)
	}

	idx.saveTwoBit = true

	return nil
}

// K returns the K value
func (idx *Index) K() int {
	return int(idx.k)
}

// Threads is the maximum concurrency number for Insert() and ParallelizedSearch().
var Threads = runtime.NumCPU()

// GenomeInfo is a struct to store some basic information of a ref seq
type RefSeqInfo struct {
	GenomeSize int   // bases of all sequences
	Len        int   // length of contatenated sequences
	NumSeqs    int   // number of sequences
	SeqSizes   []int // sizes of sequences
}

// RefSeq represents a reference sequence to insert.
type RefSeq struct {
	ID  []byte
	Seq []byte

	RefSeqSize int   // genome size
	SeqSizes   []int // lengths of each sequences
}

// PoolRefSeq is the object pool of RefSeq.
var PoolRefSeq = &sync.Pool{New: func() interface{} {
	return &RefSeq{
		ID:         make([]byte, 0, 128),
		Seq:        make([]byte, 0, 10<<20),
		RefSeqSize: 0,
		SeqSizes:   make([]int, 0, 128),
	}
}}

// Reset resets RefSeq.
func (r *RefSeq) Reset() {
	r.ID = r.ID[:0]
	r.Seq = r.Seq[:0]
	r.RefSeqSize = 0
	r.SeqSizes = r.SeqSizes[:0]
}

// _MaskResult represents a mask result, it's only used in BatchInsert.
type _MaskResult struct {
	ID     []byte
	Kmers  *[]uint64
	Locses *[][]int

	RefSeqSize int   // genome size
	SeqSizes   []int // lengths of each sequences

	TwoBit *[]byte // back-packed sequence
}

var poolMaskResult = &sync.Pool{New: func() interface{} {
	return &_MaskResult{
		ID:       make([]byte, 0, 128),
		SeqSizes: make([]int, 0, 128),
	}
}}

func (r *_MaskResult) Reset() {
	r.ID = r.ID[:0]
	r.SeqSizes = r.SeqSizes[:0]
}

// BatchInsert inserts reference sequences in parallel.
// It returns:
//
//	chan RefSeq, for sending sequence.
//	sync.WaitGroup, for wait all masks being computed.
//	chan int, for waiting all the insertions to be done.
//
// Example:
//
//	input, done := BatchInsert()
//
//	refseq := index.PoolRefSeq.Get().(*index.RefSeq)
//	refseq.Reset()
//
//	// record is a fastx.Record//
//	refseq.ID = append(refseq.ID, record.ID...)
//	refseq.Seq = append(refseq.Seq, record.Seq.Seq...)
//	refseq.SeqSizes = append(refseq.SeqSizes, len(record.Seq.Seq))
//	refseq.RefSeqSize = len(record.Seq.Seq)
//
//	input <- refseq
//
//	close(input)
//	<- done
//
// Multiple sequences can also be concatenated with (K-1) N's for being a single sequence.
// In this case, k-mers around the (K-1) N's regions will be ignored.
func (idx *Index) BatchInsert() (chan *RefSeq, chan int) {
	if idx.batchInsert {
		panic(ErrKConcurrentInsert)
	}
	idx.batchInsert = true

	input := make(chan *RefSeq, Threads)
	doneAll := make(chan int)

	go func() {
		ch := make(chan *_MaskResult, Threads)
		doneInsert := make(chan int)
		saveTwoBit := idx.saveTwoBit

		// insert to tree
		go func() {
			var wg sync.WaitGroup
			tokens := make(chan int, Threads)
			trees := idx.Trees
			var nMasks int
			var n int
			var j, start, end int
			var refIdx uint32
			var k int = idx.lh.K
			var sumLen int

			for m := range ch {
				nMasks = len(*(m.Kmers))
				n = nMasks/Threads + 1
				refIdx = idx.i
				sumLen = m.RefSeqSize + (len(m.SeqSizes)-1)*(k-1)

				if saveTwoBit {
					wg.Add(1)
					go func() {
						// here we need to write the 2bit seq into file
						idx.twobitWriter.Write2Bit(*m.TwoBit, sumLen)
						twobit.RecycleTwoBit(m.TwoBit)
						wg.Done()
					}()
				}

				for j = 0; j <= Threads; j++ {
					start, end = j*n, (j+1)*n
					if end > nMasks {
						end = nMasks
					}

					wg.Add(1)
					tokens <- 1
					go func(start, end int) {
						var kmer uint64
						var loc int
						var refpos uint64
						for i := start; i < end; i++ {
							kmer = (*m.Kmers)[i]
							for _, loc = range (*m.Locses)[i] {
								//  ref idx: 26 bits
								//  pos:     36 bits
								//  strand:   2 bits
								// here, the location already contain the strand information from Mask().
								refpos = uint64(refIdx)<<38 | uint64(loc)
								trees[i].Insert(kmer, refpos)
							}
						}
						wg.Done()
						<-tokens
					}(start, end)
				}

				wg.Wait()

				idx.IDs = append(idx.IDs, []byte(string(m.ID)))
				idx.RefSeqInfos = append(idx.RefSeqInfos, RefSeqInfo{
					GenomeSize: m.RefSeqSize,
					Len:        sumLen,
					NumSeqs:    len(m.SeqSizes),
					SeqSizes:   m.SeqSizes,
				})
				idx.i++

				idx.lh.RecycleMaskResult(m.Kmers, m.Locses)
				poolMaskResult.Put(m)

			}
			if saveTwoBit {
				idx.twobitWriter.Close()
			}
			doneInsert <- 1
		}()

		// compute mask
		var wg sync.WaitGroup
		tokens := make(chan int, Threads)

		for ref := range input {
			tokens <- 1
			wg.Add(1)
			go func(ref *RefSeq) {
				// compute regions to skip
				var skipRegions [][2]int
				if len(ref.SeqSizes) > 1 {
					k := idx.K()
					skipRegions = make([][2]int, len(ref.SeqSizes)-1)
					var n int // len of concatenated seqs
					for i, s := range ref.SeqSizes {
						if i > 0 {
							skipRegions[i-1] = [2]int{n, n + k - 1}
							n += k - 1
						}
						n += s
					}
				}

				// capture k-mers
				_kmers, locses, err := idx.lh.Mask(ref.Seq, skipRegions)
				if err != nil {
					panic(err)
				}

				m := poolMaskResult.Get().(*_MaskResult)
				m.Reset()
				m.Kmers = _kmers
				m.Locses = locses
				m.ID = append(m.ID, ref.ID...)
				m.RefSeqSize = ref.RefSeqSize
				m.SeqSizes = append(m.SeqSizes, ref.SeqSizes...)

				// here we convert the sequence to 2bit format
				// and save into MaskResult.
				if saveTwoBit {
					m.TwoBit = twobit.Seq2TwoBit(ref.Seq)
				}

				ch <- m

				PoolRefSeq.Put(ref)

				wg.Done()
				<-tokens
			}(ref)
		}

		wg.Wait()
		close(ch)
		<-doneInsert

		doneAll <- 1
	}()

	// compute

	return input, doneAll
}

// --------------------------------------------------------------

var poolPathResult = &sync.Pool{New: func() interface{} {
	return &Path{}
}}

var poolPathResults = &sync.Pool{New: func() interface{} {
	paths := make([]*Path, 0, 16)
	return &paths
}}

// RecyclePathResult recycles the node list.
func (idx *Index) RecyclePathResult(paths *[]*Path) {
	for _, p := range *paths {
		idx.Trees[p.TreeIdx].RecyclePathResult(p.Nodes)
		poolPathResult.Put(p)
	}
	poolPathResults.Put(paths)
}

// Path represents the path of query in a tree.
type Path struct {
	TreeIdx int
	Nodes   *[]string
	Bases   uint8
}

// Paths returned the paths in all trees.
// Do not forget to call RecyclePathResult after using the results.
func (idx *Index) Paths(key uint64, k uint8, minPrefix uint8) *[]*Path {
	var bases uint8
	// paths := make([]Path, 0, 8)
	paths := poolPathResults.Get().(*[]*Path)
	*paths = (*paths)[:0]
	for i, tree := range idx.Trees {
		var nodes *[]string
		// nodes, bases = tree.Path(key, uint8(k), minPrefix)
		nodes, bases = tree.Path(key, minPrefix)
		if bases >= minPrefix {
			// path := Path{TreeIdx: i, Nodes: nodes, Bases: bases}
			path := poolPathResult.Get().(*Path)
			path.TreeIdx = i
			path.Nodes = nodes
			path.Bases = bases
			*paths = append(*paths, path)
		}
	}
	return paths
}
