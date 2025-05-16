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
	"cmp"
	"fmt"
	"io"
	"io/fs"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"slices"
	"strconv"
	"sync"
	"time"
	"unsafe"

	"github.com/dustin/go-humanize"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/kv"
	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
	"github.com/shenwei356/bio/taxdump"
	"github.com/shenwei356/kmers"
	"github.com/shenwei356/lexichash"
	"github.com/shenwei356/util/pathutil"
	"github.com/shenwei356/wfa"
	"github.com/vbauerster/mpb/v8"
	"github.com/vbauerster/mpb/v8/decor"
	"github.com/zeebo/wyhash"
)

// IndexSearchingOptions contains all options for searching
type IndexSearchingOptions struct {
	// general
	NumCPUs      int
	Verbose      bool // show log
	Log2File     bool // log file
	MaxOpenFiles int  // maximum opened files, used in merging indexes

	// seed searching
	// MaxSeedingConcurrency int

	InMemorySearch bool  // load the seed/kv data into memory
	MinPrefix      uint8 // minimum prefix length, e.g., 15
	// MaxMismatch     int   // maximum mismatch, e.g., 3
	MinSinglePrefix uint8 // minimum prefix length of the single seed, e.g., 20
	// MinMatchedBases uint8 // the total matched bases
	TopN int // keep the topN scores, e.g, 10

	// seeds chaining
	MaxGap      float64 // e.g., 5000
	MaxDistance float64 // e.g., 20k

	// alignment
	ExtendLength  int // the length of extra sequence on the flanking of seeds.
	ExtendLength2 int // the length of extra sequence on the flanking of pseudo-alignment region.

	// alignment filtering
	MinQueryAlignedFractionInAGenome float64 // minimum query aligned fraction in the target genome
	MaxEvalue                        float64

	// Output
	OutputSeq bool

	// debug
	Debug bool

	// filter results by taxid
	TaxdumpDir              string
	Genome2TaxIdFile        string
	TaxIds                  []uint32
	NegativeTaxIds          []uint32
	KeepGenomesWithoutTaxId bool
}

func CheckIndexSearchingOptions(opt *IndexSearchingOptions) error {
	if opt.NumCPUs < 1 {
		return fmt.Errorf("invalid number of CPUs: %d, should be >= 1", opt.NumCPUs)
	}
	if opt.MaxOpenFiles < 2 {
		return fmt.Errorf("invalid max open files: %d, should be >= 2", opt.MaxOpenFiles)
	}

	// ------------------------
	if opt.MinPrefix < 3 || opt.MinPrefix > 32 {
		return fmt.Errorf("invalid MinPrefix: %d, valid range: [3, 32]", opt.MinPrefix)
	}

	return nil
}

var DefaultIndexSearchingOptions = IndexSearchingOptions{
	NumCPUs:      runtime.NumCPU(),
	MaxOpenFiles: 512,

	// MaxSeedingConcurrency: 8,

	MinPrefix: 15,
	// MaxMismatch:     -1,
	MinSinglePrefix: 17,
	// MinMatchedBases: 20,
	TopN: 500,

	MaxGap:      5000,
	MaxDistance: 10000,

	ExtendLength:                     2000,
	ExtendLength2:                    50,
	MinQueryAlignedFractionInAGenome: 70,
	MaxEvalue:                        10,
}

// Index creates a LexicMap index from a path
// and supports searching with query sequences.
type Index struct {
	path string

	openFileTokens chan int // control the max open files

	// lexichash
	lh *lexichash.LexicHash
	k  int
	k8 uint8

	// k-mer-value searchers
	Searchers         []*kv.Searcher
	InMemorySearchers []*kv.InMemorySearcher
	searcherTokens    []chan int // make sure one seachers is only used by one query
	poolKmers         *sync.Pool // for suffix index
	poolLocses        *sync.Pool // for suffix index

	// general options, and some for seed searching
	opt *IndexSearchingOptions

	// for seed chaining
	chainingOptions *ChainingOptions
	poolChainers    *sync.Pool

	// for sequence comparing
	contigInterval    int // read from info file
	seqCompareOption  *SeqComparatorOptions
	poolSeqComparator *sync.Pool
	poolChainers2     *sync.Pool

	// genome data reader
	poolGenomeRdrs []chan *genome.Reader
	hasGenomeRdrs  bool

	// genome chunks
	hasGenomeChunks              bool       // file FileGenomeChunks exists and it's not empty
	genomeChunks                 [][]uint64 // original chunk information
	poolGenomeChunksIdx2List     *sync.Pool
	poolGenomeChunksPointer2List *sync.Pool

	// totalBases
	totalBases int64

	// filter results by taxid
	filterByTaxId         bool
	filterByPositiveTaxId bool
	filterByNegativeTaxId bool
	genomeIdx2TaxId       map[uint64]uint32
	Taxonomy              *taxdump.Taxonomy
	poolTaxIDfilter       *sync.Pool
}

// SetSeqCompareOptions sets the sequence comparing options
func (idx *Index) SetSeqCompareOptions(sco *SeqComparatorOptions) {
	idx.seqCompareOption = sco
	idx.poolChainers2 = &sync.Pool{New: func() interface{} {
		return NewChainer2(&sco.Chaining2Options)
	}}
	idx.poolSeqComparator = &sync.Pool{New: func() interface{} {
		return NewSeqComparator(sco, idx.poolChainers2)
	}}
}

// NewIndexSearcher creates a new searcher
func NewIndexSearcher(outDir string, opt *IndexSearchingOptions) (*Index, error) {
	ok, err := pathutil.DirExists(outDir)
	if err != nil {
		return nil, err
	}
	if !ok {
		return nil, fmt.Errorf("index path not found: %s", outDir)
	}

	idx := &Index{path: outDir, opt: opt}

	// -----------------------------------------------------
	// info file
	fileInfo := filepath.Join(outDir, FileInfo)
	info, err := readIndexInfo(fileInfo)
	if err != nil {
		return nil, fmt.Errorf("failed to read info file: %s", err)
	}
	if info.MainVersion != MainVersion {
		checkError(fmt.Errorf("index main versions do not match: %d (index) != %d (tool). please re-create the index", info.MainVersion, MainVersion))
	}
	if info.InputBases == 0 {
		// checkError(fmt.Errorf(`please run "lexicmap utils recount-bases -d %s"`, outDir))
		startTime := time.Now()
		if opt.Verbose {
			log.Info("  counting total bases for this index (run only once) ...")
		}
		totalBases, err := updateInputBases(info, outDir, opt.NumCPUs)
		checkError(err)
		if opt.Verbose {
			log.Infof("  done counting total bases (%s) in %s", humanize.Comma(totalBases), time.Since(startTime))
		}
		info.InputBases = totalBases
	}
	idx.totalBases = info.InputBases

	if idx.opt.MaxOpenFiles < info.Chunks+2 {
		return nil, fmt.Errorf("max open files (%d) should not be < chunks (%d) + 2",
			idx.opt.MaxOpenFiles, info.Chunks)
	}

	idx.contigInterval = info.ContigInterval

	// -----------------------------------------------------
	// taxid-related files

	var wgT sync.WaitGroup
	if len(idx.opt.TaxIds)+len(idx.opt.NegativeTaxIds) > 0 {
		idx.filterByTaxId = true
		idx.filterByPositiveTaxId = len(idx.opt.TaxIds) > 0
		idx.filterByNegativeTaxId = len(idx.opt.NegativeTaxIds) > 0

		wgT.Add(1)
		go func() {
			defer wgT.Done()

			idx.Taxonomy, err = taxdump.NewTaxonomyFromNCBI(filepath.Join(idx.opt.TaxdumpDir, "nodes.dmp"))
			if err != nil {
				checkError(fmt.Errorf("  failed to load taxonomy data: %s", idx.opt.TaxdumpDir))
			}
			idx.Taxonomy.CacheLCA()

			if opt.Verbose || opt.Log2File {
				log.Infof("  taxonomy data loaded from: %s", idx.opt.TaxdumpDir)
			}
		}()

		wgT.Add(1)
		go func() {
			defer wgT.Done()

			idx.genomeIdx2TaxId = make(map[uint64]uint32, info.Genomes)

			// genome2taxid
			genome2taxids, err := readKVsUint32(idx.opt.Genome2TaxIdFile, false)
			if err != nil {
				checkError(fmt.Errorf("  failed to read genome2taxid file: %s", idx.opt.Genome2TaxIdFile))
			}
			if opt.Verbose || opt.Log2File {
				log.Infof("  %d genome2taxid records loaded", len(genome2taxids))
			}

			// genomes.map.bin
			fh, err := os.Open(filepath.Join(idx.path, FileGenomeIndex))
			if err != nil {
				checkError(fmt.Errorf("  failed to read genome index mapping file: %s", err))
			}
			defer fh.Close()

			r := bufio.NewReader(fh)

			buf := make([]byte, 8)
			var n, lenID int
			var batchIDAndRefID uint64
			var ok bool
			var taxid uint32

			for {
				n, err = io.ReadFull(r, buf[:2])
				if err != nil {
					if err == io.EOF {
						break
					}
					checkError(fmt.Errorf("failed to read genome index mapping file: %s", err))
				}
				if n < 2 {
					checkError(fmt.Errorf("broken genome map file"))
				}
				lenID = int(be.Uint16(buf[:2]))
				genomeId := make([]byte, lenID)

				n, err = io.ReadFull(r, genomeId)
				if err != nil {
					checkError(fmt.Errorf("broken genome map file"))
				}
				if n < lenID {
					checkError(fmt.Errorf("broken genome map file"))
				}

				n, err = io.ReadFull(r, buf)
				if err != nil {
					checkError(fmt.Errorf("broken genome map file"))
				}
				if n < 8 {
					checkError(fmt.Errorf("broken genome map file"))
				}

				batchIDAndRefID = be.Uint64(buf)

				if taxid, ok = genome2taxids[string(genomeId)]; ok {
					idx.genomeIdx2TaxId[batchIDAndRefID] = taxid
				} else {
					if opt.Verbose || opt.Log2File {
						log.Warningf("  taxid of %s is not given in the genome2taxid file: %s", genomeId, idx.opt.Genome2TaxIdFile)
					}
				}
			}

			if opt.Verbose || opt.Log2File {
				log.Infof("  taxid information loaded")
			}
		}()

		idx.poolTaxIDfilter = &sync.Pool{New: func() interface{} {
			tmp := make(map[uint64]bool, max(1024, len(idx.genomeIdx2TaxId)/10))
			return &tmp
		}}
	}

	// -----------------------------------------------------
	// read masks
	fileMask := filepath.Join(outDir, FileMasks)
	if opt.Verbose || opt.Log2File {
		log.Infof("  reading masks...")
	}
	idx.lh, err = lexichash.NewFromFile(fileMask)
	if err != nil {
		return nil, err
	}

	// create a lookup table for faster masking
	lenPrefix := 1
	for 1<<(lenPrefix<<1) <= len(idx.lh.Masks) {
		lenPrefix++
	}
	lenPrefix--
	err = idx.lh.IndexMasks(lenPrefix)
	if err != nil {
		return nil, err
	}
	err = idx.lh.IndexMasksWithDistinctPrefixes(lenPrefix + 1)
	if err != nil {
		return nil, err
	}

	idx.k8 = uint8(idx.lh.K)
	idx.k = idx.lh.K

	if opt.MinPrefix > idx.k8 { // check again
		return nil, fmt.Errorf("MinPrefix (%d) should not be <= k (%d)", opt.MinPrefix, idx.k8)
	}

	// -----------------------------------------------------
	// read genome chunks data if existed
	fileGenomeChunks := filepath.Join(outDir, FileGenomeChunks)
	idx.genomeChunks, err = readGenomeChunksLists(fileGenomeChunks)
	if err != nil {
		return nil, err
	}
	if len(idx.genomeChunks) > 0 {
		idx.hasGenomeChunks = true

		// used to stored slice indexes of genome chunks of the same genome
		idx.poolGenomeChunksIdx2List = &sync.Pool{New: func() interface{} {
			data := make(map[uint64]*[]int)
			for _, idxs := range idx.genomeChunks {
				li := make([]int, 0, len(idxs))
				for _, idx := range idxs {
					data[idx] = &li
				}
			}
			return &data
		}}

		idx.poolGenomeChunksPointer2List = &sync.Pool{New: func() interface{} {
			data := make(map[uintptr]*[]int, 1024)
			return &data
		}}
	}
	// -----------------------------------------------------
	// read index of seeds

	inMemorySearch := idx.opt.InMemorySearch

	threads := opt.NumCPUs
	dirSeeds := filepath.Join(outDir, DirSeeds)
	fileSeeds := make([]string, 0, 64)
	fs.WalkDir(os.DirFS(dirSeeds), ".", func(p string, d fs.DirEntry, err error) error {
		if filepath.Ext(p) == ExtSeeds {
			fileSeeds = append(fileSeeds, filepath.Join(dirSeeds, p))
		}
		return nil
	})

	if len(fileSeeds) == 0 {
		return nil, fmt.Errorf("seeds file not found in: %s", dirSeeds)
	}
	if inMemorySearch {
		idx.InMemorySearchers = make([]*kv.InMemorySearcher, 0, len(fileSeeds))
	} else {
		idx.Searchers = make([]*kv.Searcher, 0, len(fileSeeds))
	}
	idx.searcherTokens = make([]chan int, len(fileSeeds))
	for i := range idx.searcherTokens {
		idx.searcherTokens[i] = make(chan int, 1)
	}
	idx.poolKmers = &sync.Pool{New: func() interface{} {
		tmp := make([]*[]uint64, len(idx.lh.Masks))
		for i := range tmp {
			tmp[i] = &[]uint64{}
		}
		return &tmp
	}}

	idx.poolLocses = &sync.Pool{New: func() interface{} {
		tmp := make([]*[]int, len(idx.lh.Masks))
		for i := range tmp {
			tmp[i] = &[]int{}
		}
		return &tmp
	}}

	// check options again
	if opt.MaxOpenFiles < len(fileSeeds) {
		return nil, fmt.Errorf("MaxOpenFiles (%d) should be > number of seeds files (%d), or even bigger", opt.MaxOpenFiles, len(fileSeeds))
	}
	idx.openFileTokens = make(chan int, opt.MaxOpenFiles) // tokens

	// read indexes

	if opt.Verbose || opt.Log2File {
		if inMemorySearch {
			log.Infof("  reading seeds (k-mer-value) data into memory...")
		} else {
			log.Infof("  reading indexes of seeds (k-mer-value) data...")
		}
	}
	done := make(chan int)
	var ch chan *kv.Searcher
	var chIM chan *kv.InMemorySearcher

	if inMemorySearch {
		chIM = make(chan *kv.InMemorySearcher, threads)
		go func() {
			for scr := range chIM {
				idx.InMemorySearchers = append(idx.InMemorySearchers, scr)
			}
			done <- 1
		}()
	} else {
		ch = make(chan *kv.Searcher, threads)
		go func() {
			for scr := range ch {
				idx.Searchers = append(idx.Searchers, scr)

				idx.openFileTokens <- 1 // increase the number of open files
			}
			done <- 1
		}()
	}
	var wg sync.WaitGroup
	tokens := make(chan int, threads)
	for _, file := range fileSeeds {
		wg.Add(1)
		tokens <- 1
		go func(file string) {
			if inMemorySearch { // read all the k-mer-value data into memory
				scr, err := kv.NewInMemomrySearcher(file)
				if err != nil {
					checkError(fmt.Errorf("failed to create a in-memory searcher from file: %s: %s", file, err))
				}

				chIM <- scr
			} else { // just read the index data
				scr, err := kv.NewSearcher(file)
				if err != nil {
					checkError(fmt.Errorf("failed to create a searcher from file: %s: %s", file, err))
				}

				ch <- scr
			}

			wg.Done()
			<-tokens
		}(file)
	}
	wg.Wait()
	if inMemorySearch {
		close(chIM)
	} else {
		close(ch)
	}
	<-done

	// we can create genome reader pools
	n := (idx.opt.MaxOpenFiles - len(fileSeeds) - 1) / info.GenomeBatches // 1 is for the output file
	if n < 2 {
	} else {
		if n > opt.NumCPUs {
			n = opt.NumCPUs
		}
		if opt.Verbose || opt.Log2File {
			log.Infof("  creating reader pools for %d genome batches, each with %d readers...", info.GenomeBatches, n)
		}
		idx.poolGenomeRdrs = make([]chan *genome.Reader, info.GenomeBatches)
		for i := 0; i < info.GenomeBatches; i++ {
			idx.poolGenomeRdrs[i] = make(chan *genome.Reader, n)
		}

		// parallelize it
		var wg sync.WaitGroup
		tokens := make(chan int, opt.NumCPUs)
		for i := 0; i < info.GenomeBatches; i++ {
			for j := 0; j < n; j++ {
				tokens <- 1
				wg.Add(1)
				go func(i int) {
					fileGenomes := filepath.Join(outDir, DirGenomes, batchDir(i), FileGenomes)
					rdr, err := genome.NewReader(fileGenomes)
					if err != nil {
						checkError(fmt.Errorf("failed to create genome reader: %s", err))
					}
					idx.poolGenomeRdrs[i] <- rdr

					idx.openFileTokens <- 1 // genome file

					wg.Done()
					<-tokens
				}(i)
			}
		}
		wg.Wait()

		idx.hasGenomeRdrs = true
	}

	// other resources
	co := &ChainingOptions{
		MaxGap:   opt.MaxGap,
		MinLen:   opt.MinSinglePrefix,
		MinScore: seedWeight(float64(opt.MinSinglePrefix)),
		// MinScore:    seedWeight(float64(opt.MinMatchedBases)),
		MaxDistance: opt.MaxDistance,
	}
	idx.chainingOptions = co
	idx.poolChainers = &sync.Pool{New: func() interface{} {
		return NewChainer(co)
	}}

	if idx.filterByTaxId {
		wgT.Wait()
	}

	return idx, nil
}

// Close closes the searcher.
func (idx *Index) Close() error {
	var _err error

	// seed data
	if idx.opt.InMemorySearch {
		for _, scr := range idx.InMemorySearchers {
			err := scr.Close()
			if err != nil {
				_err = err
			}
		}
	} else {
		for _, scr := range idx.Searchers {
			err := scr.Close()
			if err != nil {
				_err = err
			}
		}
	}

	// genome reader
	if idx.hasGenomeRdrs {
		var wg sync.WaitGroup
		for _, pool := range idx.poolGenomeRdrs {
			wg.Add(1)
			go func(pool chan *genome.Reader) {
				close(pool)
				for rdr := range pool {
					err := rdr.Close()
					if err != nil {
						_err = err
					}
				}
				wg.Done()
			}(pool)
		}
		wg.Wait()
	}
	return _err
}

// --------------------------------------------------------------------------
// structs for seeding results

// SubstrPair represents a pair of found substrings/seeds, it's also called an anchor.
type SubstrPair struct {
	QBegin int32 // start position of the substring (0-based) in query
	TBegin int32 // start position of the substring (0-based) in reference

	Len uint8 // prefix length
	// Mismatch uint8 // number of mismatches

	TRC bool // is the substring from the reference seq on the negative strand.
	QRC bool // is the substring from the query seq on the negative strand.

	// QCode uint64 // k-mer, for computing matched k-mers
	// TCode uint64
}

// func (s SubstrPair) Matches(k uint8) uint8 {
// 	return s.Len + util.MustSharingPrefixKmersSuffixMatches(s.QCode, s.TCode, k, s.Len)
// }

func (s SubstrPair) String() string {
	s1 := "+"
	s2 := "+"
	if s.QRC {
		s1 = "-"
	}
	if s.TRC {
		s2 = "-"
	}
	// return fmt.Sprintf("%3d-%3d (%s) vs %3d-%3d (%s), len:%2d, mismatches:%d",
	// 	s.QBegin+1, s.QBegin+int32(s.Len), s1, s.TBegin+1, s.TBegin+int32(s.Len), s2, s.Len, s.Mismatch)
	return fmt.Sprintf("%3d-%3d (%s) vs %3d-%3d (%s), len:%2d",
		s.QBegin+1, s.QBegin+int32(s.Len), s1, s.TBegin+1, s.TBegin+int32(s.Len), s2, s.Len)
}

var poolSub = &sync.Pool{New: func() interface{} {
	return &SubstrPair{}
}}

var poolSubs = &sync.Pool{New: func() interface{} {
	tmp := make([]*SubstrPair, 0, 64)
	return &tmp
}}

// RecycleSubstrPairs recycles a list of SubstrPairs
func RecycleSubstrPairs(poolSub *sync.Pool, poolSubs *sync.Pool, subs *[]*SubstrPair) {
	for _, sub := range *subs {
		poolSub.Put(sub)
	}
	*subs = (*subs)[:0]
	poolSubs.Put(subs)
}

var poolSubsLong = &sync.Pool{New: func() interface{} {
	tmp := make([]*SubstrPair, 0, 102400)
	return &tmp
}}

// ClearSubstrPairs removes nested/embedded and same anchors. k is the largest k-mer size.
func ClearSubstrPairs(poolSub *sync.Pool, subs *[]*SubstrPair, k int) {
	if len(*subs) < 2 {
		return
	}

	// sort substrings/seeds in ascending order based on the starting position
	// and in descending order based on the ending position.
	// sort.Slice(*subs, func(i, j int) bool {
	// 	a := (*subs)[i]
	// 	b := (*subs)[j]
	// 	if a.QBegin == b.QBegin {
	// 		// return a.QBegin+int32(a.Len) >= b.QBegin+int32(b.Len)
	// 		if a.QBegin+int32(a.Len) == b.QBegin+int32(b.Len) {
	// 			return a.TBegin <= b.TBegin
	// 		}
	// 		return a.QBegin+int32(a.Len) > b.QBegin+int32(b.Len)
	// 	}
	// 	return a.QBegin < b.QBegin
	// })
	slices.SortFunc(*subs, func(a, b *SubstrPair) int {
		if a.QBegin == b.QBegin {
			if a.QBegin+int32(a.Len) == b.QBegin+int32(b.Len) {
				return int(a.TBegin - b.TBegin)
			}
			return int(b.QBegin) + int(b.Len) - (int(a.QBegin) + int(a.Len))
		}
		return int(a.QBegin - b.QBegin)
	})

	var p *SubstrPair
	var upbound, vQEnd, vTEnd int32
	var j int
	markers := poolBoolList.Get().(*[]bool)
	*markers = (*markers)[:0]
	for range *subs {
		*markers = append(*markers, false)
	}
	for i, v := range (*subs)[1:] {
		vQEnd = int32(v.QBegin) + int32(v.Len)
		upbound = int32(vQEnd) - int32(k)
		vTEnd = int32(v.TBegin) + int32(v.Len)
		j = i
		for j >= 0 { // have to check previous N seeds
			p = (*subs)[j]
			if p.QBegin < upbound { // no need to check
				break
			}

			// same or nested region
			if vQEnd <= p.QBegin+int32(p.Len) &&
				v.TBegin >= p.TBegin && vTEnd <= p.TBegin+int32(p.Len) {
				poolSub.Put(v)         // do not forget to recycle the object
				(*markers)[i+1] = true // because of: range (*subs)[1:]
				break
			}

			j--
		}
	}

	j = 0
	for i, embedded := range *markers {
		if !embedded {
			(*subs)[j] = (*subs)[i]
			j++
		}
	}
	if j > 0 {
		*subs = (*subs)[:j]
	}

	poolBoolList.Put(markers)
}

var poolBoolList = &sync.Pool{New: func() interface{} {
	m := make([]bool, 0, 1024)
	return &m
}}

// --------------------------------------------------------------------------
// structs for searching result

var poolSimilarityDetail = &sync.Pool{New: func() interface{} {
	return &SimilarityDetail{
		SeqID: make([]byte, 0, 128),
	}
}}

var poolSimilarityDetails = &sync.Pool{New: func() interface{} {
	tmp := make([]*SimilarityDetail, 0, 8)
	return &tmp
}}

var poolSearchResult = &sync.Pool{New: func() interface{} {
	return &SearchResult{
		ID: make([]byte, 0, 128),
	}
}}

var poolSearchResults = &sync.Pool{New: func() interface{} {
	tmp := make([]*SearchResult, 0, 16)
	return &tmp
}}

// SearchResult stores a search result in a genome for the given query sequence.
type SearchResult struct {
	BatchGenomeIndex uint64 // just for finding genome chunks of the same genome

	GenomeBatch int
	GenomeIndex int
	ID          []byte
	GenomeSize  int

	Subs *[]*SubstrPair // matched substring pairs (query,target)

	Score  float64 //  score for soring
	Chains *[]*[]int

	// more about the alignment detail
	SimilarityDetails *[]*SimilarityDetail // sequence comparing
	AlignedFraction   float64              // query coverage per genome
}

func (sr *SearchResult) SortBySeqID() {
	if len(*sr.SimilarityDetails) <= 1 {
		return
	}

	sds0 := sr.SimilarityDetails

	sds := poolSimilarityDetails.Get().(*[]*SimilarityDetail) // HSPs in a reference

	// first one
	*sds = append(*sds, (*sds0)[0])
	seqid := (*sds0)[0].SeqID
	(*sds0)[0] = nil

	i := 1
	var nextSeqId []byte
	var sds1 []*SimilarityDetail
	var first bool
	var j, nextI int
	var sd *SimilarityDetail
	n := len(*sds0)
	for {
		sds1 = (*sds0)[i:]
		first = true
		for j, sd = range sds1 {
			if sd == nil {
				continue
			}

			if bytes.Equal(sd.SeqID, seqid) { // has the same seqid
				*sds = append(*sds, sd)
				(*sds0)[i+j] = nil
			} else if first { // record another seqid
				first = false
				nextSeqId = sd.SeqID
				nextI = i + j
			}
		}
		if first { // no other seqids
			break
		}

		*sds = append(*sds, (*sds0)[nextI])
		(*sds0)[nextI] = nil
		seqid = nextSeqId
		i = nextI + 1
		if i >= n { // it's alreadly the last one
			break
		}
	}

	*sr.SimilarityDetails = (*sr.SimilarityDetails)[:0]
	poolSimilarityDetails.Put(sr.SimilarityDetails)
	sr.SimilarityDetails = sds
}

// SimilarityDetail is the similarity detail of one reference sequence
type SimilarityDetail struct {
	// QBegin int
	// QEnd   int
	// TBegin int
	// TEnd   int
	RC bool

	SimilarityScore float64
	Similarity      *SeqComparatorResult
	// Chain           *[]int
	NSeeds int

	// sequence details
	SeqLen int
	SeqID  []byte // seqid of the region
}

func (r *SearchResult) Reset() {
	r.GenomeBatch = -1
	r.GenomeIndex = -1
	r.ID = r.ID[:0]
	r.GenomeSize = 0
	r.Subs = nil
	r.Score = 0
	r.Chains = nil
	r.SimilarityDetails = nil
	r.AlignedFraction = 0
}

// RecycleSearchResults recycles a search result object
func (idx *Index) RecycleSearchResult(r *SearchResult) {
	if r.Subs != nil {
		RecycleSubstrPairs(poolSub, poolSubs, r.Subs)
		r.Subs = nil
	}

	if r.Chains != nil {
		RecycleChainingResult(r.Chains)
		r.Chains = nil
	}

	// yes, it might be nil for some failed in chaining
	if r.SimilarityDetails != nil {
		idx.RecycleSimilarityDetails(r.SimilarityDetails)

		r.SimilarityDetails = nil
	}

	poolSearchResult.Put(r)
}

// RecycleSimilarityDetails recycles a list of SimilarityDetails
func (idx *Index) RecycleSimilarityDetails(sds *[]*SimilarityDetail) {
	for _, sd := range *sds {
		RecycleSeqComparatorResult(sd.Similarity)
		poolSimilarityDetail.Put(sd)
	}
	*sds = (*sds)[:0]
	poolSimilarityDetails.Put(sds)
}

// RecycleSearchResults recycles search results objects
func (idx *Index) RecycleSearchResults(sr *[]*SearchResult) {
	if sr == nil {
		return
	}

	for _, r := range *sr {
		idx.RecycleSearchResult(r)
	}
	poolSearchResults.Put(sr)
}

var poolSearchResultsMap = &sync.Pool{New: func() interface{} {
	m := make(map[int]*SearchResult, 1024)
	return &m
}}

// --------------------------------------------------------------------------
// searching

// Search queries the index with a sequence.
// After using the result, do not forget to call RecycleSearchResult().
func (idx *Index) Search(query *Query) (*[]*SearchResult, error) {
	var startTime time.Time
	debug := idx.opt.Debug

	if debug {
		log.Debugf("%s (%d bp): start to search", query.seqID, len(query.seq))
		startTime = time.Now()
	}

	s := query.seq

	// ----------------------------------------------------------------
	// 1) mask the query sequence

	// _kmers, _locses, err := idx.lh.Mask(s, nil)
	// _kmers, _locses, err := idx.lh.MaskKnownPrefixes(s, nil)
	_kmers, _locses, err := idx.lh.MaskKnownDistinctPrefixes(s, nil, true)
	if err != nil {
		return nil, err
	}
	defer idx.lh.RecycleMaskResult(_kmers, _locses)

	// remove low-complexity k-mers
	k8 := idx.k8
	ccc := util.Ns(0b01, k8)
	ggg := util.Ns(0b10, k8)
	ttt := (uint64(1) << (k8 << 1)) - 1

	for i, kmer := range *_kmers {
		if kmer == 0 || kmer == ccc || kmer == ggg || kmer == ttt ||
			util.IsLowComplexityDust(kmer, k8) {
			// if kmer != 0 {
			// 	fmt.Printf("low-complexity k-mer #%d: %s\n", i, lexichash.MustDecode(kmer, k8))
			// }

			(*_kmers)[i] = 0
			// (*_locses)[i] = (*_locses)[i][:0]
			continue
		}
	}

	// ----------------------------------------------------------------
	// 2) matching the captured k-mers in databases

	// a map for collecting matches for each reference: IdIdx -> result
	m := poolSearchResultsMap.Get().(*map[int]*SearchResult)
	clear(*m) // requires go >= v1.21

	inMemorySearch := idx.opt.InMemorySearch

	var searchers []*kv.Searcher
	var searchersIM []*kv.InMemorySearcher
	var nSearchers int

	if inMemorySearch {
		searchersIM = idx.InMemorySearchers
		nSearchers = len(searchersIM)
	} else {
		searchers = idx.Searchers
		nSearchers = len(searchers)
	}

	minPrefix := idx.opt.MinPrefix
	// maxMismatch := idx.opt.MaxMismatch

	ch := make(chan *[]*kv.SearchResult, nSearchers)
	done := make(chan int) // later, we will reuse this
	var wg sync.WaitGroup
	var beginM, endM int // range of mask of a chunk

	// -----------------------
	// reverse k-mers
	_kmersR := idx.poolKmers.Get().(*[]*[]uint64)
	_locsesR := idx.poolLocses.Get().(*[]*[]int)

	chR := make(chan [3]uint64, nSearchers)
	doneR := make(chan int)
	go func() {
		var v *[]uint64
		for _, v = range *_kmersR {
			*v = (*v)[:0]
		}

		var vl *[]int
		for _, vl = range *_locsesR {
			*vl = (*vl)[:0]
		}

		var _kmer, _v, newMask, oldMask uint64
		var existed bool

		for i2k := range chR {
			// multiple oldMask might points to the same newMask
			newMask, _kmer, oldMask = i2k[0], i2k[1], i2k[2]
			v = (*_kmersR)[newMask]
			vl = (*_locsesR)[newMask]

			existed = false
			for _, _v = range *v {
				if _kmer == _v {
					existed = true
					break
				}
			}
			if !existed {
				*v = append(*v, _kmer)
				*vl = append(*vl, int(oldMask))
			}
		}

		doneR <- 1
	}()

	for iS := 0; iS < nSearchers; iS++ {
		if inMemorySearch {
			beginM = searchersIM[iS].ChunkIndex
			endM = searchersIM[iS].ChunkIndex + searchersIM[iS].ChunkSize
		} else {
			beginM = searchers[iS].ChunkIndex
			endM = searchers[iS].ChunkIndex + searchers[iS].ChunkSize
		}

		wg.Add(1)
		go func(iS, beginM, endM int) {
			var iMasks *[]int
			var j int
			var minj int
			var mask, h, minh uint64
			k := idx.lh.K
			lh := idx.lh

			for i, kmer := range (*_kmers)[beginM:endM] {
				if kmer == 0 {
					continue
				}
				// fmt.Printf("mask: %d, %s\n", i+1, kmers.MustDecode(kmer, k))
				kmer = kmers.MustReverse(kmer, k) // reverse the k-mer
				iMasks = lh.MaskKmer(kmer)
				minh = math.MaxUint64
				for _, j = range *iMasks {
					mask = lh.Masks[j]
					h = mask ^ kmer
					if h < minh {
						minj, minh = j, h
					}
				}
				// fmt.Printf("mask: %d, kmer: %s, locs: %d,  new mask: %d\n",
				// 	beginM+i, kmers.MustDecode(kmer, k), (*_locses)[beginM+i], minj)

				// multiple beginM + i will points to the same minj
				chR <- [3]uint64{uint64(minj), kmer, uint64(beginM + i)}
				lh.RecycleMaskKmerResult(iMasks)
			}

			wg.Done()
		}(iS, beginM, endM)
	}
	wg.Wait()
	close(chR)
	<-doneR
	// -----------------------

	// 2.2) collect search results, they will be kept in RAM.
	// For quries with a lot of hits, the memory would be high.
	// And it's inevitable currently, but if we do want to decrease the memory usage,
	// we can write these matches in temporal files.
	go func() {
		var refpos uint64

		// query substring
		var posQ int
		var beginQ int
		var rcQ bool

		// var qCode, tCode uint64
		var kPrefix int
		var refBatchAndIdx, posT, beginT int
		// var mismatch uint8
		var rcT bool
		var rvT bool

		K := idx.k
		// K8 := idx.k8
		var locs []int
		var sr *kv.SearchResult
		var ok bool

		var refBatchAndIdxUint64 uint64
		var filter *map[uint64]bool
		filterByTaxId := idx.filterByTaxId
		filterByPositiveTaxId := idx.filterByPositiveTaxId
		filterByNegativeTaxId := idx.filterByNegativeTaxId
		if filterByTaxId {
			filter = idx.poolTaxIDfilter.Get().(*map[uint64]bool)
		}
		var keepGenome, matchOne bool
		taxon := idx.Taxonomy
		genomeIdx2TaxId := idx.genomeIdx2TaxId
		keepGenomesWithoutTaxId := idx.opt.KeepGenomesWithoutTaxId
		var _taxid, taxid uint32
		taxids := idx.opt.TaxIds
		negativeTaxids := idx.opt.NegativeTaxIds

		for srs := range ch {
			// different k-mers in subjects,
			// most of cases, there are more than one
			for _, sr = range *srs {
				// matched length
				kPrefix = int(sr.Len)
				// mismatch = sr.Mismatch
				// qCode = (*_kmers)[sr.IQuery]

				// locations in the query
				// multiple locations for each QUERY k-mer,
				// but most of cases, there's only one.
				if !sr.IsSuffix {
					locs = (*_locses)[sr.IQuery] // the mask is unknown
				} else {
					locs = (*_locses)[(*(*_locsesR)[sr.IQuery])[sr.IQuery2]] // the mask is unknown
					// fmt.Println(sr.IQuery, (*_locsesR)[sr.IQuery], locs)
				}
				for _, posQ = range locs {
					// query k-mers do not have the reverse flag !!!!
					rcQ = posQ&BITS_STRAND > 0 // if on the reverse complement sequence
					posQ >>= BITS_STRAND

					// matched
					// code = util.KmerPrefix(sr.Kmer, K8, sr.LenPrefix)
					// tCode = sr.Kmer

					// multiple locations for each MATCHED k-mer
					// but most of cases, there's only one.
					for _, refpos = range sr.Values {
						// refBatchAndIdx = int(refpos >> 30) // batch+refIdx
						refBatchAndIdxUint64 = refpos >> BITS_NONE_IDX // batch+refIdx
						refBatchAndIdx = int(refBatchAndIdxUint64)

						// filter by taxid
						if filterByTaxId {
							if keepGenome, ok = (*filter)[refBatchAndIdxUint64]; ok {
								if !keepGenome {
									continue
								}
							} else {

								if taxid, ok = genomeIdx2TaxId[refBatchAndIdxUint64]; ok {
									// black list
									if filterByNegativeTaxId {
										matchOne = false
										for _, _taxid = range negativeTaxids {
											if taxon.LCA(taxid, _taxid) == _taxid {
												matchOne = true
												break
											}
										}
										if matchOne {
											(*filter)[refBatchAndIdxUint64] = false
											continue
										} else if !filterByPositiveTaxId {
											(*filter)[refBatchAndIdxUint64] = true
											continue
										}
									}

									// white list
									if filterByPositiveTaxId {
										matchOne = false
										for _, _taxid = range taxids {
											if taxon.LCA(taxid, _taxid) == _taxid {
												matchOne = true
												break
											}
										}
										if matchOne {
											(*filter)[refBatchAndIdxUint64] = true
										} else {
											(*filter)[refBatchAndIdxUint64] = false
											continue
										}
									}
								} else if !keepGenomesWithoutTaxId {
									(*filter)[refBatchAndIdxUint64] = false
									continue
								}
								(*filter)[refBatchAndIdxUint64] = true
							}
						}

						// posT = int(refpos << 34 >> 35)
						posT = int(refpos << BITS_IDX >> BITS_IDX_FLAGS)
						rvT = refpos&BITS_REVERSE > 0
						rcT = refpos>>BITS_REVERSE&BITS_REVERSE > 0

						if !rvT {
							// query location
							if rcQ { // on the negative strand
								beginQ = posQ + K - kPrefix
							} else {
								beginQ = posQ
							}

							// subject location
							if rcT {
								beginT = posT + K - kPrefix
							} else {
								beginT = posT
							}
						} else {
							// query location
							if rcQ { // on the negative strand
								beginQ = posQ
							} else {
								beginQ = posQ + K - kPrefix
							}

							// subject location
							if rcT {
								beginT = posT
							} else {
								beginT = posT + K - kPrefix
							}
						}

						_sub2 := poolSub.Get().(*SubstrPair)
						_sub2.QBegin = int32(beginQ)
						_sub2.TBegin = int32(beginT)
						// _sub2.QCode = qCode
						// _sub2.TCode = tCode
						_sub2.Len = uint8(kPrefix)
						// _sub2.Mismatch = mismatch
						_sub2.QRC = rcQ
						_sub2.TRC = rcT

						var r *SearchResult
						if r, ok = (*m)[refBatchAndIdx]; !ok {
							subs := poolSubs.Get().(*[]*SubstrPair)

							r = poolSearchResult.Get().(*SearchResult)
							r.BatchGenomeIndex = uint64(refBatchAndIdx)
							r.GenomeBatch = refBatchAndIdx >> BITS_GENOME_IDX
							r.GenomeIndex = refBatchAndIdx & MASK_GENOME_IDX
							r.ID = r.ID[:0] // extract it from genome file later
							r.GenomeSize = 0
							r.Subs = subs
							r.Score = 0
							r.Chains = nil            // important
							r.SimilarityDetails = nil // important
							r.AlignedFraction = 0

							(*m)[refBatchAndIdx] = r
						}

						*r.Subs = append(*r.Subs, _sub2)
					}
				}
			}

			kv.RecycleSearchResults(srs)
		}

		if filterByTaxId {
			clear(*filter)
			idx.poolTaxIDfilter.Put(filter)
		}
		done <- 1
	}()

	// 2.1) search with multiple searchers
	// tokensS := make(chan int, idx.opt.MaxSeedingConcurrency)
	for iS := 0; iS < nSearchers; iS++ {
		if inMemorySearch {
			beginM = searchersIM[iS].ChunkIndex
			endM = searchersIM[iS].ChunkIndex + searchersIM[iS].ChunkSize
		} else {
			beginM = searchers[iS].ChunkIndex
			endM = searchers[iS].ChunkIndex + searchers[iS].ChunkSize
		}

		wg.Add(1)
		// tokensS <- 1
		go func(iS, beginM, endM int) {
			idx.searcherTokens[iS] <- 1 // get the access to the searcher
			var srs *[]*kv.SearchResult
			var srs2 *[]*kv.SearchResult
			var err error
			if inMemorySearch {
				// prefix search
				// srs, err = searchersIM[iS].Search((*_kmers)[beginM:endM], minPrefix, maxMismatch)
				srs, err = searchersIM[iS].Search((*_kmers)[beginM:endM], minPrefix, true, false)
				if err != nil {
					checkError(err)
				}

				// suffix search
				srs2, err = searchersIM[iS].Search2((*_kmersR)[beginM:endM], minPrefix, true, true)
				if err != nil {
					checkError(err)
				}
				if len(*srs2) > 0 {
					*srs = append(*srs, (*srs2)...)
					*srs2 = (*srs2)[:0]
				}
				kv.RecycleSearchResults(srs2)
			} else {
				// prefix search
				// srs, err = searchers[iS].Search((*_kmers)[beginM:endM], minPrefix, maxMismatch)
				srs, err = searchers[iS].Search((*_kmers)[beginM:endM], minPrefix, true, false)
				if err != nil {
					checkError(err)
				}

				// suffix search
				srs2, err = searchers[iS].Search2((*_kmersR)[beginM:endM], minPrefix, true, true)
				if err != nil {
					checkError(err)
				}
				if len(*srs2) > 0 {
					*srs = append(*srs, (*srs2)...)
					*srs2 = (*srs2)[:0]
				}
				kv.RecycleSearchResults(srs2)
			}
			if err != nil {
				checkError(err)
			}

			if len(*srs) == 0 { // no matcheds
				kv.RecycleSearchResults(srs)
			} else {
				ch <- srs // send result
			}

			<-idx.searcherTokens[iS] // return the access
			wg.Done()
			// <-tokensS
		}(iS, beginM, endM)
	}
	wg.Wait()
	close(ch)
	<-done

	idx.poolKmers.Put(_kmersR)
	idx.poolLocses.Put(_locsesR)

	if debug {
		if idx.filterByTaxId {
			log.Debugf("%s (%d bp): finished seed-matching with filtering by TaxId (%d genome hits) in %s", query.seqID, len(query.seq), len(*m), time.Since(startTime))
		} else {
			log.Debugf("%s (%d bp): finished seed-matching (%d genome hits) in %s", query.seqID, len(query.seq), len(*m), time.Since(startTime))
		}

		startTime = time.Now()
	}

	if len(*m) == 0 { // no results
		poolSearchResultsMap.Put(m)
		return nil, nil
	}

	// ----------------------------------------------------------------
	// 3) chaining matches for all reference genomes, and alignment

	// minMatchedBases := idx.opt.MinMatchedBases

	// 3.1) preprocess substring matches and chaining for each reference genome
	rs := poolSearchResults.Get().(*[]*SearchResult)
	*rs = (*rs)[:0]

	K := idx.k
	// k8 := idx.k8
	// var matches uint8
	// checkMismatch := maxMismatch >= 0 && maxMismatch < K-int(idx.opt.MinPrefix)

	ch1 := make(chan *SearchResult, idx.opt.NumCPUs)
	tokens := make(chan int, idx.opt.NumCPUs)

	// collect chaining result
	go func() {
		for r := range ch1 {
			*rs = append(*rs, r)
		}

		done <- 1
	}()

	minScore := idx.chainingOptions.MinScore
	for _, r := range *m {
		tokens <- 1
		wg.Add(1)

		go func(r *SearchResult) {
			ClearSubstrPairs(poolSub, r.Subs, K) // remove duplicates and nested anchors

			// -----------------------------------------------------
			// chaining

			chainer := idx.poolChainers.Get().(*Chainer)
			r.Chains, r.Score = chainer.Chain(r.Subs)

			defer func() {
				idx.poolChainers.Put(chainer)
				<-tokens
				wg.Done()
			}()

			if r.Score < minScore {
				idx.RecycleSearchResult(r) // do not forget to recycle unused objects
				return
			}

			ch1 <- r
		}(r)
	}

	wg.Wait()
	close(ch1)
	<-done

	poolSearchResultsMap.Put(m)

	// 3.2) only keep the top N targets
	topN := idx.opt.TopN
	if topN > 0 && len(*rs) > topN {
		// sort subjects in descending order based on the score
		// just use the standard library for a few seed pairs.
		// sort.Slice(*rs, func(i, j int) bool {
		// 	return (*rs)[i].Score > (*rs)[j].Score
		// })
		slices.SortFunc(*rs, func(a, b *SearchResult) int {
			return cmp.Compare[float64](b.Score, a.Score)
		})

		var r *SearchResult
		for i := topN; i < len(*rs); i++ {
			r = (*rs)[i]

			// do not forget to recycle the filtered result
			idx.RecycleSearchResult(r)
		}
		*rs = (*rs)[:topN]
	}

	if debug {
		log.Debugf("%s (%d bp): finished chaining (%d genome hits) in %s", query.seqID, len(query.seq), len(*rs), time.Since(startTime))
		startTime = time.Now()
	}

	if len(*rs) == 0 { // It happens when there's only one anchor which is shorter than MinSinglePrefix.
		poolSearchResults.Put(rs)
		return nil, nil
	}

	// 3.3) alignment

	rs2 := poolSearchResults.Get().(*[]*SearchResult)
	*rs2 = (*rs2)[:0]

	ch2 := make(chan *SearchResult, idx.opt.NumCPUs)

	// collect hits with good alignment
	go func() {
		for r := range ch2 {
			*rs2 = append(*rs2, r)
		}

		done <- 1
	}()

	cpr := idx.poolSeqComparator.Get().(*SeqComparator)
	// recycle the previous tree data
	cpr.RecycleIndex()
	err = cpr.Index(s) // index the query sequence
	if err != nil {
		checkError(err)
	}

	alignOption := &wfa.Options{GlobalAlignment: true}

	// do not do this, as multiple genome readers would compete in reading sequences from the same file.
	// sort.Slice(*rs, func(i, j int) bool { return (*rs)[i].BatchGenomeIndex < (*rs)[j].BatchGenomeIndex })

	// sort genomes according to the index in a genome data file for slightly faster file seeking
	if len(*rs) > 1 { // GenomeIndex
		// sort.Slice(*rs, func(i, j int) bool { return (*rs)[i].GenomeIndex < (*rs)[j].GenomeIndex })
		slices.SortFunc(*rs, func(a, b *SearchResult) int {
			return a.GenomeIndex - b.GenomeIndex
		})
	}

	// process bar
	var pbs *mpb.Progress
	var bar *mpb.Bar
	var chDuration chan time.Duration
	var doneDuration chan int
	if debug {
		pbs = mpb.New(mpb.WithWidth(40), mpb.WithOutput(os.Stderr))
		bar = pbs.AddBar(int64(len(*rs)),
			mpb.PrependDecorators(
				decor.Name("checked genomes: ", decor.WC{W: len("checked genomes: "), C: decor.DindentRight}),
				decor.Name("", decor.WCSyncSpaceR),
				decor.CountersNoUnit("%d / %d", decor.WCSyncWidth),
			),
			mpb.AppendDecorators(
				decor.Name("ETA: ", decor.WC{W: len("ETA: ")}),
				decor.EwmaETA(decor.ET_STYLE_GO, 1024),
				decor.OnComplete(decor.Name(""), ". done"),
			),
		)

		chDuration = make(chan time.Duration, idx.opt.NumCPUs)
		doneDuration = make(chan int)
		go func() {
			for t := range chDuration {
				bar.EwmaIncrBy(1, t)
			}
			doneDuration <- 1
		}()
	}

	fcpus := float64(idx.opt.NumCPUs)
	for _, r := range *rs { // multiple references
		tokens <- 1
		wg.Add(1)

		go func(r *SearchResult) { // for a reference genome
			timeStart := time.Now()
			defer func() {
				<-tokens
				if debug {
					chDuration <- time.Duration(float64(time.Since(timeStart)) / fcpus)
				}
				wg.Done()
			}()

			// -----------------------------------------------------
			// alignment

			minQcovGnm := idx.opt.MinQueryAlignedFractionInAGenome
			minQcovHSP := idx.seqCompareOption.MinAlignedFraction
			minPIdent := idx.seqCompareOption.MinIdentity
			maxEvalue := idx.opt.MaxEvalue
			extLen := idx.opt.ExtendLength
			extLen2 := idx.opt.ExtendLength2
			contigInterval := idx.contigInterval
			outSeq := idx.opt.OutputSeq

			algn := wfa.New(wfa.DefaultPenalties, alignOption)
			algn.AdaptiveReduction(wfa.DefaultAdaptiveOption)
			// algn.AdaptiveReduction(&wfa.AdaptiveReductionOption{
			// 	MinWFLen:    10,
			// 	MaxDistDiff: 50,
			// })

			// from blastn_values_2_3 in ncbi-blast-2.15.0+-src/c++/src/algo/blast/core/blast_stat.c
			fScoreAndEvalue := scoreAndEvalue(2, -3, 5, 2, int(idx.totalBases), 0.625, 0.41)

			var _qseq, _tseq []byte
			var cigar *wfa.AlignmentResult
			// var op *wfa.CIGARRecord
			var op uint64
			var Q, A, T *[]byte

			// -----------------------------------------------------

			refBatch := r.GenomeBatch
			refID := r.GenomeIndex

			var rdr *genome.Reader
			// sequence reader
			if idx.hasGenomeRdrs {
				rdr = <-idx.poolGenomeRdrs[refBatch]
			} else {
				idx.openFileTokens <- 1 // genome file
				fileGenome := filepath.Join(idx.path, DirGenomes, batchDir(refBatch), FileGenomes)
				rdr, err = genome.NewReader(fileGenome)
				if err != nil {
					checkError(fmt.Errorf("failed to read genome data file: %s", err))
				}
			}

			var sub *SubstrPair
			qlen := len(s)
			var rc bool
			var qb, qe, tb, te, tBegin, tEnd, qBegin, qEnd int
			var l, iSeq, iSeqPre, tPosOffsetBegin, tPosOffsetEnd int
			var _begin, _end int

			sds := poolSimilarityDetails.Get().(*[]*SimilarityDetail) // HSPs in a reference

			// fragments of a HSP.
			// Since HSP fragments in a HSP might come from different contigs.
			// Multiple contigs are concatenated, remember?
			// So we need to create seperate HPSs for these fragments.
			var crChains2 *[]*Chain2Result

			// for remove duplicated alignments
			var duplicated bool
			hashes := poolHashes.Get().(*map[uint64]interface{})
			clear(*hashes)
			var hash uint64

			var tSeq *genome.Genome

			// sort chains according to coordinates for faster file seeking
			if len(*r.Chains) > 1 {
				// sort.Slice(*r.Chains, func(i, j int) bool {
				// 	return (*r.Subs)[(*(*r.Chains)[i])[0]].TBegin < (*r.Subs)[(*(*r.Chains)[j])[0]].TBegin
				// })
				slices.SortFunc(*r.Chains, func(a, b *[]int) int {
					return int((*r.Subs)[(*a)[0]].TBegin - (*r.Subs)[(*b)[0]].TBegin)
				})
			}

			// check sequences from all chains
			var nSeeds int
			for i, chain := range *r.Chains { // for each lexichash chain
				// ------------------------------------------------------------------------
				// extract subsequence from the refseq for comparing

				// fmt.Printf("\n----------------- [ lexichash chain %d ] --------------\n", i+1)
				// for _i, _c := range *chain {
				// 	fmt.Printf("  %d, %s\n", _i, (*r.Subs)[_c])
				// }

				nSeeds = len(*chain)

				// the first seed pair
				sub = (*r.Subs)[(*chain)[0]]
				// fmt.Printf("  first: %s\n", sub)
				qb = int(sub.QBegin)
				tb = int(sub.TBegin)

				// the last seed pair
				sub = (*r.Subs)[(*chain)[nSeeds-1]]
				// fmt.Printf("  last: %s\n", sub)
				qe = int(sub.QBegin) + int(sub.Len) - 1
				te = int(sub.TBegin) + int(sub.Len) - 1
				// fmt.Printf("  (%d, %d) vs (%d, %d) rc:%v\n", qb, qe, tb, te, rc)

				if nSeeds == 1 { // if there's only one seed, need to check the strand information
					rc = sub.QRC != sub.TRC
				} else { // check the strand according to coordinates of seeds
					rc = tb > int(sub.TBegin)
				}
				// fmt.Printf("  rc: %v\n", rc)

				// recycle chain ASAP
				*chain = (*chain)[:0]
				poolChain.Put(chain)
				(*r.Chains)[i] = nil

				// extend the locations in the reference
				if rc { // reverse complement
					// tBegin = int(sub.TBegin) - min(qlen-qe-1, extLen)
					tBegin = int(sub.TBegin) - extLen
					if tBegin < 0 {
						tBegin = 0
					}
					// tEnd = tb + int(sub.Len) - 1 + min(qb, extLen)
					tEnd = tb + int(sub.Len) - 1 + extLen
				} else {
					// tBegin = tb - min(qb, extLen)
					tBegin = tb - extLen
					if tBegin < 0 {
						tBegin = 0
					}
					// tEnd = te + min(qlen-qe-1, extLen)
					tEnd = te + extLen
				}

				// extend the locations in the query
				qBegin = qb - min(qb, extLen)
				qEnd = qe + min(qlen-qe-1, extLen)

				// extract target sequence for comparison.
				// Right now, we fetch seq from disk for each seq,
				// In the future, we might buffer frequently accessed references for improving speed.
				tSeq, err = rdr.SubSeq3(refID, tBegin, tEnd, tSeq)
				if err != nil {
					checkError(err)
				}
				// this happens when the matched sequene is the last one in the gneome
				if len(tSeq.Seq) < tEnd-tBegin+1 {
					tEnd -= tEnd - tBegin + 1 - len(tSeq.Seq)
				}

				if rc { // reverse complement
					RC(tSeq.Seq)
				}

				// fmt.Printf("---------\nchain:%d, query:%d-%d, subject:%d.%d:%d-%d(len:%d), rc:%v, genome:%s, seq:%s\n",
				// 	i+1, qBegin+1, qEnd+1, refBatch, refID, tBegin+1, tEnd+1, tEnd-tBegin+1, rc, r.ID, tSeq.ID)
				// fmt.Printf("%s\n", tSeq.Seq)

				// ------------------------------------------------------------------------
				// comparing the two sequences with pseudo-alignment

				// fmt.Printf("qBegin: %d, qEnd: %d, len(tseq): %d\n", qBegin, qEnd, len(tSeq.Seq))
				cr, err := cpr.Compare(uint32(qBegin), uint32(qEnd), tSeq.Seq, qlen)
				if err != nil {
					checkError(err)
				}
				if cr == nil {
					// recycle target sequence
					// genome.RecycleGenome(tSeq)

					continue
				}

				if len(r.ID) == 0 { // record genome information, do it once
					r.ID = append(r.ID, tSeq.ID...)
					// if debug {
					// 	log.Debugf("  checking genome: %s", r.ID)
					// }
					r.GenomeSize = tSeq.GenomeSize
				}

				iSeqPre = -1 // the index of previous sequence in this HSP

				crChains2 = poolChains2.Get().(*[]*Chain2Result)

				for _i, c := range *cr.Chains { // for each pseudo alignment chain
					qb, qe, tb, te = c.QBegin, c.QEnd, c.TBegin, c.TEnd
					// fmt.Printf("\n--------------------------------------------\n")
					// fmt.Printf("--- %s: lexichash chain: %d, pseudo alignment chain: %d ---\n", r.ID, i+1, _i+1)
					// fmt.Printf("q: %d-%d, t: %d-%d\n", qb, qe, tb, te)

					// ------------------------------------------------------------
					// get the index of target seq according to the position

					iSeq = 0
					tPosOffsetBegin = 0 // start position of current sequence
					tPosOffsetEnd = 0   // end pososition of current sequence
					var j int
					if tSeq.NumSeqs > 1 { // just for genomes with multiple contigs
						iSeq = -1
						// ===========aaaaaaa================================aaaaaaa=======
						//                   | tPosOffsetBegin              | tPosOffsetEnd
						//                     tb ---------------te (matched region, substring region)

						// fmt.Printf("genome: %s, nSeqs: %d\n", tSeq.ID, tSeq.NumSeqs)
						// fmt.Printf("  qb: %d, qe: %d\n", qb, qe)
						// fmt.Printf("  tBegin: %d, tEnd: %d, tb: %d, te: %d, rc: %v\n", tBegin, tEnd, tb, te, rc)

						// minusing K is because the interval A's might be matched.
						if rc {
							_begin, _end = tEnd-te+K, tEnd-tb-K
						} else {
							_begin, _end = tBegin+tb+K, tBegin+te-K
						}

						if _begin >= _end { // sequences shorter than 2*k
							// fmt.Printf("  _begin (%d) >= _end (%d)\n", _begin, _end)
							if rc {
								_begin, _end = tEnd-te, tEnd-tb
							} else {
								_begin, _end = tBegin+tb, tBegin+te
							}
						}

						// fmt.Printf("  try %d: %d-%d\n", j, _begin, _end)

						for j, l = range tSeq.SeqSizes {
							// end position of current contig
							tPosOffsetEnd += l - 1 // length sum of 0..j

							// fmt.Printf("  seq %d: %d-%d\n", j, tPosOffsetBegin, tPosOffsetEnd)

							// +K -K is because chained region might have little overlaps with contig intervals
							if _begin+K >= tPosOffsetBegin && _end-K <= tPosOffsetEnd {
								iSeq = j
								// fmt.Printf("    iSeq: %d, tPosOffsetBegin: %d, tPosOffsetEnd: %d, seqlen: %d\n",
								// 	iSeq, tPosOffsetBegin, tPosOffsetEnd, l)
								break
							} else if _end < tPosOffsetBegin { // no need to find
								iSeq = -1
								break
							}

							tPosOffsetEnd += contigInterval + 1
							tPosOffsetBegin = tPosOffsetEnd // begin position of the next contig
						}

						// it will not happen now.
						if iSeq < 0 { // this means the aligned sequence crosses two sequences.
							// fmt.Printf("  invalid fragment: seqid: %s, aligned: %d, %d-%d, rc:%v, %d-%d\n",
							// 	tSeq.ID, cr.AlignedBases, tBegin, tEnd, rc, _begin, _end)

							poolChain2.Put(c)
							(*cr.Chains)[_i] = nil

							continue
						}

						if iSeqPre >= 0 && iSeq != iSeqPre { // two HSP fragments belong to different sequences ~~~~~
							// fmt.Printf("  %d != %d\n", iSeq, iSeqPre)

							iSeq0 := iSeq // used to restore the value later

							iSeq = iSeqPre // do not want to change iSeq below to iSeqPre

							// ------------------------------------------------------------
							// convert the positions

							// fmt.Printf("  aligned: (%d, %d) vs (%d, %d) rc:%v\n", qb, qe, tb, te, rc)
							c.QBegin = qb
							c.QEnd = qe
							// alignments might belong to different seqs, so we have to store it for later use
							c.tPosOffsetBegin = tPosOffsetBegin
							if rc {
								c.TBegin = tBegin - tPosOffsetBegin + (len(tSeq.Seq) - te - 1)
								if c.TBegin < 0 { // position in the interval
									c.QEnd += c.TBegin
									c.AlignedBasesQ += c.TBegin
									c.TBegin = 0
								}
								c.TEnd = tBegin - tPosOffsetBegin + (len(tSeq.Seq) - tb - 1)
								if c.TEnd > tSeq.SeqSizes[iSeq]-1 {
									c.QBegin += c.TEnd - (tSeq.SeqSizes[iSeq] - 1)
									c.TEnd = tSeq.SeqSizes[iSeq] - 1
								}
							} else {
								// fmt.Printf("tBegin: %d, tPosOffsetBegin: %d, tPosOffsetEnd: %d\n",
								// 	tBegin, tPosOffsetBegin, tPosOffsetEnd)
								// fmt.Printf("tb: %d, te: %d, tBegin+tb: %d, tBegin+te: %d\n", tb, te, tBegin+tb, tBegin+te)
								c.TBegin = tBegin - tPosOffsetBegin + tb
								if c.TBegin < 0 { // position in the interval
									c.QBegin -= c.TBegin
									c.AlignedBasesQ += c.TBegin
									c.TBegin = 0
								}
								c.TEnd = tBegin - tPosOffsetBegin + te
								// fmt.Printf("tmp: t: %d-%d, seqlen: %d \n", c.TBegin, c.TEnd, tSeq.SeqSizes[iSeq])
								if c.TEnd > tSeq.SeqSizes[iSeq]-1 {
									c.QEnd -= c.TEnd - (tSeq.SeqSizes[iSeq] - 1)
									c.TEnd = tSeq.SeqSizes[iSeq] - 1
								}
							}
							c.MaxExtLen = tSeq.SeqSizes[iSeq] - 1 - c.TEnd
							// fmt.Printf("  adjusted: (%d, %d) vs (%d, %d) rc:%v, iSeq: %d, %s\n", c.QBegin, c.QEnd, c.TBegin, c.TEnd, rc, iSeq, *tSeq.SeqIDs[iSeq])

							// ------------------------------------------------------------

							// fmt.Printf("  add previous one: %d fragments, aligned-bases: %d\n", len(*crChains2), (*crChains2)[0].AlignedBases)

							if len(*crChains2) > 0 { // it might be empty after duplicated results are removed
								// only include valid chains
								r2 := poolSeqComparatorResult.Get().(*SeqComparatorResult)
								r2.Update2(crChains2, cr.QueryLen) // r2.Chains = crChains2
								// there's no need
								// sort.Slice(*r2.Chains, func(i, j int) bool {
								// 	return (*r2.Chains)[i].AlignedBasesQ >= (*r2.Chains)[j].AlignedBasesQ
								// })

								hasResult := false
								// j := 0
								var start, end int
								var _s1, _e1, _s2, _e2 int // extend length
								var similarityScore, maxSimilarityScore float64
								for i, c := range *r2.Chains {
									if c.QBegin >= c.QEnd+1 { // rare case when the contig interval is two small
										poolChain2.Put(c)
										(*r2.Chains)[i] = nil
										continue
									}

									if rc {
										start, end = tEnd-c.TEnd-c.tPosOffsetBegin, tEnd-c.TBegin-c.tPosOffsetBegin+1
									} else {
										start, end = c.tPosOffsetBegin+c.TBegin-tBegin, c.tPosOffsetBegin+c.TEnd-tBegin+1

									}
									if start >= end {
										// fmt.Println(string(*tSeq.SeqIDs[0]))
										poolChain2.Put(c)
										(*r2.Chains)[i] = nil
										continue
									}

									// _qseq = s[c.QBegin : c.QEnd+1]
									// _tseq = tSeq.Seq[start:end]
									_qseq, _tseq, _s1, _e1, _s2, _e2, err = extendMatch(s, tSeq.Seq, c.QBegin, c.QEnd+1, start, end, extLen2, c.TBegin, c.MaxExtLen, rc)
									if err != nil {
										checkError(fmt.Errorf("fail to extend aligned region"))
									}

									// fmt.Printf("q: %s\nt: %s\n", _qseq, _tseq)
									cigar, err = algn.Align(_qseq, _tseq)
									if err != nil {
										checkError(fmt.Errorf("fail to align sequence"))
									}

									// score and e-value
									c.Score, c.BitScore, c.Evalue = fScoreAndEvalue(len(_qseq), cigar)
									if c.Evalue > maxEvalue {
										poolChain2.Put(c)
										(*r2.Chains)[i] = nil
										continue
									}

									// update sequence regions
									c.QBegin -= _s1
									c.QEnd += _e1

									c.QBegin = c.QBegin + cigar.QBegin - 1
									c.QEnd = c.QEnd - (len(_qseq) - cigar.QEnd)
									if rc {
										c.TBegin -= _e2
										c.TEnd += _s2

										c.TBegin = c.TBegin + (len(_tseq) - cigar.TEnd)
										c.TEnd = c.TEnd - cigar.TBegin - 1
									} else {
										c.TBegin -= _s2
										c.TEnd += _e2

										c.TBegin = c.TBegin + cigar.TBegin - 1
										c.TEnd = c.TEnd - (len(_tseq) - cigar.TEnd)
									}

									c.AlignedBasesQ = c.QEnd - c.QBegin + 1
									c.AlignedLength = int(cigar.AlignLen)
									c.MatchedBases = int(cigar.Matches)
									c.Gaps = int(cigar.Gaps)
									c.AlignedFraction = float64(c.AlignedBasesQ) / float64(cr.QueryLen) * 100
									if c.AlignedFraction > 100 {
										c.AlignedFraction = 100
									}
									c.PIdent = float64(c.MatchedBases) / float64(cigar.AlignLen) * 100

									if !outSeq {
										wfa.RecycleAlignmentResult(cigar)
									}

									if c.AlignedFraction < minQcovHSP || c.PIdent < minPIdent {
										poolChain2.Put(c)
										(*r2.Chains)[i] = nil
										continue
									}

									if outSeq {
										if c.CIGAR == nil {
											c.CIGAR = make([]byte, 0, 128)
											c.QSeq = make([]byte, 0, cigar.AlignLen)
											c.TSeq = make([]byte, 0, cigar.AlignLen)
											c.Alignment = make([]byte, 0, cigar.AlignLen)
										} else {
											c.CIGAR = c.CIGAR[:0]
											c.QSeq = c.QSeq[:0]
											c.TSeq = c.TSeq[:0]
											c.Alignment = c.Alignment[:0]
										}

										for _, op = range cigar.Ops {
											// c.CIGAR = append(c.CIGAR, []byte(strconv.Itoa(int(op.N)))...)
											c.CIGAR = append(c.CIGAR, []byte(strconv.Itoa(int(op&4294967295)))...)
											// c.CIGAR = append(c.CIGAR, op.Op)
											c.CIGAR = append(c.CIGAR, byte(op>>32))
										}

										Q, A, T = cigar.AlignmentText(&_qseq, &_tseq, true)

										c.QSeq = append(c.QSeq, *Q...)
										c.TSeq = append(c.TSeq, *T...)
										c.Alignment = append(c.Alignment, *A...)

										wfa.RecycleAlignmentText(Q, A, T)
										wfa.RecycleAlignmentResult(cigar)
									}

									similarityScore = float64(c.BitScore) * c.PIdent
									if similarityScore > maxSimilarityScore {
										maxSimilarityScore = similarityScore
									}
									hasResult = true
								}

								if hasResult {
									sd := poolSimilarityDetail.Get().(*SimilarityDetail)
									sd.RC = rc
									// sd.Chain = (*r.Chains)[i]
									sd.NSeeds = nSeeds
									sd.Similarity = r2
									// sd.SimilarityScore = float64(r2.AlignedBases) * (*r2.Chains)[j].PIdent // chain's aligned base * pident of 1st hsp.
									sd.SimilarityScore = maxSimilarityScore
									sd.SeqID = sd.SeqID[:0]
									// fmt.Printf("target seq a: iSeq:%d, %s, pident:%f\n", iSeq, *tSeq.SeqIDs[iSeq], (*r2.Chains)[j].PIdent)
									sd.SeqID = append(sd.SeqID, (*tSeq.SeqIDs[iSeq])...)
									sd.SeqLen = tSeq.SeqSizes[iSeq]

									*sds = append(*sds, sd)
								} else {
									RecycleChaining2Result(r2.Chains)
								}
							}

							// ----------

							// create anther HSP
							iSeqPre = -1
							crChains2 = poolChains2.Get().(*[]*Chain2Result)

							// ------------------------------------------------------------
							// remove duplicated alignments

							hash = wyhash.HashString(fmt.Sprintf("[%d, %d] vs [%d, %d], %v, %d", c.QBegin, c.QEnd, c.TBegin, c.TEnd, rc, iSeq), 0)
							if _, duplicated = (*hashes)[hash]; duplicated {
								// fmt.Printf("  duplicated: (%d, %d) vs (%d, %d) rc:%v\n", c.QBegin, c.QEnd, c.TBegin, c.TEnd, rc)
								poolChain2.Put(c)
								(*cr.Chains)[_i] = nil
							} else {
								*crChains2 = append(*crChains2, c)
								(*hashes)[hash] = struct{}{}
							}

							iSeq = iSeq0
							continue
						}
					}
					iSeqPre = iSeq

					// ------------------------------------------------------------
					// convert the positions

					// fmt.Printf("  aligned: (%d, %d) vs (%d, %d) rc:%v\n", qb, qe, tb, te, rc)
					c.QBegin = qb
					c.QEnd = qe
					// alignments might belong to different seqs, so we have to store it for later use
					c.tPosOffsetBegin = tPosOffsetBegin
					if rc {
						c.TBegin = tBegin - tPosOffsetBegin + (len(tSeq.Seq) - te - 1)
						// fmt.Printf("c.TBegin: %d\n", c.TBegin)
						if c.TBegin < 0 { // position in the interval
							c.QEnd += c.TBegin
							c.AlignedBasesQ += c.TBegin
							c.TBegin = 0
						}
						c.TEnd = tBegin - tPosOffsetBegin + (len(tSeq.Seq) - tb - 1)
						// fmt.Printf("c.TEnd: %d\n", c.TEnd)
						if c.TEnd > tSeq.SeqSizes[iSeq]-1 {
							c.QBegin += c.TEnd - (tSeq.SeqSizes[iSeq] - 1)
							c.TEnd = tSeq.SeqSizes[iSeq] - 1
						}
					} else {
						c.TBegin = tBegin - tPosOffsetBegin + tb
						if c.TBegin < 0 { // position in the interval
							c.QBegin -= c.TBegin
							c.AlignedBasesQ += c.TBegin
							c.TBegin = 0
						}
						c.TEnd = tBegin - tPosOffsetBegin + te
						if c.TEnd > tSeq.SeqSizes[iSeq]-1 {
							c.QEnd -= c.TEnd - (tSeq.SeqSizes[iSeq] - 1)
							c.TEnd = tSeq.SeqSizes[iSeq] - 1
						}
					}
					c.MaxExtLen = tSeq.SeqSizes[iSeq] - 1 - c.TEnd
					// fmt.Printf("  adjusted: (%d, %d) vs (%d, %d) rc:%v, %s\n", c.QBegin, c.QEnd, c.TBegin, c.TEnd, rc, *tSeq.SeqIDs[iSeq])

					// ------------------------------------------------------------
					// remove duplicated alignments

					hash = wyhash.HashString(fmt.Sprintf("[%d, %d] vs [%d, %d], %v, %d", c.QBegin, c.QEnd, c.TBegin, c.TEnd, rc, iSeq), 0)
					if _, duplicated = (*hashes)[hash]; duplicated {
						// fmt.Printf("  duplicated: (%d, %d) vs (%d, %d) rc:%v\n", c.QBegin, c.QEnd, c.TBegin, c.TEnd, rc)
						poolChain2.Put(c)
						(*cr.Chains)[_i] = nil
					} else {
						*crChains2 = append(*crChains2, c)
						(*hashes)[hash] = struct{}{}
					}
				}

				// fmt.Printf("  add current one: %d fragments, aligned-bases: %d\n", len(*crChains2), (*crChains2)[0].AlignedBases)

				if iSeq >= 0 {
					if len(*crChains2) > 0 { // it might be empty after duplicated results are removed
						// only include valid chains
						r2 := poolSeqComparatorResult.Get().(*SeqComparatorResult)
						r2.Update2(crChains2, cr.QueryLen) // r2.Chains = crChains2
						// there's no need
						// sort.Slice(*r2.Chains, func(i, j int) bool {
						// 	return (*r2.Chains)[i].AlignedBasesQ >= (*r2.Chains)[j].AlignedBasesQ
						// })

						hasResult := false
						// j := 0
						var start, end int
						var _s1, _e1, _s2, _e2 int // extend length
						var similarityScore, maxSimilarityScore float64
						for i, c := range *r2.Chains {
							if c.QBegin >= c.QEnd+1 { // rare case when the contig interval is two small
								poolChain2.Put(c)
								(*r2.Chains)[i] = nil
								continue
							}

							// Attention, it's different from previous code
							if rc {
								start, end = tEnd-c.TEnd-c.tPosOffsetBegin, tEnd-c.TBegin-c.tPosOffsetBegin+1
							} else {
								start, end = c.tPosOffsetBegin+c.TBegin-tBegin, c.tPosOffsetBegin+c.TEnd-tBegin+1
							}
							if start >= end {
								// fmt.Println(string(*tSeq.SeqIDs[0]))
								poolChain2.Put(c)
								(*r2.Chains)[i] = nil
								continue
							}

							// _qseq = s[c.QBegin : c.QEnd+1]
							// _tseq = tSeq.Seq[start:end]
							_qseq, _tseq, _s1, _e1, _s2, _e2, err = extendMatch(s, tSeq.Seq, c.QBegin, c.QEnd+1, start, end, extLen2, c.TBegin, c.MaxExtLen, rc)
							if err != nil {
								checkError(fmt.Errorf("fail to extend aligned region"))
							}

							// fmt.Printf("q: %s\nt: %s\n", _qseq, _tseq)
							cigar, err = algn.Align(_qseq, _tseq)
							if err != nil {
								checkError(fmt.Errorf("fail to align sequence"))
							}

							// score and e-value
							c.Score, c.BitScore, c.Evalue = fScoreAndEvalue(len(_qseq), cigar)
							if c.Evalue > maxEvalue {
								poolChain2.Put(c)
								(*r2.Chains)[i] = nil
								continue
							}

							// update sequence regions
							c.QBegin -= _s1
							c.QEnd += _e1

							c.QBegin = c.QBegin + (cigar.QBegin - 1)
							c.QEnd = c.QEnd - (len(_qseq) - cigar.QEnd)
							if rc {
								c.TBegin -= _e2
								c.TEnd += _s2

								c.TBegin = c.TBegin + (len(_tseq) - cigar.TEnd)
								c.TEnd = c.TEnd - (cigar.TBegin - 1)
							} else {
								c.TBegin -= _s2
								c.TEnd += _e2

								c.TBegin = c.TBegin + (cigar.TBegin - 1)
								c.TEnd = c.TEnd - (len(_tseq) - cigar.TEnd)
							}

							// fmt.Println(c.QBegin, c.QEnd, c.TBegin, c.TEnd)

							c.AlignedBasesQ = c.QEnd - c.QBegin + 1
							c.AlignedLength = int(cigar.AlignLen)
							c.MatchedBases = int(cigar.Matches)
							c.Gaps = int(cigar.Gaps)
							c.AlignedFraction = float64(c.AlignedBasesQ) / float64(cr.QueryLen) * 100
							if c.AlignedFraction > 100 {
								c.AlignedFraction = 100
							}
							c.PIdent = float64(c.MatchedBases) / float64(cigar.AlignLen) * 100

							if !outSeq {
								wfa.RecycleAlignmentResult(cigar)
							}

							if c.AlignedFraction < minQcovHSP || c.PIdent < minPIdent {
								poolChain2.Put(c)
								(*r2.Chains)[i] = nil
								continue
							}

							if outSeq {
								if c.CIGAR == nil {
									c.CIGAR = make([]byte, 0, 128)
									c.QSeq = make([]byte, 0, cigar.AlignLen)
									c.TSeq = make([]byte, 0, cigar.AlignLen)
									c.Alignment = make([]byte, 0, cigar.AlignLen)
								} else {
									c.CIGAR = c.CIGAR[:0]
									c.QSeq = c.QSeq[:0]
									c.TSeq = c.TSeq[:0]
									c.Alignment = c.Alignment[:0]
								}

								for _, op = range cigar.Ops {
									// c.CIGAR = append(c.CIGAR, []byte(strconv.Itoa(int(op.N)))...)
									c.CIGAR = append(c.CIGAR, []byte(strconv.Itoa(int(op&4294967295)))...)
									// c.CIGAR = append(c.CIGAR, op.Op)
									c.CIGAR = append(c.CIGAR, byte(op>>32))
								}

								Q, A, T = cigar.AlignmentText(&_qseq, &_tseq, true)

								c.QSeq = append(c.QSeq, *Q...)
								c.TSeq = append(c.TSeq, *T...)
								c.Alignment = append(c.Alignment, *A...)

								wfa.RecycleAlignmentText(Q, A, T)
								wfa.RecycleAlignmentResult(cigar)
							}

							similarityScore = float64(c.BitScore) * c.PIdent
							if similarityScore > maxSimilarityScore {
								maxSimilarityScore = similarityScore
							}
							hasResult = true
						}

						if hasResult {
							sd := poolSimilarityDetail.Get().(*SimilarityDetail)
							sd.RC = rc
							sd.NSeeds = nSeeds
							sd.Similarity = r2
							// sd.SimilarityScore = float64(r2.AlignedBases) * (*r2.Chains)[j].PIdent // chain's aligned base * pident of 1st hsp.
							sd.SimilarityScore = maxSimilarityScore
							sd.SeqID = sd.SeqID[:0]
							// fmt.Printf("target seq b: iSeq:%d, %s, pident:%f\n", iSeq, *tSeq.SeqIDs[iSeq], (*r2.Chains)[j].PIdent)
							sd.SeqID = append(sd.SeqID, (*tSeq.SeqIDs[iSeq])...)
							sd.SeqLen = tSeq.SeqSizes[iSeq]

							*sds = append(*sds, sd)
						} else {
							RecycleChaining2Result(r2.Chains)
						}
					}
				}

				// recycle target sequence

				*cr.Chains = (*cr.Chains)[:0]
				poolChains2.Put(cr.Chains)
				// genome.RecycleGenome(tSeq)
			}

			// recyle chains ASAP
			RecycleChainingResult(r.Chains)
			r.Chains = nil

			RecycleSubstrPairs(poolSub, poolSubs, r.Subs)
			r.Subs = nil

			genome.RecycleGenome(tSeq)

			poolHashes.Put(hashes)
			wfa.RecycleAligner(algn)

			if len(*sds) == 0 { // no valid alignments
				idx.RecycleSimilarityDetails(sds)
				idx.RecycleSearchResult(r) // do not forget to recycle unused objects

				if idx.hasGenomeRdrs {
					idx.poolGenomeRdrs[refBatch] <- rdr
				} else {
					err = rdr.Close()
					if err != nil {
						checkError(fmt.Errorf("failed to close genome data file: %s", err))
					}
					<-idx.openFileTokens
				}

				return
			}

			if !idx.hasGenomeChunks { // if hasGenomeChunks, do not filter results now
				// compute aligned bases per genome
				var alignedBasesGenome int
				regions := poolRegions.Get().(*[]*[2]int)
				*regions = (*regions)[:0]
				for _, sd := range *sds {
					for _, c := range *sd.Similarity.Chains {
						if c != nil {
							region := poolRegion.Get().(*[2]int)
							region[0], region[1] = c.QBegin, c.QEnd
							*regions = append(*regions, region)
						}
					}
				}
				alignedBasesGenome = coverageLen(regions)
				recycleRegions(regions)

				// filter by query coverage per genome
				r.AlignedFraction = float64(alignedBasesGenome) / float64(len(s)) * 100
				if r.AlignedFraction > 100 {
					r.AlignedFraction = 100
				}
				if r.AlignedFraction < minQcovGnm { // no valid alignments
					idx.RecycleSimilarityDetails(sds)
					idx.RecycleSearchResult(r) // do not forget to recycle unused objects

					if idx.hasGenomeRdrs {
						idx.poolGenomeRdrs[refBatch] <- rdr
					} else {
						err = rdr.Close()
						if err != nil {
							checkError(fmt.Errorf("failed to close genome data file: %s", err))
						}
						<-idx.openFileTokens
					}
					return
				}
			}

			// Within each subject genome, alignments (HSP) are sorted by qcovHSP*pident
			// r.AlignResults = ars
			// sort.Slice(*sds, func(i, j int) bool {
			// 	return (*sds)[i].SimilarityScore > (*sds)[j].SimilarityScore
			// })
			slices.SortFunc(*sds, func(a, b *SimilarityDetail) int {
				return cmp.Compare[float64](b.SimilarityScore, a.SimilarityScore)
			})

			r.SimilarityDetails = sds

			// recycle genome reader
			if idx.hasGenomeRdrs {
				idx.poolGenomeRdrs[refBatch] <- rdr
			} else {
				err = rdr.Close()
				if err != nil {
					checkError(fmt.Errorf("failed to close genome data file: %s", err))
				}
				<-idx.openFileTokens
			}

			ch2 <- r
		}(r)
	}

	wg.Wait()
	close(ch2)
	<-done
	// process bar
	if debug {
		close(chDuration)
		<-doneDuration
		pbs.Wait()
	}
	poolSearchResults.Put(rs)

	// recycle this comparator
	idx.poolSeqComparator.Put(cpr)

	if debug {
		log.Debugf("%s (%d bp): finished alignment (%d genome hits) in %s", query.seqID, len(query.seq), len(*rs2), time.Since(startTime))
		startTime = time.Now()
	}

	if len(*rs2) == 0 {
		poolSearchResults.Put(rs2)
		return nil, nil
	}

	// merge search result from genome chunks, if has split genome
	if idx.hasGenomeChunks {
		var r, rp *SearchResult
		var i, j int
		var a uint64 // SearchResult.BatchGenomeIndex
		var ok bool
		var li *[]int
		var ptr uintptr
		gcIdx2List := idx.poolGenomeChunksIdx2List.Get().(*map[uint64]*[]int)
		gcPtr2List := idx.poolGenomeChunksPointer2List.Get().(*map[uintptr]*[]int)

		// collect "i"s belonging to the same genome
		for i, r = range *rs2 {
			a = r.BatchGenomeIndex
			if li, ok = (*gcIdx2List)[a]; !ok { // not a chunk of some genome
				continue
			}

			*li = append(*li, i)

			ptr = uintptr(unsafe.Pointer(li))
			if _, ok = (*gcPtr2List)[ptr]; !ok {
				(*gcPtr2List)[ptr] = li
			}
		}

		// merge alignments of the same genome
		for _, li = range *gcPtr2List {
			if len(*li) == 1 { // there's only one genome chunk
				// reset the list
				(*li) = (*li)[:0]
				continue
			}
			i = (*li)[0]
			rp = (*rs2)[i]
			for _, j = range (*li)[1:] { // merge j -> i
				r = (*rs2)[j]

				// only need to update SimilarityDetails, AlignedFraction
				*rp.SimilarityDetails = append(*rp.SimilarityDetails, *r.SimilarityDetails...)

				poolSearchResult.Put(r)
				(*rs2)[j] = nil
			}

			// reset the list
			(*li) = (*li)[:0]
		}

		// recycle datastructure
		clear(*gcPtr2List)
		idx.poolGenomeChunksPointer2List.Put(gcPtr2List)
		idx.poolGenomeChunksIdx2List.Put(gcIdx2List)

		if debug {
			log.Debugf("%s (%d bp): finished merging alignment results (%d genome hits) in %s", query.seqID, len(query.seq), len(*rs2), time.Since(startTime))
			startTime = time.Now()
		}

		// recompute query coverage per genome
		var alignedBasesGenome int
		minQcovGnm := idx.opt.MinQueryAlignedFractionInAGenome
		j = 0
		for _, r = range *rs2 {
			if r == nil {
				continue
			}

			// compute aligned bases per genome
			regions := poolRegions.Get().(*[]*[2]int)
			*regions = (*regions)[:0]
			for _, sd := range *r.SimilarityDetails {
				for _, c := range *sd.Similarity.Chains {
					if c != nil {
						region := poolRegion.Get().(*[2]int)
						region[0], region[1] = c.QBegin, c.QEnd
						*regions = append(*regions, region)
					}
				}
			}
			alignedBasesGenome = coverageLen(regions)
			recycleRegions(regions)

			// filter by query coverage per genome
			r.AlignedFraction = float64(alignedBasesGenome) / float64(len(s)) * 100
			if r.AlignedFraction > 100 {
				r.AlignedFraction = 100
			}
			if r.AlignedFraction < minQcovGnm { // no valid alignments
				idx.RecycleSearchResult(r) // do not forget to recycle unused objects

				continue
			}

			// Within each subject genome, alignments (HSP) are sorted by qcovHSP*pident
			// r.AlignResults = ars
			// sort.Slice(*r.SimilarityDetails, func(i, j int) bool {
			// 	return (*r.SimilarityDetails)[i].SimilarityScore > (*r.SimilarityDetails)[j].SimilarityScore
			// })
			slices.SortFunc(*r.SimilarityDetails, func(a, b *SimilarityDetail) int {
				return cmp.Compare[float64](b.SimilarityScore, a.SimilarityScore)
			})

			(*rs2)[j] = r
			j++
		}
		*rs2 = (*rs2)[:j]

		if debug {
			log.Debugf("%s (%d bp): finished filtering merged alignment results (%d genome hits) in %s", query.seqID, len(query.seq), len(*rs2), time.Since(startTime))
			startTime = time.Now()
		}
	}

	// sort all genomes, by qcovHSP*pident of the best alignment.
	// sort.Slice(*rs2, func(i, j int) bool {
	// 	return (*(*rs2)[i].SimilarityDetails)[0].SimilarityScore > (*(*rs2)[j].SimilarityDetails)[0].SimilarityScore
	// })
	slices.SortFunc(*rs2, func(a, b *SearchResult) int {
		return cmp.Compare[float64]((*b.SimilarityDetails)[0].SimilarityScore, (*a.SimilarityDetails)[0].SimilarityScore)
	})

	// ----------------------------------
	// In each target genome, sort alignments by target seq first
	// query
	//    result in a genome         SearchResult
	//        result in a sequence   SimilarityDetai
	//            alignments         Chain2Result
	//
	for _, r := range *rs2 {
		r.SortBySeqID()
	}

	if debug {
		log.Debugf("%s (%d bp): finished sorting alignment results (%d genome hits) in %s", query.seqID, len(query.seq), len(*rs2), time.Since(startTime))
	}

	return rs2, nil
}

// RC computes the reverse complement sequence
func RC(s []byte) []byte {
	n := len(s)
	for i := 0; i < n; i++ {
		s[i] = rcTable[s[i]]
	}
	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		s[i], s[j] = s[j], s[i]
	}
	return s
}

var rcTable = [256]byte{
	0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
	16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
	32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
	48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
	64, 84, 86, 71, 72, 69, 70, 67, 68, 73, 74, 77, 76, 75, 78, 79,
	80, 81, 89, 83, 65, 85, 66, 87, 88, 82, 90, 91, 92, 93, 94, 95,
	96, 116, 118, 103, 104, 101, 102, 99, 100, 105, 106, 109, 108, 107, 110, 111,
	112, 113, 121, 115, 97, 117, 98, 119, 120, 114, 122, 123, 124, 125, 126, 127,
	128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
	144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
	160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
	176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
	192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
	208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
	224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
	240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255,
}

var poolBounds = &sync.Pool{New: func() interface{} {
	tmp := make([]int, 128)
	return &tmp
}}

var poolHashes = &sync.Pool{New: func() interface{} {
	tmp := make(map[uint64]interface{}, 128)
	return &tmp
}}

// func parseKmerValue(v uint64) (int, int, int, int) {
// 	return int(v >> 47), int(v << 17 >> 47), int(v << 34 >> 35), int(v & 1)
// }

// func kmerValueString(v uint64) string {
// 	return fmt.Sprintf("batchIdx: %d, genomeIdx: %d, pos: %d, rc: %v",
// 		int(v>>47), int(v<<17>>47), int(v<<34>>35), v&1 > 0)
// }
