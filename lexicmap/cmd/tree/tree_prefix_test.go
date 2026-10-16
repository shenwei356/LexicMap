package tree

import (
	"fmt"
	"math/rand"
	"reflect"
	"slices"
	"sync"
	"testing"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
)

// Comparing target pointers proves that both APIs visit exactly the same
// subtree, including compressed edges, absent prefixes, and p clamping.
func TestPrefixDirectoryMatchesRoot(t *testing.T) {
	for _, tc := range []struct {
		k          uint8
		n          int
		compressed bool
		depth      uint8
	}{
		{5, 1024, false, 0}, {6, 4096, false, 6},
		{31, 4095, false, 0}, {31, 4096, false, 6},
		{31, 16384, false, 7}, {31, 65536, false, 8},
		{31, 4096, true, 6}, {32, 4096, false, 6},
	} {
		t.Run(fmt.Sprintf("k=%d/n=%d/compressed=%v", tc.k, tc.n, tc.compressed), func(t *testing.T) {
			tr := NewTree(tc.k)
			defer RecycleTree(tr)
			rng := rand.New(rand.NewSource(55))
			mask := ^uint64(0) >> (64 - 2*tc.k)
			entries := make([]BatchEntry, 0, tc.n*2)
			queries := []uint64{0, mask, mask >> 1}
			for i := 0; i < tc.n; i++ {
				// Distinct suffixes guarantee the intended directory size.
				key := (rng.Uint64() & mask &^ uint64(0xffff)) | uint64(i)
				if tc.k <= 6 {
					key = uint64(i)
				}
				if tc.compressed {
					key = uint64(0x1b)<<54 | uint64(i%16)<<50 | uint64(0x12345)<<16 | uint64(i)
				}
				entries = append(entries, BatchEntry{key, uint32(i*2 + 1)}, BatchEntry{key, uint32(i * 2)})
				if i < 128 {
					queries = append(queries, key, key^(uint64(1)<<(2*tc.k-9)), key^1)
				}
			}
			tr.InsertBatch(entries)
			tr.BuildPrefixDirectory()
			if tr.prefixBases != tc.depth {
				t.Fatalf("depth %d, want %d", tr.prefixBases, tc.depth)
			}
			for i := 0; i < 128; i++ {
				queries = append(queries, rng.Uint64()&mask)
			}
			for _, q := range queries {
				for p := uint8(0); p <= tc.k+1; p++ {
					got, want := tr.searchNode(q, p), tr.searchNodeFromRoot(q, p)
					if got != want {
						t.Fatalf("q=%x p=%d: target differs", q, p)
					}
				}
				if tr.searchNode(q, 255) != tr.searchNodeFromRoot(q, 255) {
					t.Fatal("p=255 differs")
				}
			}
			// Values include both strand bits and duplicates in their existing
			// insertion order. Match the public visitor and allocating API too.
			for _, q := range queries[:20] {
				var want, got []SearchResult
				p := uint8(11)
				n := tr.searchNodeFromRoot(q, p)
				if n != nil {
					recursiveWalk(n, func(key uint64, vals []uint32) bool {
						prefix := uint8(0)
						for j := tc.k; j > 0; j-- {
							if (q>>(2*(j-1)))&3 != (key>>(2*(j-1)))&3 {
								break
							}
							prefix++
						}
						want = append(want, SearchResult{key, prefix, slices.Clone(vals)})
						return false
					})
				}
				found := tr.SearchEach(q, p, func(key uint64, prefix uint8, vals []uint32) {
					got = append(got, SearchResult{key, prefix, slices.Clone(vals)})
				})
				if found != (n != nil) || !reflect.DeepEqual(got, want) {
					t.Fatal("visitor values/order differ")
				}
				rs, found := tr.Search(q, p)
				got = nil
				if found {
					for _, r := range *rs {
						got = append(got, *r)
					}
					if !reflect.DeepEqual(got, want) {
						t.Fatal("Search values/order differ")
					}
					tr.RecycleSearchResult(rs)
				}
			}
			// Lookup never writes shared directory state.
			var wg sync.WaitGroup
			for i := 0; i < 4; i++ {
				wg.Add(1)
				go func() {
					defer wg.Done()
					for _, q := range queries {
						if tr.searchNode(q, 11) != tr.searchNodeFromRoot(q, 11) {
							t.Error("parallel lookup differs")
						}
					}
				}()
			}
			wg.Wait()
		})
	}
}

func TestPrefixDirectoryInvalidation(t *testing.T) {
	tr := NewTree(31)
	for i := uint64(0); i < 4096; i++ {
		tr.Insert(i<<20, uint32(i))
	}
	tr.BuildPrefixDirectory()
	if tr.prefixBases == 0 {
		t.Fatal("missing directory")
	}
	tr.Insert(1, 999)
	if tr.prefixBases != 0 || tr.prefixStarts != nil {
		t.Fatal("Insert retained stale directory")
	}
	tr.BuildPrefixDirectory()
	if tr.searchNode(1, 31) != tr.searchNodeFromRoot(1, 31) {
		t.Fatal("inserted key missing")
	}
	RecycleTree(tr)
	if tr.prefixBases != 0 || tr.prefixStarts != nil {
		t.Fatal("recycled directory retained")
	}
	small := NewTree(31)
	defer RecycleTree(small)
	small.Insert(7, 3)
	small.BuildPrefixDirectory()
	if small.prefixBases != 0 || small.prefixStarts != nil {
		t.Fatal("small tree retained directory")
	}
	empty := NewTree(31)
	defer RecycleTree(empty)
	empty.BuildPrefixDirectory()
	if empty.prefixBases != 0 {
		t.Fatal("empty tree indexed")
	}
}

// Frozen pre-directory traversal is the semantic oracle.
func (t *Tree) searchNodeFromRoot(key uint64, p uint8) *node {
	if p < 1 {
		p = 1
	}
	k := t.k
	if p > k {
		p = k
	}
	var target *node
	n := t.root
	search := key
	var lenPrefix uint8
	var atleast uint8
	for {
		// Check for key exhaution
		if k == 0 {
			break
		}

		// Look for a child
		n = n.children[util.KmerBaseAt(search, k, 0)]
		if n == nil {
			break
		}

		// Consume the search prefix
		// if KmerHasPrefix(search, n.prefix, k, n.k) {
		// if search>>((k-n.k)<<1) == n.prefix { // manually inline code
		// this line is slow, because of RAM access of node information (cache miss)
		if util.MustKmerHasPrefix(search, n.prefix, k, n.k) {
			lenPrefix += n.k
			// already matched at least m bases
			// we can output all leaves below n
			if lenPrefix >= p {
				target = n
				break
			}

			search = util.KmerSuffix(search, k, n.k)
			k = k - n.k
		} else {
			// also check the prefix, because the prefix of some nodes
			// might be very long. Only checking prefix will ignore theses.
			// For example, the strings below shared 4 bases,
			// the third node would be this case.
			//   A C AGCT
			//   A C AGGC
			// lenPrefix += KmerLongestPrefix(search, n.prefix, k, n.k)
			// because k >= n.k
			//
			// lenPrefix += MustKmerLongestPrefix(search, n.prefix, k, n.k)
			// if lenPrefix >= m {
			// 	target = n
			// 	break
			// }

			atleast = p - lenPrefix // just check if prefixes of m - lenPrefix bases are equal
			if search>>((k-atleast)<<1) == n.prefix>>((n.k-atleast)<<1) {
				target = n
			}

			break
		}
	}

	return target
}
