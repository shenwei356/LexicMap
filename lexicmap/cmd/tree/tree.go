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

package tree

import (
	"math/bits"
	"sync"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
	"github.com/twotwotwo/sorts"
)

// leafNode is used to represent a value.
type leafNode struct {
	key uint64   // ALL the bases in the node, the k-mer
	val []uint32 // yes, multiple values
}

// Reset resets the leafNode object.
func (n *leafNode) Reset() {
	n.val = n.val[:0]
}

// node represents a node in the tree, it might be the root, inner or leaf node.
type node struct {
	prefix uint64 // prefix of the current node
	k      uint8  // bases length of the prefix

	numChildren uint8
	children    [4]*node // just use an array

	leaf *leafNode // optional
}

// Reset resets the node object.
func (n *node) Reset() {
	n.numChildren = 0
	for i := 0; i < 4; i++ {
		n.children[i] = nil
	}
	n.leaf = nil
}

// Tree is a radix tree for storing bit-packed k-mer information.
type Tree struct {
	k    uint8 // use a global K
	root *node // root node

	nodes  nodeArena
	leaves leafArena
}

// Each tree owns its own chunked arenas for nodes and leaves so that all
// allocations for a single tree are physically contiguous, drastically
// improving cache locality during Insert / Walk / Search. Chunks are fixed
// size, so the address of any element is stable for the tree's lifetime
// (and across Tree recycling, since chunks are reused in place).
const (
	nodeChunkSize = 1024
	leafChunkSize = 1024
)

type nodeArena struct {
	chunks [][]node
	cur    int
	pos    int
}

func (a *nodeArena) alloc() *node {
	if a.pos >= nodeChunkSize {
		a.cur++
		a.pos = 0
	}
	if a.cur >= len(a.chunks) {
		a.chunks = append(a.chunks, make([]node, nodeChunkSize))
	}
	n := &a.chunks[a.cur][a.pos]
	a.pos++
	n.Reset()
	return n
}

func (a *nodeArena) reset() {
	a.cur = 0
	a.pos = 0
}

type leafArena struct {
	chunks [][]leafNode
	cur    int
	pos    int
}

func (a *leafArena) alloc() *leafNode {
	if a.pos >= leafChunkSize {
		a.cur++
		a.pos = 0
	}
	if a.cur >= len(a.chunks) {
		a.chunks = append(a.chunks, make([]leafNode, leafChunkSize))
	}
	n := &a.chunks[a.cur][a.pos]
	a.pos++
	n.Reset()
	return n
}

func (a *leafArena) reset() {
	a.cur = 0
	a.pos = 0
}

var poolTree = &sync.Pool{New: func() interface{} {
	return &Tree{}
}}

// NewTree implements a radix tree for k-mer extractly or prefix querying.
func NewTree(k uint8) *Tree {
	t := poolTree.Get().(*Tree)
	t.k = k
	t.nodes.reset()
	t.leaves.reset()
	t.root = t.nodes.alloc()
	return t
}

// RecycleTree recycles the tree object. Node and leaf memory stays with
// the tree's arenas for reuse on the next NewTree from the pool.
func RecycleTree(t *Tree) {
	poolTree.Put(t)
}

// K returns the K value of k-mers.
func (t *Tree) K() int {
	return int(t.k)
}

// Insert is used to add a new entry or update
// an existing entry. Returns true if an existing record is updated.
func (t *Tree) Insert(key uint64, v uint32) bool {
	key0 := key // will save it into the leaf node
	k := t.k

	var parent *node
	n := t.root
	search := key // current key
	for {
		// Handle key exhaustion
		if k == 0 {
			if n.leaf != nil {
				n.leaf.val = append(n.leaf.val, v)
				return true
			}

			// n is not a leaf node, that means
			// the current key is a prefix of some other keys.
			leaf := t.leaves.alloc()
			leaf.key = key0
			leaf.val = append(leaf.val, v)
			n.leaf = leaf

			return false
		}

		// Look for the child
		parent = n
		firstBase := util.KmerBaseAt(search, k, 0)
		n = n.children[firstBase]

		// No child, create one
		if n == nil {
			n = t.nodes.alloc()
			n.prefix = search
			n.k = k
			leaf := t.leaves.alloc()
			leaf.key = key0
			leaf.val = append(leaf.val, v)
			n.leaf = leaf

			parent.children[firstBase] = n
			parent.numChildren++

			return false
		}

		// has a child -- exists a path

		// Fast path: the common case is that search fully contains n.prefix
		// (commonPrefix == n.k), in which case we only need an equality test,
		// not a precise LCP. Slow path computes the LCP only when we must split.
		shifted := search >> ((k - n.k) << 1)
		if shifted == n.prefix {
			search = util.KmerSuffix(search, k, n.k)
			k = k - n.k
			continue
		}

		// search and n.prefix diverge before n.k bases: need to split n.
		commonPrefix := uint8(bits.LeadingZeros64(shifted^n.prefix)>>1) + n.k - 32
		child := t.nodes.alloc()
		// o---<=8, here the prefix of one of the 8 is ---,
		child.prefix = util.KmerPrefix(search, k, commonPrefix)
		child.k = commonPrefix
		parent.children[firstBase] = child // change from n to c

		// child points to n now
		child.children[util.KmerBaseAt(n.prefix, n.k, commonPrefix)] = n
		child.numChildren++
		n.prefix = util.KmerSuffix(n.prefix, n.k, commonPrefix)
		n.k = n.k - commonPrefix

		// Create a new leaf node for the new key
		leaf := t.leaves.alloc()
		leaf.key = key0
		leaf.val = append(leaf.val, v)

		// the new key is a prefix of the old n, add the leaf node to this node. len(new) = len(prefix)
		search = util.KmerSuffix(search, k, commonPrefix)
		k = k - commonPrefix
		if k == 0 {
			child.leaf = leaf
			return false
		}

		// the new key and the key of node n share a prefix shorter than both of them
		// Create a new child node for the node
		n = t.nodes.alloc()
		n.prefix = search
		n.k = k
		n.leaf = leaf
		child.children[util.KmerBaseAt(search, k, 0)] = n

		child.numChildren++
		return false
	}
}

// BatchEntry is a (key, value) pair for InsertBatch.
type BatchEntry struct {
	Key uint64
	Val uint32
}

// BatchEntries is a sortable slice of BatchEntry, ordered by Key ascending.
// It implements both sort.Interface and sorts.Uint64Interface, so callers
// may pass it to sort.Sort, sorts.Quicksort, or sorts.ByUint64 (parallel
// radix sort). Since Key is bit-packed in lexicographic order, sorting by
// uint64 Key is equivalent to lexicographic order on the k-mer.
type BatchEntries []BatchEntry

func (s BatchEntries) Len() int           { return len(s) }
func (s BatchEntries) Less(i, j int) bool { return s[i].Key < s[j].Key }
func (s BatchEntries) Swap(i, j int)      { s[i], s[j] = s[j], s[i] }
func (s BatchEntries) Key(i int) uint64   { return s[i].Key }

// InsertBatch builds the tree from the given entries in O(n log n) time
// (dominated by sort) without per-key tree descent or split operations.
// The tree must be empty (freshly created via NewTree). entries is sorted
// in place by Key.
func (t *Tree) InsertBatch(entries []BatchEntry) {
	if len(entries) == 0 {
		return
	}
	K := t.k

	// slowest:
	// sort.Slice(entries, func(i, j int) bool {
	// 	return entries[i].Key < entries[j].Key
	// })
	//
	// much slower:
	// slices.SortFunc(entries, func(a, b BatchEntry) int {
	// 	if a.Key < b.Key {
	// 		return -1
	// 	} else if a.Key > b.Key {
	// 		return 1
	// 	}
	// 	return 0
	// })
	//
	// slightly slower:
	// sorts.Quicksort(BatchEntries(entries))
	//
	sorts.ByUint64(BatchEntries(entries))

	// Path stack: each entry is a node on the current rightmost path from
	// root to the most recently inserted leaf. depth[i] is the cumulative
	// k-mer base depth at that node (i.e. number of bases from root to the
	// END of nodes[i].prefix). depth[0] == 0 for the root.
	stack := make([]frame, 0, int(K)+1)
	stack = append(stack, frame{n: t.root, depth: 0})

	shift := int(K) - 32

	// Insert the very first leaf as a child of root.
	first := entries[0]
	insertFreshLeaf(t, &stack, first.Key, first.Val, K)

	prev := first.Key
	for i := 1; i < len(entries); i++ {
		cur := entries[i]
		key := cur.Key

		if key == prev {
			// Same k-mer as previous: append to the current leaf, which is
			// the leaf attached to the deepest node on the stack.
			top := stack[len(stack)-1].n
			top.leaf.val = append(top.leaf.val, cur.Val)
			continue
		}

		// LCP between prev and cur, in bases.
		lcp := uint8(bits.LeadingZeros64(prev^key)>>1 + shift)

		// Pop frames whose START depth is >= lcp; they cannot contain the
		// branch point. The branch point lies inside the frame whose
		// (start_depth <= lcp) and (end_depth > lcp), or exactly between
		// frames if end_depth == lcp.
		for len(stack) > 1 {
			top := stack[len(stack)-1]
			startDepth := top.depth - top.n.k
			if startDepth >= lcp {
				stack = stack[:len(stack)-1]
				continue
			}
			break
		}

		top := stack[len(stack)-1]
		if top.depth > lcp {
			// Branch point lies strictly inside top.n.prefix: split it.
			// top.n keeps the suffix; insert a new parent that holds the
			// shared prefix, then the new branch hangs off the new parent.
			parent := stack[len(stack)-2].n
			startDepth := top.depth - top.n.k
			splitAt := lcp - startDepth // # bases of top.n.prefix to keep above

			oldN := top.n
			oldFirstBase := util.KmerBaseAt(oldN.prefix, oldN.k, 0)

			mid := t.nodes.alloc()
			mid.prefix = util.KmerPrefix(oldN.prefix, oldN.k, splitAt)
			mid.k = splitAt
			parent.children[oldFirstBase] = mid
			// parent.numChildren unchanged: replaced one child with one child.

			oldN.prefix = util.KmerSuffix(oldN.prefix, oldN.k, splitAt)
			oldN.k = oldN.k - splitAt
			mid.children[util.KmerBaseAt(oldN.prefix, oldN.k, 0)] = oldN
			mid.numChildren++

			// Replace top of stack with mid (oldN is no longer on the
			// rightmost path; the new leaf will be).
			stack[len(stack)-1] = frame{n: mid, depth: lcp}
		}

		// Now stack top's end_depth == lcp. Hang a new branch holding the
		// suffix of cur from there.
		insertFreshLeaf(t, &stack, key, cur.Val, K)
		prev = key
	}
}

// insertFreshLeaf creates a new node holding key's suffix below the current
// stack top (whose end depth must be < K) and pushes it onto the stack.
// It also creates the leaf and attaches it to the new node.
func insertFreshLeaf(t *Tree, stack *[]frame, key uint64, v uint32, K uint8) {
	top := (*stack)[len(*stack)-1]
	parent := top.n
	parentDepth := top.depth

	suffix := util.KmerSuffix(key, K, parentDepth)
	suffixK := K - parentDepth

	n := t.nodes.alloc()
	n.prefix = suffix
	n.k = suffixK

	leaf := t.leaves.alloc()
	leaf.key = key
	leaf.val = append(leaf.val, v)
	n.leaf = leaf

	parent.children[util.KmerBaseAt(suffix, suffixK, 0)] = n
	parent.numChildren++

	*stack = append(*stack, frame{n: n, depth: K})
}

type frame struct {
	n     *node
	depth uint8
}

// SearchResult records information of a search result.
type SearchResult struct {
	Kmer      uint64   // searched kmer
	LenPrefix uint8    // length of common prefix between the query and this k-mer
	Values    []uint32 // value of this key
}

var poolSearchResults = &sync.Pool{New: func() interface{} {
	tmp := make([]*SearchResult, 0, 128)
	return &tmp
}}

var poolSearchResult = &sync.Pool{New: func() interface{} {
	return &SearchResult{}
}}

// RecycleSearchResult recycles search results objects.
func (idx *Tree) RecycleSearchResult(sr *[]*SearchResult) {
	for _, r := range *sr {
		poolSearchResult.Put(r)
	}
	*sr = (*sr)[:0]
	poolSearchResults.Put(sr)
}

// Search finds keys that shared prefixes at least p bases.
// We assume the k values of the query k-mer and k-mers in the tree are the same.
// After using the result, do not forget to call RecycleSearchResult().
func (t *Tree) Search(key uint64, p uint8) (*[]*SearchResult, bool) {
	if p < 1 {
		p = 1
	}
	k := t.k
	if p > k {
		p = k
	}
	key0, k0 := key, k
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

	if target == nil {
		return nil, false
	}

	// output all leaves below n
	// results := make([]SearchResult, 0, 8)
	results := poolSearchResults.Get().(*[]*SearchResult)

	var shift int = int(k0 - 32) // pre calculate it, a little bit faster
	recursiveWalk(target, func(key uint64, v []uint32) bool {
		r := poolSearchResult.Get().(*SearchResult)
		r.Kmer = key
		r.LenPrefix = uint8(bits.LeadingZeros64(key0^key)>>1 + shift)
		r.Values = v

		*results = append(*results, r)
		return false
	})

	return results, true
}

// WalkFn is used for walking the tree. Takes a
// key and value, returning true if iteration should
// be terminated.
// type WalkFn func(key uint64, k uint8, v []uint64) bool
type WalkFn func(key uint64, v []uint32) bool

// Walk is used to walk the whole tree.
func (t *Tree) Walk(fn WalkFn) {
	recursiveWalk(t.root, fn)
}

// recursiveWalk is used to do a pre-order walk of a node
// recursively. Returns true if the walk should be aborted.
// The walked k-mers are in lexicographic order,
// so sorting is unneeded.
func recursiveWalk(n *node, fn WalkFn) bool {
	if n.leaf != nil && fn(n.leaf.key, n.leaf.val) {
		return true
	}

	for _, child := range n.children {
		if child != nil && recursiveWalk(child, fn) {
			return true
		}
	}

	return false
}

// WalkGroupFn is called once per maximal group of k-mers that share a prefix
// of at least lenPrefix bases. The keys and vals slices are reused across
// calls; copy them if you need to retain across invocations.
// Return true to abort iteration.
type WalkGroupFn func(keys []uint64, vals [][]uint32, lenPrefix uint8) bool

// WalkGroups visits each maximal subtree whose leaves share a prefix of at
// least p bases, calling fn once per group with all k-mers in that subtree.
// lenPrefix is the guaranteed minimum LCP within the group; specific pairs
// may share more. Singleton groups (one k-mer) are also reported.
// Returns true if iteration was aborted.
func (t *Tree) WalkGroups(p uint8, fn WalkGroupFn) bool {
	if p < 1 {
		p = 1
	}
	if p > t.k {
		p = t.k
	}
	keys := make([]uint64, 0, 1024)
	vals := make([][]uint32, 0, 1024)
	for _, child := range t.root.children {
		if child != nil && recursiveWalkGroups(child, 0, p, &keys, &vals, fn) {
			return true
		}
	}
	return false
}

func recursiveWalkGroups(n *node, parentDepth, p uint8, keys *[]uint64, vals *[][]uint32, fn WalkGroupFn) bool {
	depth := parentDepth + n.k
	if depth >= p {
		*keys = (*keys)[:0]
		*vals = (*vals)[:0]
		collectLeaves(n, keys, vals)
		return fn(*keys, *vals, depth)
	}
	for _, child := range n.children {
		if child != nil && recursiveWalkGroups(child, depth, p, keys, vals, fn) {
			return true
		}
	}
	return false
}

func collectLeaves(n *node, keys *[]uint64, vals *[][]uint32) {
	if n.leaf != nil {
		*keys = append(*keys, n.leaf.key)
		*vals = append(*vals, n.leaf.val)
	}
	for _, child := range n.children {
		if child != nil {
			collectLeaves(child, keys, vals)
		}
	}
}

// WalkPairFn is called for each unordered pair of k-mers (a, b) whose LCP is
// at least p bases; lenPrefix is the exact LCP of this pair.
// Return true to abort iteration.
type WalkPairFn func(keyA, keyB uint64, valsA, valsB []uint32, lenPrefix uint8) bool

// WalkPairs visits every unordered pair of distinct k-mers in the tree
// whose LCP is at least p bases. Returns true if iteration was aborted.
func (t *Tree) WalkPairs(p uint8, fn WalkPairFn) bool {
	shift := int(t.k) - 32
	var stop bool
	t.WalkGroups(p, func(keys []uint64, vals [][]uint32, _ uint8) bool {
		n := len(keys)
		for i := 0; i < n; i++ {
			for j := i + 1; j < n; j++ {
				lcp := uint8(bits.LeadingZeros64(keys[i]^keys[j])>>1 + shift)
				if fn(keys[i], keys[j], vals[i], vals[j], lcp) {
					stop = true
					return true
				}
			}
		}
		return false
	})
	return stop
}

// KmerPrefix returns the first n bases. n needs to be > 0.
// The length of the prefix is n.
func KmerPrefix(code uint64, k uint8, n uint8) uint64 {
	return code >> ((k - n) << 1)
}
