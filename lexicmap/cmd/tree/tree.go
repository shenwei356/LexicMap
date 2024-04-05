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

package tree

import (
	"math/bits"
	"sync"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
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
}

var poolLeafNode = &sync.Pool{New: func() interface{} {
	return &leafNode{val: make([]uint32, 0, 1)}
}}

var poolNode = &sync.Pool{New: func() interface{} {
	return &node{}
}}

var poolTree = &sync.Pool{New: func() interface{} {
	return &Tree{}
}}

// NewTree implements a radix tree for k-mer extractly or prefix querying.
func NewTree(k uint8) *Tree {
	t := poolTree.Get().(*Tree)
	t.k = k
	n := poolNode.Get().(*node)
	n.Reset()
	t.root = n
	return t
}

// RecycleTree recycles the tree object.
func RecycleTree(t *Tree) {
	recursiveRecycle(t.root)
	poolTree.Put(t)
}

// recursiveRecycle recycle all nodes, including leaf nodes.
func recursiveRecycle(n *node) {
	if n.leaf != nil {
		poolLeafNode.Put(n.leaf)
	}

	for _, child := range n.children {
		if child == nil {
			continue
		}
		recursiveRecycle(child)
	}

	poolNode.Put(n)
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
			leaf := poolLeafNode.Get().(*leafNode)
			leaf.Reset()
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
			n = poolNode.Get().(*node)
			n.Reset()
			n.prefix = search
			n.k = k
			leaf := poolLeafNode.Get().(*leafNode)
			leaf.Reset()
			leaf.key = key0
			leaf.val = append(leaf.val, v)
			n.leaf = leaf

			parent.children[firstBase] = n
			parent.numChildren++

			return false
		}

		// has a child -- exists a path

		// Determine longest prefix of the search key on match
		// commonPrefix := KmerLongestPrefix(search, n.prefix, k, n.k)
		// because k >= n.k
		commonPrefix := util.MustKmerLongestPrefix(search, n.prefix, k, n.k)
		// the new key is longer than key of n, continue to search. len(prefix) = len(n)
		if commonPrefix == n.k {
			search = util.KmerSuffix(search, k, commonPrefix) // left bases
			k = k - commonPrefix                              // need to update it
			continue
		}

		// the new key and the key of node n share a prefix, len(prefix) < len(n)
		// Split the node n
		child := poolNode.Get().(*node)
		child.Reset()
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
		leaf := poolLeafNode.Get().(*leafNode)
		leaf.Reset()
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
		n = poolNode.Get().(*node)
		n.Reset()
		n.prefix = search
		n.k = k
		n.leaf = leaf
		child.children[util.KmerBaseAt(search, k, 0)] = n

		child.numChildren++
		return false
	}
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
	*results = (*results)[:0]

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
// key and value, returning if iteration should
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

// KmerPrefix returns the first n bases. n needs to be > 0.
// The length of the prefix is n.
func KmerPrefix(code uint64, k uint8, n uint8) uint64 {
	return code >> ((k - n) << 1)
}
