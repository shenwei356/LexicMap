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

	"github.com/shenwei356/LexicMap/lexicmap/util"
	"github.com/shenwei356/lexichash"
)

// leafNode is used to represent a value.
type leafNode struct {
	key uint64   // ALL the bases in the node, the k-mer
	val []uint64 // yes, multiple values
}

// node represents a node in the tree, it might be the root, inner or leaf node.
type node struct {
	prefix uint64 // prefix of the current node
	k      uint8  // bases length of the prefix

	numChildren uint8
	children    [4]*node // just use an array

	leaf *leafNode // optional
}

// Tree is a radix tree for storing bit-packed k-mer information.
type Tree struct {
	k    uint8 // use a global K
	root *node // root node

	numNodes     int // the number of nodes, including leaf nodes
	numLeafNodes int // the number of leaf nodes

	poolPath *sync.Pool // pool for path
}

// Tree implements a radix tree for k-mer extractly or prefix querying.
func New(k uint8) *Tree {
	t := &Tree{k: k, root: &node{}}
	t.poolPath = &sync.Pool{New: func() interface{} {
		tmp := make([]string, k)
		return &tmp
	}}
	return t
}

// K returns the K value of k-mers.
func (t *Tree) K() int {
	return int(t.k)
}

// NumNodes returns the number of nodes, including leaf nodes.
func (t *Tree) NumNodes() int {
	return t.numNodes
}

// NumLeafNodes returns the number of leaf nodes.
func (t *Tree) NumLeafNodes() int {
	return t.numLeafNodes
}

// Insert is used to add a new entry or update
// an existing entry. Returns true if an existing record is updated.
func (t *Tree) Insert(key uint64, v uint64) bool {
	key0 := key // will save it into the leaf node
	k := t.k

	var parent *node
	n := t.root
	search := key // current key
	for {
		// Handle key exhaustion
		if k == 0 {
			if n.leaf != nil {
				if n.leaf.val == nil {
					n.leaf.val = []uint64{v}
				} else {
					n.leaf.val = append(n.leaf.val, v)
				}
				return true
			}

			// n is not a leaf node, that means
			// the current key is a prefix of some other keys.
			n.leaf = &leafNode{
				key: key0,
				val: []uint64{v},
			}
			t.numLeafNodes++

			return false
		}

		// Look for the child
		parent = n
		firstBase := util.KmerBaseAt(search, k, 0)
		n = n.children[firstBase]

		// No child, create one
		if n == nil {
			parent.children[firstBase] = &node{
				leaf: &leafNode{
					key: key0,
					val: []uint64{v},
				},
				prefix: search,
				k:      k,
			}
			parent.numChildren++

			t.numNodes++
			t.numLeafNodes++
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
		child := &node{
			// o---<=8, here the prefix of one of the 8 is ---,
			prefix: util.KmerPrefix(search, k, commonPrefix),
			k:      commonPrefix,
		}
		t.numNodes++
		parent.children[firstBase] = child // change from n to c

		// child points to n now
		child.children[util.KmerBaseAt(n.prefix, n.k, commonPrefix)] = n
		child.numChildren++
		n.prefix = util.KmerSuffix(n.prefix, n.k, commonPrefix)
		n.k = n.k - commonPrefix

		// Create a new leaf node for the new key
		leaf := &leafNode{
			key: key0,
			val: []uint64{v},
		}
		t.numLeafNodes++

		// the new key is a prefix of the old n, add the leaf node to this node. len(new) = len(prefix)
		search = util.KmerSuffix(search, k, commonPrefix)
		k = k - commonPrefix
		if k == 0 {
			child.leaf = leaf
			return false
		}

		// the new key and the key of node n share a prefix shorter than both of them
		// Create a new child node for the node
		child.children[util.KmerBaseAt(search, k, 0)] = &node{
			leaf:   leaf,
			prefix: search,
			k:      k,
		}
		child.numChildren++
		t.numNodes++
		return false
	}
}

// Get is used to lookup a specific key, returning
// the value and if it was found.
func (t *Tree) Get(key uint64) ([]uint64, bool) {
	n := t.root
	search := key
	k := t.k
	for {
		// Check for key exhaution
		if k == 0 {
			if n.leaf != nil {
				return n.leaf.val, true
			}
			break
		}

		// Look for a child
		n = n.children[util.KmerBaseAt(search, k, 0)]
		if n == nil { // not found
			break
		}

		// Consume the search prefix
		if util.KmerHasPrefix(search, n.prefix, k, n.k) {
			search = util.KmerSuffix(search, k, n.k)
			k = k - n.k
		} else {
			break
		}
	}
	return nil, false
}

// RecyclePathResult recycles the node list.
func (t *Tree) RecyclePathResult(nodes *[]string) {
	t.poolPath.Put(nodes)
}

// Path returns the path of a key, i.e., the nodes list.
// and the number of visited/matched bases.
// Do not forgot to call RecyclePathResult() after using the results.
func (t *Tree) Path(key uint64, minPrefix uint8) (*[]string, uint8) {
	n := t.root
	search := key
	k := t.k
	//	nodes := make([]string, 0, k)
	nodes := t.poolPath.Get().(*[]string)
	*nodes = (*nodes)[:0]
	var matched uint8
	for {
		// Check for key exhaution
		if k == 0 {
			if n.leaf != nil {
				// nodes = append(nodes, string(lexichash.MustDecode(n.leaf.key, t.k)))
				return nodes, matched
			}
			break
		}

		// Look for a child
		n = n.children[util.KmerBaseAt(search, k, 0)]
		if n == nil { // not found
			break
		}

		// Consume the search prefix
		if util.KmerHasPrefix(search, n.prefix, k, n.k) {
			matched += n.k
			*nodes = append(*nodes, string(lexichash.MustDecode(n.prefix, n.k)))
			search = util.KmerSuffix(search, k, n.k)
			k = k - n.k
		} else {
			// also check the prefix, because the prefix of some nodes
			// might be very long. Only checking prefix will ignore theses.
			// For example, the strings below shared 4 bases,
			// the third node would be this case.
			//   A C AGCT
			//   A C AGGC
			// matched += KmerLongestPrefix(search, n.prefix, k, n.k)
			// because k >= n.k
			matched += util.MustKmerLongestPrefix(search, n.prefix, k, n.k)
			if matched >= minPrefix {
				*nodes = append(*nodes, string(lexichash.MustDecode(n.prefix, n.k)))
				break
			}

			break
		}
	}
	return nodes, matched
}

// SearchResult records information of a search result.
type SearchResult struct {
	Kmer      uint64   // searched kmer
	LenPrefix uint8    // length of common prefix between the query and this k-mer
	Values    []uint64 // value of this key
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

// Search finds keys that shared prefixes at least m bases.
// We assume the k values of the query k-mer and k-mers in the tree are the same.
// After using the result, do not forget to call RecycleSearchResult().
func (t *Tree) Search(key uint64, m uint8) (*[]*SearchResult, bool) {
	if m < 1 {
		m = 1
	}
	k := t.k
	if m > k {
		m = k
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
			if lenPrefix >= m {
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

			atleast = m - lenPrefix // just check if prefixes of m - lenPrefix bases are equal
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
	recursiveWalk(target, func(key uint64, v []uint64) bool {
		r := poolSearchResult.Get().(*SearchResult)
		r.Kmer = key
		r.LenPrefix = uint8(bits.LeadingZeros64(key0^key)>>1 + shift)
		r.Values = v

		*results = append(*results, r)
		return false
	})

	return results, true
}

// LongestPrefix is like Get, but instead of an
// exact match, it will return the longest prefix match.
func (t *Tree) LongestPrefix(key uint64) (uint64, []uint64, bool) {
	var last *leafNode
	n := t.root
	search := key
	k := t.k
	for {
		// Look for a leaf node
		if n.leaf != nil {
			last = n.leaf
		}

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
		if util.KmerHasPrefix(search, n.prefix, k, n.k) {
			search = util.KmerSuffix(search, k, n.k)
			k = k - n.k
		} else {
			break
		}
	}
	if last != nil {
		return last.key, last.val, true
	}
	return 0, nil, false
}

// WalkFn is used for walking the tree. Takes a
// key and value, returning if iteration should
// be terminated.
// type WalkFn func(key uint64, k uint8, v []uint64) bool
type WalkFn func(key uint64, v []uint64) bool

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
