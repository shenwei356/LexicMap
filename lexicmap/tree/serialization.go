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
	"bytes"
	"encoding/binary"
	"errors"
	"io"
	"sync"

	"github.com/shenwei356/LexicMap/lexicmap/util"
	"github.com/shenwei356/xopen"
)

var be = binary.BigEndian

// Magic number for checking file format
var Magic = [8]byte{'k', 'm', 'e', 'r', 't', 'r', 'e', 'e'}

// MainVersion is use for checking compatibility
var MainVersion uint8 = 0

// MinorVersion is less important
var MinorVersion uint8 = 1

// ErrInvalidFileFormat means invalid file format.
var ErrInvalidFileFormat = errors.New("k-mer tree: invalid binary format")

// ErrBrokenFile means the file is not complete.
var ErrBrokenFile = errors.New("k-mer tree: broken file")

// ErrKOverflow means K < 1 or K > 32.
var ErrKOverflow = errors.New("k-mer tree: k-mer size [1, 32] overflow")

// ErrVersionMismatch means version mismatch between files and program
var ErrVersionMismatch = errors.New("k-mer tree: version mismatch")

// NewFromFile creates a Tree from a file.
func NewFromFile(file string) (*Tree, error) {
	fh, err := xopen.Ropen(file)
	if err != nil {
		return nil, err
	}
	defer fh.Close()

	return Read(fh)
}

// WriteToFile writes a tree to a file, optional with file extension of .gz, .xz, .zst, .bz2.
// But using uncompressed format is good enough, and faster and memory efficient
// when writing a lot of trees.
func (t *Tree) WriteToFile(file string) (int, error) {
	outfh, err := xopen.Wopen(file)
	if err != nil {
		return 0, err
	}
	defer outfh.Close()

	return t.Write(outfh)
}

var poolBuf = &sync.Pool{New: func() interface{} {
	return &bytes.Buffer{}
}}

// Write writes the tree to a writer.
// Currently, it just writes all the leaves, i.e., k-mers and values.
//
// Header (24 bytes):
//
//	Magic number, 8 bytes, kmertree
//	Main and minor versions, 2 bytes
//	K, 1 byte
//	Blank, 5 bytes
//	Number of keys: 8 bytes
//
// In practice, all-kmers of a tree is very similar, sharing common prefixes.
// And the results of Walk() are already sorted.
// So we use varints for saving two k-mers and numbers of values.
// It would save a lot of space and it is also easy to read and write.
//
// Data: k-mer-values pairs.
//
//	Control byte for 2 k-mers, 1 byte
//	Delta values of the 2 k-mers, 2-16 bytes
//	Control byte for numbers of values, 1 byte
//	Numbers of values of the 2 k-mers, 2-16 bytes, 2 bytes for most cases.
//	Values of the 2 k-mers, 8*n bytes, 16 bytes for most cases.
func (t *Tree) Write(w io.Writer) (int, error) {
	var N int // the number of bytes.
	var err error

	// 8-byte magic number
	err = binary.Write(w, be, Magic)
	if err != nil {
		return N, err
	}
	N += 8

	// 8-byte meta info
	err = binary.Write(w, be, [8]uint8{MainVersion, MinorVersion, t.k})
	if err != nil {
		return N, err
	}
	N += 8

	// 8-byte the number of k-mers
	err = binary.Write(w, be, uint64(t.numLeafNodes))
	if err != nil {
		return N, err
	}
	N += 8

	var hasPrev bool
	var preKey uint64
	var preVal []uint64
	var offset uint64
	var ctrlByteKey, ctrlByteVal byte
	var nBytesKey, nBytesVal, n int
	bufVar := make([]byte, 16) // needs at most 8+8=16
	buf := make([]byte, 36)    // needs at most 1+16+1+16=34
	buf8 := make([]byte, 8)    // for writing uint8
	bufVals := poolBuf.Get().(*bytes.Buffer)
	defer poolBuf.Put(bufVals)
	var _v uint64

	t.Walk(func(key uint64, v []uint64) bool {
		if !hasPrev { // write it later
			preKey = key
			preVal = v
			hasPrev = true

			return false
		}

		// 2 k-mers and numbers of values

		// only save key2 - key1, which is small so it could be saved in few bytes
		ctrlByteKey, nBytesKey = util.PutUint64s(bufVar, preKey-offset, key-preKey)
		buf[0] = ctrlByteKey
		copy(buf[1:nBytesKey+1], bufVar[:nBytesKey])
		n = nBytesKey + 1

		// save lengths of values
		ctrlByteVal, nBytesVal = util.PutUint64s(bufVar, uint64(len(preVal)), uint64(len(v)))
		buf[n] = ctrlByteVal
		copy(buf[n+1:n+nBytesVal+1], bufVar[:nBytesVal])
		n += nBytesVal + 1

		_, err = w.Write(buf[:n])
		if err != nil {
			return true // stop the walk
		}
		N += n

		// values

		bufVals.Reset()
		for _, _v = range preVal {
			be.PutUint64(buf8, _v)
			bufVals.Write(buf8)
		}
		for _, _v = range v {
			be.PutUint64(buf8, _v)
			bufVals.Write(buf8)
		}

		_, err = w.Write(bufVals.Bytes())
		if err != nil {
			return true // stop the walk
		}
		N += bufVals.Len()

		// update

		offset = key
		hasPrev = false

		return false
	})

	if err != nil {
		return N, err
	}

	if hasPrev { // the last single one
		// save key2 - key1 and length of the values
		ctrlByteKey, nBytesKey = util.PutUint64s(bufVar, preKey-offset, uint64(len(preVal)))
		buf[0] = ctrlByteKey // set the first bit to one, as a flag of last record

		copy(buf[1:nBytesKey+1], bufVar[:nBytesKey])
		n = nBytesKey + 1

		_, err = w.Write(buf[:n])
		if err != nil {
			return N, err
		}
		N += n

		bufVals.Reset()
		for _, _v = range preVal {
			be.PutUint64(buf8, _v)
			bufVals.Write(buf8)
		}

		_, err = w.Write(bufVals.Bytes())
		if err != nil {
			return N, err
		}
		N += bufVals.Len()
	}

	return N, nil
}

// Read reads a tree from an io.Reader.
func Read(r io.Reader) (*Tree, error) {
	buf := make([]byte, 64)

	var err error
	var n int

	// check the magic number
	n, err = io.ReadFull(r, buf[:8])
	if err != nil {
		return nil, err
	}
	if n < 8 {
		return nil, ErrBrokenFile
	}
	same := true
	for i := 0; i < 8; i++ {
		if Magic[i] != buf[i] {
			same = false
			break
		}
	}
	if !same {
		return nil, ErrInvalidFileFormat
	}

	// read metadata
	n, err = io.ReadFull(r, buf[:8])
	if err != nil {
		return nil, err
	}
	if n < 8 {
		return nil, ErrBrokenFile
	}
	// check compatibility
	if MainVersion != buf[0] {
		return nil, ErrVersionMismatch
	}
	// check k-mer size
	if buf[2] > 32 {
		return nil, ErrKOverflow
	}

	t := &Tree{root: &node{}}
	t.k = buf[2]

	// the number of k-mers
	_, err = io.ReadFull(r, buf[:8])
	if err != nil {
		return nil, err
	}

	nKmers := be.Uint64(buf[:8])
	nReadRounds := nKmers / 2
	hasSingle := nKmers%2 > 0
	var offset uint64
	var ctrlByte byte
	var bytes [2]uint8
	var nBytes int
	var nReaded, nDecoded int
	var decodedVals [2]uint64
	var kmer1, kmer2 uint64
	var lenVal1, lenVal2 uint64
	var j uint64
	buf8 := make([]byte, 8)
	for i := 0; i < int(nReadRounds); i++ {
		// ------------------ k-mers -------------------

		// read the control byte
		_, err = io.ReadFull(r, buf[:1])
		if err != nil {
			return nil, err
		}
		ctrlByte = buf[0]

		// parse the control byte
		bytes = util.CtrlByte2ByteLengths[ctrlByte]
		nBytes = int(bytes[0] + bytes[1])

		// read encoded bytes
		nReaded, err = io.ReadFull(r, buf[:nBytes])
		if err != nil {
			return nil, err
		}
		if nReaded < nBytes {
			return nil, ErrBrokenFile
		}

		decodedVals, nDecoded = util.Uint64s(ctrlByte, buf[:nBytes])
		if nDecoded == 0 {
			return nil, ErrBrokenFile
		}

		kmer1 = decodedVals[0] + offset
		kmer2 = kmer1 + decodedVals[1]

		offset = kmer2

		// ------------------ lengths of values -------------------

		// read the control byte
		_, err = io.ReadFull(r, buf[:1])
		if err != nil {
			return nil, err
		}
		ctrlByte = buf[0]

		// parse the control byte
		bytes = util.CtrlByte2ByteLengths[ctrlByte]
		nBytes = int(bytes[0] + bytes[1])

		// read encoded bytes
		nReaded, err = io.ReadFull(r, buf[:nBytes])
		if err != nil {
			return nil, err
		}
		if nReaded < nBytes {
			return nil, ErrBrokenFile
		}

		decodedVals, nDecoded = util.Uint64s(ctrlByte, buf[:nBytes])
		if nDecoded == 0 {
			return nil, ErrBrokenFile
		}

		lenVal1 = decodedVals[0]
		lenVal2 = decodedVals[1]

		// ------------------ values -------------------

		v1 := make([]uint64, lenVal1)
		for j = 0; j < lenVal1; j++ {
			nReaded, err = io.ReadFull(r, buf8)
			if err != nil {
				return nil, err
			}
			if nReaded < 8 {
				return nil, ErrBrokenFile
			}

			v1[j] = be.Uint64(buf8)
		}

		v2 := make([]uint64, lenVal2)
		for j = 0; j < lenVal2; j++ {
			nReaded, err = io.ReadFull(r, buf8)
			if err != nil {
				return nil, err
			}
			if nReaded < 8 {
				return nil, ErrBrokenFile
			}

			v2[j] = be.Uint64(buf8)
		}

		// ------------------ insert into the tree -------------------

		t.insertKeyVals(kmer1, v1)
		t.insertKeyVals(kmer2, v2)
	}

	if hasSingle { // read the last one
		// read the control byte
		_, err = io.ReadFull(r, buf[:1])
		if err != nil {
			return nil, err
		}
		ctrlByte = buf[0]

		// parse the control byte
		bytes = util.CtrlByte2ByteLengths[ctrlByte]
		nBytes = int(bytes[0] + bytes[1])

		// read encoded bytes
		nReaded, err = io.ReadFull(r, buf[:nBytes])
		if err != nil {
			return nil, err
		}
		if nReaded < nBytes {
			return nil, ErrBrokenFile
		}

		decodedVals, nDecoded = util.Uint64s(ctrlByte, buf[:nBytes])
		if nDecoded == 0 {
			return nil, ErrBrokenFile
		}

		kmer1 = decodedVals[0] + offset
		lenVal1 = decodedVals[1]

		v1 := make([]uint64, lenVal1)
		for j = 0; j < lenVal1; j++ {
			nReaded, err = io.ReadFull(r, buf8)
			if err != nil {
				return nil, err
			}
			if nReaded < 8 {
				return nil, ErrBrokenFile
			}

			v1[j] = be.Uint64(buf8)
		}

		t.insertKeyVals(kmer1, v1)
	}

	t.poolPath = &sync.Pool{New: func() interface{} {
		tmp := make([]string, t.k)
		return &tmp
	}}

	return t, nil
}

// Different from Insert(), this method adds the values,
// rather than a single value.
func (t *Tree) insertKeyVals(key uint64, v []uint64) bool {
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
					n.leaf.val = v
				}
				return true
			}

			// n is not a leaf node, that means
			// the current key is a prefix of some other keys.
			n.leaf = &leafNode{
				key: key0,
				val: v,
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
					val: v,
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
			val: v,
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
