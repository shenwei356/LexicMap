// Copyright © 2026-2026 Wei Shen <shenwei356@gmail.com>
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
package rangeindex

import (
	"math/rand"
	"testing"
	"time"
)

func TestBasicQuery(t *testing.T) {
	ri := NewRangeIndex()
	defer ri.Release()

	ri.Add(100, 1)
	ri.Add(105, 2)
	ri.Add(110, 3)

	results := ri.Query(100, 110)

	if len(results) != 3 {
		t.Fatalf("expected 3 results, got %d", len(results))
	}

	for _, v := range results {
		loc, val := Unpack(v)
		if loc < 100 || loc > 110 {
			t.Fatalf("unexpected loc: %d", loc)
		}
		if val < 1 || val > 3 {
			t.Fatalf("unexpected value: %d", val)
		}
	}
}

func TestDuplicateLocations(t *testing.T) {
	ri := NewRangeIndex()
	defer ri.Release()

	ri.Add(100, 1)
	ri.Add(100, 2)
	ri.Add(100, 3)

	results := ri.Query(100, 100)

	if len(results) != 3 {
		t.Fatalf("expected 3 results, got %d", len(results))
	}
}

func TestEmptyIndex(t *testing.T) {
	ri := NewRangeIndex()
	defer ri.Release()

	results := ri.Query(0, 100)

	if len(results) != 0 {
		t.Fatalf("expected empty result")
	}
}

func TestOutOfRange(t *testing.T) {
	ri := NewRangeIndex()
	defer ri.Release()

	ri.Add(100, 1)
	ri.Add(200, 2)

	results := ri.Query(300, 400)

	if len(results) != 0 {
		t.Fatalf("expected empty result")
	}
}

func TestUnpack(t *testing.T) {
	loc := uint32(12345)
	val := uint32(67890)

	entry := (uint64(loc) << 32) | uint64(val)
	l, v := Unpack(entry)

	if l != loc || v != val {
		t.Fatalf("unpack failed: got (%d,%d)", l, v)
	}
}

func TestMultipleQueries(t *testing.T) {
	ri := NewRangeIndex()
	defer ri.Release()

	for i := uint32(0); i < 1000; i++ {
		ri.Add(i, i+1000)
	}

	for i := uint32(100); i < 200; i++ {
		results := ri.Query(i, i)
		if len(results) != 1 {
			t.Fatalf("expected 1 result at %d", i)
		}
	}
}

func TestRandomDataCorrectness(t *testing.T) {
	r := rand.New(rand.NewSource(time.Now().UnixNano()))

	ri := NewRangeIndex()
	defer ri.Release()

	var reference []uint64

	for i := 0; i < 10000; i++ {
		loc := uint32(r.Intn(10000))
		val := uint32(i)
		ri.Add(loc, val)
		reference = append(reference, (uint64(loc)<<32)|uint64(val))
	}

	left := uint32(2000)
	right := uint32(4000)

	results := ri.Query(left, right)

	// brute force reference check
	var expected []uint64
	for _, v := range reference {
		loc := uint32(v >> 32)
		if loc >= left && loc <= right {
			expected = append(expected, v)
		}
	}

	if len(results) != len(expected) {
		t.Fatalf("mismatch: expected %d got %d", len(expected), len(results))
	}
}

func TestPoolReuse(t *testing.T) {
	ri := NewRangeIndex()

	ri.Add(1, 1)
	ri.Release()

	ri2 := NewRangeIndex()
	defer ri2.Release()

	if len(ri2.data) != 0 {
		t.Fatalf("expected empty after reuse")
	}
}
