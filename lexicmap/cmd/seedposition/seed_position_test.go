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

package seedposition

import (
	"math/rand"
	"os"
	"testing"
)

func TestSeedPositions(t *testing.T) {
	tests := [][]uint32{
		{},
		{1},
		{1, 15},
		{1, 15, 300},
		{1, 15, 300, 301},
		{1, 15, 300, 301, 2500},
		{1, 15, 300, 301, 2500, 3100},
		{1, 15, 300, 301, 2500, 3100, 3111},
		{1, 15, 300, 301, 2500, 3100, 3111, 5000},
		{1, 15, 300, 301, 2500, 3100, 3111, 5000, 10000},
	}

	file := "test.bin"

	// ---------------------------------------

	wtr, err := NewWriter(file, 0)
	if err != nil {
		t.Error(err)
		return
	}
	for i, test := range tests {
		err = wtr.Write(test)
		if err != nil {
			t.Errorf("write #%d data: %s", i+1, err)
			return
		}
	}
	err = wtr.Close()
	if err != nil {
		t.Error(err)
		return
	}

	idxs := make([]int, len(tests))
	for i := range tests {
		idxs[i] = i
	}
	rand.Shuffle(len(tests), func(i, j int) { idxs[i], idxs[j] = idxs[j], idxs[i] })

	// ---------------------------------------

	rdr, err := NewReader(file)
	if err != nil {
		t.Error(err)
		return
	}

	locs := make([]uint32, 64)
	var test []uint32
	var v uint32
	var j int

	for _, i := range idxs {
		test = tests[i]
		err = rdr.SeedPositions(i, &locs)
		if err != nil {
			t.Errorf("read #%d data: %s", i, err)
			return
		}

		if len(locs) != len(test) {
			t.Errorf("[#%d] unequal of position numbers, expected: %d, returned %d",
				i, len(test), len(locs))
			return
		}

		for j, v = range locs {
			if v != test[j] {
				t.Errorf("[#%d] unequal of positions, expected: %d, returned %d", i, test[j], v)
				return
			}
		}
	}

	// clean up

	err = os.RemoveAll(file)
	if err != nil {
		t.Error(err)
		return
	}

	err = os.RemoveAll(file + PositionsIndexFileExt)
	if err != nil {
		t.Error(err)
		return
	}
}
