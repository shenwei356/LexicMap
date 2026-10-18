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

package cmd

import (
	"encoding/json"
	"os"
	"slices"
	"testing"
)

func TestSampleQueryFragmentMatchesPointerBasedSampler(t *testing.T) {
	// Generated with the sampler at ffc3b1f before changing buffer ownership.
	data, err := os.ReadFile("testdata/query-seeds.json")
	if err != nil {
		t.Fatal(err)
	}
	var cases []struct {
		Fragment string
		Seeds    []uint64
	}
	if err := json.Unmarshal(data, &cases); err != nil {
		t.Fatal(err)
	}
	var buffer []uint64
	for i, tt := range cases {
		buffer, err = sampleQueryFragment([]byte(tt.Fragment), buffer)
		if err != nil || !slices.Equal(buffer, tt.Seeds) {
			t.Fatalf("case %d: sampled k-mers changed, err=%v", i, err)
		}
	}
}

func TestRecycleQuerySeedsBoundsUnusedSlots(t *testing.T) {
	seeds := [][]uint64{make([]uint64, 1, 8), make([]uint64, 1, thresholdNSubsLong+1)}
	backing := seeds
	bufferStart := &seeds[0][0]
	seeds = seeds[:1]
	recycleQuerySeeds(&seeds)
	if len(seeds) != 0 || len(backing[0]) != 0 || cap(backing[0]) != 8 || backing[1] != nil {
		t.Fatal("recycling lost a small buffer or retained an oversized unused buffer")
	}
	if &backing[0][:1][0] != bufferStart {
		t.Fatal("recycling replaced the reusable sampling buffer")
	}
}
