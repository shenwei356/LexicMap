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
	posting := &seeds[0][0]
	seeds = seeds[:1]
	recycleQuerySeeds(&seeds)
	if len(seeds) != 0 || len(backing[0]) != 0 || cap(backing[0]) != 8 || backing[1] != nil {
		t.Fatal("recycling lost a small buffer or retained an oversized unused buffer")
	}
	if &backing[0][:1][0] != posting {
		t.Fatal("recycling replaced the reusable sampling buffer")
	}
}
