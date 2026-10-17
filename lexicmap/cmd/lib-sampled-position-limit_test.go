package cmd

import (
	"fmt"
	"reflect"
	"testing"
)

// referenceAddSampledPosition is the uncapped sketch builder. Keep its original
// packed list metadata so the allocation benchmark includes the count's cost.
func referenceAddSampledPosition(kmers map[uint64]uint32, repeated map[uint64]uint64, positions *[]repeatedKmerPosition, kmer uint64, pos uint32) {
	if _, ok := kmers[kmer]; !ok {
		kmers[kmer] = pos + 1
		return
	}
	list := repeated[kmer]
	head, tail := uint32(list>>32), uint32(list)
	*positions = append(*positions, repeatedKmerPosition{position: pos})
	index := uint32(len(*positions))
	if head == 0 {
		head = index
	} else {
		(*positions)[tail-1].next = index
	}
	repeated[kmer] = uint64(head)<<32 | uint64(index)
}

func sampledPositionPrefix(first, head uint32, positions []repeatedKmerPosition) []uint32 {
	result := []uint32{first - 1}
	for head != 0 && len(result) < 100 { // the original consumer's limit
		p := positions[head-1]
		result = append(result, p.position)
		head = p.next
	}
	return result
}

func TestSampledPositionLimitPreservesOriginalPrefix(t *testing.T) {
	for _, n := range []int{1, 2, 99, 100, 101, 1000} {
		t.Run(fmt.Sprint(n), func(t *testing.T) {
			kmers := make(map[uint64]uint32)
			repeated := make(map[uint64]repeatedKmerList)
			oldKmers := make(map[uint64]uint32)
			oldRepeated := make(map[uint64]uint64)
			var positions, oldPositions []repeatedKmerPosition
			// Interleave keys and both encoded strands. Each key's retained
			// prefix must match even when arena nodes are not contiguous.
			for i := range n {
				for key := uint64(1); key <= 3; key++ {
					pos := uint32((i*31+int(key))<<1 | i&1)
					addSampledKmerPosition(&kmers, &repeated, &positions, key, pos)
					referenceAddSampledPosition(oldKmers, oldRepeated, &oldPositions, key, pos)
				}
			}
			if len(positions) != 3*(min(n, 100)-1) {
				t.Fatalf("retained %d repeated positions", len(positions))
			}
			for key := uint64(1); key <= 3; key++ {
				got := sampledPositionPrefix(kmers[key], repeated[key].head, positions)
				want := sampledPositionPrefix(oldKmers[key], uint32(oldRepeated[key]>>32), oldPositions)
				if !reflect.DeepEqual(got, want) {
					t.Fatalf("key %d: got %v, want %v", key, got, want)
				}
				if n > 1 && repeated[key].count != uint32(min(n, 100)) {
					t.Fatalf("key %d: wrong retained count %d", key, repeated[key].count)
				}
				if tail := repeated[key].tail; tail != 0 && positions[tail-1].next != 0 {
					t.Fatal("retained list does not terminate")
				}
			}
		})
	}
}

func BenchmarkSampledPositionLimit(b *testing.B) {
	for _, copies := range []int{2, 1000} {
		for _, limited := range []bool{false, true} {
			b.Run(fmt.Sprintf("copies=%d/limited=%t", copies, limited), func(b *testing.B) {
				b.ReportAllocs()
				for range b.N {
					kmers := make(map[uint64]uint32, 1000)
					var positions []repeatedKmerPosition
					if limited {
						repeated := make(map[uint64]repeatedKmerList, 1000)
						for i := range copies {
							for key := uint64(0); key < 1000; key++ {
								addSampledKmerPosition(&kmers, &repeated, &positions, key, uint32(i))
							}
						}
					} else {
						repeated := make(map[uint64]uint64, 1000)
						for i := range copies {
							for key := uint64(0); key < 1000; key++ {
								referenceAddSampledPosition(kmers, repeated, &positions, key, uint32(i))
							}
						}
					}
				}
			})
		}
	}
}
