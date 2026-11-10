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
	"bytes"
	"fmt"
	"math"
	"math/rand"
	"os"
	"path/filepath"
	"reflect"
	"strconv"
	"strings"
	"testing"
)

// TestGenomeScreenRepeatedSeeds checks that splitting long posting lists into
// borrowed batches preserves per-mask maxima, cutoff ties, and query isolation.
// All references are generated here; no local database is needed.
func TestGenomeScreenRepeatedSeeds(t *testing.T) {
	dir := t.TempDir()
	rng := rand.New(rand.NewSource(1))
	motif := make([]byte, 512)
	for i := range motif {
		motif[i] = "ACGT"[rng.Intn(4)]
	}
	// Each captured k-mer has more occurrences than a streaming decoder batch.
	sequence := bytes.Repeat(motif, 300)
	var refs []string
	for i := range 4 {
		file := filepath.Join(dir, fmt.Sprintf("ref%d.fa", i))
		if err := os.WriteFile(file, append([]byte(">seq\n"), sequence...), 0600); err != nil {
			t.Fatal(err)
		}
		refs = append(refs, file)
	}
	query := filepath.Join(dir, "partial.fa")
	if err := os.WriteFile(query, append([]byte(">seq\n"), sequence[:4096]...), 0600); err != nil {
		t.Fatal(err)
	}
	db := filepath.Join(dir, "index.lmi")
	args := []string{"index", "--quiet", "-j", "2", "-b", "2", "-c", "2", "-m", "1024", "-k", "16", "--partitions", "4", "-O", db}
	runSeedIndexCommand(t, append(args, refs...)...)
	for _, masks := range []int{0, 256} {
		t.Run(fmt.Sprintf("masks=%d", masks), func(t *testing.T) {
			for _, memory := range []bool{false, true} {
				t.Run(fmt.Sprintf("rolling-memory=%t", memory), func(t *testing.T) {
					opt := DefaultIndexSearchingOptions
					opt.NumCPUs, opt.MaxSeedSearchingConcurrency = 2, 2
					opt.MinPrefix, opt.TopN, opt.SearchMaskCount = 12, 1, masks
					opt.InMemorySearch = memory
					idx, err := NewIndexSearcher(db, &opt)
					if err != nil {
						t.Fatal(err)
					}
					defer idx.Close()
					gr := NewGenomeReader(16, nil)
					q, err := gr.Read(query, true, false)
					if err != nil {
						t.Fatal(err)
					}
					defer gr.Recycle(q)
					subset := make(map[int]struct{})
					for mask := 0; mask < 1024; mask += 3 {
						subset[mask] = struct{}{}
					}
					for _, selection := range []map[int]struct{}{nil, subset, nil} {
						aIDs, a, err := idx.GSearchScreen(q, 2, selection)
						if err != nil {
							t.Fatal(err)
						}
						bIDs, b, err := idx.gsearchScreen(q, 2, selection, false)
						if err != nil {
							t.Fatal(err)
						}
						if a == nil || b == nil || len(*a) != 4 || len(*b) != 4 || !reflect.DeepEqual(aIDs, bIDs) {
							t.Fatal("rolling screening changed candidates or cutoff ties")
						}
						for i, r := range *a {
							s := (*b)[i]
							if r.SumPrefix != s.SumPrefix || r.HitMasks != s.HitMasks || !reflect.DeepEqual(r.BatchGenomeIndex, s.BatchGenomeIndex) {
								t.Fatal("rolling screening changed scores or candidate order")
							}
							if len(s.lastMatches) == 0 || len(s.LongestMatches) != 0 {
								t.Fatal("score-only screening retained per-mask details")
							}
						}
						_, count := compactMaskSelection(1024, idx.maskSelection, selection)
						fraction := float64((*a)[0].HitMasks) / float64(count)
						if fraction <= 0 || fraction >= 1 {
							t.Fatalf("fixture must match a proper subset of masks, got %g", fraction)
						}
						idx.RecycleGSearchScreenResult(aIDs)
						idx.RecycleGSearchScreenResult(bIDs)
						idx.RecycleGSearchScreenDetailResults(a)
						idx.RecycleGSearchScreenDetailResults(b)
						for _, details := range []bool{true, false} {
							for _, threshold := range []float64{fraction, math.Nextafter(fraction, 1), 0} {
								opt.MinMaskFrac = threshold
								ids, results, err := idx.gsearchScreen(q, 2, selection, details)
								if err != nil {
									t.Fatal(err)
								}
								if threshold > fraction {
									if ids != nil || results != nil {
										t.Fatal("filtered candidates reached alignment")
									}
									continue
								}
								if ids == nil || results == nil || len(*results) != 4 {
									t.Fatal("selected-mask fraction dropped boundary ties")
								}
								idx.RecycleGSearchScreenResult(ids)
								idx.RecycleGSearchScreenDetailResults(results)
								opt.MinMaskFrac = 0
							}
						}
					}
				})
			}
			var want []byte
			for _, mode := range [][]string{{"-j", "1"}, {"-j", "4"}, {"-j", "4", "-w"}} {
				args := []string{"genome", "search", "--quiet", "-d", db, "-S", "--extra", "-m", strconv.Itoa(masks), "-p", "12", "-n", "1", "--windows", "2"}
				args = append(args, mode...)
				// Repeated partial queries also exercise concurrent use and recycling.
				got := runSeedIndexCommand(t, append(args, query, query)...)
				lines := strings.Split(strings.TrimSpace(string(got)), "\n")
				if len(lines) != 9 { // header + both queries' four tied references
					t.Fatalf("expected all cutoff ties for both queries, got %s", got)
				}
				for _, line := range lines[1:] {
					fields := strings.Split(line, "\t")
					hits, err := strconv.Atoi(fields[4])
					if err != nil || hits == 0 {
						t.Fatalf("invalid hit count: %s", line)
					}
					score, err := strconv.Atoi(fields[5])
					if err != nil || score != hits*16 {
						t.Fatalf("repeated exact seeds must contribute once per mask: %s", line)
					}
					if fields[8] != strings.TrimSuffix(strings.Repeat("16,", hits), ",") {
						t.Fatalf("incorrect per-mask maxima: %s", line)
					}
				}
				if want == nil {
					want = got
				} else if !bytes.Equal(got, want) {
					t.Fatalf("screening changed with mode %v:\ngot %s\nwant %s", mode, got, want)
				}
			}
		})
	}
}
