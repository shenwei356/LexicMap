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
	"path/filepath"
	"reflect"
	"regexp"
	"slices"
	"strconv"
	"strings"
	"testing"
	"time"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/util"
	"github.com/shenwei356/lexichash"
)

// TestQueryDesertsMatchIndex checks the independent index implementation,
// including repeat positions, low complexity, N gaps, and both orientations.
func TestQueryDesertsMatchIndex(t *testing.T) {
	for _, reverse := range []bool{false, true} {
		t.Run(fmt.Sprint(reverse), func(t *testing.T) {
			sequence := []byte(indexTestDNA(60001))
			// Short N runs also change global winners when the sequence is reversed.
			for p := 100; p < len(sequence); p += 503 {
				sequence[p] = 'N'
			}
			copy(sequence[5000:6500], sequence[8000:9500])
			copy(sequence[20000:21000], bytes.Repeat([]byte{'N'}, 1000))
			copy(sequence[40000:41000], bytes.Repeat([]byte{'A'}, 1000))
			if reverse {
				RC(sequence)
			}
			dir := t.TempDir()
			ref := writeTestInput(t, filepath.Join(dir, "ref.fa"), ">s\n"+string(sequence)+"\n")
			db := filepath.Join(dir, "db")
			runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-c", "2", "--partitions", "4", "-O", db, ref)
			lh, err := lexichash.NewFromFile(filepath.Join(db, FileMasks))
			if err != nil {
				t.Fatal(err)
			}
			if err = lh.IndexMasks(3); err != nil {
				t.Fatal(err)
			}
			if err = lh.IndexMasksWithDistinctPrefixes(4); err != nil {
				t.Fatal(err)
			}
			seeds, locations, err := lh.MaskKnownDistinctPrefixes(sequence, queryMaskSkipRegions(sequence), true)
			if err != nil {
				t.Fatal(err)
			}
			defer lh.RecycleMaskResult(seeds, locations)
			for i, code := range *seeds {
				if util.IsLowComplexityDust(code, uint8(lh.K)) {
					(*seeds)[i] = 0
				}
			}
			for _, length := range []int{49999, 50000} {
				extra, locs, err := fillQueryDeserts(lh, sequence[:length], *seeds, *locations, defaultQueryDesertThreshold)
				if err != nil || extra != nil || !reflect.DeepEqual(locs, *locations) {
					t.Fatalf("filled short query %d: %v", length, err)
				}
			}
			// 50,001 must take the filling path; use matching mask locations.
			shortSeeds, shortLocs, err := lh.MaskKnownDistinctPrefixes(sequence[:50001], nil, true)
			if err != nil {
				t.Fatal(err)
			}
			extra, _, err := fillQueryDeserts(lh, sequence[:50001], *shortSeeds, *shortLocs, defaultQueryDesertThreshold)
			lh.RecycleMaskResult(shortSeeds, shortLocs)
			if err != nil || extra == nil {
				t.Fatalf("did not fill 50,001 bases: %v", err)
			}
			extra, locs, err := fillQueryDeserts(lh, sequence, *seeds, *locations, defaultQueryDesertThreshold)
			if err != nil {
				t.Fatal(err)
			}
			if len(locs) <= len(*locations) {
				t.Fatal("no extra seeds")
			}
			if len(*locations) != len(lh.Masks) {
				t.Fatal("changed pooled slice length")
			}
			var got []string
			add := func(mask int, code uint64, positions []int) {
				for _, loc := range positions {
					got = append(got, fmt.Sprintf("%d\t%s\t%d\t%c", mask+1, lexichash.MustDecode(code, uint8(lh.K)), (loc>>1)+1, lexichash.Strands[loc&1]))
				}
			}
			for mask, code := range *seeds {
				if code != 0 {
					add(mask, code, locs[mask])
				}
			}
			for mask, codes := range extra.kmers {
				for i, code := range codes {
					add(mask, code, locs[extra.locIndexes[mask][i]])
				}
			}
			// Independently build both grids on both orientations. The query
			// must contain their exact union, with no duplicate seed occurrences.
			wantSet := make(map[string]bool)
			for _, cfg := range [][2]string{{"100", "50"}, {"60", "30"}} {
				for _, rc := range []bool{false, true} {
					indexRef := ref
					if rc {
						reversed := bytes.Clone(sequence)
						RC(reversed)
						indexRef = writeTestInput(t, filepath.Join(dir, fmt.Sprintf("rc-%s.fa", cfg[0])), ">rc\n"+string(reversed)+"\n")
					}
					indexDir := filepath.Join(dir, fmt.Sprintf("expected-%s-%s-%v", cfg[0], cfg[1], rc))
					runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-c", "2", "--partitions", "4", "-D", cfg[0], "-d", cfg[1], "-O", indexDir, indexRef)

					output := runSeedIndexCommand(t, "utils", "kmers", "--quiet", "-d", indexDir, "-m", "0", "--only-forward")
					for _, line := range strings.Split(strings.TrimSpace(string(output)), "\n")[1:] {
						fields := strings.Split(line, "\t")
						if rc {
							position, err := strconv.Atoi(fields[5])
							if err != nil {
								t.Fatal(err)
							}
							fields[5] = strconv.Itoa(len(sequence) - lh.K - position + 2)
							if fields[6] == "+" {
								fields[6] = "-"
							} else {
								fields[6] = "+"
							}
						}
						wantSet[strings.Join([]string{fields[0], fields[1], fields[5], fields[6]}, "\t")] = true
					}
				}
			}

			want := make([]string, 0, len(wantSet))
			for seed := range wantSet {
				want = append(want, seed)
			}
			slices.Sort(got)
			slices.Sort(want)
			if !slices.Equal(got, want) {
				for i := 0; i < min(len(got), len(want)); i++ {
					if got[i] != want[i] {
						t.Fatalf("seed %d: query %s; index %s (counts %d/%d)", i, got[i], want[i], len(got), len(want))
					}
				}
				t.Fatalf("seed counts: query %d; index %d", len(got), len(want))
			}
		})
	}
}

// TestSearchQueryDeserts covers serial/parallel collection and disk/memory/spill
// readers with identical full-length self hits on both query strands.
func TestSearchQueryDeserts(t *testing.T) {
	dir := t.TempDir()
	sequence := []byte(indexTestDNA(60001))
	ref := writeTestInput(t, filepath.Join(dir, "ref.fa"), ">s\n"+string(sequence)+"\n")
	RC(sequence)
	query := writeTestInput(t, filepath.Join(dir, "query.fa"), ">reverse\n"+string(sequence)+"\n")
	db := filepath.Join(dir, "db")
	// Separate genome batches force the parallel collector path at -j 2.
	other := writeTestInput(t, filepath.Join(dir, "other.fa"), ">other\n"+string(sequence)+"\n")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-m", "64", "-c", "2", "--partitions", "4", "-b", "1", "-O", db, ref, other)
	var want []string
	for _, threads := range []string{"1", "2"} {
		for _, mode := range [][]string{nil, {"--load-whole-seeds"}, {"--max-seed-memory", "64K"}, {"--load-whole-seeds", "--max-seed-memory", "64K"}} {
			t.Run(threads+strings.Join(mode, "_"), func(t *testing.T) {
				args := []string{"search", "--quiet", "-j", threads, "-d", db, "-n", "0", "-N", "0", "-a", ref, query}
				output := runSeedIndexCommand(t, append(args, mode...)...)
				rows := strings.Split(strings.TrimSpace(string(output)), "\n")[1:]
				full := make(map[string]bool)
				for _, row := range rows {
					f := strings.Split(row, "\t")
					if f[3] != "ref" && f[3] != "other" {
						continue
					}
					if f[9] != strconv.Itoa(len(sequence)) || f[10] != "100.000" || f[11] != "0" || f[12] != "1" || f[13] != "60001" || f[14] != "1" || f[15] != "60001" {
						t.Fatalf("incomplete self hit: %s", row)
					}
					strand := "+"
					if (f[0] == "reverse") != (f[3] == "other") {
						strand = "-"
					}
					if f[16] != strand {
						t.Fatalf("wrong strand: %s", row)
					}
					full[f[0]+"/"+f[3]] = true
				}
				if len(full) != 4 {
					t.Fatalf("missing forward/reverse self hit: %v", full)
				}
				slices.Sort(rows) // Only concurrent query completion order is unspecified.
				if want == nil {
					want = rows
				} else if !slices.Equal(rows, want) {
					t.Fatal("output differs across search paths")
				}
			})
		}
	}
	// Timing messages belong between query masking and seed matching, and must
	// distinguish the disabled/exclusive-boundary path from actual filling.
	for _, tc := range []struct {
		name   string
		flags  []string
		filled bool
	}{
		{"filled", nil, true},
		{"disabled", []string{"--fill-desert-threshold", "0"}, false},
		{"boundary", []string{"--fill-desert-threshold", "60001"}, false},
	} {
		t.Run("masking-log/"+tc.name, func(t *testing.T) {
			args := []string{"search", "--debug", "-j", "1", "-d", db, "-n", "0", "-N", "0", ref}
			_, stderr, err := runAuditCommand(t, "", append(args, tc.flags...)...)
			if err != nil {
				t.Fatalf("search failed: %v\n%s", err, stderr)
			}
			logs := string(stderr)
			start := strings.Index(logs, "start to search")
			masked := strings.Index(logs, "finished query masking")
			matched := strings.Index(logs, "finished seed-matching")
			if start < 0 || masked <= start || matched <= masked {
				t.Fatalf("incorrect log order:\n%s", logs)
			}
			pattern := `finished query masking without desert filling in ([^\s]+)`
			if tc.filled {
				pattern = `finished query masking with desert filling in ([^\s]+) \(masking: ([^,]+), desert filling: ([^)]+)\)`
			}
			fields := regexp.MustCompile(pattern).FindStringSubmatch(logs)
			if fields == nil {
				t.Fatalf("missing masking timing:\n%s", logs)
			}
			var durations []time.Duration
			for _, field := range fields[1:] {
				duration, err := time.ParseDuration(field)
				if err != nil || duration < 0 {
					t.Fatalf("invalid duration %q: %v", field, err)
				}
				durations = append(durations, duration)
			}
			if tc.filled && durations[0] != durations[1]+durations[2] {
				t.Fatalf("masking times do not sum to total: %v", durations)
			}
		})
	}

}

// TestQueryDesertThreshold checks exclusive custom boundaries and the disabled
// path, including queries shorter than the default and option validation.
func TestQueryDesertThreshold(t *testing.T) {
	lh, err := lexichash.New(31, 64, 1)
	if err != nil {
		t.Fatal(err)
	}
	if err = lh.IndexMasks(3); err != nil {
		t.Fatal(err)
	}
	if err = lh.IndexMasksWithDistinctPrefixes(4); err != nil {
		t.Fatal(err)
	}
	sequence := []byte(indexTestDNA(30001))
	seeds, locations, err := lh.MaskKnownDistinctPrefixes(sequence, nil, true)
	if err != nil {
		t.Fatal(err)
	}
	defer lh.RecycleMaskResult(seeds, locations)
	for i, code := range *seeds {
		if util.IsLowComplexityDust(code, uint8(lh.K)) {
			(*seeds)[i] = 0
		}
	}
	for _, tc := range []struct {
		threshold int
		enabled   bool
	}{
		{0, false}, {30000, true}, {30001, false}, {50000, false},
	} {
		t.Run(strconv.Itoa(tc.threshold), func(t *testing.T) {
			extra, locs, err := fillQueryDeserts(lh, sequence, *seeds, *locations, tc.threshold)
			if err != nil {
				t.Fatal(err)
			}
			if (extra != nil) != tc.enabled {
				t.Fatalf("threshold %d: enabled=%v, want %v", tc.threshold, extra != nil, tc.enabled)
			}
			if tc.enabled {
				if len(locs) <= len(*locations) {
					t.Fatal("custom threshold added no seeds")
				}
			} else if !reflect.DeepEqual(locs, *locations) {
				t.Fatal("disabled filling changed locations")
			}
		})
	}
	opt := DefaultIndexSearchingOptions
	if opt.QueryDesertThreshold != 50000 {
		t.Fatalf("unexpected default threshold: %d", opt.QueryDesertThreshold)
	}
	for _, threshold := range []int{0, 1, 50000, -1} {
		opt.QueryDesertThreshold = threshold
		err := CheckIndexSearchingOptions(&opt)
		if (err != nil) != (threshold < 0) {
			t.Fatalf("threshold %d: %v", threshold, err)
		}
	}
	if _, err := NewIndexSearcher("", &opt); err == nil || !strings.Contains(err.Error(), "QueryDesertThreshold") {
		t.Fatalf("constructor accepted a negative threshold or returned the wrong error: %v", err)
	}
	_, stderr, err := runAuditCommand(t, "", "search", "--quiet", "-d", "unused", "--fill-desert-threshold", "-1")
	if err == nil || !strings.Contains(string(stderr), "fill-desert-threshold") {
		t.Fatalf("negative CLI threshold accepted or wrong error: %v, %s", err, stderr)
	}
}
