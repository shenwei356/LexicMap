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

package kv

import (
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"testing"
)

// Skipped postings must report truncation even in the last pair, where no
// subsequent read would otherwise expose the error. Every exit returns its kit.
func TestSearcherTruncatedSkippedPositions(t *testing.T) {
	const k uint8 = 16
	const key uint64 = 1 << 28
	for _, threeBytes := range []bool{false, true} {
		for _, tc := range []struct {
			name      string
			keys      int
			query     uint64
			checkFlag bool
		}{
			{"first/strand", 1, key, true},
			{"second/strand", 2, key + 1, true},
			{"first/below-query", 1, key + 2, false},
			{"second/below-query", 2, key + 2, false},
		} {
			t.Run(fmt.Sprintf("threeBytes=%t/%s", threeBytes, tc.name), func(t *testing.T) {
				file := filepath.Join(t.TempDir(), "seeds")
				writer, err := NewWriter(k, 0, 1, file, 1, 1, threeBytes)
				if err != nil {
					t.Fatal(err)
				}
				data := make(map[uint64]*[]uint64)
				for i := range tc.keys {
					values := []uint64{0, 4} // Both postings have the forward flag.
					data[key+uint64(i)] = &values
				}
				writeErr := writer.WriteDataOfAMask(data)
				if err := errors.Join(writeErr, writer.Close()); err != nil {
					t.Fatal(err)
				}
				original, err := os.ReadFile(file)
				if err != nil {
					t.Fatal(err)
				}
				searcher, err := NewSearcher(file, 1)
				if err != nil {
					t.Fatal(err)
				}
				defer func() {
					if err := searcher.Close(); err != nil {
						t.Error(err)
					}
				}()
				for _, multi := range []bool{false, true} {
					for _, stream := range []bool{false, true} {
						t.Run(fmt.Sprintf("multi=%t/stream=%t", multi, stream), func(t *testing.T) {
							search := func() error {
								var results *[]SearchResult
								var err error
								consume := func(SearchResult) error {
									t.Error("returned a skipped key")
									return nil
								}
								switch {
								case multi && stream:
									err = searcher.Search2Stream([][]uint64{{tc.query}}, k, tc.checkFlag, true, consume)
								case multi:
									results, err = searcher.Search2([][]uint64{{tc.query}}, k, tc.checkFlag, true)
								case stream:
									err = searcher.SearchStream([]uint64{tc.query}, k, tc.checkFlag, true, consume)
								default:
									results, err = searcher.Search([]uint64{tc.query}, k, tc.checkFlag, true)
								}
								if results != nil {
									if err != nil || len(*results) != 0 {
										t.Error("returned results after a skip or read failure")
									}
									RecycleSearchResults(results)
								}
								if len(searcher.searchKits) != 1 {
									t.Fatal("search did not return its decoding kit")
								}
								return err
							}
							if err := search(); err != nil {
								t.Fatalf("intact input: %v", err)
							}
							positionBytes := 8
							if threeBytes {
								positionBytes = 7
							}
							if err := os.Truncate(file, int64(len(original)-positionBytes)); err != nil {
								t.Fatal(err)
							}
							if err := search(); !errors.Is(err, ErrBrokenFile) || !errors.Is(err, io.EOF) {
								t.Errorf("truncated skipped positions: got %v, want broken file with EOF", err)
							}
							if err := os.WriteFile(file, original, 0600); err != nil {
								t.Fatal(err)
							}
							if err := search(); err != nil {
								t.Fatalf("reuse after read failure: %v", err)
							}
						})
					}
				}
			})
		}
	}
}
