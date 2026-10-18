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

import "testing"

func TestCoverageLen(t *testing.T) {
	for _, tt := range []struct {
		name    string
		regions [][2]int
		want    int
	}{
		{name: "empty"},
		{name: "single base", regions: [][2]int{{4, 4}}, want: 1},
		{name: "inclusive endpoints", regions: [][2]int{{0, 4}}, want: 5},
		{name: "duplicates", regions: [][2]int{{0, 4}, {0, 4}}, want: 5},
		{name: "nested", regions: [][2]int{{2, 4}, {0, 8}, {0, 3}}, want: 9},
		{name: "disjoint", regions: [][2]int{{8, 10}, {1, 3}}, want: 6},
		{name: "adjacent", regions: [][2]int{{3, 5}, {0, 2}}, want: 6},
		{name: "bridged", regions: [][2]int{{7, 9}, {1, 3}, {3, 7}}, want: 9},
		{name: "shared endpoint", regions: [][2]int{{0, 3}, {3, 5}}, want: 6},
	} {
		t.Run(tt.name, func(t *testing.T) {
			if got := coverageLen(&tt.regions); got != tt.want {
				t.Fatalf("coverage = %d, want %d", got, tt.want)
			}
		})
	}
}
