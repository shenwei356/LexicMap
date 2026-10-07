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
