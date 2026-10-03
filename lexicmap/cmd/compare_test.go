package cmd

import (
	"os"
	"path/filepath"
	"reflect"
	"testing"
)

func TestCombinationCount(t *testing.T) {
	if got, err := combinationCount(3, "genome IDs"); err != nil || got != 3 {
		t.Fatalf("count: got %d, err=%v, want 3", got, err)
	}
	if _, err := combinationCount(1, "genome IDs"); err == nil {
		t.Fatal("expected an error for fewer than two genome IDs")
	}
}

func TestReadPairFilesSkipsHeaderInEachFile(t *testing.T) {
	dir := t.TempDir()
	file1 := filepath.Join(dir, "pairs1.tsv")
	file2 := filepath.Join(dir, "pairs2.tsv")
	if err := os.WriteFile(file1, []byte("genome1\tgenome2\nA\tB\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(file2, []byte("genome1\tgenome2\nC\tD\textra\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	got, err := readPairFiles([]string{file1, file2}, true)
	if err != nil {
		t.Fatal(err)
	}
	want := []string{"A", "B", "C", "D"}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("pairs: got %v, want %v", got, want)
	}
}
