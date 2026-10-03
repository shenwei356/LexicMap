package cmd

import (
	"bytes"
	"fmt"
	"os"
	"path/filepath"
	"testing"

	"github.com/shenwei356/lexichash"
)

func TestLoadIndexMasksUsesCustomMaskCount(t *testing.T) {
	lh, err := lexichash.NewWithSeed(31, 64, 1, 0)
	if err != nil {
		t.Fatal(err)
	}

	var data bytes.Buffer
	for _, mask := range lh.Masks {
		if _, err := fmt.Fprintf(&data, "%s\n", lexichash.MustDecode(mask, uint8(lh.K))); err != nil {
			t.Fatal(err)
		}
	}

	file := filepath.Join(t.TempDir(), "masks.txt")
	if err := os.WriteFile(file, data.Bytes(), 0o644); err != nil {
		t.Fatal(err)
	}

	opt := &IndexBuildingOptions{MaskFile: file, K: 21, Masks: 20_000, Chunks: 16}
	got, err := loadIndexMasks(opt)
	if err != nil {
		t.Fatal(err)
	}
	if got.K != 31 || opt.K != 31 {
		t.Fatalf("k-mer size: loaded=%d option=%d, want 31", got.K, opt.K)
	}
	if len(got.Masks) != 64 || opt.Masks != 64 {
		t.Fatalf("mask count: loaded=%d option=%d, want 64", len(got.Masks), opt.Masks)
	}
}
