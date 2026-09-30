package cmd

import (
	"fmt"
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/shenwei356/LexicMap/lexicmap/cmd/genome"
)

func TestSubseqGenomeReadersWithoutPool(t *testing.T) {
	dbDir := t.TempDir()
	for batch := 0; batch < 2; batch++ {
		writeSubseqTestGenome(t, dbDir, batch)
	}

	readers, err := newSubseqGenomeReaders(dbDir, 2, 1, 2)
	if err != nil {
		t.Fatal(err)
	}
	defer func() {
		if err := readers.close(); err != nil {
			t.Errorf("close readers: %v", err)
		}
	}()
	if readers.pooled() {
		t.Fatal("expected on-demand readers when the file budget cannot create one reader per batch")
	}

	errs := make(chan error, 2)
	for batch := 0; batch < 2; batch++ {
		go func(batch int) {
			rdr, err := readers.acquire(batch)
			if err != nil {
				errs <- err
				return
			}
			g, err := rdr.SubSeq(0, 0, 3)
			if err == nil {
				genome.RecycleGenome(g)
			}
			if releaseErr := readers.release(batch, rdr); err == nil {
				err = releaseErr
			}
			errs <- err
		}(batch)
	}

	for range 2 {
		select {
		case err := <-errs:
			if err != nil {
				t.Fatal(err)
			}
		case <-time.After(2 * time.Second):
			t.Fatal("timed out waiting for an on-demand genome reader")
		}
	}
}

func writeSubseqTestGenome(t *testing.T, dbDir string, batch int) {
	t.Helper()
	dir := filepath.Join(dbDir, DirGenomes, batchDir(batch))
	if err := os.MkdirAll(dir, 0o755); err != nil {
		t.Fatal(err)
	}
	wtr, err := genome.NewWriter(filepath.Join(dir, FileGenomes), uint32(batch))
	if err != nil {
		t.Fatal(err)
	}

	g := genome.PoolGenome.Get().(*genome.Genome)
	g.Reset()
	g.ID = append(g.ID, fmt.Sprintf("genome-%d", batch)...)
	g.Seq = append(g.Seq, "ACGTACGT"...)
	g.GenomeSize = len(g.Seq)
	g.Len = len(g.Seq)
	g.NumSeqs = 1
	g.SeqSizes = append(g.SeqSizes, len(g.Seq))
	seqID := []byte("seq")
	g.SeqIDs = append(g.SeqIDs, &seqID)
	if err := wtr.Write(g); err != nil {
		genome.RecycleGenome(g)
		_ = wtr.Close()
		t.Fatal(err)
	}
	genome.RecycleGenome(g)
	if err := wtr.Close(); err != nil {
		t.Fatal(err)
	}
}
