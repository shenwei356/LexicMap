package cmd

import (
	"bytes"
	"compress/gzip"
	"io"
	"os"
	"path/filepath"
	"testing"
)

// TestForceGzipOutput checks inherited flags and complete gzip streams through the CLI.
func TestForceGzipOutput(t *testing.T) {
	args := []string{"utils", "masks", "--quiet", "-j", "1", "-k", "16", "-m", "64"}
	want := runSeedIndexCommand(t, args...)
	if len(want) == 0 {
		t.Fatal("empty reference output")
	}
	for _, tt := range []struct {
		name     string
		file     string
		flag     string
		wantGzip bool
	}{
		{name: "plain stdout", file: "-"},
		{name: "plain file", file: "output.tsv"},
		{name: "gzip suffix", file: "output.tsv.gz", wantGzip: true},
		{name: "forced stdout", file: "-", flag: "--force-gzip-output", wantGzip: true},
		{name: "forced file", file: "output.tsv", flag: "--force-gzip-output", wantGzip: true},
		{name: "forced gzip suffix", file: "output.tsv.gz", flag: "--force-gzip-output", wantGzip: true},
		{name: "disabled stdout", file: "-", flag: "--force-gzip-output=false"},
		{name: "disabled file", file: "output.tsv", flag: "--force-gzip-output=false"},
		{name: "disabled with gzip suffix", file: "output.tsv.gz", flag: "--force-gzip-output=false", wantGzip: true},
	} {
		t.Run(tt.name, func(t *testing.T) {
			file := tt.file
			if file != "-" {
				file = filepath.Join(t.TempDir(), file)
			}
			commandArgs := append(append([]string(nil), args...), "-o", file)
			if tt.flag != "" {
				commandArgs = append(commandArgs, tt.flag)
			}
			got := runSeedIndexCommand(t, commandArgs...)
			if file != "-" {
				var err error
				got, err = os.ReadFile(file)
				if err != nil {
					t.Fatal(err)
				}
			}
			if tt.wantGzip {
				reader, err := gzip.NewReader(bytes.NewReader(got))
				if err != nil {
					t.Fatal(err)
				}
				// Read through EOF to validate the gzip trailer as well as the payload.
				got, err = io.ReadAll(reader)
				if err != nil {
					t.Fatal(err)
				}
				if err := reader.Close(); err != nil {
					t.Fatal(err)
				}
			}
			if !bytes.Equal(got, want) {
				t.Fatalf("output differs from uncompressed reference: %q", got)
			}
		})
	}
}
