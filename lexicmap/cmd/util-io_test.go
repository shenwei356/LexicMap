package cmd

import (
	"bytes"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"testing"
)

// A streaming flush must expose complete records without closing the gzip stream.
func TestStreamingGzipFlush(t *testing.T) {
	for _, level := range []int{gzip.HuffmanOnly, gzip.NoCompression, gzip.BestSpeed, 5, gzip.BestCompression} {
		t.Run(fmt.Sprint(level), func(t *testing.T) {
			file := filepath.Join(t.TempDir(), "records.gz")
			out, compressed, fh, err := openOutStream(file, true, level, true)
			if err != nil {
				t.Fatal(err)
			}
			defer func() {
				if err := closeOutStream(out, compressed, fh); err != nil {
					t.Error(err)
				}
			}()
			var want []byte
			for _, record := range []string{"first record\n", "second record\n"} {
				if _, err := out.WriteString(record); err != nil {
					t.Fatal(err)
				}
				if err := flushOutStream(out, compressed); err != nil {
					t.Fatal(err)
				}
				want = append(want, record...)
				data, err := os.ReadFile(file)
				if err != nil {
					t.Fatal(err)
				}
				r, err := gzip.NewReader(bytes.NewReader(data))
				if err != nil {
					t.Fatal(err)
				}
				got := make([]byte, len(want))
				_, readErr := io.ReadFull(r, got)
				closeErr := r.Close()
				if readErr != nil || closeErr != nil || !bytes.Equal(got, want) {
					t.Fatalf("record remained buffered: read=%v, close=%v, got=%q", readErr, closeErr, got)
				}
			}
		})
	}
}
