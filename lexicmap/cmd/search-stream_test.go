package cmd

import (
	"bufio"
	"bytes"
	"io"
	"os"
	"testing"
)

// streamTestResult builds a minimal SearchResult with the fields emit()
// reads: SimilarityDetails[0].SimilarityScore (ordering key) and SeqID
// (SortBySeqID). Each seqids entry becomes one SimilarityDetail sharing
// the given score.
func streamTestResult(score float64, seqids ...string) *SearchResult {
	sds := make([]*SimilarityDetail, len(seqids))
	for i, id := range seqids {
		sds[i] = &SimilarityDetail{SimilarityScore: score, SeqID: []byte(id)}
	}
	return &SearchResult{SimilarityDetails: &sds}
}

// streamTestRows returns a writeRows callback that writes a fixed row
// payload regardless of the result.
func streamTestRows(rows string) func(io.Writer, *SearchResult) {
	return func(w io.Writer, _ *SearchResult) {
		io.WriteString(w, rows)
	}
}

func drainToString(t *testing.T, s *resultStreamer, rowPrefix string) string {
	t.Helper()
	var out bytes.Buffer
	bw := bufio.NewWriter(&out)
	s.drain(bw, rowPrefix)
	if err := bw.Flush(); err != nil {
		t.Fatal(err)
	}
	return out.String()
}

// Drain must emit results ordered by descending SimilarityScore, matching
// the buffered path's final sort.
func TestResultStreamerDrainOrdersByScore(t *testing.T) {
	s := &resultStreamer{dir: t.TempDir()}
	defer s.close()

	s.emit(streamTestResult(10, "a"), streamTestRows("g10\n"))
	s.emit(streamTestResult(30, "b"), streamTestRows("g30\n"))
	s.emit(streamTestResult(20, "c"), streamTestRows("g20\n"))

	got := drainToString(t, s, "q\t100\t3\t")
	want := "q\t100\t3\tg30\n" +
		"q\t100\t3\tg20\n" +
		"q\t100\t3\tg10\n"
	if got != want {
		t.Fatalf("drain order/prefix:\n got %q\nwant %q", got, want)
	}
}

// A genome hit can produce multiple rows (one per HSP); every physical
// line must get the query/qlen/hits prefix, not just the first.
// Regression test for the multi-HSP output corruption.
func TestResultStreamerDrainPrefixesEveryLine(t *testing.T) {
	s := &resultStreamer{dir: t.TempDir()}
	defer s.close()

	s.emit(streamTestResult(1, "a"), streamTestRows("hsp1\nhsp2\n"))
	s.emit(streamTestResult(2, "b"), streamTestRows("only\n"))

	got := drainToString(t, s, "p\t")
	want := "p\tonly\np\thsp1\np\thsp2\n"
	if got != want {
		t.Fatalf("per-line prefix:\n got %q\nwant %q", got, want)
	}
}

// The ordering score is SimilarityDetails[0].SimilarityScore *before*
// SortBySeqID reorders the details. Capturing it after the sort would
// order this record by 9 (the "a" detail) instead of 5.
func TestResultStreamerEmitCapturesScoreBeforeSeqIDSort(t *testing.T) {
	s := &resultStreamer{dir: t.TempDir()}
	defer s.close()

	sds := []*SimilarityDetail{
		{SimilarityScore: 5, SeqID: []byte("z")},
		{SimilarityScore: 9, SeqID: []byte("a")},
	}
	s.emit(&SearchResult{SimilarityDetails: &sds}, streamTestRows("five\n"))
	s.emit(streamTestResult(7, "m"), streamTestRows("seven\n"))

	got := drainToString(t, s, "p\t")
	want := "p\tseven\np\tfive\n"
	if got != want {
		t.Fatalf("ordering used post-sort score:\n got %q\nwant %q", got, want)
	}
}

// emit() appends to the same file using a reused scratch buffer; a second
// record must not inherit trailing bytes of the first.
func TestResultStreamerBufferReuse(t *testing.T) {
	s := &resultStreamer{dir: t.TempDir()}
	defer s.close()

	s.emit(streamTestResult(1, "a"), streamTestRows("a-very-long-row-payload\n"))
	s.emit(streamTestResult(2, "b"), streamTestRows("x\n"))

	got := drainToString(t, s, "p\t")
	want := "p\tx\np\ta-very-long-row-payload\n"
	if got != want {
		t.Fatalf("buffer reuse corrupted rows:\n got %q\nwant %q", got, want)
	}
}

// The temp file is unlinked on open (POSIX), so it must not be visible in
// the directory while still in use, and close() must leave nothing behind.
func TestResultStreamerTempFileUnlinkedOnOpen(t *testing.T) {
	dir := t.TempDir()
	s := &resultStreamer{dir: dir}

	s.emit(streamTestResult(1, "a"), streamTestRows("g\n"))
	if s.fh == nil {
		t.Fatal("emit did not create the temp file")
	}
	if !s.unlinked {
		t.Skip("filesystem does not support unlink-while-open (e.g. Windows)")
	}
	if _, err := os.Stat(s.name); !os.IsNotExist(err) {
		t.Fatalf("temp file %s still linked: %v", s.name, err)
	}

	s.close()
	entries, err := os.ReadDir(dir)
	if err != nil {
		t.Fatal(err)
	}
	if len(entries) != 0 {
		t.Fatalf("temp dir not clean after close: %v", entries)
	}
}

// A query with no hits creates no file and drains nothing.
func TestResultStreamerEmpty(t *testing.T) {
	s := &resultStreamer{dir: t.TempDir()}
	defer s.close()

	if got := drainToString(t, s, "q\t0\t0\t"); got != "" {
		t.Fatalf("empty streamer produced output: %q", got)
	}
	if s.fh != nil {
		t.Fatal("temp file created although nothing was emitted")
	}
	// close() on a never-used streamer must not panic
}
