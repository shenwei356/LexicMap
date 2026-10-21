package cmd

import (
	"bytes"
	"fmt"
	"os"
	"os/exec"
	"os/signal"
	"path/filepath"
	"slices"
	"strings"
	"syscall"
	"testing"
)

// Parsing can hold an eligibility reader beside the FASTA reader; seed writing
// holds two files per chunk, or three with idx15. Both phases share the cap.
func TestIndexConstructionFileBudget(t *testing.T) {
	dir := t.TempDir()
	dna := indexTestDNA(6000)
	fasta := ">s0\n" + dna[:2000] + "\n>s1\n" + dna[2000:4000] + "\n>s2\n" + dna[4000:] + "\n"
	files := make([]string, 24)
	for i := range files {
		files[i] = writeTestInput(t, filepath.Join(dir, fmt.Sprintf("g%02d.fa", i)), fasta)
	}
	query := writeTestInput(t, filepath.Join(dir, "query.fa"), ">q\n"+dna[:800]+"\n")
	build := []string{"index", "--quiet", "-j", "24", "-m", "64", "-c", "16", "-k", "16", "--partitions", "4",
		"--max-genome", "4000", "--log", filepath.Join(dir, "index.log"), "--big-genomes", filepath.Join(dir, "skipped.txt")}
	baseline := filepath.Join(dir, "baseline")
	runSeedIndexCommand(t, append(append(slices.Clone(build), "-O", baseline), files...)...)
	search := func(db string) []string {
		out := runSeedIndexCommand(t, "search", "--quiet", "-j", "1", "-d", db, "-p", "12", "-P", "12", "-n", "0", "-N", "0", "-a", query)
		rows := strings.Split(strings.TrimRight(string(out), "\r\n"), "\n")
		if len(rows) != len(files)+1 {
			t.Fatalf("expected one match per genome, got %d rows", len(rows))
		}
		slices.Sort(rows[1:]) // batch record IDs and equal-score genome order may differ
		return rows
	}
	want := search(baseline)
	for _, adaptive := range []bool{false, true} {
		for _, budget := range []int{13, 14, 24} {
			db := filepath.Join(dir, fmt.Sprintf("budget%d-adaptive%v", budget, adaptive))
			args := append(slices.Clone(build), "-O", db, "--max-open-files", fmt.Sprint(budget))
			if budget >= 14 {
				args = append(args, "--save-seed-pos")
			}
			if adaptive {
				args = append(args, "-2")
			}
			_, stderr, err := runResourceLimitedCommand(t, "open-files="+fmt.Sprint(budget), append(args, files...)...)
			if err != nil {
				log, _ := os.ReadFile(filepath.Join(dir, "index.log"))
				t.Fatalf("budget=%d adaptive=%v: %v, %s\n%s", budget, adaptive, err, stderr, log)
			}
			if !slices.Equal(search(db), want) {
				t.Fatalf("budget=%d adaptive=%v: search output changed", budget, adaptive)
			}
		}
	}
}

func TestReplayableGenomeInputPipe(t *testing.T) {
	dir := t.TempDir()
	input := filepath.Join(dir, "input.fifo")
	if err := syscall.Mkfifo(input, 0600); err != nil {
		t.Fatal(err)
	}
	want := []byte(">sequence\nACGT\n")
	done := make(chan error, 1)
	go func() {
		done <- os.WriteFile(input, want, 0600)
	}()
	file, err := replayableGenomeInput(input, dir)
	if err != nil {
		t.Fatal(err)
	}
	if err := <-done; err != nil {
		t.Fatal(err)
	}
	got, err := os.ReadFile(file)
	if err != nil || !bytes.Equal(got, want) || file == input {
		t.Fatalf("pipe staging: %q, %v", got, err)
	}
}

// TestWriteIndexInfoWriteFailure forces a real short write in a subprocess.
// The file-size limit must not affect the rest of the test suite.
func TestWriteIndexInfoWriteFailure(t *testing.T) {
	if file := os.Getenv("LEXICMAP_INFO_WRITE_FAILURE"); file != "" {
		signal.Ignore(syscall.SIGXFSZ)
		if err := syscall.Setrlimit(syscall.RLIMIT_FSIZE, &syscall.Rlimit{Cur: 128, Max: 128}); err != nil {
			t.Fatal(err)
		}
		if err := writeIndexInfo(file, &IndexInfo{K: 31, Genomes: 9}); err == nil {
			t.Fatal("ignored a partial metadata write")
		}
		return
	}
	dir := t.TempDir()
	file := writeTestInput(t, filepath.Join(dir, FileInfo), "original metadata\n")
	command := exec.Command(os.Args[0], "-test.run=^TestWriteIndexInfoWriteFailure$")
	command.Env = append(os.Environ(), "LEXICMAP_INFO_WRITE_FAILURE="+file)
	if out, err := command.CombinedOutput(); err != nil {
		t.Fatalf("write-failure subprocess: %v\n%s", err, out)
	}
	got, err := os.ReadFile(file)
	if err != nil || !bytes.Equal(got, []byte("original metadata\n")) {
		t.Fatalf("failed write changed old metadata: %q, %v", got, err)
	}
	entries, err := os.ReadDir(dir)
	if err != nil || len(entries) != 1 {
		t.Fatalf("temporary metadata remains: %v, %v", entries, err)
	}
}

func TestWriteIndexInfoUmask(t *testing.T) {
	if file := os.Getenv("LEXICMAP_INFO_UMASK"); file != "" {
		syscall.Umask(0077)
		if err := writeIndexInfo(file, &IndexInfo{K: 31}); err != nil {
			t.Fatal(err)
		}
		stat, err := os.Stat(file)
		if err != nil || stat.Mode().Perm() != 0600 {
			t.Fatalf("new metadata ignored umask: %v, %v", stat, err)
		}
		if err := os.Chmod(file, 0640); err != nil {
			t.Fatal(err)
		}
		if err := writeIndexInfo(file, &IndexInfo{K: 31, Genomes: 1}); err != nil {
			t.Fatal(err)
		}
		stat, err = os.Stat(file)
		if err != nil || stat.Mode().Perm() != 0640 {
			t.Fatalf("update changed existing permissions: %v, %v", stat, err)
		}
		return
	}
	command := exec.Command(os.Args[0], "-test.run=^TestWriteIndexInfoUmask$")
	command.Env = append(os.Environ(), "LEXICMAP_INFO_UMASK="+filepath.Join(t.TempDir(), FileInfo))
	if out, err := command.CombinedOutput(); err != nil {
		t.Fatalf("umask subprocess: %v\n%s", err, out)
	}
}
