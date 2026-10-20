package cmd

import (
	"bytes"
	"os"
	"os/exec"
	"os/signal"
	"path/filepath"
	"syscall"
	"testing"
)

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
