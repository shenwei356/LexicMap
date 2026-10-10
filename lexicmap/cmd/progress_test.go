//go:build linux

package cmd

import (
	"bytes"
	"context"
	"errors"
	"fmt"
	"io"
	"math/rand"
	"os"
	"os/exec"
	"path/filepath"
	"testing"
	"time"

	"golang.org/x/sys/unix"
)

// progressTestTerminal opens a Linux pseudo-terminal so mpb renders its bars;
// writing stderr to a pipe would disable rendering regardless of quiet mode.
func progressTestTerminal(t *testing.T) (*os.File, *os.File) {
	t.Helper()
	master, err := os.OpenFile("/dev/ptmx", os.O_RDWR|unix.O_NOCTTY, 0)
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { master.Close() })
	if err := unix.IoctlSetPointerInt(int(master.Fd()), unix.TIOCSPTLCK, 0); err != nil {
		t.Fatal(err)
	}
	number, err := unix.IoctlGetInt(int(master.Fd()), unix.TIOCGPTN)
	if err != nil {
		t.Fatal(err)
	}
	slave, err := os.OpenFile(fmt.Sprintf("/dev/pts/%d", number), os.O_RDWR|unix.O_NOCTTY, 0)
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { slave.Close() })
	if err := unix.IoctlSetWinsize(int(slave.Fd()), unix.TIOCSWINSZ, &unix.Winsize{Row: 24, Col: 100}); err != nil {
		t.Fatal(err)
	}
	return master, slave
}

// TestDebugProgressRespectsQuiet exercises progress creation, updates, and shutdown
// through each public command, including both genome-search alignment modes.
func TestDebugProgressRespectsQuiet(t *testing.T) {
	dir := t.TempDir()
	rng := rand.New(rand.NewSource(1))
	sequence := make([]byte, 4080)
	for i := range sequence {
		sequence[i] = "ACGT"[rng.Intn(4)]
	}
	query := filepath.Join(dir, "query.fa")
	reference := filepath.Join(dir, "reference.fa")
	for _, file := range []string{query, reference} {
		if err := os.WriteFile(file, append([]byte(">seq\n"), sequence...), 0600); err != nil {
			t.Fatal(err)
		}
	}
	db := filepath.Join(dir, "index.lmi")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-c", "2", "-m", "64", "-k", "16", "--partitions", "4", "-O", db, reference)

	for _, command := range []struct {
		name  string
		args  []string
		label string
	}{
		{"search", []string{"search", "-d", db, "-p", "12", "-P", "12", query}, "checked genomes:"},
		{"genome-search", []string{"genome", "search", "-d", db, "-m", "0", "-p", "12", query}, "checked subject genomes:"},
		{"genome-search-orthoani", []string{"genome", "search", "--OrthoANI", "-d", db, "-m", "0", "-p", "12", query}, "checked subject genomes:"},
		{"genome-compare", []string{"genome", "compare", query, reference}, "compared genome pairs:"},
	} {
		t.Run(command.name, func(t *testing.T) {
			var baseline []byte
			for _, mode := range []struct {
				name  string
				debug bool
				quiet bool
			}{
				{"normal", false, false},
				{"quiet", false, true},
				{"debug", true, false},
				{"debug-quiet", true, true},
			} {
				t.Run(mode.name, func(t *testing.T) {
					args := append([]string{}, command.args...)
					args = append(args, "-j", "2")
					if mode.debug {
						args = append(args, "--debug")
					}
					if mode.quiet {
						args = append(args, "--quiet")
					}
					// A log destination must not re-enable terminal progress in quiet mode.
					args = append(args, "--log", filepath.Join(t.TempDir(), "run.log"))
					ctx, cancel := context.WithTimeout(context.Background(), 30*time.Second)
					defer cancel()
					process := exec.CommandContext(ctx, os.Args[0], append([]string{"-test.run=^TestSeedIndexCommand$", "--"}, args...)...)
					process.Env = append(os.Environ(), "LEXICMAP_SEED_INDEX_TEST_HELPER=1")
					var stdout, stderr bytes.Buffer
					master, slave := progressTestTerminal(t)
					process.Stdout, process.Stderr = &stdout, slave
					done := make(chan error, 1)
					go func() {
						_, err := io.Copy(&stderr, master)
						done <- err
					}()
					err := process.Run()
					slave.Close()
					// Linux PTY reads end with EIO when the last slave is closed.
					if readErr := <-done; readErr != nil && !errors.Is(readErr, unix.EIO) {
						t.Fatal(readErr)
					}
					if err != nil {
						t.Fatalf("command %v failed: %v\n%s", args, err, stderr.String())
					}
					if got, want := bytes.Contains(stderr.Bytes(), []byte(command.label)), mode.debug && !mode.quiet; got != want {
						t.Fatalf("progress visible=%v, want %v; stderr:\n%s", got, want, stderr.String())
					}
					if bytes.Count(stdout.Bytes(), []byte("\n")) < 2 {
						t.Fatalf("fixture produced no results: %s", stdout.String())
					}
					if baseline == nil {
						baseline = bytes.Clone(stdout.Bytes())
					} else if !bytes.Equal(stdout.Bytes(), baseline) {
						t.Fatalf("results changed with logging mode:\ngot: %s\nwant: %s", stdout.Bytes(), baseline)
					}
				})
			}
		})
	}
}
