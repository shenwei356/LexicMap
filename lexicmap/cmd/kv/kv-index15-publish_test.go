// Copyright © 2023-2026 Wei Shen <shenwei356@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

package kv

import (
	"bytes"
	"errors"
	"os"
	"path/filepath"
	"reflect"
	"testing"
)

// index15PublicationFixture stages a different valid generation while leaving
// the previous pair installed, just before the publication step under test.
type index15PublicationFixture struct {
	file, primaryTmp, sidecarTmp string
	primary, sidecar             []byte
	results                      [][]index15ResultSnapshot
}

func stagedIndex15Publication(t *testing.T, adaptive bool) index15PublicationFixture {
	t.Helper()
	f := index15PublicationFixture{file: reindexFixture(t)}
	if !adaptive {
		if err := CreateKVIndex(f.file, 4096); err != nil {
			t.Fatal(err)
		}
	}
	var err error
	f.primary, err = os.ReadFile(f.file + KVIndexFileExt)
	if err != nil {
		t.Fatal(err)
	}
	if adaptive {
		f.sidecar, err = os.ReadFile(f.file + KVIndex15FileExt)
		if err != nil {
			t.Fatal(err)
		}
	}
	f.results = publicationSearchResults(t, f.file)
	if _, err = CreateKVIndex15(f.file, 1024, MinIndex15Threshold); err != nil {
		t.Fatal(err)
	}
	f.primaryTmp = filepath.Join(filepath.Dir(f.file), "candidate.idx")
	f.sidecarTmp = filepath.Join(filepath.Dir(f.file), "candidate.idx15")
	for _, pair := range [][2]string{{f.file + KVIndexFileExt, f.primaryTmp}, {f.file + KVIndex15FileExt, f.sidecarTmp}} {
		if err = os.Rename(pair[0], pair[1]); err != nil {
			t.Fatal(err)
		}
	}
	if err = os.WriteFile(f.file+KVIndexFileExt, f.primary, 0600); err != nil {
		t.Fatal(err)
	}
	if adaptive {
		if err = os.WriteFile(f.file+KVIndex15FileExt, f.sidecar, 0600); err != nil {
			t.Fatal(err)
		}
	}
	return f
}

func publicationSearchResults(t *testing.T, file string) [][]index15ResultSnapshot {
	t.Helper()
	r, err := NewSearcher(file, 1)
	if err != nil {
		t.Fatal(err)
	}
	results := collectIndex15Searches(t, r, []uint64{uint64(54) << 18})
	if err = r.Close(); err != nil {
		t.Fatal(err)
	}
	return results
}

func TestIndex15PublicationRollback(t *testing.T) {
	for _, adaptive := range []bool{false, true} {
		for _, stage := range []string{"secondary", "primary", "rollback"} {
			t.Run(stage+map[bool]string{false: "/legacy", true: "/adaptive"}[adaptive], func(t *testing.T) {
				f := stagedIndex15Publication(t, adaptive)
				injected := errors.New("injected publication failure")
				backup := f.file + KVIndex15FileExt + index15BackupSuffix
				rename := func(from, to string) error {
					if (stage == "secondary" && from == f.sidecarTmp) ||
						(stage != "secondary" && from == f.primaryTmp) ||
						(stage == "rollback" && from == backup) {
						return injected
					}
					return os.Rename(from, to)
				}
				if err := publishIndex15Files(f.file, f.primaryTmp, f.sidecarTmp, rename); !errors.Is(err, injected) {
					t.Fatalf("missing injected error: %v", err)
				}
				primary, err := os.ReadFile(f.file + KVIndexFileExt)
				if err != nil || !bytes.Equal(primary, f.primary) {
					t.Fatalf("failed publication changed primary: %v", err)
				}
				if stage != "rollback" || !adaptive {
					sidecar, err := os.ReadFile(f.file + KVIndex15FileExt)
					if adaptive && (err != nil || !bytes.Equal(sidecar, f.sidecar)) {
						t.Fatalf("rollback changed old sidecar: %v", err)
					}
					if !adaptive && !os.IsNotExist(err) {
						t.Fatalf("failed first upgrade retained a sidecar: %v", err)
					}
				}
				if got := publicationSearchResults(t, f.file); !reflect.DeepEqual(got, f.results) {
					t.Fatal("search results changed after publication failure")
				}
				// A retry repairs any interrupted rollback and removes its backup.
				if _, err = CreateKVIndex15(f.file, 256, MinIndex15Threshold); err != nil {
					t.Fatal(err)
				}
				if _, err = os.Stat(backup); !os.IsNotExist(err) {
					t.Fatalf("backup remains after retry: %v", err)
				}
			})
		}
	}
}

// Stop immediately after either rename, without rollback or cleanup, to model
// process interruption. A search must work before any recovery writes occur.
func TestIndex15InterruptedPublication(t *testing.T) {
	for _, stopAfterPrimary := range []bool{false, true} {
		t.Run(map[bool]string{false: "after secondary", true: "after primary"}[stopAfterPrimary], func(t *testing.T) {
			f := stagedIndex15Publication(t, true)
			stopped := errors.New("simulated process interruption")
			func() {
				defer func() {
					if got := recover(); got != stopped {
						t.Fatalf("unexpected stop: %v", got)
					}
				}()
				err := publishIndex15Files(f.file, f.primaryTmp, f.sidecarTmp, func(from, to string) error {
					if err := os.Rename(from, to); err != nil {
						return err
					}
					if (!stopAfterPrimary && from == f.sidecarTmp) || (stopAfterPrimary && from == f.primaryTmp) {
						panic(stopped)
					}
					return nil
				})
				t.Fatalf("publication did not stop: %v", err)
			}()
			paths := []string{f.file + KVIndexFileExt, f.file + KVIndex15FileExt, f.file + KVIndex15FileExt + index15BackupSuffix}
			before := make([][]byte, len(paths))
			for i, path := range paths {
				var err error
				before[i], err = os.ReadFile(path)
				if err != nil {
					t.Fatal(err)
				}
			}
			if got := publicationSearchResults(t, f.file); !reflect.DeepEqual(got, f.results) {
				t.Fatal("interruption changed search results")
			}
			for i, path := range paths {
				after, err := os.ReadFile(path)
				if err != nil || !bytes.Equal(after, before[i]) {
					t.Fatalf("read-only fallback modified %s: %v", path, err)
				}
			}
			if _, err := CreateKVIndex15(f.file, 256, MinIndex15Threshold); err != nil {
				t.Fatal(err)
			}
			if _, err := os.Stat(paths[2]); !os.IsNotExist(err) {
				t.Fatalf("backup remains after recovery: %v", err)
			}
		})
	}
}

func TestIndex15FallbackRejectsUnmatchedBackup(t *testing.T) {
	f := stagedIndex15Publication(t, true)
	if err := publishIndex15Files(f.file, f.primaryTmp, f.sidecarTmp, os.Rename); err != nil {
		t.Fatal(err)
	}
	backup := f.file + KVIndex15FileExt + index15BackupSuffix
	if err := os.WriteFile(backup, f.sidecar, 0600); err != nil {
		t.Fatal(err)
	}
	if err := os.Remove(f.file + KVIndex15FileExt); err != nil {
		t.Fatal(err)
	}
	if r, err := NewSearcher(f.file, 1); err == nil {
		r.Close()
		t.Fatal("accepted a backup from another primary generation")
	}
	// A single-level rebuild safely discards both obsolete secondary files.
	if err := CreateKVIndex(f.file, 4096); err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(backup); !os.IsNotExist(err) {
		t.Fatalf("obsolete backup remains: %v", err)
	}
}
