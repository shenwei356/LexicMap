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
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
)

// recoverIndex15Backup resolves a previous interrupted publication before a
// new rebuild. The primary is the commit point: keep its matching sidecar.
// Searchers perform the same validation without modifying files.
func recoverIndex15Backup(file string) error {
	sidecar := filepath.Clean(file) + KVIndex15FileExt
	backup := sidecar + index15BackupSuffix
	if _, err := os.Lstat(backup); err != nil {
		if os.IsNotExist(err) {
			return nil
		}
		return err
	}
	primary := filepath.Clean(file) + KVIndexFileExt
	k, chunk, masks, maskPrefix, anchorPrefix, err := ReadKVIndexInfo(primary)
	if err != nil {
		return err
	}
	id, err := readTaggedIndexBuildID(primary)
	if errors.Is(err, ErrInvalidFileFormat) {
		// A committed single-level index no longer references any sidecar.
		return os.Remove(backup)
	}
	if err != nil {
		return err
	}
	info, err := os.Stat(file)
	if err != nil {
		return err
	}
	current, currentErr := openIndex15File(sidecar, k, maskPrefix, anchorPrefix, uint64(chunk), uint64(masks), uint64(info.Size()), id)
	if currentErr == nil {
		if err = current.close(); err != nil {
			return err
		}
		return os.Remove(backup) // publication committed; only cleanup was interrupted
	}
	old, err := openIndex15File(backup, k, maskPrefix, anchorPrefix, uint64(chunk), uint64(masks), uint64(info.Size()), id)
	if err != nil {
		return fmt.Errorf("cannot recover secondary seed index: %w", errors.Join(currentErr, err))
	}
	if err = old.close(); err != nil {
		return err
	}
	return os.Rename(backup, sidecar)
}

// backupIndex15File preserves the old inode cheaply where hard links work.
// Other filesystems use a synced temporary copy, renamed only when complete.
func backupIndex15File(sidecar, backup string) (bool, error) {
	if _, err := os.Lstat(backup); err == nil {
		return false, fmt.Errorf("secondary seed index backup already exists: %s", backup)
	} else if !os.IsNotExist(err) {
		return false, err
	}
	info, err := os.Stat(sidecar)
	if os.IsNotExist(err) {
		return false, nil // first upgrade from a single-level index
	}
	if err != nil {
		return false, err
	}
	if !info.Mode().IsRegular() {
		return false, fmt.Errorf("secondary seed index is not a regular file: %s", sidecar)
	}
	if err = os.Link(sidecar, backup); err == nil {
		return true, nil
	}
	src, err := os.Open(sidecar)
	if err != nil {
		return false, err
	}
	defer src.Close()
	dst, err := os.CreateTemp(filepath.Dir(backup), "."+filepath.Base(backup)+"-*")
	if err != nil {
		return false, err
	}
	defer dst.Close()
	defer os.Remove(dst.Name())
	if _, err = io.Copy(dst, src); err != nil {
		return false, err
	}
	if err = dst.Chmod(info.Mode().Perm()); err != nil {
		return false, err
	}
	if err = dst.Sync(); err != nil {
		return false, err
	}
	if err = dst.Close(); err != nil {
		return false, err
	}
	if err = os.Rename(dst.Name(), backup); err != nil {
		return false, err
	}
	return true, nil
}

// publishIndex15Files commits the primary last, keeping the previous sidecar
// readable until then. rename is passed explicitly so tests can stop at either
// publication step without changing global file operations or production data.
func publishIndex15Files(file, primaryTmp, sidecarTmp string, rename func(string, string) error) error {
	sidecar := filepath.Clean(file) + KVIndex15FileExt
	backup := sidecar + index15BackupSuffix
	hadSidecar, err := backupIndex15File(sidecar, backup)
	if err != nil {
		return err
	}
	if err = rename(sidecarTmp, sidecar); err != nil {
		if hadSidecar {
			err = errors.Join(err, os.Remove(backup))
		}
		return err
	}
	if err = rename(primaryTmp, filepath.Clean(file)+KVIndexFileExt); err != nil {
		var rollbackErr error
		if hadSidecar {
			rollbackErr = rename(backup, sidecar)
		} else {
			rollbackErr = os.Remove(sidecar)
		}
		// If rollback fails, preserve the backup for validated read-only fallback
		// and for recovery by the next reindex-seeds2 invocation.
		return errors.Join(err, rollbackErr)
	}
	if hadSidecar {
		return os.Remove(backup)
	}
	return nil
}
