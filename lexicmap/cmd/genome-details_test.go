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

package cmd

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// genomeDetailsTestCache includes empty contigs/IDs, split genomes, and multiple records.
func genomeDetailsTestCache(t testing.TB, saveIDs bool) []byte {
	t.Helper()
	var b bytes.Buffer
	put := func(v any) {
		t.Helper()
		if err := binary.Write(&b, binary.BigEndian, v); err != nil {
			t.Fatal(err)
		}
	}
	var flags uint64
	if saveIDs {
		flags = FLAG_SAVE_SEQIDS
	}
	put(flags)
	for record, chunks := range []int{2, 1} {
		id := fmt.Sprintf("g%d", record)
		put(uint16(len(id)))
		put(uint32(chunks))
		b.WriteString(id)
		for chunk := 0; chunk < chunks; chunk++ {
			put(uint64(record*2 + chunk))
			put(uint32(5))
			put(uint32(3))
			for _, size := range []uint32{2, 0, 3} {
				put(size)
			}
			if saveIDs {
				put(uint32(8)) // Three lengths and two one-byte IDs.
				for _, id := range []string{"a", "", "b"} {
					put(uint16(len(id)))
					b.WriteString(id)
				}
			}
		}
	}
	return b.Bytes()
}

// TestGenomeDetailsOutput checks the complete default/extra output and every truncation.
func TestGenomeDetailsOutput(t *testing.T) {
	for _, ids := range []bool{false, true} {
		for _, extra := range []bool{false, true} {
			t.Run(fmt.Sprintf("ids=%t/extra=%t", ids, extra), func(t *testing.T) {
				data := genomeDetailsTestCache(t, ids)
				file := filepath.Join(t.TempDir(), "details.bin")
				for length := 0; length <= len(data); length++ {
					if err := os.WriteFile(file, data[:length], 0600); err != nil {
						t.Fatal(err)
					}
					var out bytes.Buffer
					w := bufio.NewWriter(&out)
					err := readGenomeDetails(file, w, extra, 3)
					if length < len(data) {
						if !errors.Is(err, ErrBrokenFile) {
							t.Fatalf("truncation at %d: %v", length, err)
						}
						continue
					}
					if err != nil {
						t.Fatal(err)
					}
					if err = w.Flush(); err != nil {
						t.Fatal(err)
					}
					want := "ref\tgenome_size\tchunks\tchunk\tcidx\tgidx\tchunk_size\tseqs"
					suffix := ""
					if extra {
						want += "\tseqsizes\tseqids"
						suffix = "\t2,0,3\t"
						if ids {
							suffix += "a,,b"
						}
					}
					want += "\n" + "g0\t10\t2\t1\t0\t0\t5\t3" + suffix + "\n" +
						"g0\t10\t2\t2\t0\t1\t5\t3" + suffix + "\n" +
						"g1\t5\t1\t1\t0\t2\t5\t3" + suffix + "\n"
					if out.String() != want {
						t.Fatalf("output changed:\n%s\nwant:\n%s", out.String(), want)
					}
				}
			})
		}
	}
}

// TestGenomeDetailsLongID exercises buffered and unread bytes, including a truncated skipped ID.
func TestGenomeDetailsLongID(t *testing.T) {
	id := strings.Repeat("x", 65535)
	data := make([]byte, 41+len(id))
	binary.BigEndian.PutUint64(data[:8], FLAG_SAVE_SEQIDS)
	binary.BigEndian.PutUint16(data[8:10], 1)
	binary.BigEndian.PutUint32(data[10:14], 1)
	data[14] = 'g'
	binary.BigEndian.PutUint32(data[23:27], 1)
	binary.BigEndian.PutUint32(data[27:31], 1)
	binary.BigEndian.PutUint32(data[31:35], 1)
	binary.BigEndian.PutUint32(data[35:39], uint32(2+len(id)))
	binary.BigEndian.PutUint16(data[39:41], uint16(len(id)))
	copy(data[41:], id)
	file := filepath.Join(t.TempDir(), "details.bin")
	for _, extra := range []bool{false, true} {
		for _, length := range []int{IOBufferSize - 1, IOBufferSize, len(data) - 1, len(data)} {
			if err := os.WriteFile(file, data[:length], 0600); err != nil {
				t.Fatal(err)
			}
			var out bytes.Buffer
			w := bufio.NewWriter(&out)
			err := readGenomeDetails(file, w, extra, 1)
			if length < len(data) {
				if !errors.Is(err, ErrBrokenFile) {
					t.Fatalf("accepted truncated ID: extra=%t length=%d err=%v", extra, length, err)
				}
				continue
			}
			if err != nil {
				t.Fatal(err)
			}
			if err = w.Flush(); err != nil {
				t.Fatal(err)
			}
			want := "g\t1\t1\t1\t0\t0\t1\t1"
			if extra {
				want += "\t1\t" + id
			}
			rows := strings.Split(strings.TrimSuffix(out.String(), "\n"), "\n")
			if len(rows) != 2 || rows[1] != want {
				t.Fatal("long-ID output changed")
			}
		}
	}
}
