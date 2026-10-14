package cmd

import (
	"bytes"
	"fmt"
	"math/rand"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
)

// TestGenomeScreenRepeatedSeeds checks that splitting long posting lists into
// borrowed batches preserves per-mask maxima, cutoff ties, and query isolation.
// All references are generated here; no local database is needed.
func TestGenomeScreenRepeatedSeeds(t *testing.T) {
	dir := t.TempDir()
	rng := rand.New(rand.NewSource(1))
	motif := make([]byte, 512)
	for i := range motif {
		motif[i] = "ACGT"[rng.Intn(4)]
	}
	// Each captured k-mer has more occurrences than a streaming decoder batch.
	sequence := bytes.Repeat(motif, 300)
	var refs []string
	for i := range 4 {
		file := filepath.Join(dir, fmt.Sprintf("ref%d.fa", i))
		if err := os.WriteFile(file, append([]byte(">seq\n"), sequence...), 0600); err != nil {
			t.Fatal(err)
		}
		refs = append(refs, file)
	}
	query := filepath.Join(dir, "partial.fa")
	if err := os.WriteFile(query, append([]byte(">seq\n"), sequence[:4096]...), 0600); err != nil {
		t.Fatal(err)
	}
	db := filepath.Join(dir, "index.lmi")
	args := []string{"index", "--quiet", "-j", "2", "-b", "2", "-c", "2", "-m", "1024", "-k", "16", "--partitions", "4", "-O", db}
	runSeedIndexCommand(t, append(args, refs...)...)
	for _, masks := range []int{0, 256} {
		t.Run(fmt.Sprintf("masks=%d", masks), func(t *testing.T) {
			var want []byte
			for _, mode := range [][]string{{"-j", "1"}, {"-j", "4"}, {"-j", "4", "-w"}} {
				args := []string{"genome", "search", "--quiet", "-d", db, "-S", "--extra", "-m", strconv.Itoa(masks), "-p", "12", "-n", "1", "--windows", "2"}
				args = append(args, mode...)
				// Repeated partial queries also exercise concurrent use and recycling.
				got := runSeedIndexCommand(t, append(args, query, query)...)
				lines := strings.Split(strings.TrimSpace(string(got)), "\n")
				if len(lines) != 9 { // header + both queries' four tied references
					t.Fatalf("expected all cutoff ties for both queries, got %s", got)
				}
				for _, line := range lines[1:] {
					fields := strings.Split(line, "\t")
					hits, err := strconv.Atoi(fields[4])
					if err != nil || hits == 0 {
						t.Fatalf("invalid hit count: %s", line)
					}
					score, err := strconv.Atoi(fields[5])
					if err != nil || score != hits*16 {
						t.Fatalf("repeated exact seeds must contribute once per mask: %s", line)
					}
					if fields[8] != strings.TrimSuffix(strings.Repeat("16,", hits), ",") {
						t.Fatalf("incorrect per-mask maxima: %s", line)
					}
				}
				if want == nil {
					want = got
				} else if !bytes.Equal(got, want) {
					t.Fatalf("screening changed with mode %v:\ngot %s\nwant %s", mode, got, want)
				}
			}
		})
	}
}
