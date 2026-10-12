package cmd

import (
	"math/rand"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// TestGenomeANIOutputPrecision checks ANI and AF precision through the public
// search and compare commands, including OrthoANI and quality-adjusted ANI.
func TestGenomeANIOutputPrecision(t *testing.T) {
	dir := t.TempDir()
	rng := rand.New(rand.NewSource(1))
	sequence := make([]byte, 4080)
	for i := range sequence {
		sequence[i] = "ACGT"[rng.Intn(4)]
	}
	query := filepath.Join(dir, "query.fa")
	reference := filepath.Join(dir, "reference.fa")
	fastq := filepath.Join(dir, "query.fq")
	if err := os.WriteFile(query, append([]byte(">seq\n"), sequence...), 0600); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(fastq, []byte("@seq\n"+string(sequence)+"\n+\n"+strings.Repeat("I", len(sequence))+"\n"), 0600); err != nil {
		t.Fatal(err)
	}
	// A substitution gives a non-integer ANI so the test also covers rounding.
	sequence[1040] = "ACGT"[(strings.IndexByte("ACGT", sequence[1040])+1)%4]
	if err := os.WriteFile(reference, append([]byte(">seq\n"), sequence...), 0600); err != nil {
		t.Fatal(err)
	}
	db := filepath.Join(dir, "index.lmi")
	runSeedIndexCommand(t, "index", "--quiet", "-j", "1", "-c", "2", "-m", "64", "-k", "16", "--partitions", "4", "-O", db, reference)
	for _, command := range []struct {
		name string
		args []string
	}{
		{"search", []string{"genome", "search", "-d", db, "-m", "0", "-p", "12", query}},
		{"search OrthoANI", []string{"genome", "search", "-O", "-d", db, "-m", "0", "-p", "12", query}},
		{"search adjusted ANI", []string{"genome", "search", "--adjust-ani-by-quality", "-d", db, "-m", "0", "-p", "12", fastq}},
		{"compare", []string{"genome", "compare", query, reference}},
		{"compare OrthoANI", []string{"genome", "compare", "-O", query, reference}},
	} {
		t.Run(command.name, func(t *testing.T) {
			output := runSeedIndexCommand(t, append(command.args, "--quiet", "-j", "1")...)
			lines := strings.Split(strings.TrimSuffix(string(output), "\n"), "\n")
			if len(lines) < 2 {
				t.Fatalf("expected genome matches, got %s", output)
			}
			header := strings.Split(lines[0], "\t")
			for _, line := range lines[1:] {
				fields := strings.Split(line, "\t")
				if len(fields) != len(header) {
					t.Fatalf("unexpected column count: %s", line)
				}
				for i, name := range header {
					decimals := 0
					if strings.Contains(name, "ANI") {
						decimals = 4
					} else if strings.Contains(name, "AF") {
						decimals = 3
					} else {
						continue
					}
					parts := strings.Split(fields[i], ".")
					if len(parts) != 2 || len(parts[1]) != decimals {
						t.Fatalf("%s: got %q, want %d decimal places", name, fields[i], decimals)
					}
				}
			}
		})
	}
}
