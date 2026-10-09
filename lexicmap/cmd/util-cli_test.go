package cmd

import (
	"bytes"
	"strings"
	"testing"
	"time"

	"github.com/shenwei356/go-logging"
)

// TestLogStartupInfoDate checks that the local date follows the command line.
func TestLogStartupInfoDate(t *testing.T) {
	var output bytes.Buffer
	previousLog := log
	log = logging.MustGetLogger("startup-test")
	backend := logging.NewBackendFormatter(logging.NewLogBackend(&output, "", 0),
		logging.MustStringFormatter("%{message}"))
	log.SetBackend(logging.AddModuleLevel(backend))
	t.Cleanup(func() { log = previousLog })

	before := time.Now().Format("2006-01-02")
	logStartupInfo()
	after := time.Now().Format("2006-01-02")
	lines := strings.Split(output.String(), "\n")
	for i, line := range lines {
		if !strings.HasPrefix(line, " CMD: ") {
			continue
		}
		if i+1 >= len(lines) {
			t.Fatal("missing date after CMD")
		}
		// Accept either date if logging crosses local midnight.
		if date := lines[i+1]; date != "DATE: "+before && date != "DATE: "+after {
			t.Fatalf("unexpected line after CMD: %q", date)
		}
		return
	}
	t.Fatal("missing CMD in startup log")
}
