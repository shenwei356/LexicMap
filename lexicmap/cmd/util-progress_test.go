package cmd

import (
	"strings"
	"testing"
	"testing/synctest"
	"time"

	"github.com/vbauerster/mpb/v8/decor"
)

// progressTimeText renders the complete suffix without requiring a terminal.
func progressTimeText(ds []decor.Decorator, s decor.Statistics) string {
	var text strings.Builder
	for _, d := range ds {
		value, _ := d.Decor(s)
		text.WriteString(value)
	}
	return text.String()
}

func TestProgressTimeDecorators(t *testing.T) {
	t.Run("running-and-completed", func(t *testing.T) {
		synctest.Test(t, func(t *testing.T) {
			ds := progressTimeDecorators(10)
			check := func(s decor.Statistics, want string) {
				t.Helper()
				if got := progressTimeText(ds, s); got != want {
					t.Fatalf("progress suffix = %q, want %q", got, want)
				}
			}
			check(decor.Statistics{}, "Elapsed: 0s, ETA: 0s")
			time.Sleep(135*time.Second + 500*time.Millisecond)
			check(decor.Statistics{}, "Elapsed: 2m15s, ETA: 0s")
			// Completion must include time since the last running refresh.
			time.Sleep(45 * time.Second)
			check(decor.Statistics{Completed: true}, "done in 3m0s")
			time.Sleep(time.Hour)
			check(decor.Statistics{Completed: true}, "done in 3m0s")
		})
	})
	t.Run("completed-before-first-refresh", func(t *testing.T) {
		synctest.Test(t, func(t *testing.T) {
			for _, duration := range []time.Duration{0, 500 * time.Millisecond, time.Hour + time.Second} {
				ds := progressTimeDecorators(10)
				time.Sleep(duration)
				want := "done in " + duration.Truncate(time.Second).String()
				if got := progressTimeText(ds, decor.Statistics{Completed: true}); got != want {
					t.Fatalf("first refresh after %s = %q, want %q", duration, got, want)
				}
			}
		})
	})
	t.Run("independent-timers", func(t *testing.T) {
		synctest.Test(t, func(t *testing.T) {
			first := progressTimeDecorators(3)
			time.Sleep(time.Minute)
			second := progressTimeDecorators(1024)
			time.Sleep(time.Minute)
			if got := progressTimeText(first, decor.Statistics{Completed: true}); got != "done in 2m0s" {
				t.Fatal(got)
			}
			if got := progressTimeText(second, decor.Statistics{Completed: true}); got != "done in 1m0s" {
				t.Fatal(got)
			}
		})
	})
	t.Run("ewma-eta", func(t *testing.T) {
		for _, age := range []float64{3, 5, 10, 64, 1024} {
			ds := progressTimeDecorators(age)
			eta := ds[2].(decor.Wrapper).Unwrap()
			baseline := decor.EwmaETA(decor.ET_STYLE_GO, age)
			for i := 1; i <= 64; i++ {
				duration := time.Duration(i) * time.Second
				eta.(decor.EwmaDecorator).EwmaUpdate(1, duration)
				baseline.(decor.EwmaDecorator).EwmaUpdate(1, duration)
			}
			stats := decor.Statistics{Total: 100, Current: 64}
			want, _ := baseline.Decor(stats)
			if got, _ := ds[2].Decor(stats); got != want || got == "0s" {
				t.Fatalf("ETA with age %g = %q, want %q", age, got, want)
			}
		}
	})
}
