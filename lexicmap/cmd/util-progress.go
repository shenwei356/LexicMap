package cmd

import (
	"time"

	"github.com/vbauerster/mpb/v8/decor"
)

// progressTimeDecorators shows elapsed time and ETA, then only the final duration.
// Each bar gets its own timer; mpb serializes calls to its decorators.
func progressTimeDecorators(age float64) []decor.Decorator {
	started := time.Now()
	var completedElapsed string
	return []decor.Decorator{
		decor.Any(func(s decor.Statistics) string {
			elapsed := time.Since(started).Truncate(time.Second).String()
			if s.Completed {
				// Capture the final duration once, including for short tasks.
				if completedElapsed == "" {
					completedElapsed = elapsed
				}
				return "done in " + completedElapsed
			}
			return "Elapsed: " + elapsed
		}),
		decor.OnComplete(decor.Name(", ETA: "), ""),
		decor.OnComplete(decor.EwmaETA(decor.ET_STYLE_GO, age), ""),
	}
}
