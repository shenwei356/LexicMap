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
