// Copyright © 2023-2024 Wei Shen <shenwei356@gmail.com>
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
	"math"
	"sync"
)

// ChainingOptions contains all options in chaining.
type ChainingOptions struct {
	MaxGap      float64
	MinScore    float64
	MaxDistance float64
}

// DefaultChainingOptions is the defalt vaule of ChainingOption.
var DefaultChainingOptions = ChainingOptions{
	MaxGap:      5000,
	MinScore:    40,
	MaxDistance: 10000,
}

// Chainer is an object for chaining the lexichash substrings between query and reference sequences.
type Chainer struct {
	options *ChainingOptions

	// scores        []float64 // actually, it's not necessary
	maxscores     []float64
	maxscoresIdxs []int
	visited       []bool
}

// NewChainer creates a new chainer.
func NewChainer(options *ChainingOptions) *Chainer {
	c := &Chainer{
		options: options,

		// scores:        make([]float64, 0, 128),
		maxscores:     make([]float64, 0, 10240),
		maxscoresIdxs: make([]int, 0, 10240),
		visited:       make([]bool, 0, 10240),
	}
	return c
}

// RecycleChainingResult reycles the chaining results.
// Please remember to call this after using the results.
func RecycleChainingResult(chains *[]*[]int) {
	for _, chain := range *chains {
		poolChain.Put(chain)
	}
	poolChains.Put(chains)
}

var poolChains = &sync.Pool{New: func() interface{} {
	tmp := make([]*[]int, 0, 8)
	return &tmp
}}

var poolChain = &sync.Pool{New: func() interface{} {
	tmp := make([]int, 0, 32)
	return &tmp
}}

// Chain finds the possible seed paths.
// Please remember to call RecycleChainingResult after using the results.
func (ce *Chainer) Chain(subs *[]*SubstrPair) (*[]*[]int, float64) {
	n := len(*subs)

	var sumMaxScore float64

	if n == 1 { // for one seed, just check the seed weight
		paths := poolChains.Get().(*[]*[]int)
		*paths = (*paths)[:0]

		w := seedWeight(float64((*subs)[0].Len))
		if w >= ce.options.MinScore {
			path := poolChain.Get().(*[]int)
			*path = (*path)[:0]

			*path = append(*path, 0)

			*paths = append(*paths, path)
		}

		return paths, w
	}

	var i, j, mj int

	// a list for storing triangular score matrix, the size is n*(n+1)>>1
	// scores := ce.scores[:0]
	// _n := n * (n + 1) >> 1
	// for k = 0; k < _n; k++ {
	// 	scores = append(scores, 0)
	// }
	// the maximum score for each seed, the size is n
	maxscores := ce.maxscores[:0]
	// index of previous seed, the size is n
	maxscoresIdxs := ce.maxscoresIdxs[:0]
	// for i = 0; i < n; i++ {
	// 	maxscores = append(maxscores, 0)
	// 	maxscoresIdxs = append(maxscoresIdxs, 0)
	// }
	// initialize
	// for i, b := range *subs { // j == i, means a path starting from this seed
	// 	j0 = i * (i + 1) >> 1
	// 	k = j0 + i
	// 	//scores[k] = seedWeight(float64(b.Len))
	// }
	maxscores = append(maxscores, seedWeight(float64((*subs)[0].Len)))
	maxscoresIdxs = append(maxscoresIdxs, 0)

	// compute scores
	var s, m, d, g float64
	var a, b *SubstrPair
	maxGap := ce.options.MaxGap
	maxDistance := ce.options.MaxDistance
	for i = 1; i < n; i++ {
		a = (*subs)[i]

		// just initialize the max score, which comes from the current seed
		// m = scores[k]
		m, mj = seedWeight(float64(a.Len)), i

		for j = 0; j < i; j++ { // try all previous seeds, no bound
			b = (*subs)[j]

			d = distance(a, b)
			if d > maxDistance {
				continue
			}

			g = gap(a, b)
			if g > maxGap {
				continue
			}

			s = maxscores[j] + seedWeight(float64(b.Len)) - distanceScore(d) - gapScore(g)
			// scores[k] = s

			if s >= m { // update the max score
				m = s
				mj = j
			}
		}
		maxscores = append(maxscores, m) // save the max score
		maxscoresIdxs = append(maxscoresIdxs, mj)
	}
	// print the score matrix
	// fmt.Printf("i\tiMax\tscores\n")
	// for i = 0; i < n; i++ {
	// 	fmt.Printf("%d\t%s\t%d", i, (*subs)[i], maxscoresIdxs[i])
	// 	for j = 0; j <= i; j++ {
	// 		k = i*(i+1)/2 + j
	// 		if scores[k] > 0 {
	// 			fmt.Printf("\t%d:%6.2f", j, scores[k])
	// 		}
	// 	}
	// 	fmt.Printf("\n")
	// }

	// backtrack
	visited := ce.visited[:0]
	for i = 0; i < n; i++ {
		visited = append(visited, false)
	}
	paths := poolChains.Get().(*[]*[]int)
	*paths = (*paths)[:0]
	var first bool
	minScore := ce.options.MinScore

	path := poolChain.Get().(*[]int)
	*path = (*path)[:0]

	var M float64
	var Mi int

	for {
		// find the next highest score
		M = 0
		for i, m = range maxscores {
			if visited[i] {
				continue
			}
			if m > M {
				M, Mi = m, i
			}
		}
		if M < minScore { // no valid anchors
			break
		}

		i = Mi
		first = true
		for {
			j = maxscoresIdxs[i] // previous anchor
			if visited[j] {      // current anchor is abandoned
				if len(*path) > 0 { // but don't forget already added path
					// reverseInts(*path)
					// *paths = append(*paths, path)
					// // fmt.Printf("stop at %d, %s\n", i, (*subs)[i])

					// path = poolChain.Get().(*[]int)
					// *path = (*path)[:0]
					*path = (*path)[:0]
				}
				visited[i] = true // do not check it again
				break
			}

			*path = append(*path, i) // record the anchor
			visited[i] = true        // mark as visited
			if first {
				// fmt.Printf("start from %d, %s\n", i, (*subs)[i])
				sumMaxScore += maxscores[i]
				first = false
			}
			// else {
			// 	fmt.Printf("  add %d, %s\n", i, (*subs)[i])
			// }
			if i != j {
				i = j
			} else { // the path starts here
				reverseInts(*path)
				*paths = append(*paths, path)
				// fmt.Printf("stop at %d, %s\n", i, (*subs)[i])

				path = poolChain.Get().(*[]int)
				*path = (*path)[:0]

				break
			}
		}
	}

	return paths, sumMaxScore
}

func seedWeight(l float64) float64 {
	return 0.1 * l * l
}

func distance(a, b *SubstrPair) float64 {
	return math.Max(math.Abs(float64(a.QBegin-b.QBegin)), math.Abs(float64(a.TBegin-b.TBegin)))
}

func distanceScore(d float64) float64 {
	return 0.01 * d
}

func gap(a, b *SubstrPair) float64 {
	return math.Abs(math.Abs(float64(a.QBegin-b.QBegin)) - math.Abs(float64(a.TBegin-b.TBegin)))
}

func gapScore(gap float64) float64 {
	if gap == 0 {
		return 0
	}
	return 0.1*gap + 0.5*math.Log2(gap)
}

func reverseInts(s []int) {
	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		s[i], s[j] = s[j], s[i]
	}
}
