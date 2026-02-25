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
	"cmp"
	"math"
	"slices"
	"sync"
)

// ChainingOptions contains all options in chaining.
type ChainingOptions struct {
	MaxGap      float64
	MinLen      uint8
	MinScore    float64
	MaxDistance float64
	TopChains   int // only keep the top N chains
}

// DefaultChainingOptions is the defalt vaule of ChainingOption.
var DefaultChainingOptions = ChainingOptions{
	MaxGap:      5000,
	MinScore:    40,
	MaxDistance: 10000,
	TopChains:   -1,
}

// Chainer is an object for chaining the lexichash substrings between query and reference sequences.
type Chainer struct {
	options *ChainingOptions

	// scores        []float64 // actually, it's not necessary
	prevQIdx      []int
	maxscores     []float64
	maxscoresIdxs []int
	directions    []float64
	visited       []bool

	score2idx [][2]float64

	topChains int
}

// NewChainer creates a new chainer.
func NewChainer(options *ChainingOptions) *Chainer {
	c := &Chainer{
		options: options,

		// scores:        make([]float64, 0, 128),
		prevQIdx:      make([]int, 0, 10240),
		maxscores:     make([]float64, 0, 10240),
		maxscoresIdxs: make([]int, 0, 10240),
		directions:    make([]float64, 0, 10240),
		visited:       make([]bool, 0, 10240),

		score2idx: make([][2]float64, 0, 10240),

		topChains: options.TopChains,
	}
	return c
}

// RecycleChainingResult reycles the chaining results.
// Please remember to call this after using the results.
func RecycleChainingResult(chains *[]*[]int32) {
	for _, chain := range *chains {
		if chain != nil {
			*chain = (*chain)[:0]
			poolChain.Put(chain)
		}
	}

	*chains = (*chains)[:0]
	poolChains.Put(chains)
}

var poolChains = &sync.Pool{New: func() interface{} {
	tmp := make([]*[]int32, 0, 8)
	return &tmp
}}

var poolChain = &sync.Pool{New: func() interface{} {
	tmp := make([]int32, 0, 8)
	return &tmp
}}

// Chain finds the possible seed paths.
// Please remember to call RecycleChainingResult after using the results.
func (ce *Chainer) Chain(subs *[]*SubstrPair) (*[]*[]int32, float64) {
	n := len(*subs)

	if n == 1 { // for one seed, just check the seed weight
		paths := poolChains.Get().(*[]*[]int32)

		w := seedWeight(float64((*subs)[0].Len))
		if w >= ce.options.MinScore {
			path := poolChain.Get().(*[]int32)

			*path = append(*path, 0)

			*paths = append(*paths, path)
		}

		return paths, w
	}

	// minLen := ce.options.MinLen
	minScore := ce.options.MinScore

	var i, j, mj int

	// a list for storing triangular score matrix, the size is n*(n+1)>>1
	// scores := ce.scores[:0]
	// _n := n * (n + 1) >> 1
	// for k = 0; k < _n; k++ {
	// 	scores = append(scores, 0)
	// }
	// the maximum score for each seed, the size is n
	maxscores := &ce.maxscores
	*maxscores = (*maxscores)[:0]
	// index of previous seed, the size is n
	maxscoresIdxs := &ce.maxscoresIdxs
	*maxscoresIdxs = (*maxscoresIdxs)[:0]
	directions := &ce.directions
	*directions = (*directions)[:0]
	score2idx := ce.score2idx
	score2idx = score2idx[:0]
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
	*maxscores = append(*maxscores, seedWeight(float64((*subs)[0].Len)))
	*maxscoresIdxs = append(*maxscoresIdxs, 0)
	*directions = append(*directions, 0)
	score2idx = append(score2idx, [2]float64{(*maxscores)[0], 0})

	// compute scores
	var length int32 // substr/anchor length
	var s, m, w, g float64
	// var d float64
	var d int32
	var dir, mdir float64
	var a, b *SubstrPair
	maxGap := ce.options.MaxGap
	maxDistance := ce.options.MaxDistance
	maxDistanceInt32 := int32(maxDistance)
	maxDistanceUint32 := uint32(maxDistance)

	// store i of 'a' with a different QBegin
	prevQIdx := ce.prevQIdx[:0]
	lastDiffGroup := -1
	for i = 0; i < n; i++ {
		if i > 0 && (*subs)[i].QBegin != (*subs)[i-1].QBegin {
			lastDiffGroup = i - 1
		}
		prevQIdx = append(prevQIdx, lastDiffGroup)
	}

	var minJ, rightBound int
	var targetQ int32
	for i = 1; i < n; i++ {
		a = (*subs)[i]

		// fmt.Printf("i:%d/%d, a: %s\n", i, n, a)

		// just initialize the max score, which comes from the current seed
		// m = scores[k]
		// w = seedWeight(float64(a.Len))
		// m, mj, mdir = w, i, 0
		m, mj, mdir = seedWeight(float64(a.Len)), i, 0

		// for j = 0; j < i; j++ { // try all previous seeds, no bound
		//
		// j = i
		// for {
		// 	j--
		// 	if j < 0 {
		// 		break
		// 	}
		//
		rightBound = prevQIdx[i]
		if rightBound >= 0 {
			targetQ = a.QBegin - maxDistanceInt32

			// leftBound of j
			minJ, _ = slices.BinarySearchFunc((*subs)[:rightBound+1], targetQ, func(e *SubstrPair, t int32) int {
				if e.QBegin < t {
					return -1
				}
				if e.QBegin > t {
					return 1
				}
				return 0
			})

			for j = rightBound; j >= minJ; j-- {
				b = (*subs)[j]

				// fmt.Printf("  j:%d, b: %s\n", j, b)

				// avoided by using prevQIdx
				//
				// if a.QBegin == b.QBegin || a.TBegin == b.TBegin {
				// 	// fmt.Printf("   skip\n")
				// 	continue
				// }

				// avoided by using minJ
				//
				// if a.QBegin-b.QBegin > maxDistanceInt32 {
				// 	fmt.Printf("   distant in target too long: %d > %d\n", a.QBegin-b.QBegin, maxDistanceInt32)
				// 	break
				// }

				// d = distance(a, b)
				// if d > maxDistance {
				// fmt.Printf("   distant too long: %f > %f\n", d, maxDistance)
				//
				// if a.TBegin-b.TBegin > maxDistanceInt32 || b.TBegin-a.TBegin > maxDistanceInt32 {
				d = a.TBegin - b.TBegin
				// if d > maxDistanceInt32 || d < -maxDistanceInt32 || d == 0 {
				if d < 0 {
					d = -d
				}
				if uint32(d)-1 >= maxDistanceUint32 { // including d == 0
					// fmt.Printf("   distant in target too long: %d > %d\n", a.TBegin-b.TBegin, maxDistanceInt32)
					continue // must not be break
				}

				g = gap(a, b)
				if g > maxGap {
					// fmt.Printf("   gap too big: %f > %f\n", g, maxGap)
					continue // must not be break
				}

				// effective seed length
				if a.QBegin > b.QBegin+int32(b.Len) { // no overlap
					// b -----
					// a       ------
					length = int32(a.Len)
					w = seedWeight(float64(length))
				} else if g == 0 { // merge them into a longer anchor
					// b -----
					// a    ------
					length = a.QBegin + int32(a.Len) - int32(b.QBegin)
					w = -seedWeight(float64(b.Len)) + seedWeight(float64(length))
				} else {
					// b -----
					// a    ------
					length = a.QBegin + int32(a.Len) - (b.QBegin + int32(b.Len))
					w = seedWeight(float64(length))
				}

				// s = (*maxscores)[j] + seedWeight(float64(b.Len)) - distanceScore(d) - gapScore(g)

				dir = direction(a, b)

				if (*directions)[j] == 0 || (*directions)[j] == dir {
					s = (*maxscores)[j] + w - gapScore(g)
				} else {
					// fmt.Printf("   different directions: pre: %f, dir: %f\n", (*directions)[j], dir)
					// 	continue

					// the previous chain might be wrong in repetitive region.
					s = seedWeight(float64(b.Len)) + w - gapScore(g)
				}

				// fmt.Printf("    max: %d-%f-%f, s:%f-%f\n", mj, m, mdir, s, dir)

				if s >= minScore && s > m { //
					// fmt.Printf("   update score: %d-%f-%f -> %d-%f-%f\n", mj, m, mdir, j, s, dir)
					m = s
					mj = j
					mdir = dir
				}
			}
		}
		*maxscores = append(*maxscores, m) // save the max score
		*maxscoresIdxs = append(*maxscoresIdxs, mj)
		*directions = append(*directions, mdir)
		score2idx = append(score2idx, [2]float64{m, float64(i)})
	}
	// print the score matrix
	// fmt.Printf("i\tpair-i\tiMax\tj:scores\n")
	// for i = 0; i < n; i++ {
	// 	fmt.Printf("%d\t%s\t%d:%.0f:%.3f\n", i, (*subs)[i], (*maxscoresIdxs)[i], (*directions)[i], (*maxscores)[i])
	// }

	// backtrack
	visited := &ce.visited
	*visited = (*visited)[:0]
	for i = 0; i < n; i++ {
		*visited = append(*visited, false)
	}
	paths := poolChains.Get().(*[]*[]int32)

	var M float64
	var Mi int

	// memory inefficient
	// sort.Slice(score2idx, func(i, j int) bool {
	// 	return score2idx[i][0] > score2idx[j][0]
	// })
	slices.SortFunc(score2idx, func(a, b [2]float64) int {
		return cmp.Compare[float64](b[0], a[0])
	})

	var iMaxScore int

	// for computing the score for sorting
	var maxScore float64
	first := true // best chain
	var changeDirection bool
	var nChecked int
	for {
		nChecked++
		if ce.topChains > 0 && nChecked > ce.topChains {
			break
		}

		// find the next highest score

		// M = 0
		// for i, m = range *maxscores {
		// 	if (*visited)[i] {
		// 		continue
		// 	}
		// 	if m > M {
		// 		M, Mi = m, i
		// 	}
		// }

		M = 0
		for iMaxScore < n {
			M = 0
			Mi = int(score2idx[iMaxScore][1])
			if !(*visited)[Mi] {
				M = score2idx[iMaxScore][0]

				iMaxScore++
				break
			}

			iMaxScore++
		}

		if M < minScore { // no valid anchors
			break
		}

		path := poolChain.Get().(*[]int32)

		// fmt.Printf("max: Mi:%d(%d), %f\n", Mi, n, M)

		i = Mi
		if first {
			maxScore = M
			first = false
		}

		for {
			j = (*maxscoresIdxs)[i] // previous anchor

			changeDirection = (i != j && (*directions)[j] != 0 && (*directions)[i] != (*directions)[j])

			// fmt.Printf(" i:%d, visited:%v; j:%d, visited:%v\n", i, (*visited)[i], j, (*visited)[j])
			if (*visited)[j] && !changeDirection { // current anchor is abandoned
				// if len(*path) == 0 && !(*visited)[i] && (*subs)[i].Len >= minLen {
				// 	*path = append(*path, i) // record the anchor
				// 	// fmt.Printf(" orphan from %d, %s\n", i, (*subs)[i])
				// }

				// if len(*path) > 0 {
				// 	// but don't forget already added path
				// 	reverseInts(*path)
				// 	*paths = append(*paths, path)
				// 	// fmt.Printf("  stop at %d, %s\n", i, (*subs)[i])

				// 	path = poolChain.Get().(*[]int)
				// }

				*path = (*path)[:0]
				(*visited)[i] = true // do not check it again

				break
			}

			*path = append(*path, int32(i)) // record the anchor
			(*visited)[i] = true            // mark as visited

			// if firstAnchor {
			// fmt.Printf(" start from %d, %s\n", i, (*subs)[i])
			// firstAnchor = false
			// }
			// else {
			// fmt.Printf("  add %d, %s\n", i, (*subs)[i])
			// }
			if i == j || changeDirection { // the path starts here
				if changeDirection {
					*path = append(*path, int32(j))
				}

				reverseInt32s(*path)
				*paths = append(*paths, path)
				// fmt.Printf("  stop at %d, %s\n", i, (*subs)[i])

				break
			} else {
				i = j
			}
		}

	}

	// fmt.Println(maxScore)
	return paths, maxScore
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

func direction(a, b *SubstrPair) float64 {
	if a.TBegin >= b.TBegin {
		return 1
	}
	return -1
}

// different in plus:plus and plus:minus match
func gap(a, b *SubstrPair) float64 {
	if a.TBegin >= b.TBegin {
		return math.Abs(math.Abs(float64(a.QBegin-b.QBegin)) - math.Abs(float64(a.TBegin-b.TBegin)))
	}
	return math.Abs(math.Abs(float64(a.QBegin-b.QBegin)) - math.Abs(float64(a.TBegin+int32(a.Len)-b.TBegin-int32(b.Len))))
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

func reverseInt32s(s []int32) {
	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		s[i], s[j] = s[j], s[i]
	}
}
