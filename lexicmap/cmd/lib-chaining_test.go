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
	"testing"
)

func TestChaining(t *testing.T) {
	/* command to prepare seeds from a certain query
	 cat t.txt | csvtk grep -t -f target -p GCA_013693855.1  \
	 	| csvtk cut -t -f qstart,tstart | sed 1d \
		| awk '{print "{QBegin: "$1", TBegin: "$2", Len: 31},"}'
	*/
	subs := []*SubstrPair{
		// two sequences on different strands
		// {QBegin: 18, TBegin: 3453, Len: 31},
		// {QBegin: 18, TBegin: 3640464, Len: 31},
		// {QBegin: 1924, TBegin: 1547, Len: 31},
		// {QBegin: 1924, TBegin: 3638544, Len: 31},

		// not perfect in this case, there are two chains: 0,1 and 2., while it should be one.
		{QBegin: 552, TBegin: 3798905, Len: 17},
		{QBegin: 667, TBegin: 3799019, Len: 15},
		{QBegin: 1332, TBegin: 3799686, Len: 31},

		// a kmer has multiple matches
		{QBegin: 1384, TBegin: 628584, Len: 31},
		{QBegin: 1490, TBegin: 628690, Len: 31},
		{QBegin: 1879, TBegin: 900465, Len: 31},
		{QBegin: 1879, TBegin: 629079, Len: 31},
		{QBegin: 1879, TBegin: 627005, Len: 31},
		{QBegin: 1910, TBegin: 6123921, Len: 23},

		// same strands

		{QBegin: 182, TBegin: 1282695, Len: 26},
		{QBegin: 182, TBegin: 1769573, Len: 26},
		{QBegin: 315, TBegin: 1282830, Len: 15},
		{QBegin: 315, TBegin: 1769708, Len: 15},
		{QBegin: 343, TBegin: 1769724, Len: 27},

		{QBegin: 10, TBegin: 314159, Len: 20},

		// this case is kept in the chainning step,
		// because we can not simply limit
		// the minimum distance between two anchors.
		{QBegin: 60, TBegin: 14234, Len: 15},
		{QBegin: 61, TBegin: 14235, Len: 15},

		{QBegin: 60, TBegin: 3395374, Len: 15},
		{QBegin: 70, TBegin: 3395384, Len: 15},

		{QBegin: 50, TBegin: 950, Len: 31},
		{QBegin: 79, TBegin: 3637976, Len: 31},
		{QBegin: 100, TBegin: 3637997, Len: 31},
		{QBegin: 519, TBegin: 1419, Len: 31},
		{QBegin: 550, TBegin: 3638447, Len: 31},
		{QBegin: 647, TBegin: 3638544, Len: 31},

		{QBegin: 111, TBegin: 1146311, Len: 31},
		{QBegin: 136, TBegin: 1146336, Len: 31},
		{QBegin: 138, TBegin: 1146338, Len: 31},
		{QBegin: 139, TBegin: 1146339, Len: 31},
		{QBegin: 264, TBegin: 1146464, Len: 31},
		{QBegin: 1479, TBegin: 1147679, Len: 31},
		{QBegin: 1484, TBegin: 1147684, Len: 31},
		{QBegin: 1543, TBegin: 1147743, Len: 31},
		{QBegin: 1566, TBegin: 1147766, Len: 31},
		{QBegin: 1919, TBegin: 1148119, Len: 31},
	}
	tmp := []*SearchResult{
		{
			Subs: &subs,
		},
	}
	rs := &tmp

	cf := &DefaultChainingOptions

	chainer := NewChainer(cf)
	for _, r := range *rs {
		paths, sumMaxScore := chainer.Chain(r.Subs)

		t.Logf("sum score: %f, paths:\n", sumMaxScore)
		for _, p := range *paths {
			t.Logf("  %d\n", *p)
		}

		RecycleChainingResult(paths)
	}
}
