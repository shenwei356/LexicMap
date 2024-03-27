// Copyright © 2018-2021 Wei Shen <shenwei356@gmail.com>
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

package util

import (
	"math/rand"
	"testing"
)

var testsUint64 [][2]uint64
var testsUint32 [][4]uint32

func init() {
	ntests := 10000
	testsUint64 = make([][2]uint64, ntests)
	testsUint32 = make([][4]uint32, ntests)
	var i int
	for ; i < ntests/4; i++ {
		testsUint64[i] = [2]uint64{rand.Uint64(), rand.Uint64()}
		testsUint32[i] = [4]uint32{rand.Uint32(), rand.Uint32(), rand.Uint32(), rand.Uint32()}
	}
	for ; i < ntests/2; i++ {
		testsUint64[i] = [2]uint64{uint64(rand.Uint32()), uint64(rand.Uint32())}
		testsUint32[i] = [4]uint32{rand.Uint32(), rand.Uint32(), rand.Uint32(), rand.Uint32()}
	}
	for ; i < ntests*3/4; i++ {
		testsUint64[i] = [2]uint64{uint64(rand.Intn(65536)), uint64(rand.Intn(256))}
		testsUint32[i] = [4]uint32{uint32(rand.Intn(65536)), uint32(rand.Intn(256)), uint32(rand.Intn(65536)), uint32(rand.Intn(256))}
	}
	for ; i < ntests; i++ {
		testsUint64[i] = [2]uint64{uint64(rand.Intn(256)), uint64(rand.Intn(256))}
		testsUint32[i] = [4]uint32{uint32(rand.Intn(256)), uint32(rand.Intn(256)), uint32(rand.Intn(256)), uint32(rand.Intn(256))}
	}
}

func TestStreamVByte64(t *testing.T) {
	buf := make([]byte, 16)
	var ctrl byte
	var n, n2 int
	var v1, v2 uint64
	for i, test := range testsUint64 {
		ctrl, n = PutUint64s(buf, test[0], test[1])
		if CtrlByte2ByteLengthsUint64(ctrl) != n {
			t.Errorf("#%d, wrong byte length", i)
		}

		v1, v2, n2 = Uint64s(ctrl, buf[0:n])
		if n2 == 0 {
			t.Errorf("#%d, wrong decoded number", i)
		}

		if v1 != test[0] || v2 != test[1] {
			t.Errorf("#%d, wrong decoded result: %d, %d, answer: %d, %d", i, v1, v2, test[0], test[1])
		}
		// fmt.Printf("%d, %d => n=%d, buf=%v\n", test[0], test[1], n, buf[0:n])
	}
}

func TestStreamVByte32(t *testing.T) {
	buf := make([]byte, 16)
	var ctrl byte
	var n, n2 int
	var v1, v2, v3, v4 uint32
	for i, test := range testsUint32 {
		ctrl, n = PutUint32s(buf, test[0], test[1], test[2], test[3])
		if CtrlByte2ByteLengthsUint32(ctrl) != n {
			t.Errorf("#%d, wrong byte length", i)
		}

		v1, v2, v3, v4, n2 = Uint32s(ctrl, buf[0:n])
		if n2 == 0 {
			t.Errorf("#%d, wrong decoded number", i)
		}

		if v1 != test[0] || v2 != test[1] || v3 != test[2] || v4 != test[3] {
			t.Errorf("#%d, wrong decoded result: %d, %d, %d, %d, answer: %d, %d, %d, %d", i, v1, v2, v3, v4, test[0], test[1], test[2], test[3])
		}
		// fmt.Printf("%d, %d => n=%d, buf=%v\n", test[0], test[1], n, buf[0:n])
	}
}

var _v1, _v2 uint64

func BenchmarkUint64s(b *testing.B) {
	buf := make([]byte, 16)
	var ctrl byte
	var n, n2 int
	var v1, v2 uint64
	for i := 0; i < b.N; i++ {
		for i, test := range testsUint64 {
			ctrl, n = PutUint64s(buf, test[0], test[1])

			v1, v2, n2 = Uint64s(ctrl, buf[0:n])
			if n2 == 0 {
				b.Errorf("#%d, wrong decoded number", i)
			}

			if v1 != test[0] || v2 != test[1] {
				b.Errorf("#%d, wrong decoded result: %d, %d, answer: %d, %d", i, v1, v2, test[0], test[1])
			}
			// fmt.Printf("%d, %d => n=%d, buf=%v\n", test[0], test[1], n, buf[0:n])
		}
	}
	_v1, _v2 = v1, v2
}

var __v1, __v2, __v3, __v4 uint32

func BenchmarkUint32s(b *testing.B) {
	buf := make([]byte, 16)
	var ctrl byte
	var n, n2 int
	var v1, v2, v3, v4 uint32
	for i := 0; i < b.N; i++ {
		for i, test := range testsUint32 {
			ctrl, n = PutUint32s(buf, test[0], test[1], test[2], test[3])

			v1, v2, v3, v4, n2 = Uint32s(ctrl, buf[0:n])
			if n2 == 0 {
				b.Errorf("#%d, wrong decoded number", i)
			}

			if v1 != test[0] || v2 != test[1] || v3 != test[2] || v4 != test[3] {
				b.Errorf("#%d, wrong decoded result: %d, %d, %d, %d, answer: %d, %d, %d, %d", i, v1, v2, v3, v4, test[0], test[1], test[2], test[3])
			}
			// fmt.Printf("%d, %d => n=%d, buf=%v\n", test[0], test[1], n, buf[0:n])
		}
	}
	__v1, __v2, __v3, __v4 = v1, v2, v3, v4
}
