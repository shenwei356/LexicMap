// Copyright © 2018-2024 Wei Shen <shenwei356@gmail.com>
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

var offsetsUint64 = []uint8{56, 48, 40, 32, 24, 16, 8, 0}

// PutUint64s encodes two uint64s into 2-16 bytes, and returns control byte
// and encoded byte length.
func PutUint64s(buf []byte, v1, v2 uint64) (ctrl byte, n int) {
	blen := ByteLengthUint64(v1)
	ctrl |= byte(blen - 1)
	for _, offset := range offsetsUint64[8-blen:] {
		buf[n] = byte((v1 >> offset) & 0xff)
		n++
	}

	ctrl <<= 3
	blen = ByteLengthUint64(v2)
	ctrl |= byte(blen - 1)
	for _, offset := range offsetsUint64[8-blen:] {
		buf[n] = byte((v2 >> offset) & 0xff)
		n++
	}
	return
}

// Uint64s decodes encoded bytes.
func Uint64s(ctrl byte, buf []byte) (values [2]uint64, n int) {
	blens := CtrlByte2ByteLengthsUint64[ctrl]
	if len(buf) < int(blens[0]+blens[1]) {
		return values, 0
	}
	for i := 0; i < 2; i++ {
		for j := uint8(0); j < blens[i]; j++ {
			values[i] <<= 8
			values[i] |= uint64(buf[n])
			n++
		}
	}

	return
}

// ByteLengthUint64 returns the minimum number of bytes to store a integer.
func ByteLengthUint64(n uint64) uint8 {
	if n < 256 {
		return 1
	}
	if n < 65536 {
		return 2
	}
	if n < 16777216 {
		return 3
	}
	if n < 4294967296 {
		return 4
	}
	if n < 1099511627776 {
		return 5
	}
	if n < 281474976710656 {
		return 6
	}
	if n < 72057594037927936 {
		return 7
	}
	return 8
}

// CtrlByte2ByteLengthsUint64 is a table for query byte lenghts from the control byte.
var CtrlByte2ByteLengthsUint64 = [64][2]uint8{
	{1, 1}, // 0, 0b000000
	{1, 2},
	{1, 3},
	{1, 4},
	{1, 5},
	{1, 6},
	{1, 7},
	{1, 8},
	{2, 1}, // 8, 0b001000
	{2, 2},
	{2, 3},
	{2, 4},
	{2, 5},
	{2, 6},
	{2, 7},
	{2, 8},
	{3, 1}, // 16, 0b010000
	{3, 2},
	{3, 3},
	{3, 4},
	{3, 5},
	{3, 6},
	{3, 7},
	{3, 8},
	{4, 1}, // 24, 0b011000
	{4, 2},
	{4, 3},
	{4, 4},
	{4, 5},
	{4, 6},
	{4, 7},
	{4, 8},
	{5, 1}, // 32, 0b100000
	{5, 2},
	{5, 3},
	{5, 4},
	{5, 5},
	{5, 6},
	{5, 7},
	{5, 8},
	{6, 1}, // 40, 0b101000
	{6, 2},
	{6, 3},
	{6, 4},
	{6, 5},
	{6, 6},
	{6, 7},
	{6, 8},
	{7, 1}, // 48, 0b110000
	{7, 2},
	{7, 3},
	{7, 4},
	{7, 5},
	{7, 6},
	{7, 7},
	{7, 8},
	{8, 1}, // 56, 0b111000
	{8, 2},
	{8, 3},
	{8, 4},
	{8, 5},
	{8, 6},
	{8, 7},
	{8, 8},
}
