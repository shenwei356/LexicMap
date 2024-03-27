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
var offsetsUint32 = []uint8{24, 16, 8, 0}

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

// Uint64sOld decodes encoded bytes. Using a lookup list is slower.
func Uint64sOld(ctrl byte, buf []byte) (values [2]uint64, n int) {
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

// PutUint32s encodes four uint32s into 4-16 bytes, and returns control byte
// and encoded byte length.
func PutUint32s(buf []byte, v1, v2, v3, v4 uint32) (ctrl byte, n int) {
	blen := ByteLengthUint32(v1)
	ctrl |= byte(blen - 1)
	for _, offset := range offsetsUint32[4-blen:] {
		buf[n] = byte((v1 >> offset) & 0xff)
		n++
	}

	ctrl <<= 2
	blen = ByteLengthUint32(v2)
	ctrl |= byte(blen - 1)
	for _, offset := range offsetsUint32[4-blen:] {
		buf[n] = byte((v2 >> offset) & 0xff)
		n++
	}

	ctrl <<= 4
	blen = ByteLengthUint32(v3)
	ctrl |= byte(blen - 1)
	for _, offset := range offsetsUint32[4-blen:] {
		buf[n] = byte((v3 >> offset) & 0xff)
		n++
	}

	ctrl <<= 6
	blen = ByteLengthUint32(v4)
	ctrl |= byte(blen - 1)
	for _, offset := range offsetsUint32[4-blen:] {
		buf[n] = byte((v4 >> offset) & 0xff)
		n++
	}

	return
}

// Uint64s decodes encoded bytes.
func Uint64s(ctrl byte, buf []byte) (values [2]uint64, n int) {
	blen1 := int((ctrl>>3)&7) + 1
	blen2 := int(ctrl&7) + 1
	if len(buf) < blen1+blen2 {
		return values, 0
	}

	var j int

	for j = 0; j < blen1; j++ {
		values[0] <<= 8
		values[0] |= uint64(buf[n])
		n++
	}

	for j = 0; j < blen2; j++ {
		values[1] <<= 8
		values[1] |= uint64(buf[n])
		n++
	}

	return
}

// Uint64s decodes encoded bytes.
func Uint64s2(ctrl byte, buf []byte) (v1, v2 uint64, n int) {
	blen1 := int((ctrl>>3)&7) + 1
	blen2 := int(ctrl&7) + 1
	if len(buf) < blen1+blen2 {
		return 0, 0, 0
	}

	var j int

	for j = 0; j < blen1; j++ {
		v1 <<= 8
		v1 |= uint64(buf[n])
		n++
	}

	for j = 0; j < blen2; j++ {
		v2 <<= 8
		v2 |= uint64(buf[n])
		n++
	}

	return
}

// Uint32s decodes encoded bytes.
func Uint32s(ctrl byte, buf []byte) (values [4]uint32, n int) {
	blen1 := int((ctrl>>6)&3) + 1
	blen2 := int((ctrl>>4)&3) + 1
	blen3 := int((ctrl>>2)&3) + 1
	blen4 := int(ctrl&3) + 1
	if len(buf) < blen1+blen2+blen3+blen4 {
		return values, 0
	}

	var j int

	for j = 0; j < blen1; j++ {
		values[0] <<= 8
		values[0] |= uint32(buf[n])
		n++
	}

	for j = 0; j < blen2; j++ {
		values[1] <<= 8
		values[1] |= uint32(buf[n])
		n++
	}

	for j = 0; j < blen3; j++ {
		values[2] <<= 8
		values[2] |= uint32(buf[n])
		n++
	}

	for j = 0; j < blen4; j++ {
		values[3] <<= 8
		values[3] |= uint32(buf[n])
		n++
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

// ByteLengthUint32 returns the minimum number of bytes to store a integer.
func ByteLengthUint32(n uint32) uint8 {
	if n < 256 {
		return 1
	}
	if n < 65536 {
		return 2
	}
	if n < 16777216 {
		return 3
	}
	return 4
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
