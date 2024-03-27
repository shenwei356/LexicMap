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

	ctrl <<= 2
	blen = ByteLengthUint32(v3)
	ctrl |= byte(blen - 1)
	for _, offset := range offsetsUint32[4-blen:] {
		buf[n] = byte((v3 >> offset) & 0xff)
		n++
	}

	ctrl <<= 2
	blen = ByteLengthUint32(v4)
	ctrl |= byte(blen - 1)
	for _, offset := range offsetsUint32[4-blen:] {
		buf[n] = byte((v4 >> offset) & 0xff)
		n++
	}

	return
}

// Uint64s decodes encoded bytes.
func Uint64s(ctrl byte, buf []byte) (v1, v2 uint64, n int) {
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
func Uint32s(ctrl byte, buf []byte) (v1, v2, v3, v4 uint32, n int) {
	blen1 := int((ctrl>>6)&3) + 1
	blen2 := int((ctrl>>4)&3) + 1
	blen3 := int((ctrl>>2)&3) + 1
	blen4 := int(ctrl&3) + 1
	if len(buf) < blen1+blen2+blen3+blen4 {
		return 0, 0, 0, 0, 0
	}

	var j int

	for j = 0; j < blen1; j++ {
		v1 <<= 8
		v1 |= uint32(buf[n])
		n++
	}

	for j = 0; j < blen2; j++ {
		v2 <<= 8
		v2 |= uint32(buf[n])
		n++
	}

	for j = 0; j < blen3; j++ {
		v3 <<= 8
		v3 |= uint32(buf[n])
		n++
	}

	for j = 0; j < blen4; j++ {
		v4 <<= 8
		v4 |= uint32(buf[n])
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

// CtrlByte2ByteLengthsUint64 returns the byte length for a given control byte.
func CtrlByte2ByteLengthsUint64(ctrl byte) int {
	return int(ctrl>>3&7+ctrl&7) + 2
}

// CtrlByte2ByteLengthsUint32 returns the byte length for a given control byte.
func CtrlByte2ByteLengthsUint32(ctrl byte) int {
	return int(ctrl>>6&3+ctrl>>4&3+ctrl>>2&3+ctrl&3) + 4
}
