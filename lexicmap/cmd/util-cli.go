// Copyright © 2023-2025 Wei Shen <shenwei356@gmail.com>
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
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/pkg/errors"
	"github.com/shenwei356/util/stringutil"
	"github.com/spf13/cobra"
)

func checkError(err error) {
	if err != nil {
		log.Error(err)
		os.Exit(-1)
	}
}

func isStdin(file string) bool {
	return file == "-"
}

func isStdout(file string) bool {
	return file == "-"
}

func getFlagInt(cmd *cobra.Command, flag string) int {
	value, err := cmd.Flags().GetInt(flag)
	checkError(err)
	return value
}

func getFlagIntSlice(cmd *cobra.Command, flag string) []int {
	value, err := cmd.Flags().GetIntSlice(flag)
	checkError(err)
	return value
}

func getFlagUint8(cmd *cobra.Command, flag string) uint8 {
	value, err := cmd.Flags().GetUint8(flag)
	checkError(err)
	return value
}

func getFlagUint32(cmd *cobra.Command, flag string) uint32 {
	value, err := cmd.Flags().GetUint32(flag)
	checkError(err)
	return value
}

func getFlagUint64(cmd *cobra.Command, flag string) uint64 {
	value, err := cmd.Flags().GetUint64(flag)
	checkError(err)
	return value
}

func getFlagPositiveInt(cmd *cobra.Command, flag string) int {
	value, err := cmd.Flags().GetInt(flag)
	checkError(err)
	if value <= 0 {
		checkError(fmt.Errorf("value of flag --%s should be greater than 0", flag))
	}
	return value
}

func getFlagPositiveFloat64(cmd *cobra.Command, flag string) float64 {
	value, err := cmd.Flags().GetFloat64(flag)
	checkError(err)
	if value <= 0 {
		checkError(fmt.Errorf("value of flag --%s should be greater than 0", flag))
	}
	return value
}

func getFlagNonNegativeInt(cmd *cobra.Command, flag string) int {
	value, err := cmd.Flags().GetInt(flag)
	checkError(err)
	if value < 0 {
		checkError(fmt.Errorf("value of flag --%s should be greater than or equal to 0", flag))
	}
	return value
}

func getFlagNonNegativeFloat64(cmd *cobra.Command, flag string) float64 {
	value, err := cmd.Flags().GetFloat64(flag)
	checkError(err)
	if value < 0 {
		checkError(fmt.Errorf("value of flag --%s should be greater than or equal to ", flag))
	}
	return value
}

func getFlagBool(cmd *cobra.Command, flag string) bool {
	value, err := cmd.Flags().GetBool(flag)
	checkError(err)
	return value
}

func getFlagString(cmd *cobra.Command, flag string) string {
	value, err := cmd.Flags().GetString(flag)
	checkError(err)
	return value
}

func getFlagNonEmptyString(cmd *cobra.Command, flag string) string {
	value, err := cmd.Flags().GetString(flag)
	checkError(err)
	if value == "" {
		checkError(fmt.Errorf("flag --%s needed", flag))
	}
	return value
}

func getFlagCommaSeparatedStrings(cmd *cobra.Command, flag string) []string {
	value, err := cmd.Flags().GetString(flag)
	checkError(err)
	return stringutil.Split(value, ",")
}

func getFlagSemicolonSeparatedStrings(cmd *cobra.Command, flag string) []string {
	value, err := cmd.Flags().GetString(flag)
	checkError(err)
	return stringutil.Split(value, ";")
}

func getFlagCommaSeparatedInts(cmd *cobra.Command, flag string) []int {
	filedsStrList := getFlagCommaSeparatedStrings(cmd, flag)
	fields := make([]int, len(filedsStrList))
	for i, value := range filedsStrList {
		v, err := strconv.Atoi(value)
		if err != nil {
			checkError(fmt.Errorf("value of flag --%s should be comma separated integers", flag))
		}
		fields[i] = v
	}
	return fields
}

func getFlagRune(cmd *cobra.Command, flag string) rune {
	value, err := cmd.Flags().GetString(flag)
	checkError(err)
	if len(value) > 1 {
		checkError(fmt.Errorf("value of flag --%s should has length of 1", flag))
	}
	var v rune
	for _, r := range value {
		v = r
		break
	}
	return v
}

func getFlagFloat64(cmd *cobra.Command, flag string) float64 {
	value, err := cmd.Flags().GetFloat64(flag)
	checkError(err)
	return value
}

func getFlagInt64(cmd *cobra.Command, flag string) int64 {
	value, err := cmd.Flags().GetInt64(flag)
	checkError(err)
	return value
}

func getFlagStringSlice(cmd *cobra.Command, flag string) []string {
	value, err := cmd.Flags().GetStringSlice(flag)
	checkError(err)
	return value
}

func getFileList(args []string, checkFile bool) []string {
	files := make([]string, 0, 1024)
	if len(args) == 0 {
		files = append(files, "-")
	} else {
		for _, file := range args {
			if isStdin(file) {
				continue
			}
			if !checkFile {
				continue
			}
			if _, err := os.Stat(file); os.IsNotExist(err) {
				checkError(errors.Wrap(err, file))
			}
		}
		files = args
	}
	return files
}

func getFileListFromFile(file string, checkFile bool) ([]string, error) {
	var fh *os.File
	var err error
	if file == "-" {
		fh = os.Stdin
	} else {
		fh, err = os.Open(file)
		if err != nil {
			return nil, fmt.Errorf("read file list from '%s': %s", file, err)
		}
	}

	var _file string
	lists := make([]string, 0, 1024)
	scanner := bufio.NewScanner(fh)
	for scanner.Scan() {
		_file = scanner.Text()
		if strings.TrimSpace(_file) == "" {
			continue
		}
		lists = append(lists, _file)
	}
	if err = scanner.Err(); err != nil {
		return nil, fmt.Errorf("read file list from '%s': %s", file, err)
	}

	if !checkFile {
		return lists, fh.Close()
	}

	for _, _file = range lists {
		if !isStdin(_file) {
			if _, err = os.Stat(_file); os.IsNotExist(err) {
				return lists, fmt.Errorf("check file '%s': %s", _file, err)
			}
		}
	}

	return lists, fh.Close()
}

func getFileListFromArgsAndFile(cmd *cobra.Command, args []string, checkFileFromArgs bool, flag string, checkFileFromFile bool) []string {
	infileList := getFlagString(cmd, flag)
	files := getFileList(args, checkFileFromArgs)
	if infileList != "" {
		_files, err := getFileListFromFile(infileList, checkFileFromFile)
		checkError(err)
		if len(_files) == 0 {
			log.Warningf("no files found in file list: %s", infileList)
			return files
		}

		if len(files) == 1 && isStdin(files[0]) {
			return _files
		}
		files = append(files, _files...)
	}
	return files
}

// ParseByteSize parses byte size from string
func ParseByteSize(val string) (int64, error) {
	val = strings.Trim(val, " \t\r\n")
	if val == "" {
		return 0, nil
	}
	var u int64
	var noUnit bool
	switch val[len(val)-1] {
	case 'B', 'b':
		u = 1
	case 'K', 'k':
		u = 1 << 10
	case 'M', 'm':
		u = 1 << 20
	case 'G', 'g':
		u = 1 << 30
	case 'T', 't':
		u = 1 << 40
	default:
		noUnit = true
		u = 1
	}
	var size float64
	var err error
	if noUnit {
		size, err = strconv.ParseFloat(val, 10)
		if err != nil {
			return 0, fmt.Errorf("invalid byte size: %s", val)
		}
		if size < 0 {
			size = 0
		}
		return int64(size), nil
	}

	if len(val) == 1 { // no value
		return 0, nil
	}

	size, err = strconv.ParseFloat(strings.Trim(val[0:len(val)-1], " \t\r\n"), 10)
	if err != nil {
		return 0, fmt.Errorf("invalid byte size: %s", val)
	}
	if size < 0 {
		size = 0
	}
	return int64(size * float64(u)), nil
}

func roundup32(x uint32) uint32 {
	if x == 0 {
		return 1
	}
	x--
	x |= x >> 1
	x |= x >> 2
	x |= x >> 4
	x |= x >> 8
	x |= x >> 16
	return (x | x>>32) + 1
}
