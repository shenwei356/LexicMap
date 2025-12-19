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
	"fmt"
	"net/http"
	"strings"

	"github.com/shenwei356/util/cliutil"
	"github.com/spf13/cobra"
)

// VERSION is the version
var VERSION = "0.8.1"

// COMMIT is the last commit
// var COMMIT = func() string {
// 	if info, ok := debug.ReadBuildInfo(); ok {
// 		for _, setting := range info.Settings {
// 			if setting.Key == "vcs.revision" {
// 				return setting.Value[:7]
// 			}
// 		}
// 	}
// 	return ""
// }()

// can pass from from command line:
// commit=$(git rev-parse --short HEAD)
// go build -trimpath -o=lexicmap -ldflags="-s -w -X github.com/shenwei356/LexicMap/lexicmap/cmd.COMMIT=$commit" -tags netgo
var COMMIT = ""

// versionCmd represents the version command
var versionCmd = &cobra.Command{
	Use:   "version",
	Short: "Print version information and check for update",
	Long: `Print version information and check for update

`,
	Run: func(cmd *cobra.Command, args []string) {
		app := "LexicMap"
		if COMMIT == "" {
			fmt.Printf("%s v%s\n", app, VERSION)
		} else {
			fmt.Printf("%s v%s (%s)\n", app, VERSION, COMMIT)
		}

		if !cliutil.GetFlagBool(cmd, "check-update") {
			return
		}

		fmt.Println("\nChecking new version...")

		resp, err := http.Get(fmt.Sprintf("https://github.com/shenwei356/%s/releases/latest", app))
		if err != nil {
			checkError(fmt.Errorf("network error"))
		}
		items := strings.Split(resp.Request.URL.String(), "/")
		version := ""
		if items[len(items)-1] == "" {
			version = items[len(items)-2]
		} else {
			version = items[len(items)-1]
		}
		if version == "v"+VERSION {
			fmt.Printf("You are using the latest version of %s\n", app)
		} else {
			fmt.Printf("New version available: %s %s at %s\n", app, version, resp.Request.URL.String())
		}
	},
}

func init() {
	RootCmd.AddCommand(versionCmd)

	versionCmd.Flags().BoolP("check-update", "u", false, `check update`)
}
