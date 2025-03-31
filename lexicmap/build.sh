#!/bin/sh

commit=$(git rev-parse --short HEAD)

go build -trimpath -o=lexicmap -ldflags="-s -w -X github.com/shenwei356/LexicMap/lexicmap/cmd.COMMIT=$commit" -tags netgo

./lexicmap version
