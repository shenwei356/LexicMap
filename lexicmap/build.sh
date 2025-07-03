#!/bin/sh

commit=$(git rev-parse --short HEAD)

# export GOARCH=arm64

export GOEXPERIMENT=greenteagc # for go1.25
CGO_ENABLED=0 go build -trimpath -o=lexicmap -ldflags="-s -w -X github.com/shenwei356/LexicMap/lexicmap/cmd.COMMIT=$commit" -tags netgo

./lexicmap version
