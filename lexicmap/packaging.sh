#!/usr/bin/env sh

commit=""

if [ $# -gt 0 ]; then
    commit=" -X github.com/shenwei356/LexicMap/lexicmap/cmd.COMMIT=$(git rev-parse --short HEAD)"
fi

CGO_ENABLED=0 gox -os="windows darwin linux freebsd" -arch="amd64 arm64" -tags netgo -ldflags "-w -s $commit" -asmflags '-trimpath' \
    -output "lexicmap_{{.OS}}_{{.Arch}}"

dir=binaries
mkdir -p $dir;
rm -rf $dir/$f;

for f in lexicmap_*; do
    mkdir -p $dir/$f;
    mv $f $dir/$f;
    cd $dir/$f;
    mv $f $(echo $f | perl -pe 's/_[^\.]+//g');
    tar -zcf $f.tar.gz lexicmap*;
    mv *.tar.gz ../;
    cd ..;
    rm -rf $f;
    cd ..;
done;

ls binaries/*.tar.gz | rush 'cd {/}; md5sum {%} > {%}.md5.txt'
