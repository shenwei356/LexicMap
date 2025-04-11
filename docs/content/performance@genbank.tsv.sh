#!/bin/sh

cat performance@genbank.tsv \
    | csvtk replace -t -f RAM -p ' .+' \
    | csvtk rename -t -f RAM -n 'RAM(GB)' \
    | csvtk replace -t -f Query -p 'marker ' \
    | csvtk replace -t -f Time -p '(\d+h:\d+m):\d+s' -r '$1' \
    | csvtk cut -t -f 1,2,6,7 \
    | csvtk comma -t -f 2 \
    | csvtk csv2md -t -a l,r,r,r

echo

cat performance@genbank.tsv \
    | csvtk comma -t -f 2-5 \
    | csvtk csv2md -t -a l,r,r,r,r,r,r
