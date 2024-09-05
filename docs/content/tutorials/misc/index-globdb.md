---
title: Indexing GlobDB
weight: 20
---


Info:

- [GlobDB](https://globdb.org/) , a dereplicated dataset of the species reps of the GTDB, GEM, SPIRE and SMAG datasets a lot.
- https://x.com/daanspeth/status/1822964436950192218


Steps:

    # download data
    wget https://fileshare.lisc.univie.ac.at/globdb/globdb_r220/globdb_r220_genome_fasta.tar.gz

    tar -zxf globdb_r220_genome_fasta.tar.gz

    # file list
    find globdb_r220_genome_fasta/ -name "*.fa.gz" > files.txt

    # index with lexicmap
    # elapsed time: 3h:40m:38s
    # peak rss: 87.15 GB
    lexicmap index -S -X files.txt -O globdb_r220.lmi --log globdb_r220.lmi -g 50000000

