---
title:
geekdocNav: false
geekdocAlign: center
geekdocAnchor: false
---
# <img src="logo.svg" width="50"/> LexicMap
<!-- markdownlint-capture -->
<!-- markdownlint-disable MD033 -->
<!-- markdownlint-restore -->

[![Latest Version](https://img.shields.io/github/release/shenwei356/LexicMap.svg?style=flat?maxAge=86400)](https://github.com/shenwei356/LexicMap/releases)
[![Github Releases](https://img.shields.io/github/downloads/shenwei356/LexicMap/latest/total.svg?maxAge=3600)](http://bioinf.shenwei.me/LexicMap/download/)
[![Anaconda Cloud](https://anaconda.org/bioconda/lexicmap/badges/version.svg)](https://anaconda.org/bioconda/lexicmap)
[![Cross-platform](https://img.shields.io/badge/platform-any-ec2eb4.svg?style=flat)](http://bioinf.shenwei.me/LexicMap/download/)


## LexicMap is a sequence alignment tool aiming to query gene or plasmid sequences efficiently against up to millions of prokaryotic genomes.


{{< button size="large" relref="quick-start" >}}Quick start{{< /button >}}



## Feature overview

{{< columns >}}

### Easy to install

    conda install -c bioconda lexicmap

Or [download](https://github.com/shenwei356/lexicmap/releases) the binary files for popular patforms.


{{< button size="small" relref="installation" >}}Installation{{< /button >}}
{{< button size="small" relref="releases" >}}Releases{{< /button >}}

<--->

### Easy to use

    # indexing
    lexicmap index -I genomes/ -O db.lmi

    # searching
    lexicmap search -d db.lmi q.fasta -o r.tsv

{{< button size="small" relref="introduction" >}}More details{{< /button >}}
{{< button size="small" relref="tutorials" >}}Tutorials{{< /button >}}
{{< button size="small" relref="usage" >}}Usages{{< /button >}}
{{< button size="small" relref="faqs" >}}FAQs{{< /button >}}

<--->

### Efficient search

Querying a 51.5-kb plasmid in all 2,340,672 Genbank+Refseq prokaryotic genomes takes only 3 minutes and 32 seconds with 15.7 GB RAM and 48 CPUs, with 19,265 genome hits returned.

{{< button size="small" relref="introduction/#performance" >}}Performance{{< /button >}}


{{< /columns >}}

