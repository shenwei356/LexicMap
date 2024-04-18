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
[![Anaconda Cloud](https://anaconda.org/bioconda/lexicmap/badges/version.svg)](https://anaconda.org/bioconda/lexicmap)
[![Cross-platform](https://img.shields.io/badge/platform-any-ec2eb4.svg?style=flat)](http://bioinf.shenwei.me/LexicMap/installation/)


LexicMap is a **sequence alignment** tool aiming to efficiently query **gene/plasmid/virus/long-read sequences** against up to **millions** of prokaryotic genomes.


{{< button size="large" relref="quick-start" >}}Quick start{{< /button >}}



## Feature overview

{{< columns >}}

### Easy to install

Both x86 and ARM CPUs are supported. Just install it by

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

{{< button size="small" relref="introduction" >}}Introduction{{< /button >}}
{{< button size="small" relref="tutorials/index" >}}Tutorials{{< /button >}}
{{< button size="small" relref="usage/lexicmap" >}}Usages{{< /button >}}
{{< button size="small" relref="faqs" >}}FAQs{{< /button >}}
{{< button size="small" relref="notes/motivation" >}}Notes{{< /button >}}

<--->

### Efficient search

Querying a **51.5-kb plasmid** in **all <ins>2,340,672</ins> Genbank+Refseq prokaryotic genomes** takes only <ins>**3 minutes and 32 seconds with 15.7 GB RAM**</ins> and 48 CPUs, with <ins>**19,265 genome hits**</ins> returned.

{{< button size="small" relref="introduction/#performance" >}}Performance{{< /button >}}


{{< /columns >}}

