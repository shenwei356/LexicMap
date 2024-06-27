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



<font size=5rem>LexicMap is a **nucleotide sequence alignment** tool for efficiently querying **gene, plasmid, virus, or long-read sequences** against up to **millions** of **prokaryotic genomes**.</font>


{{< button size="medium" relref="introduction" >}}Introduction{{< /button >}}



## Feature overview

{{< columns >}}

### Easy to install

Linux, Windows, MacOS and more OS are supported.

Both x86 and ARM CPUs are supported.

Just [download](https://github.com/shenwei356/lexicmap/releases) the binary files and run!


Or install it by

    conda install -c bioconda lexicmap


{{< button size="small" relref="installation" >}}Installation{{< /button >}}
{{< button size="small" relref="releases" >}}Releases{{< /button >}}

<--->

### Easy to use

Step 1: indexing

    lexicmap index -I genomes/ -O db.lmi

Step 2: searching

    lexicmap search -d db.lmi q.fasta -o r.tsv

{{< button size="small" relref="tutorials/index" >}}Tutorials{{< /button >}}
{{< button size="small" relref="usage/lexicmap" >}}Usages{{< /button >}}
{{< button size="small" relref="faqs" >}}FAQs{{< /button >}}
{{< button size="small" relref="notes/motivation" >}}Notes{{< /button >}}

<--->

### Efficient search

Using LexicMap to searching in the whole **2,340,672** Genbank+Refseq prokaryotic genomes with 48 CPUs.

|Query                   |Genome hits|Time    |RAM    |
|:-----------------------|----------:|-------:|------:|
|One 1.3-kb  gene        |16,832     |6.1 s   |1.5 GB |
|One 52.8-kb plasmid     |508,230    |6 m 50 s|17.6 GB|
|One 16S rRNA gene       |1,923,014  |7 m 32 s|12.3 GB|
|1003 AMR genes          |18,181,903 |1 h 28 m|20.0 GB|


***Blastn** is unable to run with the same dataset on common servers as it requires >2000 GB RAM*.

{{< button size="small" relref="introduction/#performance" >}}Performance{{< /button >}}


{{< /columns >}}

