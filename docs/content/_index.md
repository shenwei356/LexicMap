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
[![license](https://img.shields.io/github/license/shenwei356/taxonkit.svg?maxAge=2592000)](https://github.com/shenwei356/taxonkit/blob/master/LICENSE)



<font size=5rem>LexicMap is a **nucleotide sequence alignment** tool for efficiently querying **gene, plasmid, virus, or long-read sequences (>150 bp)** against up to **millions** of **prokaryotic genomes**.</font>


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

[10 utility commands](https://bioinf.shenwei.me/LexicMap/usage/utils/) are available
to explore the index data, merge search results, extract matched subsequences and more.

{{< button size="small" relref="tutorials/index" >}}Tutorials{{< /button >}}
{{< button size="small" relref="usage/lexicmap" >}}Usages{{< /button >}}
{{< button size="small" relref="faqs" >}}FAQs{{< /button >}}

<--->

### Accurate and efficient alignment

Using LexicMap v0.7.0 to align against the whole **2,340,672** Genbank+Refseq prokaryotic genomes with 48 CPUs.

|Query            |Genome hits|Time   |RAM(GB)|
|:----------------|----------:|------:|------:|
|A 1.3-kb gene    |41,718     |3m:06s |3.97   |
|A 1.5-kb 16S rRNA|1,955,167  |32m:59s|11.09  |
|A 52.8-kb plasmid|560,330    |52m:22s|14.48  |
|1003 AMR genes   |30,967,882 |15h:52m|24.86  |


***Blastn** is unable to run with the same dataset on common servers as it requires >2000 GB RAM*.


{{< /columns >}}

