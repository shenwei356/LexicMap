---
title: Motivation
weight: 0
---

1. BLASTN cannot scale to millions of bacterial genomes: it is slow and has high memory usage.
   For example, it requires >2000 GB to align a 2-kb gene sequence against all 2.34 million prokaryotic genomes in Genbank and RefSeq.

2. [Large-scale sequence searching tools](https://kamimrcht.github.io/webpage/set_kmer_sets2.html) only return which genomes a query matches (color), but they can't return positional information.
