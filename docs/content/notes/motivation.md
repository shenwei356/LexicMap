---
title: Motivation
weight: 0
---

1. BLASTN is not able to scale to millions of bacterial genomes, it's slow and has a high memory occupation.
   For example, it requires >2000 GB for alignment a 2-kp gene sequence against all the 2.34 millions of prokaryotics genomes in Genbank and RefSeq.

2. [Large-scale sequence searching tools](https://kamimrcht.github.io/webpage/set_kmer_sets2.html) only return which genomes a query matches (color), but they can't return location information.
