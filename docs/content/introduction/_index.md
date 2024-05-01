---
title: Introduction
weight: 10
---


LexicMap is a nucleotide sequence pseudo-alignment tool for efficiently querying gene, plasmid, viral, or long-read sequences against up to millions of prokaryotic genomes.

**Motivation**: Alignment against a database of genomes is a fundamental operation in bioinformatics, popularised by BLAST. However, given the increasing rate at which genomes are sequenced, existing tools struggle to scale. Current tools either attempt full alignment but face challenges of high memory consumption and slow speeds, or they fall back on k-mer indexing, without information of where matches occur in the genome.

**Results**: In LexicMap, a [modified version](https://github.com/shenwei356/lexichash) of the sequence sketching method [LexicHash](https://doi.org/10.1093/bioinformatics/btad652) is adopted to compute alignment seeds.
A multi-level index enables fast and low-memory variable-length seed matching and pseudo-alignment on a single server
at the scale of millions of genomes (See [algorithm overview](#algorithm-overview)),
successfully indexing and searching both RefSeq+GenBank, and the [AllTheBacteria](https://www.biorxiv.org/content/10.1101/2024.03.08.584059v1) datasets (2.3 and 1.9 million genomes respectively).
Running at this scale has previously only been achieved by [Phylign](https://github.com/karel-brinda/Phylign) (previously called mof-search).

For example, **querying a 52.8-kb plasmid in all 2,340,672 Genbank+Refseq prokaryotic genomes takes only 4 minutes and 8 seconds with 15 GB RAM and 48 CPUs, with 494,860 genome hits returned**.
In contrast, BLASTN is unable to run with the same dataset on common servers because it requires >2000 GB RAM. See [performance](#performance).

LexicMap is easy to [install](http://bioinf.shenwei.me/LexicMap/installation/) (a binary file with no dependencies) and use ([tutorials](http://bioinf.shenwei.me/LexicMap/tutorials/index/) and [usages](http://bioinf.shenwei.me/LexicMap/usage/lexicmap/)).

More documents: http://bioinf.shenwei.me/LexicMap.

## Quick start

Building an index (see the tutorial of [building an index](http://bioinf.shenwei.me/LexicMap/tutorials/index/)).

    # From a directory with multiple genome files
    lexicmap index -I genomes/ -O db.lmi

    # From a file list with one file per line
    lexicmap index -X files.txt -O db.lmi

Querying (see the tutorial of [searching](http://bioinf.shenwei.me/LexicMap/tutorials/search/)).

    # For short queries like genes or long reads, returning top N hits.
    lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
        --min-match-pident 50 --min-qcov-per-genome 70 --min-qcov-per-hsp 70 --top-n-genomes 1000

    # For longer queries like plasmids, returning all hits.
    lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
        --min-match-pident 60 --min-qcov-per-genome 0  --min-qcov-per-hsp 0  --top-n-genomes 0


Sample output (queries are a few Nanopore Q20 reads). See [output format details](http://bioinf.shenwei.me/LexicMap/tutorials/search/#output).

```plain
query                qlen   hits   sgenome           sseqid                 qcovGnm   hsp   qcovHSP   alenHSP   pident    qstart   qend   sstart    send      sstr   slen      seeds   species
------------------   ----   ----   ---------------   --------------------   -------   ---   -------   -------   -------   ------   ----   -------   -------   ----   -------   -----   ----------------------------------
ERR5396170.1000006   796    4      GCF_013394085.1   NZ_CP040910.1          96.231    1     96.231    766       97.389    31       796    1138938   1139706   +      1887974   24      Limosilactobacillus fermentum
ERR5396170.1000006   796    4      GCF_013394085.1   NZ_CP040910.1          96.231    2     95.226    758       96.834    39       796    1280352   1282817   +      1887974   21      Limosilactobacillus fermentum
ERR5396170.1000006   796    4      GCF_013394085.1   NZ_CP040910.1          96.231    3     95.226    758       97.493    39       796    32646     33406     +      1887974   24      Limosilactobacillus fermentum
ERR5396170.1000006   796    4      GCF_013394085.1   NZ_CP040910.1          96.231    4     95.226    758       97.493    39       796    134468    135228    -      1887974   24      Limosilactobacillus fermentum
ERR5396170.1000006   796    4      GCF_013394085.1   NZ_CP040910.1          96.231    5     95.226    758       97.361    39       796    1768935   1769695   +      1887974   24      Limosilactobacillus fermentum
ERR5396170.1000006   796    4      GCF_013394085.1   NZ_CP040910.1          96.231    6     95.226    758       97.230    39       796    242012    242772    -      1887974   24      Limosilactobacillus fermentum
ERR5396170.1000006   796    4      GCF_013394085.1   NZ_CP040910.1          96.231    7     95.226    758       96.834    39       796    154380    155137    -      1887974   23      Limosilactobacillus fermentum
ERR5396170.1000006   796    4      GCF_003344625.1   NZ_QPKJ02000188.1      81.910    1     81.910    652       99.540    66       717    1         653       -      826       21      Solirubrobacter sp. CPCC 204708
ERR5396170.1000006   796    4      GCF_009663775.1   NZ_RDBR01000008.1      91.332    1     91.332    727       86.657    39       796    21391     22151     -      52610     7       Lactobacillus sp. 0.1XD8-4
ERR5396170.1000006   796    4      GCF_001591685.1   NZ_BCVJ01000102.1      87.940    1     87.940    700       88.429    66       796    434       1165      +      1933      4       Ligilactobacillus murinus
ERR5396170.1000017   516    1      GCF_013394085.1   NZ_CP040910.1          94.574    1     94.574    488       100.000   27       514    293509    293996    +      1887974   11      Limosilactobacillus fermentum
ERR5396170.1000012   848    1      GCF_013394085.1   NZ_CP040910.1          95.165    1     95.165    807       93.804    22       828    190329    191136    -      1887974   21      Limosilactobacillus fermentum
ERR5396170.1000052   330    1      GCF_013394085.1   NZ_CP040910.1          90.000    1     90.000    297       100.000   27       323    1161955   1162251   -      1887974   3       Limosilactobacillus fermentum
ERR5396170.1000000   698    1      GCF_001457615.1   NZ_LN831024.1          85.673    1     85.673    598       97.157    53       650    4452083   4452685   +      6316979   10      Pseudomonas aeruginosa

Note: the column `species` is added by mapping genome ID (column `sgenome`) to taxonomic information.
```

Matched query and subject sequences can be outputted as extra two columns via the flag `-a/-all`.

Learn more [tutorials](http://bioinf.shenwei.me/LexicMap/tutorials/index/) and [usages](http://bioinf.shenwei.me/LexicMap/usage/lexicmap/).

## Performance

Indexing

|dataset          |genomes  |gzip_size|tool    |db_size|time     |RAM   |
|:----------------|--------:|--------:|:-------|------:|--------:|-----:|
|GTDB complete    |402,538  |578 GB   |LexicMap|510 GB |3 h 26 m |35 GB |
|                 |         |         |Blastn  |360 GB |3 h 11 m |718 MB|
|Genbank+RefSeq   |2,340,672|3.5 TB   |LexicMap|2.91 TB|16 h 40 m|79 GB |
|                 |         |         |Blastn  |2.15 TB|14 h 4 m |4.3 GB|
|AllTheBacteria HQ|1,858,610|3.1 TB   |LexicMap|2.32 TB|10 h 48 m|43 GB |
|                 |         |         |Blastn  |1.76 TB|14 h 3 m |2.9 GB|
|                 |         |         |Phylign |248 GB |/        |/     |

Searching.
Blastn failed to run as it requires >2000GB RAM for Genbank+RefSeq and AllTheBacteria datasets.
Phylign only has the index for AllTheBacteria HQ dataset.

|dataset          |genomes  |query          |query_len|tool           |genome_hits|time     |RAM    |
|:----------------|--------:|:--------------|--------:|:--------------|----------:|--------:|------:|
|GTDB complete    |402,538  |a marker gene  |1,299 bp |LexicMap       |3,588      |1.4 s    |598 MB |
|                 |         |               |         |Blastn         |7,121      |36 m 11 s|351 GB |
|                 |         |a 16S rRNA gene|1,542 bp |LexicMap       |294,285    |1 m 19 s |3.0 GB |
|                 |         |               |         |Blastn         |301,197    |39 m 13 s|378 GB |
|                 |         |a plasmid      |52,830 bp|LexicMap       |58,930     |40 s     |3.2 GB |
|                 |         |               |         |Blastn         |69,311     |37 m 42 s|365 GB |
|Genbank+RefSeq   |2,340,672|a marker gene  |1,299 bp |LexicMap       |16,556     |6.2 s    |1.3 GB |
|                 |         |a 16S rRNA gene|1,542 bp |LexicMap       |1,875,260  |8 m 29 s |10.8 GB|
|                 |         |a plasmid      |52,830 bp|LexicMap       |494,860    |4 m 08 s |15.0 GB|
|AllTheBacteria HQ|1,858,610|a marker gene  |1,299 bp |LexicMap       |10,837     |11.3 s   |1.1 GB |
|                 |         |               |         |Phylign_local  |7,937      |1 h 44 m |27.1 GB|
|                 |         |               |         |Phylign_cluster|7,937      |32 m 52 s|/      |
|                 |         |a 16S rRNA gene|1,542 bp |LexicMap       |1,853,846  |12 m 31 s|9.7 GB |
|                 |         |               |         |Phylign_local  |1,032,948  |8 h 06 m |28.1 GB|
|                 |         |               |         |Phylign_cluster|1,032,948  |84 m 30 s|/      |
|                 |         |a plasmid      |52,830 bp|LexicMap       |427,112    |4 m 13 s |12.9 GB|
|                 |         |               |         |Phylign_local  |1,007      |2 h 50 m |20.6 GB|
|                 |         |               |         |Phylign_cluster|1,007      |32 m 23 s|/      |


Notes:
- All files are stored on a server with HDD disks. No files are cached in memory.
- Tests are performed in a single cluster node with 48 CPU cores (Intel Xeon Gold 6336Y CPU @ 2.40 GHz).
- LexicMap index building parameters: `-k 31 -m 40000`. Genome batch size: `-b 10000` for GTDB datasets, `-b 50000` for others.
- Main searching parameters:
    - LexicMap v0.3.0: `--threads 48 --top-n-genomes 0 --min-qcov-per-genome 0 --min-qcov-per-hsp 0 --min-match-pident 50`.
    - Blastn v2.15.0+: `-num_threads 48 -max_target_seqs 10000000`.
    - Phylign (AllTheBacteria fork 9fc65e6): `cobs_kmer_thres: 0.33, nb_best_hits: 5000000, max_ram_gb: 100`; For cluster: maxinum number of slurm jobs is 100.

## Installation

LexicMap is implemented in [Go](https://go.dev/) programming language,
executable binary files **for most popular operating systems** are freely available
in [release](https://github.com/shenwei356/lexicmap/releases) page.

Or install with `conda`:

    conda install -c bioconda lexicmap

## Algorithm overview

<img src="/LexicMap/overview.svg" alt="LexicMap overview" width="900"/>

## Related projects

- High-performance [LexicHash](https://github.com/shenwei356/LexicHash) computation in Go.

## Support

Please [open an issue](https://github.com/shenwei356/LexicMap/issues) to report bugs,
propose new functions or ask for help.

## License

[MIT License](https://github.com/shenwei356/LexicMap/blob/master/LICENSE)

