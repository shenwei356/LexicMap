---
title: Introduction
weight: 10
---


[![Latest Version](https://img.shields.io/github/release/shenwei356/LexicMap.svg?style=flat?maxAge=86400)](https://github.com/shenwei356/LexicMap/releases)
[![Anaconda Cloud](https://anaconda.org/bioconda/lexicmap/badges/version.svg)](https://anaconda.org/bioconda/lexicmap)
[![Cross-platform](https://img.shields.io/badge/platform-any-ec2eb4.svg?style=flat)](http://bioinf.shenwei.me/LexicMap/installation/)
[![license](https://img.shields.io/github/license/shenwei356/taxonkit.svg?maxAge=2592000)](https://github.com/shenwei356/taxonkit/blob/master/LICENSE)

LexicMap is a **nucleotide sequence alignment** tool for efficiently querying gene, plasmid, viral, or long-read sequences against up to **millions of prokaryotic genomes**.

Preprint:

> Wei Shen and Zamin Iqbal.
> (2024) LexicMap: efficient sequence alignment against millions of prokaryotic genomes.
> bioRxiv. [https://doi.org/10.1101/2024.08.30.610459](https://doi.org/10.1101/2024.08.30.610459)

## Table of contents

{{< toc format=html >}}

## Features

1. **LexicMap is scalable to up to millions of prokaryotic genomes**.
1. **The sensitivity of LexicMap is comparable with Blastn**.
1. **The alignment is [fast and memory-efficient](https://bioinf.shenwei.me/LexicMap/introduction/#searching)**.
1. LexicMap is **easy to [install](http://bioinf.shenwei.me/LexicMap/installation/)**,
   we provide [binary files](https://github.com/shenwei356/LexicMap/releases/) with no dependencies for Linux, Windows, MacOS (x86 and arm CPUs).
2. LexicMap is **easy to use** ([tutorials](http://bioinf.shenwei.me/LexicMap/tutorials/index/) and [usages](http://bioinf.shenwei.me/LexicMap/usage/lexicmap/)). Both tabular and Blast-style output formats are available.
3. Besides, we provide [several commands](https://bioinf.shenwei.me/LexicMap/usage/utils/) to explore the index data and extract indexed subsequences.

## Introduction

**Motivation**: Alignment against a database of genomes is a fundamental operation in bioinformatics, popularised by BLAST.
However, given the increasing rate at which genomes are sequenced, **existing tools struggle to scale**.

1. Existing full alignment tools face challenges of high memory consumption and slow speeds.
1. Alignment-free large-scale sequence searching tools only return the matched genomes,
   without the vital positional information for downstream analysis.
1. Prefilter+Align strategies have the sensitivity issue in the prefiltering step.

**Methods**: ([algorithm overview](#algorithm-overview))

1. An [improved version](https://github.com/shenwei356/lexichash) of the sequence sketching method [LexicHash](https://doi.org/10.1093/bioinformatics/btad652) is adopted to compute alignment seeds accurately and efficiently.
    1. **We solved the [sketching deserts](https://www.biorxiv.org/content/10.1101/2024.01.25.577301v1) problem of LexicHash seeds to provide a [window guarantee](https://doi.org/10.1093/bioinformatics/btab790)**.
    2. **We added the support of suffix matching of seeds, making seeds much more tolerant to mutations**. Any 31-bp seed with a common ≥15 bp prefix or suffix can be matched, which means **seeds are immune to any single SNP**.
2. **A hierarchical index enables fast and low-memory variable-length seed matching** (prefix + suffix matching).
3. A pseudo alignment algorithm is used to find similar sequence regions from chaining results for alignment.
4. A [reimplemented](https://github.com/shenwei356/wfa) [Wavefront alignment algorithm](https://doi.org/10.1093/bioinformatics/btaa777) is used for base-level alignment.

**Results**:

1. LexicMap enables efficient indexing and searching of both RefSeq+GenBank and the [AllTheBacteria](https://www.biorxiv.org/content/10.1101/2024.03.08.584059v1) datasets (**2.3 and 1.9 million prokaryotic assemblies** respectively).
Running at this scale has previously only been achieved by [Phylign](https://github.com/karel-brinda/Phylign) (previously called mof-search), which compresses genomes with phylogenetic information and provides searching
(prefiltering with [COBS](https://github.com/iqbal-lab-org/cobs) and alignment with [minimap2](https://github.com/lh3/minimap2)).
1. For searching in all **2,340,672 Genbank+Refseq prokaryotic genomes**, *Blastn is unable to run with this dataset on common servers as it requires >2000 GB RAM*.  (see [performance](#performance)).

    **With LexicMap v0.4.0** (48 CPUs),

    |Query               |Genome hits|Time      |RAM    |
    |:-------------------|----------:|---------:|------:|
    |A 1.3-kb marker gene|37,164     |36 s      |4.1 GB |
    |A 1.5-kb 16S rRNA   |1,949,496  |10 m 41 s |14.1 GB|
    |A 52.8-kb plasmid   |544,619    |19 m 20 s |19.3 GB|
    |1003 AMR genes      |25,702,419 |187 m 40 s|55.4 GB|


## Quick start

Building an index (see the tutorial of [building an index](http://bioinf.shenwei.me/LexicMap/tutorials/index/)).

    # From a directory with multiple genome files
    lexicmap index -I genomes/ -O db.lmi

    # From a file list with one file per line
    lexicmap index -X files.txt -O db.lmi

Querying (see the tutorial of [searching](http://bioinf.shenwei.me/LexicMap/tutorials/search/)).

    # For short queries like genes or long reads, returning top N hits.
    lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
        --min-qcov-per-hsp 70 --min-qcov-per-genome 70  --top-n-genomes 1000

    # For longer queries like plasmids, returning all hits.
    lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
        --min-qcov-per-hsp 0  --min-qcov-per-genome 0   --top-n-genomes 0


Sample output (queries are a few Nanopore Q20 reads). See [output format details](http://bioinf.shenwei.me/LexicMap/tutorials/search/#output-format).

```text
query                qlen   hits   sgenome           sseqid          qcovGnm   hsp   qcovHSP   alenHSP   pident    gaps   qstart   qend   sstart    send      sstr   slen
------------------   ----   ----   ---------------   -------------   -------   ---   -------   -------   -------   ----   ------   ----   -------   -------   ----   -------
ERR5396170.1000016   740    1      GCF_013394085.1   NZ_CP040910.1   89.595    1     89.595    663       99.246    0      71       733    13515     14177     +      1887974
ERR5396170.1000000   698    1      GCF_001457615.1   NZ_LN831024.1   85.673    1     85.673    603       98.010    5      53       650    4452083   4452685   +      6316979
ERR5396170.1000017   516    1      GCF_013394085.1   NZ_CP040910.1   94.574    1     94.574    489       99.591    2      27       514    293509    293996    +      1887974
ERR5396170.1000012   848    1      GCF_013394085.1   NZ_CP040910.1   95.165    1     95.165    811       97.411    7      22       828    190329    191136    -      1887974
ERR5396170.1000038   1615   1      GCA_000183865.1   CM001047.1      64.706    1     60.000    973       95.889    13     365      1333   88793     89756     -      2884551
ERR5396170.1000038   1615   1      GCA_000183865.1   CM001047.1      64.706    2     4.706     76        98.684    0      266      341    89817     89892     -      2884551
ERR5396170.1000036   1159   1      GCF_013394085.1   NZ_CP040910.1   95.427    1     95.427    1107      99.729    1      32       1137   1400097   1401203   +      1887974
ERR5396170.1000031   814    4      GCF_013394085.1   NZ_CP040910.1   86.486    1     86.486    707       99.151    3      104      807    242235    242941    -      1887974
ERR5396170.1000031   814    4      GCF_013394085.1   NZ_CP040910.1   86.486    2     86.486    707       98.444    3      104      807    1138777   1139483   +      1887974
ERR5396170.1000031   814    4      GCF_013394085.1   NZ_CP040910.1   86.486    3     84.152    688       98.983    4      104      788    154620    155306    -      1887974
ERR5396170.1000031   814    4      GCF_013394085.1   NZ_CP040910.1   86.486    4     84.029    687       99.127    3      104      787    32477     33163     +      1887974
ERR5396170.1000031   814    4      GCF_013394085.1   NZ_CP040910.1   86.486    5     72.727    595       98.992    3      104      695    1280183   1280777   +      1887974
ERR5396170.1000031   814    4      GCF_013394085.1   NZ_CP040910.1   86.486    6     11.671    95        100.000   0      693      787    1282480   1282574   +      1887974
ERR5396170.1000031   814    4      GCF_013394085.1   NZ_CP040910.1   86.486    7     82.064    671       99.106    3      120      787    1768782   1769452   +      1887974
```

CIGAR string, aligned query and subject sequences can be outputted as extra columns via the flag `-a/--all`.

    # Extracting similar sequences for a query gene.

    # search matches with query coverage >= 90%
    lexicmap search -d gtdb_complete.lmi/ b.gene_E_faecalis_SecY.fasta -o results.tsv \
        --min-qcov-per-hsp 90 --all

    # extract matched sequences as FASTA format
    sed 1d results.tsv | awk -F'\t' '{print ">"$5":"$14"-"$15":"$16"\n"$20;}' \
        | seqkit seq -g > results.fasta

    seqkit head -n 1 results.fasta | head -n 3
    >NZ_JALSCK010000007.1:39224-40522:-
    TTGTTCAAGCTATTAAAGAACGCCTTTAAAGTCAAAGACATTAGATCAAAAATCTTATTT
    ACAGTTTTAATCTTGTTTGTATTTCGCCTAGGTGCGCACATTACTGTGCCCGGGGTGAAT

Export blast-style format:

```
# here, we only align <=200 bp queries and show one low-similarity result.

$ seqkit seq -g -M 200 q.long-reads.fasta.gz \
    | lexicmap search -d demo.lmi/ -a \
    | csvtk filter2 -t -f '$pident >80 && $pident < 90' \
    | csvtk head -t -n 1 \
    | lexicmap utils 2blast --kv-file-genome ass2species.map

Query = GCF_003697165.2_r40
Length = 186

[Subject genome #1/2] = GCF_002950215.1 Shigella flexneri
Query coverage per genome = 93.548%

>NZ_CP026788.1 
Length = 4659463

 HSP #1
 Query coverage per seq = 93.548%, Aligned length = 177, Identities = 88.701%, Gaps = 6
 Query range = 13-186, Subject range = 1124816-1124989, Strand = Plus/Plus

Query  13       CGGAAACTGAAACA-CCAGATTCTACGATGATTATGATGATTTA-TGCTTTCTTTACTAA  70
                |||||||||||||| |||||||||| | |||||||||||||||| |||||||||| ||||
Sbjct  1124816  CGGAAACTGAAACAACCAGATTCTATGTTGATTATGATGATTTAATGCTTTCTTTGCTAA  1124875

Query  71       AAAGTAAGCGGCCAAAAAAATGAT-AACACCTGTAATGAGTATCAGAAAAGACACGGTAA  129
                ||    |||||||||||||||||| |||||||||||||||||||||||||||||||||||
Sbjct  1124876  AA--GCAGCGGCCAAAAAAATGATTAACACCTGTAATGAGTATCAGAAAAGACACGGTAA  1124933

Query  130      GAAAACACTCTTTTGGATACCTAGAGTCTGATAAGCGATTATTCTCTCTATGTTACT  186
                 || |||||||||    |||||  |||||||||||||||||||||||| |||| |||
Sbjct  1124934  AAAGACACTCTTTGAAGTACCTGAAGTCTGATAAGCGATTATTCTCTCCATGT-ACT  1124989

```


Learn more [tutorials](http://bioinf.shenwei.me/LexicMap/tutorials/index/) and [usages](http://bioinf.shenwei.me/LexicMap/usage/lexicmap/).

## Performance

### Indexing

|dataset          |genomes  |gzip_size|tool    |db_size|time     |RAM     |
|:----------------|--------:|--------:|:-------|------:|--------:|-------:|
|GTDB complete    |402,538  |443 GB   |LexicMap|973 GB |10 h 36 m|63.3 GB |
|                 |         |         |Blastn  |387 GB |3 h 11 m |718 MB  |
|AllTheBacteria HQ|1,858,610|2.5 TB   |LexicMap|4.26 TB|48 h 08 m|88.6 GB |
|                 |         |         |Blastn  |1.93 TB|14 h 03 m|2.9 GB  |
|                 |         |         |Phylign |248 GB |/        |/       |
|Genbank+RefSeq   |2,340,672|2.7 TB   |LexicMap|5.43 TB|54 h 33 m|178.3 GB|
|                 |         |         |Blastn  |2.37 TB|14 h 04 m|4.3 GB  |

Notes:
- All files are stored on a server with HDD disks. No files are cached in memory.
- Tests are performed in a single cluster node with 48 CPU cores (Intel Xeon Gold 6336Y CPU @ 2.40 GHz).
- LexicMap index building parameters: `-k 31 -m 40000`. Genome batch size: `-b 5000` for GTDB datasets, `-b 25000` for others.

### Searching

Blastn failed to run as it requires >2000GB RAM for Genbank+RefSeq and AllTheBacteria datasets.
Phylign only has the index for AllTheBacteria HQ dataset.


GTDB complete (402,538 genomes):

|query          |query_len    |tool           |genome_hits|genome_hits(qcov>50)|time      |RAM     |
|:--------------|------------:|:--------------|----------:|-------------------:|---------:|-------:|
|a marker gene  |1,299 bp     |LexicMap       |5,170      |5,143               |17 s      |1.4 GB  |
|               |             |Blastn         |7,121      |6,177               |2,171 s   |351.2 GB|
|a 16S rRNA gene|1,542 bp     |LexicMap       |303,925    |278,141             |235 s     |4.4 GB  |
|               |             |Blastn         |301,197    |277,042             |2,353 s   |378.4 GB|
|a plasmid      |52,830 bp    |LexicMap       |63,108     |1,190               |499 s     |4.6 GB  |
|               |             |Blastn         |69,311     |2,308               |2,262 s   |364.7 GB|
|1033 AMR genes |1 kb (median)|LexicMap       |3,867,003  |2,228,339           |4,350 s   |16.3 GB |
|               |             |Blastn         |5,357,772  |2,240,766           |4,686 s   |442.1 GB|



AllTheBacteria HQ (1,858,610 genomes):


|query          |query_len    |tool           |genome_hits|genome_hits(qcov>50)|time      |RAM     |
|:--------------|------------:|:--------------|----------:|-------------------:|---------:|-------:|
|a marker gene  |1,299 bp     |LexicMap       |27,963     |27,953              |31 s      |3.4 GB  |
|               |             |Phylign_local  |7,936      |                    |30 m 48 s |77.6 GB |
|               |             |Phylign_cluster|7,936      |                    |28 m 33 s |        |
|a 16S rRNA gene|1,542 bp     |LexicMap       |1,857,761  |1,740,000           |9 m 36 s  |14.9 GB |
|               |             |Phylign_local  |1,017,765  |                    |130 m 33 s|77.0 GB |
|               |             |Phylign_cluster|1,017,765  |                    |86 m 41 s |        |
|a plasmid      |52,830 bp    |LexicMap       |468,821    |3,618               |15 m 55 s |15.7 GB |
|               |             |Phylign_local  |46,822     |                    |47 m 33 s |82.6 GB |
|               |             |Phylign_cluster|46,822     |                    |39 m 34 s |        |
|1033 AMR genes |1 kb (median)|LexicMap       |21,288,000 |12,148,642          |138 m 55 s|49.9 GB |
|               |             |Phylign_local  |1,135,215  |                    |156 m 08 s|85.9 GB |
|               |             |Phylign_cluster|1,135,215  |                    |133 m 49 s|        |



Genbank+RefSeq (2,340,672 genomes):

|query          |query_len    |tool           |genome_hits|genome_hits(qcov>50)|time      |RAM     |
|:--------------|------------:|:--------------|----------:|-------------------:|---------:|-------:|
|a marker gene  |1,299 bp     |LexicMap       |37,164     |37,082              |36 s      |4.1 GB  |
|a 16S rRNA gene|1,542 bp     |LexicMap       |1,949,496  |1,381,974           |10 m 41 s |14.1 GB |
|a plasmid      |52,830 bp    |LexicMap       |544,619    |6,563               |19 m 20 s |19.3 GB |
|1033 AMR genes |1 kb (median)|LexicMap       |25,702,419 |14,692,624          |187 m 40 s|55.4 GB |


Notes:
- All files are stored on a server with HDD disks. No files are cached in memory.
- Tests are performed in a single cluster node with 48 CPU cores (Intel Xeon Gold 6336Y CPU @ 2.40 GHz).
- Main searching parameters:
    - LexicMap v0.4.0: `--threads 48 --top-n-genomes 0 --min-qcov-per-genome 0 --min-qcov-per-hsp 0 --min-match-pident 70`.
    - Blastn v2.15.0+: `-num_threads 48 -max_target_seqs 10000000`.
    - Phylign (AllTheBacteria fork 9fc65e6): `threads: 48, cobs_kmer_thres: 0.33, minimap_preset: "asm20", nb_best_hits: 5000000, max_ram_gb: 100`; For cluster, maximum number of slurm jobs is `100`.

## Installation

LexicMap is implemented in [Go](https://go.dev/) programming language,
executable binary files **for most popular operating systems** are freely available
in [release](https://github.com/shenwei356/lexicmap/releases) page.

Or install with `conda`:

    conda install -c bioconda lexicmap

## Algorithm overview

<img src="/LexicMap/overview.svg" alt="LexicMap overview" width="900"/>

## Citation

Wei Shen and Zamin Iqbal.
(2024) LexicMap: efficient sequence alignment against millions of prokaryotic genomes.
bioRxiv. [https://doi.org/10.1101/2024.08.30.610459](https://doi.org/10.1101/2024.08.30.610459)

## Support

Please [open an issue](https://github.com/shenwei356/LexicMap/issues) to report bugs,
propose new functions or ask for help.

## License

[MIT License](https://github.com/shenwei356/LexicMap/blob/master/LICENSE)


## Related projects

- High-performance [LexicHash](https://github.com/shenwei356/LexicHash) computation in Go.
- [Wavefront alignment algorithm (WFA) in Golang](https://github.com/shenwei356/wfa).
