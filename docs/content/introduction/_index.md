---
title: Introduction
weight: 10
---


[![Latest Version](https://img.shields.io/github/release/shenwei356/LexicMap.svg?style=flat?maxAge=86400)](https://github.com/shenwei356/LexicMap/releases)
[![Anaconda Cloud](https://anaconda.org/bioconda/lexicmap/badges/version.svg)](https://anaconda.org/bioconda/lexicmap)
[![Cross-platform](https://img.shields.io/badge/platform-any-ec2eb4.svg?style=flat)](http://bioinf.shenwei.me/LexicMap/installation/)
[![license](https://img.shields.io/github/license/shenwei356/taxonkit.svg?maxAge=2592000)](https://github.com/shenwei356/taxonkit/blob/master/LICENSE)

LexicMap is a **nucleotide sequence alignment** tool for efficiently querying **gene, plasmid, viral, or long-read sequences (>150 bp)** against up to **millions of prokaryotic genomes**.


For the latest features and improvements, please download the [pre-release binaries](https://github.com/shenwei356/LexicMap/issues/10).

Preprint:

> Wei Shen, John A. Lees, Zamin Iqbal.
> (2024) LexicMap: efficient sequence alignment against millions of prokaryotic genomes.
> bioRxiv. [https://doi.org/10.1101/2024.08.30.610459](https://doi.org/10.1101/2024.08.30.610459)


## Table of contents

- [Table of contents](#table-of-contents)
- [Features](#features)
- [Introduction](#introduction)
- [Quick start](#quick-start)
- [Performance](#performance)
- [Installation](#installation)
- [Algorithm overview](#algorithm-overview)
- [Citation](#citation)
- [Limitations](#limitations)
- [Terminology differences](#terminology-differences)
- [Support](#support)
- [License](#license)
- [Related projects](#related-projects)

## Features

1. **The accuracy of LexicMap is comparable with Blastn, MMseqs2, and Minimap2**. It
    - performs **base-level alignment**, with `qcovGnm`, `qcovHSP`, `pident`, `evalue` and `bitscore` returned,
      both in TSV and pairwise alignment format ([output format](https://bioinf.shenwei.me/LexicMap/tutorials/search/#output)).
        - provides a genome-wide query coverage metric (`qcovGnm`), 
          which enables accurate interpretation of search results - particularly for [circular queries (such as plasmid, virus, and mtDNA)](https://bioinf.shenwei.me/LexicMap/tutorials/search/#searching-with-plasmids-or-other-longer-queries)
           against both complete and fragmented assemblies.
    - returns all possible matches, including multiple copies of a gene in a genome.
1. **The alignment is fast and memory-efficient, scalable to up to millions of prokaryotic genomes**.
1. LexicMap is **easy to [install](http://bioinf.shenwei.me/LexicMap/installation/),
   we provide [binary files](https://github.com/shenwei356/LexicMap/releases/)** with no dependencies for Linux, Windows, MacOS (x86 and arm CPUs).
2. LexicMap is **easy to use** (see [tutorials](http://bioinf.shenwei.me/LexicMap/tutorials/index/), [usages](http://bioinf.shenwei.me/LexicMap/usage/lexicmap/), and [FAQs](https://bioinf.shenwei.me/LexicMap/faqs/)).
    - [Database building](https://bioinf.shenwei.me/LexicMap/tutorials/index/) requires only a simple command, accepting input from files, a file list, or even a directory.
    - [Sequence searching](https://bioinf.shenwei.me/LexicMap/tutorials/search/) supports limiting search by TaxId(s), provides a progress bar.
    - [Several utility commands](https://bioinf.shenwei.me/LexicMap/usage/utils/) are available to resume unfinished indexing, and explore the index data, extract indexed subsequences.

## Introduction

**Motivation**: Alignment against a database of genomes is a fundamental operation in bioinformatics, popularised by BLAST.
However, given the increasing rate at which genomes are sequenced, **existing tools struggle to scale**.

1. Existing full alignment tools face challenges of high memory consumption and slow speeds.
1. Alignment-free large-scale sequence searching tools only return the matched genomes,
   without the vital positional information for downstream analysis.
1. Mapping tools, or those utilizing compressed full-text indexes, return only the most similar matches.
1. Prefilter+Align strategies have the sensitivity issue in the prefiltering step.

**Methods**: ([algorithm overview](#algorithm-overview))

1. A [rewritten and improved version](https://github.com/shenwei356/lexichash) of the sequence sketching method [LexicHash](https://doi.org/10.1093/bioinformatics/btad652) is adopted to compute alignment seeds accurately and efficiently.
    - **We solved the [sketching deserts](https://www.biorxiv.org/content/10.1101/2024.01.25.577301v1) problem of LexicHash seeds to provide a [window guarantee](https://doi.org/10.1093/bioinformatics/btab790)**.
    - **We added the support of suffix matching of seeds, making seeds much more tolerant to mutations**. Any 31-bp seed with a common ≥15 bp prefix or suffix can be matched.
2. **A hierarchical index enables fast and low-memory variable-length seed matching** (prefix + suffix matching).
3. A pseudo alignment algorithm is used to find similar sequence regions from chaining results for alignment.
4. A [reimplemented](https://github.com/shenwei356/wfa) [Wavefront alignment algorithm](https://doi.org/10.1093/bioinformatics/btaa777) is used for base-level alignment.

**Results**:

1. LexicMap enables efficient indexing and searching of both RefSeq+GenBank and the [AllTheBacteria](https://www.biorxiv.org/content/10.1101/2024.03.08.584059v1) datasets (**2.3 and 1.9 million prokaryotic assemblies** respectively).
1. When searching in all **2,340,672 Genbank+Refseq prokaryotic genomes**, *Blastn is unable to run with this dataset on common servers as it requires >2000 GB RAM*.  (see [performance](#performance)).
    
    **With LexicMap v0.7.0** (48 CPUs),

    |Query               |Genome hits|Genome hits<br/>(high-similarity)|Genome hits<br/>(medium-similarity)|Genome hits<br/>(low-similarity)|Time       |RAM     |
    |:-------------------|----------:|--------------------------------:|----------------------------------:|-------------------------------:|----------:|-------:|
    |A 1.3-kb marker gene|41,718     |11,746                           |115                                |29,857                          |3m:06s     |3.97 GB |
    |A 1.5-kb 16S rRNA   |1,955,167  |245,884                          |501,691                            |1,207,592                       |32m:59s    |11.09 GB|
    |A 52.8-kb plasmid   |560,330    |96                               |15,370                             |544,864                         |52m:22s    |14.48 GB|
    |1003 AMR genes      |30,967,882 |7,636,386                        |4,858,063                          |18,473,433                      |15h:52m:08s|24.86 GB|
    
    Notes:
    1. Default paramters are used, for returning all possible matches.
    1. Only the best alignment of a genome is used to evaluate alignment similarity:
        - high-similarity: (a) qcov >= 90% (genes) or 70% (plasmids), (b) pident>=90%.
        - medium-similarity: (a) not belong to high-similarity, (b) qcov >= 50% (genes) or 30% (plasmids), (c) pident>=80%.
        - low-similarity: left.
    1. The search time varies in different computing environments and mainly depends on the I/O speed.

## Quick start

Building an index (see the tutorial of [building an index](http://bioinf.shenwei.me/LexicMap/tutorials/index/)).

```plain
# From a directory with multiple genome files
lexicmap index -I genomes/ -O db.lmi

# From a file list with one file per line
lexicmap index -S -X files.txt -O db.lmi
```

Querying (see the tutorial of [searching](http://bioinf.shenwei.me/LexicMap/tutorials/search/)).

```plain
# For short queries like genes or long reads, returning top N hits.
lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
    --min-qcov-per-hsp 70 --min-qcov-per-genome 70  --top-n-genomes 10000

# For longer queries like plasmids, returning all hits.
lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
    --min-qcov-per-hsp 0  --min-qcov-per-genome 0   --top-n-genomes 0
```

Sample output (queries are a few Nanopore Q20 reads). See [output format details](http://bioinf.shenwei.me/LexicMap/tutorials/search/#output-format).

```plain
query                qlen   hits   sgenome           sseqid              qcovGnm   cls   hsp   qcovHSP   alenHSP   pident   gaps   qstart   qend   sstart    send      sstr   slen      evalue      bitscore
------------------   ----   ----   ---------------   -----------------   -------   ---   ---   -------   -------   ------   ----   ------   ----   -------   -------   ----   -------   ---------   --------
ERR5396170.1000004   190    1      GCF_000227465.1   NC_016047.1         84.211    1     1     84.211    165       89.091   5      14       173    4189372   4189536   -      4207222   1.93e-63    253     
ERR5396170.1000006   796    3      GCF_013394085.1   NZ_CP040910.1       99.623    1     1     99.623    801       97.628   9      4        796    1138907   1139706   +      1887974   0.00e+00    1431    
ERR5396170.1000006   796    3      GCF_013394085.1   NZ_CP040910.1       99.623    2     2     99.623    801       97.628   9      4        796    32607     33406     +      1887974   0.00e+00    1431    
ERR5396170.1000006   796    3      GCF_013394085.1   NZ_CP040910.1       99.623    3     3     99.623    801       97.628   9      4        796    134468    135267    -      1887974   0.00e+00    1431    
ERR5396170.1000006   796    3      GCF_013394085.1   NZ_CP040910.1       99.623    4     4     99.623    801       97.503   9      4        796    1768896   1769695   +      1887974   0.00e+00    1427    
ERR5396170.1000006   796    3      GCF_013394085.1   NZ_CP040910.1       99.623    5     5     99.623    801       97.378   9      4        796    242012    242811    -      1887974   0.00e+00    1422    
ERR5396170.1000006   796    3      GCF_013394085.1   NZ_CP040910.1       99.623    6     6     99.623    801       96.879   12     4        796    154380    155176    -      1887974   0.00e+00    1431    
ERR5396170.1000006   796    3      GCF_013394085.1   NZ_CP040910.1       99.623    7     7     57.915    469       95.736   9      4        464    1280313   1280780   +      1887974   3.71e-236   829     
ERR5396170.1000006   796    3      GCF_013394085.1   NZ_CP040910.1       99.623    8     8     42.839    341       99.120   0      456      796    1282477   1282817   +      1887974   6.91e-168   601     
ERR5396170.1000006   796    3      GCF_009663775.1   NZ_RDBR01000008.1   99.623    1     1     99.623    801       93.383   9      4        796    21391     22190     -      52610     0.00e+00    1278    
ERR5396170.1000006   796    3      GCF_003344625.1   NZ_QPKJ02000188.1   97.362    1     1     87.437    700       98.143   5      22       717    1         699       -      826       0.00e+00    1249    
ERR5396170.1000006   796    3      GCF_003344625.1   NZ_QPKJ02000423.1   97.362    2     2     27.889    222       99.550   0      575      796    1         222       +      510       3.47e-106   396     
ERR5396170.1000000   698    2      GCF_001457615.1   NZ_LN831024.1       92.264    1     1     92.264    656       96.341   13     53       696    4452083   4452737   +      6316979   0.00e+00    1169    
ERR5396170.1000000   698    2      GCF_000949385.2   NZ_JYKO02000001.1   91.977    1     1     91.977    654       78.135   13     55       696    5638788   5639440   -      5912440   2.68e-176   630     
ERR5396170.1000001   2505   3      GCF_000307025.1   NC_018584.1         67.066    1     1     67.066    1690      97.633   16     47       1726   1905511   1907194   -      2951805   0.00e+00    2985    
ERR5396170.1000001   2505   3      GCF_900187225.1   NZ_LT906436.1       65.070    1     1     65.070    1641      93.723   20     95       1724   1869503   1871134   -      2864663   0.00e+00    2626    
ERR5396170.1000001   2505   3      GCF_013394085.1   NZ_CP040910.1       30.858    1     1     30.858    780       97.692   9      1726     2498   183873    184650    +      1887974   0.00e+00    1384    
ERR5396170.1000001   2505   3      GCF_013394085.1   NZ_CP040910.1       30.858    2     2     5.030     127       87.402   1      2233     2358   1236170   1236296   +      1887974   1.73e-37    167     
ERR5396170.1000001   2505   3      GCF_013394085.1   NZ_CP040910.1       30.858    3     3     5.150     130       80.769   12     2233     2361   930381    930499    -      1887974   6.61e-43    185     
ERR5396170.1000001   2505   3      GCF_013394085.1   NZ_CP040910.1       30.858    4     4     3.713     93        93.548   0      2257     2349   1104581   1104673   -      1887974   5.09e-30    141
```

CIGAR string, aligned query and subject sequences can be outputted as extra columns via the flag `-a/--all`.

```plain
# Extracting similar sequences for a query gene.

# search matches with query coverage >= 90%
lexicmap search -d gtdb_complete.lmi/ b.gene_E_faecalis_SecY.fasta -o results.tsv \
    --min-qcov-per-hsp 90 --all

# extract matched sequences as FASTA format
sed 1d results.tsv | awk -F'\t' '{print ">"$5":"$15"-"$16":"$17"\n"$23;}' \
    | seqkit seq -g > results.fasta

seqkit head -n 1 results.fasta | head -n 3
>NZ_JALSCK010000007.1:39224-40522:-
TTGTTCAAGCTATTAAAGAACGCCTTTAAAGTCAAAGACATTAGATCAAAAATCTTATTT
ACAGTTTTAATCTTGTTTGTATTTCGCCTAGGTGCGCACATTACTGTGCCCGGGGTGAAT
```

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

 HSP cluster #1, HSP #1
 Score = 279 bits, Expect = 9.66e-75
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


Learn more: [demo](https://github.com/shenwei356/LexicMap/tree/main/demo),
[tutorials](http://bioinf.shenwei.me/LexicMap/tutorials/index/),
or [usages](http://bioinf.shenwei.me/LexicMap/usage/lexicmap/).

## Performance

See the [paper](#citation).

## Installation

LexicMap is implemented in [Go](https://go.dev/) programming language,
executable binary files **for most popular operating systems** are freely available
in [release page](https://github.com/shenwei356/lexicmap/releases).

Or install with conda or pixi:

    conda install -c bioconda lexicmap

**We also provide [pre-release binaries](https://github.com/shenwei356/LexicMap/issues/10), with new features and improvements**.

## Algorithm overview

<img src="/LexicMap/overview.svg" alt="LexicMap overview" width="800"/>

## Citation

Wei Shen, John A. Lees, Zamin Iqbal.
(2024) LexicMap: efficient sequence alignment against millions of prokaryotic genomes.
bioRxiv. [https://doi.org/10.1101/2024.08.30.610459](https://doi.org/10.1101/2024.08.30.610459)

## Limitations

- The queries need to be longer than 100 bp, though some shorter one can also be aligned.
- LexicMap is slow for >1Mb queries, and the alignment might be fragmented.
- LexicMap is slow for batch searching with more than hundreds of queries. However, there are [some ways to improve the search speed of lexicmap search](http://bioinf.shenwei.me/LexicMap/tutorials/search/#improving-searching-speed), such as keeping the top N genome matches via `-n/--top-n-genomes` or storing the index on solid state drives (SSDs).

## Terminology differences

- In the LexicMap source code and command line options, the term **"mask"** is used, following the terminology in the LexicHash paper.
- In the LexicMap manuscript, however, we use **"probe"** as it is easier to understand.
  Because these masks, which consist of thousands of k-mers and capture k-mers from sequences through prefix matching, function similarly to DNA probes in molecular biology.

## Support

Please [open an issue](https://github.com/shenwei356/LexicMap/issues) to report bugs,
propose new functions or ask for help.

## License

[MIT License](https://github.com/shenwei356/LexicMap/blob/master/LICENSE)


## Related projects

- High-performance [LexicHash](https://github.com/shenwei356/LexicHash) computation in Go.
- [Wavefront alignment algorithm (WFA) in Golang](https://github.com/shenwei356/wfa).
