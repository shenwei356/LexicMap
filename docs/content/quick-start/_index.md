---
title: Quick start
weight: 10
---

**Installing LexicMap** (see [installation](http://bioinf.shenwei.me/LexicMap/installation/)).

    conda install -c bioconda lexicmap

**Building an index** (see the tutorial of [building an index](http://bioinf.shenwei.me/LexicMap/tutorials/index/)).

    # From a directory with multiple genome files
    lexicmap index -I genomes/ -O db.lmi

    # From a file list with one file per line
    lexicmap index -X files.txt -O db.lmi

**Querying** (see the tutorial of [searching](http://bioinf.shenwei.me/LexicMap/tutorials/search/)).

    # For short queries like genes or long reads, returning top N hits.
    lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
        --min-match-pident 50 --min-qcov-per-genome 70 --min-qcov-per-hsp 70 --top-n-genomes 1000

    # For longer queries like plasmids, returning all hits.
    lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
        --min-match-pident 50 --min-qcov-per-genome 0  --min-qcov-per-hsp 0  --top-n-genomes 0


**Sample output** (queries are a few Nanopore Q20 reads). See [output format details](http://bioinf.shenwei.me/LexicMap/tutorials/search/#output).

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

Learn more [tutorials](http://bioinf.shenwei.me/LexicMap/tutorials/index/) and [usages](http://bioinf.shenwei.me/LexicMap/usage/lexicmap/).
