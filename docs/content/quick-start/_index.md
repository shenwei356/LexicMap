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
        --min-qcov-per-genome 70 --min-match-pident 70 --min-qcov-per-hsp 70 --top-n-genomes 500

    # For longer queries like plasmids, returning all hits.
    lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
        --min-qcov-per-genome 50 --min-match-pident 70 --min-qcov-per-hsp 0  --top-n-genomes 0


**Sample output** (queries are a few Nanopore Q20 reads). See [output format details](http://bioinf.shenwei.me/LexicMap/tutorials/search/#output).

```plain
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 query                qlen   qstart   qend   sgnms   sgnm              seqid               qcovGnm   hsp   qcovHSP   alenHSP   alenSeg    pident    slen      sstart    send      sstr   seeds   species
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ERR5396170.1000016   740    71       733    1       GCF_013394085.1   NZ_CP040910.1       89.595    1     89.595    663       663        98.492    1887974   13515     14177     +      19      Limosilactobacillus fermentum
 ERR5396170.1000017   516    27       514    1       GCF_013394085.1   NZ_CP040910.1       94.574    1     94.574    488       488        100.000   1887974   293509    293996    +      6       Limosilactobacillus fermentum
 ERR5396170.1000029   480    37       474    1       GCF_001027105.1   NZ_CP011526.1       91.042    1     91.042    437       437        95.423    2755072   821078    821514    +      1       Staphylococcus aureus
 ERR5396170.1000047   960    24       812    2       GCF_001027105.1   NZ_CP011526.1       91.979    1     91.979    883       803        89.015    2755072   2204718   2205520   -      7       Staphylococcus aureus
 ERR5396170.1000047   960    881      960    2       GCF_001027105.1   NZ_CP011526.1       91.979    1     91.979    883       80         89.015    2755072   2204568   2204647   -      7       Staphylococcus aureus
 ERR5396170.1000047   960    42       960    2       GCF_002902405.1   NZ_PPQS01000020.1   100.000   1     97.500    936       936        77.457    50421     25900     26835     +      3       Staphylococcus schweitzeri
 ERR5396170.1000047   960    42       950    2       GCF_002902405.1   NZ_PPQS01000020.1   100.000   2     96.458    926       926        77.214    50421     25900     26825     +      1       Staphylococcus schweitzeri
 ERR5396170.1000000   698    53       650    1       GCF_001457615.1   NZ_LN831024.1       86.390    1     86.390    603       603        96.517    6316979   4452083   4452685   +      4       Pseudomonas aeruginosa
 ERR5396170.1000005   2516   38       2510   5       GCF_000006945.2   NC_003197.2         98.490    1     98.490    2478      2478       98.951    4857450   3198806   3201283   +      14      Salmonella enterica
 ERR5396170.1000005   2516   38       2497   5       GCF_008692785.1   NZ_VXJV01000001.1   98.013    1     98.013    2466      2466       95.864    797633    423400    425865    +      8       Salmonella enterica
 ERR5396170.1000005   2516   40       2510   5       GCA_900478215.1   LS483478.1          98.450    1     98.450    2477      2477       95.317    4624613   785866    788342    -      12      Salmonella enterica
 ERR5396170.1000005   2516   1350     2497   5       GCF_008692845.1   NZ_VXJW01000004.1   87.599    1     87.599    2204      1151       91.742    366711    6705      7855      +      9       Salmonella enterica
 ERR5396170.1000005   2516   634      1309   5       GCF_008692845.1   NZ_VXJW01000004.1   87.599    1     87.599    2204      674        91.742    366711    5991      6664      +      9       Salmonella enterica
 ERR5396170.1000005   2516   387      608    5       GCF_008692845.1   NZ_VXJW01000004.1   87.599    1     87.599    2204      221        91.742    366711    5745      5965      +      9       Salmonella enterica
 ERR5396170.1000005   2516   69       205    5       GCF_008692845.1   NZ_VXJW01000004.1   87.599    1     87.599    2204      138        91.742    366711    5426      5563      +      9       Salmonella enterica
 ERR5396170.1000005   2516   306      325    5       GCF_008692845.1   NZ_VXJW01000004.1   87.599    1     87.599    2204      20         91.742    366711    5664      5683      +      9       Salmonella enterica
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Note: the column `species` is added by mapping genome ID (column `sgnm`) to taxonomic information.
```

Learn more [tutorials](http://bioinf.shenwei.me/LexicMap/tutorials/index/) and [usages](http://bioinf.shenwei.me/LexicMap/usage/lexicmap/).
