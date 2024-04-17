---
title: Searching
weight: 10
---

## Table of contents

{{< toc format=html >}}

## TL;DR

1. [Build](https://bioinf.shenwei.me/LexicMap/tutorials/index/) or download a LexicMap index.

1. Run:

    - For short queries like genes or long reads, returning top N hits.

          lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
              --min-qcov-per-genome 70 --min-match-pident 70 --min-qcov-per-hsp 70 --top-n-genomes 500

    - For longer queries like plasmids, returning all hits.

          lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
              --min-qcov-per-genome 50 --min-match-pident 70 --min-qcov-per-hsp 0  --top-n-genomes 0

## Input

{{< hint type=note >}}
**Query length**\
LexicMap only supports long (>=500 bp) reads or gene/genome/virus/plasmid/phage sequences.\
However, some short queries can also be aligned.
{{< /hint >}}

Input should be (gzipped) FASTA or FASTQ records from files or STDIN.


## Hardware requirements

LexicMap is designed to provide fast and low-memory sequence alignment against millions of prokaryotic genomes.

- **CPU:**
    - No specific requirements on CPU type and instruction sets. Both x86 and ARM chips are supported.
    - More is better as LexicMap is a CPU-intensive software. It uses all CPUs by default (`-j/--threads`).
- **RAM**
    - More RAM (> 16 GB) is preferred. The memory usage in searching is mainly related to:
        - The number of matched genomes and sequences.
        - The length of query sequences.
- **Disk**
    - Sufficient space is required to store the index size.
    - No temporary files are generated during searching.

## Algorithm

<img src="/LexicMap/searching.svg" alt="" width="900"/>

1. **Masking:**
   Query sequence is masked by the masks of the index. In other word, each mask captures the most similar k-mer and stores its posistion and strand information.
1. **Seeding:**
   For each mask, the captured k-mer is used to search seeds (captured k-mers in reference genomes) sharing prefixes of at least *p* bases.
    1. **Setting the search range**: Since the seeded k-mers are stored in lexicographic order, the k-mer matching turns into a range query.
       For example, for a query `CATGCT` requiring matching at least 4-bp prefix is equal to extract k-mers ranging from `CATGAA`, `CATGAC`, `CATGAG`, ...,  to `CATGTT`.
    2. **Finding the nearest offset**: The index file of each seed data file stores a list (default 512) of k-mers and offsets in the data file, and the index is Loaded in RAM.
       The nearest k-mer smaller than the range start k-mer (`CATGAA`) is found by binary search, i.e., `CATCAC` (blue text in the fingure), and the offset is returned as the start position in traversing the seed data file.
    3. **Retrieving seed data**: Seed k-mers are read from the file and checked one by one, and k-mers in the search range are returned, along with the k-mer information (genome batch, genome number, location, and strand).
1. **Chaining:**
    1. Seeding results, i.e., anchors (matched k-mers from the query and subject sequence), are summarized by genome, and deduplicated.
    2. Performing chaining.
1. **Alignment** for each chain.
    1. Extending the anchor region. for extracting sequences from the query and reference genome. For example, extending 2 kb in upstream and downstream of anchor region.
    2. Fast alignment of query and subject sequences.
    3. Filtering aligned segments and the whole HSPs (all alignment segments) based on user options.
       - For these HSPs that accross more than one reference sequences, splitting them into multiple HSPs.


## Parameters

**Flags in bold text** are important and frequently used.

{{< tabs "t1" >}}

{{< tab "General" >}}

|Flag                       |Value                 |Function                                              |Comment                                                                |
|:--------------------------|:---------------------|:-----------------------------------------------------|:----------------------------------------------------------------------|
|**`-w/--load-whole-seeds`**|                      |Load the whole seed data into memory for faster search|Use this if the index is not big and many queries are needed to search.|
|**`-n/--top-n-genomes`**   |Default 500, 0 for all|Keep top N genome matches for a query                 |                                                                       |

{{< /tab>}}

{{< tab "Chaining" >}}

|Flag                              |Value        |Function                                                               |Comment|
|:---------------------------------|:------------|:----------------------------------------------------------------------|:------|
|**`-p, --seed-min-prefix`**       |Default 15   |Minimum length of shared substrings (anchors).                         |       |
|**`-P, --seed-min-single-prefix`**|Default 20   |Minimum length of shared substrings (anchors) if there's only one pair.|       |
|`--seed-max-dist`                 |Default 10000|Max distance between seeds in seed chaining.                           |       |
|`--seed-max-gap`                  |Default 2000 |Max gap in seed chaining.                                              |       |
|`-m/--seed-max-mismatch`          |Default -1   |Minimum mismatch between non-prefix regions of shared substrings.      |       |

{{< /tab>}}

{{< tab "Alignment" >}}

|Flag                             |Value       |Function                                                                                                          |Comment|
|:--------------------------------|:-----------|:-----------------------------------------------------------------------------------------------------------------|:------|
|**`-Q/--min-qcov-per-genome`**   |Default 50  |Minimum query coverage (percentage) per genome.                                                                   |       |
|**`-q/--min-qcov-per-hsp`**      |Default 0   |Minimum query coverage (percentage) per HSP.                                                                      |       |
|**`-l/--align-min-match-len`**   |Default 50  |Minimum aligned length in a HSP segment.                                                                          |       |
|**`-i/--align-min-match-pident`**|Default 70  |Minimum base identity (percentage) in a HSP segment.                                                              |       |
|`--align-band`                   |Default 100 |Band size in backtracking the score matrix.                                                                       |       |
|`--align-ext-len`                |Default 2000|Extend length of upstream and downstream of seed regions, for extracting query and target sequences for alignment.|       |
|`--align-max-gap`                |Default 50  |Maximum gap in a HSP segment.                                                                                     |       |
|`--align-max-mismatch`           |Default 50  |Maximum mismatch in a HSP segment.                                                                                |       |


{{< /tab>}}

{{< /tabs >}}


## Steps


- For short queries like genes or long reads, returning top N hits.

      lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
          --min-qcov-per-genome 70 \
          --min-match-pident 70 \
          --min-qcov-per-hsp 70 \
          --top-n-genomes 500

- For longer queries like plasmids, returning all hits.

      lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
          --min-qcov-per-genome 50 \
          --min-match-pident 70 \
          --min-qcov-per-hsp 0 \
          --top-n-genomes 0


{{< expand "Click to show the log of a demo run." "..." >}}

    $ lexicmap search -d demo.lmi/  q.gene.fasta -o q.gene.fasta.lexicmap.tsv
    09:32:55.551 [INFO] LexicProf v0.3.0
    09:32:55.551 [INFO]   https://github.com/shenwei356/LexicMap
    09:32:55.551 [INFO]
    09:32:55.551 [INFO] checking input files ...
    09:32:55.551 [INFO]   1 input file(s) given
    09:32:55.551 [INFO]
    09:32:55.551 [INFO] loading index: demo.lmi/
    09:32:55.551 [INFO]   reading masks...
    09:32:55.552 [INFO]   reading indexes of seeds (k-mer-value) data...
    09:32:55.555 [INFO]   creating genome reader pools, each batch with 16 readers...
    09:32:55.555 [INFO] index loaded in 4.192051ms
    09:32:55.555 [INFO]
    09:32:55.555 [INFO] searching ...

    09:32:55.596 [INFO]
    09:32:55.596 [INFO] processed queries: 1, speed: 1467.452 queries per minute
    09:32:55.596 [INFO] 100.0000% (1/1) queries matched
    09:32:55.596 [INFO] done searching
    09:32:55.596 [INFO] search results saved to: q.gene.fasta.lexicmap.tsv
    09:32:55.596 [INFO]
    09:32:55.596 [INFO] elapsed time: 45.230604ms
    09:32:55.596 [INFO]

{{< /expand >}}


## Output

### Alignment result relationship

    Query
    ├── Subject genome                             # A query might have one or more genome hits,
        ├── Subject sequence                       # in different sequences.
            ├── High-Scoring segment Pairs (HSP)   # HSP is a cluster of alignment segments.
                ├── HSP segment                    # A local alignment with no gaps.

Here, the defination of HSP is slightly different from that in BLAST.

> A High-scoring Segment Pair (HSP) is a local alignment with no gaps that achieves one of the highest alignment scores in a given search.
> https://www.ncbi.nlm.nih.gov/books/NBK62051/


### Output format

Tab-delimited format with 18 columns. (The positions are 1-based).

    1.  query,    Query sequence ID.
    2.  qlen,     Query sequence length.
    3.  qstart,   Start of alignment in query sequence.
    4.  qend,     End of alignment in query sequence.
    5.  sgnms,    The number of Subject genomes.
    6.  sgnm,     Subject genome ID.
    7.  seqid,    Subject sequence ID.
    8.  qcovGnm,  Query coverage (percentage) per genome: $(aligned bases in the genome)/$qlen.
    9.  hsp,      Nth HSP in the genome.
    10. qcovHSP   Query coverage (percentage) per HSP: $(aligned bases in a HSP)/$qlen.
    11. alen,     Aligned length in the current HSP, a HSP might have >=1 HSP segments.
    12. alenSeg,  Aligned length in the current HSP segment.
    13. pident,   Percentage of identical matches in the current HSP segment.
    14. slen,     Subject sequence length.
    15. sstart,   Start of HSP segment in subject sequence.
    16. send,     End of HSP segment in subject sequence.
    17. sstr,     Subject strand.
    18. seeds,    Number of seeds in the current HSP.

### Examples

{{< tabs "t2" >}}

{{< tab "A 16S rRNA gene" >}}

The (part) result shows the 16S rRNA gene has 12 genome hits (column `sgnms`). And in genome `GCF_003697165.2`, it has 7 highly similar matches with query coverage per HSP (column `qcovHSP`) > 99% and percentage of identity (`pident`) > 99%. It makes sense as 16S rRNA genes might have multiple copies in a genome.

```plain
query                         qlen   qstart   qend   sgnms   sgnm              seqid           qcovGnm   hsp   qcovHSP   alenHSP   alenSeg   pident   slen      sstart    send      sstr   seeds
---------------------------   ----   ------   ----   -----   ---------------   -------------   -------   ---   -------   -------   -------   ------   -------   -------   -------   ----   -----
NC_000913.3:4166659-4168200   1542   1        1542   12      GCF_003697165.2   NZ_CP033092.2   100.000   1     99.287    1542      1542      99.287   4903501   4591684   4593225   -      18
NC_000913.3:4166659-4168200   1542   1        1542   12      GCF_003697165.2   NZ_CP033092.2   100.000   2     99.287    1542      1542      99.287   4903501   4551515   4553056   -      18
NC_000913.3:4166659-4168200   1542   1        1542   12      GCF_003697165.2   NZ_CP033092.2   100.000   3     99.287    1542      1542      99.287   4903501   3780640   3782181   -      18
NC_000913.3:4166659-4168200   1542   1        1542   12      GCF_003697165.2   NZ_CP033092.2   100.000   4     99.287    1542      1542      99.287   4903501   1285123   1286664   +      18
NC_000913.3:4166659-4168200   1542   1        1542   12      GCF_003697165.2   NZ_CP033092.2   100.000   5     99.287    1542      1542      99.287   4903501   4844587   4846128   -      18
NC_000913.3:4166659-4168200   1542   1        1542   12      GCF_003697165.2   NZ_CP033092.2   100.000   6     99.287    1542      1542      99.287   4903501   458559    460100    +      18
NC_000913.3:4166659-4168200   1542   1        1542   12      GCF_003697165.2   NZ_CP033092.2   100.000   7     99.092    1542      1542      99.092   4903501   4726193   4727734   -      18
NC_000913.3:4166659-4168200   1542   1        1542   12      GCF_002950215.1   NZ_CP026788.1   100.000   1     99.027    1542      1542      99.027   4659463   3216505   3218046   +      15
NC_000913.3:4166659-4168200   1542   1        1542   12      GCF_002950215.1   NZ_CP026788.1   100.000   2     98.962    1542      1542      98.962   4659463   3396068   3397609   +      16
```


{{< /tab>}}

{{< tab "A plasmid" >}}

This plasmid has 19,265 genome hits in the index (2.34 millions prokaryotic genome from Genbank and RefSeq).
It is highly similar to the sequence `CP058621.1` with query coverage of `99.388` and percentage of identity of `99.967`.
They also have similar sequence length `qlen` (`51466`) vs `slen` (`51480`).
Besides there are a few duplicated sequences (might be multiple-copy genes) in the query and matches,
e.g., a ~821-bp sequence (43100-43920) appears in multiple positions of sequences `CP058620.1`.

```plain
query     qlen    qstart   qend    sgnms   sgnm              seqid           qcovGnm   hsp   qcovHSP   alenHSP   alenFrag   pident    slen      sstart    send      sstr   seeds
-------   -----   ------   -----   -----   ---------------   -------------   -------   ---   -------   -------   --------   -------   -------   -------   -------   ----   -----
plasmid   51466   296      51450   19265   GCA_014826015.1   CP058621.1      100.000   1     99.388    51168     51168      99.967    51480     313       51480     -      497
plasmid   51466   27861    29179   19265   GCA_014826015.1   CP058618.1      100.000   2     2.563     1319      1319       100.000   4718403   3019498   3020816   +      9
plasmid   51466   43100    43920   19265   GCA_014826015.1   CP058620.1      100.000   3     1.673     861       821        100.000   67321     20250     21070     +      6
plasmid   51466   45163    45202   19265   GCA_014826015.1   CP058620.1      100.000   3     1.673     861       40         100.000   67321     22219     22258     +      6
plasmid   51466   43100    43920   19265   GCA_014826015.1   CP058620.1      100.000   4     1.595     821       821        100.000   67321     14895     15715     +      5
plasmid   51466   43101    43920   19265   GCA_014826015.1   CP058622.1      100.000   5     1.593     820       820        100.000   9135      3424      4243      -      5
plasmid   51466   43101    43920   19265   GCA_014826015.1   CP058620.1      100.000   6     1.593     820       820        100.000   67321     9690      10509     -      5
plasmid   51466   29791    30608   19265   GCA_014826015.1   CP058620.1      100.000   7     1.570     821       821        98.417    67321     20251     21071     +      5
plasmid   51466   29791    30608   19265   GCA_014826015.1   CP058620.1      100.000   8     1.570     821       821        98.417    67321     14896     15716     +      5
plasmid   51466   29791    30608   19265   GCA_014826015.1   CP058620.1      100.000   9     1.570     821       821        98.417    67321     9689      10509     -      5
plasmid   51466   29791    30607   19265   GCA_014826015.1   CP058622.1      100.000   10    1.568     820       820        98.415    9135      3424      4243      -      5
plasmid   51466   29791    30607   19265   GCA_014826015.1   CP058621.1      100.000   11    1.568     820       820        98.415    51480     7844      8663      -      5
plasmid   51466   43100    43865   19265   GCA_014826015.1   CP058622.1      100.000   12    1.488     766       766        100.000   9135      8370      9135      +      4
plasmid   51466   18302    18401   19265   GCA_014826015.1   CP058621.1      100.000   14    0.194     100       100        100.000   51480     27281     27380     +      2
plasmid   51466   34042    34125   19265   GCA_014826015.1   CP058621.1      100.000   15    0.146     84        84         89.286    51480     7143      7226      -      1
plasmid   51466   43866    43921   19265   GCA_014826015.1   CP058622.1      100.000   16    0.109     56        56         100.000   9135      1         56        +      2
plasmid   51466   30553    30607   19265   GCA_014826015.1   CP058622.1      100.000   17    0.107     55        55         100.000   9135      1         55        +      1
plasmid   51466   296      51450   19265   GCF_014826015.1   NZ_CP058621.1   100.000   1     99.388    51168     51168      99.967    51480     313       51480     -      497
plasmid   51466   27861    29179   19265   GCF_014826015.1   NZ_CP058618.1   100.000   2     2.563     1319      1319       100.000   4718403   3019498   3020816   +      9
```

{{< /tab>}}

{{< tab "Long reads" >}}


Queries are a few Nanopore Q20 read. The column `species` is added by mapping genome ID (column `sgnm`) to taxonomic information.

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
```

{{< /tab>}}

{{< /tabs >}}

Search results (TSV format) above are formatted with [csvtk pretty](https://github.com/shenwei356/csvtk).
