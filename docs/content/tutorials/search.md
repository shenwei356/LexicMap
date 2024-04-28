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
              --min-qcov-per-genome 70 --min-match-pident 70 --min-qcov-per-hsp 70 --top-n-genomes 1000

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
   Query sequence is masked by the masks of the index. In other word, each mask captures the most similar k-mer which shares the longest prefix with the mask, and stores its posistion and strand information.
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
    5.  hits,     The number of Subject genomes.
    6.  sgenome,  Subject genome ID.
    7.  sseqid,   Subject sequence ID.
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

The (part) result shows the 16S rRNA gene has 8 genome hits (column `sgenome`). And in genome `GCF_003697165.2`, it has 7 highly similar matches with query coverage per HSP (column `qcovHSP`) of 100% and percentage of identity (`pident`) > 99%. It makes sense as 16S rRNA genes might have multiple copies in a genome.

```plain
query                         qlen   qstart   qend   hits   sgenome           sseqid          qcovGnm   hsp   qcovHSP   alenHSP   alenSeg   pident   slen      sstart    send      sstr   seeds
---------------------------   ----   ------   ----   ----   ---------------   -------------   -------   ---   -------   -------   -------   ------   -------   -------   -------   ----   -----
NC_000913.3:4166659-4168200   1542   1        1542   8      GCF_003697165.2   NZ_CP033092.2   100.000   1     100.000   1542      1542      99.287   4903501   4844587   4846128   -      26
NC_000913.3:4166659-4168200   1542   1        1542   8      GCF_003697165.2   NZ_CP033092.2   100.000   2     100.000   1542      1542      99.287   4903501   4591684   4593225   -      26
NC_000913.3:4166659-4168200   1542   1        1542   8      GCF_003697165.2   NZ_CP033092.2   100.000   3     100.000   1542      1542      99.287   4903501   4551515   4553056   -      26
NC_000913.3:4166659-4168200   1542   1        1542   8      GCF_003697165.2   NZ_CP033092.2   100.000   4     100.000   1542      1542      99.287   4903501   3780640   3782181   -      26
NC_000913.3:4166659-4168200   1542   1        1542   8      GCF_003697165.2   NZ_CP033092.2   100.000   5     100.000   1542      1542      99.287   4903501   458559    460100    +      26
NC_000913.3:4166659-4168200   1542   1        1542   8      GCF_003697165.2   NZ_CP033092.2   100.000   6     100.000   1542      1542      99.287   4903501   1285123   1286664   +      26
NC_000913.3:4166659-4168200   1542   1        1542   8      GCF_003697165.2   NZ_CP033092.2   100.000   7     100.000   1542      1542      99.092   4903501   4726193   4727734   -      26
NC_000913.3:4166659-4168200   1542   1        1542   8      GCF_002950215.1   NZ_CP026788.1   100.000   1     100.000   1542      1542      99.027   4659463   3216505   3218046   +      25
NC_000913.3:4166659-4168200   1542   1        1542   8      GCF_002950215.1   NZ_CP026788.1   100.000   2     100.000   1542      1542      98.962   4659463   3396068   3397609   +      26
```


{{< /tab>}}

{{< tab "A plasmid" >}}

This plasmid has 17,716 genome hits in the index (2.34 millions prokaryotic genome from Genbank and RefSeq).

```plain
query     qlen    qstart   qend    hits    sgenome           sseqid                 qcovGnm   hsp   qcovHSP   alenHSP   alenSeg   pident    slen      sstart    send      sstr   seeds
-------   -----   ------   -----   -----   ---------------   --------------------   -------   ---   -------   -------   -------   -------   -------   -------   -------   ----   -----
plasmid   51466   1        50164   17716   GCA_032192075.1   JAVTRN010000003.1      100.000   1     100.000   51466     50164     100.000   51479     1303      51479     +      487
plasmid   51466   50165    51466   17716   GCA_032192075.1   JAVTRN010000003.1      100.000   1     100.000   51466     1302      100.000   51479     1         1302      +      487
plasmid   51466   18302    18401   17716   GCA_032192075.1   JAVTRN010000003.1      100.000   2     0.194     100       100       100.000   51479     25698     25797     -      2
plasmid   51466   34042    34125   17716   GCA_032192075.1   JAVTRN010000003.1      100.000   3     0.163     84        84        89.286    51479     45852     45935     +      1
plasmid   51466   1        50164   17716   GCF_032192075.1   NZ_JAVTRN010000003.1   100.000   1     100.000   51466     50164     100.000   51479     1303      51479     +      487
plasmid   51466   50165    51466   17716   GCF_032192075.1   NZ_JAVTRN010000003.1   100.000   1     100.000   51466     1302      100.000   51479     1         1302      +      487
plasmid   51466   18302    18401   17716   GCF_032192075.1   NZ_JAVTRN010000003.1   100.000   2     0.194     100       100       100.000   51479     25698     25797     -      2
plasmid   51466   34042    34125   17716   GCF_032192075.1   NZ_JAVTRN010000003.1   100.000   3     0.163     84        84        89.286    51479     45852     45935     +      1
plasmid   51466   1        50164   17716   GCA_030863645.1   CP114982.1             100.000   1     100.000   51466     50164     99.994    51479     1303      51479     +      530
plasmid   51466   50165    51466   17716   GCA_030863645.1   CP114982.1             100.000   1     100.000   51466     1302      100.000   51479     1         1302      +      530
plasmid   51466   27858    29180   17716   GCA_030863645.1   CP114979.1             100.000   2     2.571     1323      1323      100.000   4731337   1893701   1895023   -      12
plasmid   51466   43101    43921   17716   GCA_030863645.1   CP114980.1             100.000   3     1.595     821       821       100.000   165063    42818     43638     -      5
plasmid   51466   43100    43920   17716   GCA_030863645.1   CP114980.1             100.000   4     1.595     821       821       100.000   165063    50985     51805     -      5
plasmid   51466   43101    43920   17716   GCA_030863645.1   CP114981.1             100.000   5     1.593     820       820       100.000   82723     55698     56517     +      5
plasmid   51466   43101    43920   17716   GCA_030863645.1   CP114980.1             100.000   6     1.593     820       820       100.000   165063    153475    154294    -      5
plasmid   51466   43101    43920   17716   GCA_030863645.1   CP114980.1             100.000   7     1.593     820       820       100.000   165063    160027    160846    -      5
plasmid   51466   43101    43452   17716   GCA_030863645.1   CP114980.1             100.000   7     1.593     820       352       99.716    165063    158276    158627    -      5
plasmid   51466   43101    43920   17716   GCA_030863645.1   CP114980.1             100.000   8     1.593     820       820       100.000   165063    145617    146436    +      5
plasmid   51466   43101    43920   17716   GCA_030863645.1   CP114980.1             100.000   9     1.593     820       820       99.878    165063    157808    158627    -      4
plasmid   51466   43555    43920   17716   GCA_030863645.1   CP114980.1             100.000   9     1.593     820       366       100.000   165063    160027    160392    -      4
plasmid   51466   29790    30607   17716   GCA_030863645.1   CP114980.1             100.000   10    1.589     818       818       98.778    165063    42819     43639     -      5
plasmid   51466   29791    30608   17716   GCA_030863645.1   CP114980.1             100.000   11    1.589     818       818       98.778    165063    50984     51804     -      5
plasmid   51466   29791    30608   17716   GCA_030863645.1   CP114980.1             100.000   12    1.589     818       818       98.778    165063    145617    146437    +      5
plasmid   51466   29790    30607   17716   GCA_030863645.1   CP114981.1             100.000   13    1.589     818       818       98.778    82723     55697     56517     +      5
plasmid   51466   29791    30608   17716   GCA_030863645.1   CP114980.1             100.000   14    1.589     818       818       98.778    165063    160026    160846    -      5
plasmid   51466   29791    30142   17716   GCA_030863645.1   CP114980.1             100.000   14    1.589     818       352       99.148    165063    158276    158627    -      5
plasmid   51466   29791    30608   17716   GCA_030863645.1   CP114980.1             100.000   15    1.589     818       818       98.778    165063    153474    154294    -      5
plasmid   51466   29791    30609   17716   GCA_030863645.1   CP114980.1             100.000   16    1.591     819       819       98.657    165063    157806    158627    -      4
plasmid   51466   30242    30608   17716   GCA_030863645.1   CP114980.1             100.000   16    1.591     819       367       100.000   165063    160026    160392    -      4
plasmid   51466   29791    30607   17716   GCA_030863645.1   CP114982.1             100.000   17    1.587     817       817       98.776    51479     44415     45234     +      5
plasmid   51466   29791    30607   17716   GCA_030863645.1   CP114980.1             100.000   18    1.587     817       817       98.776    165063    31910     32729     +      5
plasmid   51466   18302    18401   17716   GCA_030863645.1   CP114982.1             100.000   19    0.194     100       100       100.000   51479     25698     25797     -      2
plasmid   51466   44536    44622   17716   GCA_030863645.1   CP114981.1             100.000   20    0.169     87        87        87.356    82723     49696     49782     +      1
plasmid   51466   34042    34125   17716   GCA_030863645.1   CP114982.1             100.000   21    0.163     84        84        89.286    51479     45852     45935     +      1
plasmid   51466   34042    34125   17716   GCA_030863645.1   CP114981.1             100.000   22    0.163     84        84        86.905    82723     49698     49781     +      1
```

{{< /tab>}}

{{< tab "Long reads" >}}


Queries are a few Nanopore Q20 reads from a mock metagenomic community.

```plain
query                qlen   qstart   qend   hits   sgenome           sseqid              qcovGnm   hsp   qcovHSP   alenHSP   alenSeg   pident    slen      sstart    send      sstr   seeds
------------------   ----   ------   ----   ----   ---------------   -----------------   -------   ---   -------   -------   -------   -------   -------   -------   -------   ----   -----
ERR5396170.1000017   516    27       514    1      GCF_013394085.1   NZ_CP040910.1       94.574    1     94.574    488       488       100.000   1887974   293509    293996    +      3
ERR5396170.1000047   960    24       812    1      GCF_001027105.1   NZ_CP011526.1       90.521    1     90.521    869       789       89.480    2755072   2204718   2205520   -      6
ERR5396170.1000047   960    881      960    1      GCF_001027105.1   NZ_CP011526.1       90.521    1     90.521    869       80        100.000   2755072   2204568   2204647   -      6
ERR5396170.1000016   740    71       733    1      GCF_013394085.1   NZ_CP040910.1       89.595    1     89.595    663       663       98.492    1887974   13515     14177     +      12
ERR5396170.1000000   698    53       650    1      GCF_001457615.1   NZ_LN831024.1       85.673    1     85.673    598       598       97.324    6316979   4452083   4452685   +      4
ERR5396170.1000005   2516   38       2510   5      GCF_000006945.2   NC_003197.2         98.291    1     98.291    2473      2473      99.151    4857450   3198806   3201283   +      15
ERR5396170.1000005   2516   38       2497   5      GCF_008692785.1   NZ_VXJV01000001.1   97.774    1     97.774    2460      2460      96.098    797633    423400    425865    +      14
ERR5396170.1000005   2516   40       2510   5      GCA_900478215.1   LS483478.1          98.211    1     98.211    2471      2471      95.548    4624613   785866    788342    -      13
ERR5396170.1000005   2516   1350     2497   5      GCF_008692845.1   NZ_VXJW01000004.1   86.765    1     86.765    2183      1148      95.557    366711    6705      7855      +      12
ERR5396170.1000005   2516   634      1309   5      GCF_008692845.1   NZ_VXJW01000004.1   86.765    1     86.765    2183      676       89.053    366711    5991      6664      +      12
ERR5396170.1000005   2516   387      608    5      GCF_008692845.1   NZ_VXJW01000004.1   86.765    1     86.765    2183      222       85.135    366711    5745      5965      +      12
ERR5396170.1000005   2516   69       205    5      GCF_008692845.1   NZ_VXJW01000004.1   86.765    1     86.765    2183      137       83.212    366711    5426      5563      +      12
ERR5396170.1000005   2516   1830     2263   5      GCF_000252995.1   NC_015761.1         78.378    1     78.378    1972      434       97.696    4460105   2898281   2898717   +      7
ERR5396170.1000005   2516   2307     2497   5      GCF_000252995.1   NC_015761.1         78.378    1     78.378    1972      191       76.963    4460105   2898761   2898951   +      7
ERR5396170.1000005   2516   415      938    5      GCF_000252995.1   NC_015761.1         78.378    1     78.378    1972      524       86.641    4460105   2896865   2897391   +      7
ERR5396170.1000005   2516   1113     1807   5      GCF_000252995.1   NC_015761.1         78.378    1     78.378    1972      695       71.511    4460105   2897564   2898258   +      7
ERR5396170.1000005   2516   961      1088   5      GCF_000252995.1   NC_015761.1         78.378    1     78.378    1972      128       85.156    4460105   2897414   2897541   +      7
```

{{< /tab>}}

{{< /tabs >}}

Search results (TSV format) above are formatted with [csvtk pretty](https://github.com/shenwei356/csvtk).
