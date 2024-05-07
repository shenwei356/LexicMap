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
               --min-match-pident 50 --min-qcov-per-genome 70 --min-qcov-per-hsp 70 --top-n-genomes 1000

    - For longer queries like plasmids, returning all hits.

          lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
              --min-match-pident 50 --min-qcov-per-genome 0  --min-qcov-per-hsp 0  --top-n-genomes 0

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
   Query sequence is masked by the masks of the index. In other words, each mask captures the most similar k-mer which shares the longest prefix with the mask, and stores its position and strand information.
1. **Seeding:**
   For each mask, the captured k-mer is used to search seeds (captured k-mers in reference genomes) sharing prefixes of at least *p* bases.
    1. **Setting the search range**: Since the seeded k-mers are stored in lexicographic order, the k-mer matching turns into a range query.
       For example, for a query `CATGCT` requiring matching at least 4-bp prefix is equal to extract k-mers ranging from `CATGAA`, `CATGAC`, `CATGAG`, ...,  to `CATGTT`.
    2. **Finding the nearest smaller k-mer**: The index file of each seed data file stores a list (default 512) of k-mers and offsets in the data file, and the index is Loaded in RAM.
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

|Flag                       |Value               |Function                                              |Comment                                                                |
|:--------------------------|:-------------------|:-----------------------------------------------------|:----------------------------------------------------------------------|
|**`-w/--load-whole-seeds`**|                    |Load the whole seed data into memory for faster search|Use this if the index is not big and many queries are needed to search.|
|**`-n/--top-n-genomes`**   |Default 0, 0 for all|Keep top N genome matches for a query                 |                                                                       |
|**`-a/--all`**             |                    |Output more columns, e.g., matched sequences.         |                                                                       |

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
|**`-Q/--min-qcov-per-genome`**   |Default 0   |Minimum query coverage (percentage) per genome.                                                                   |       |
|**`-q/--min-qcov-per-hsp`**      |Default 0   |Minimum query coverage (percentage) per HSP.                                                                      |       |
|**`-l/--align-min-match-len`**   |Default 50  |Minimum aligned length in a HSP segment.                                                                          |       |
|**`-i/--align-min-match-pident`**|Default 50  |Minimum base identity (percentage) in a HSP segment.                                                              |       |
|`--align-band`                   |Default 100 |Band size in backtracking the score matrix.                                                                       |       |
|`--align-ext-len`                |Default 2000|Extend length of upstream and downstream of seed regions, for extracting query and target sequences for alignment.|       |
|`--align-max-gap`                |Default 50  |Maximum gap in a HSP segment.                                                                                     |       |
|`--align-max-mismatch`           |Default 50  |Maximum mismatch in a HSP segment.                                                                                |       |


{{< /tab>}}

{{< /tabs >}}


## Steps


- For short queries like genes or long reads, returning top N hits.

      lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
          --min-match-pident 50 \
          --min-qcov-per-genome 70 \
          --min-qcov-per-hsp 70 \
          --top-n-genomes 1000

- For longer queries like plasmids, returning all hits.

      lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
          --min-match-pident 50 \
          --min-qcov-per-genome 0  \
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

- Extracting similar sequences for a query gene.

      # search matches with query cover >= 90%
      lexicmap search -d gtdb_complete.lmi/ b.gene_E_faecalis_SecY.fasta --min-qcov-per-hsp 90 --all -o results.tsv

      # extract matched sequences as FASTA format
      sed 1d results.tsv | awk '{print ">"$5":"$13"-"$14":"$15"\n"$19;}' > results.fasta

      seqkit head -n 1 results.fasta | head -n 3
      >NZ_JALSCK010000007.1:39224-40522:-
      TTGTTCAAGCTATTAAAGAACGCCTTTAAAGTCAAAGACATTAGATCAAAAATCTTATTT
      ACAGTTTTAATCTTGTTTGTATTTCGCCTAGGTGCGCACATTACTGTGCCCGGGGTGAAT

## Output

### Alignment result relationship

    Query
    ├── Subject genome                             # A query might have one or more genome hits,
        ├── Subject sequence                       # in different sequences.
            ├── High-Scoring segment Pairs (HSP)   # HSP is a cluster of alignment segments.
                ├── HSP segment (not outputted)

Here, the defination of HSP is slightly different from that in BLAST.

> A High-scoring Segment Pair (HSP) is a local alignment with no gaps that achieves one of the highest alignment scores in a given search.
> https://www.ncbi.nlm.nih.gov/books/NBK62051/


### Output format

Tab-delimited format with 17+ columns, with 1-based positions.

    1.  query,    Query sequence ID.
    2.  qlen,     Query sequence length.
    3.  hits,     Number of Subject genomes.
    4.  sgenome,  Subject genome ID.
    5.  sseqid,   Subject sequence ID.
    6.  qcovGnm,  Query coverage (percentage) per genome: $(aligned bases in the genome)/$qlen.
    7.  hsp,      Nth HSP in the genome.
    8.  qcovHSP   Query coverage (percentage) per HSP: $(aligned bases in a HSP)/$qlen.
    9.  alenHSP,  Aligned length in the current HSP.
    10. pident,   Percentage of identical matches in the current HSP.
    11. qstart,   Start of alignment in query sequence.
    12. qend,     End of alignment in query sequence.
    13. sstart,   Start of alignment in subject sequence.
    14. send,     End of alignment in subject sequence.
    15. sstr,     Subject strand.
    16. slen,     Subject sequence length.
    17. seeds,    Number of seeds in the current HSP.
    18. qseq,     Aligned part of query sequence.   (optional with -a/--all)
    19. sseq,     Aligned part of subject sequence. (optional with -a/--all)

### Examples

{{< tabs "t2" >}}

{{< tab "A single-copy gene (SecY)" >}}

```plain
query                                      qlen   hits   sgenome           sseqid                 qcovGnm   hsp   qcovHSP   alenHSP   pident    qstart   qend   sstart   send     sstr   slen     seeds
----------------------------------------   ----   ----   ---------------   --------------------   -------   ---   -------   -------   -------   ------   ----   ------   ------   ----   ------   -----
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   3588   GCF_000763355.1   NZ_JPTY01000036.1      100.000   1     100.000   1299      100.000   1        1299   39247    40545    -      136653   19
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   3588   GCF_002944655.1   NZ_PTUM01000018.1      100.000   1     100.000   1299      100.000   1        1299   51657    52955    +      62674    24
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   3588   GCF_009735285.1   NZ_WMFU01000011.1      100.000   1     100.000   1299      100.000   1        1299   63314    64612    +      103681   28
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   3588   GCF_024339425.1   NZ_JAAOEK010000012.1   100.000   1     100.000   1299      100.000   1        1299   38854    40152    -      84246    27
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   3588   GCF_021122725.1   NZ_JADMEX010000011.1   100.000   1     100.000   1299      100.000   1        1299   38881    40179    -      84306    26
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   3588   GCF_014325545.1   NZ_JAAIGO010000013.1   100.000   1     100.000   1299      100.000   1        1299   38880    40178    -      84305    27
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   3588   GCF_024160695.1   NZ_JAMYXD010000005.1   100.000   1     100.000   1299      100.000   1        1299   39412    40710    -      257638   18
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   3588   GCF_902164315.1   NZ_CABHCQ010000007.1   100.000   1     100.000   1299      100.000   1        1299   132860   134158   +      173437   22
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   3588   GCF_003319505.1   NZ_KZ845864.1          100.000   1     100.000   1299      100.000   1        1299   231655   232953   +      271872   22
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   3588   GCF_902160975.1   NZ_CABGPX010000010.1   100.000   1     100.000   1299      100.000   1        1299   39275    40573    -      96796    20
```
{{< /tab>}}

{{< tab "A 16S rRNA gene" >}}


```plain
query                         qlen   hits     sgenome           sseqid                 qcovGnm   hsp   qcovHSP   alenHSP   pident    qstart   qend   sstart    send      sstr   slen      seeds
---------------------------   ----   ------   ---------------   --------------------   -------   ---   -------   -------   -------   ------   ----   -------   -------   ----   -------   -----
NC_000913.3:4166659-4168200   1542   294285   GCF_016774415.1   NZ_CP068706.1          100.000   1     100.000   1542      100.000   1        1542   354260    355801    +      4682670   24
NC_000913.3:4166659-4168200   1542   294285   GCF_016774415.1   NZ_CP068706.1          100.000   2     100.000   1542      100.000   1        1542   1053203   1054744   +      4682670   24
NC_000913.3:4166659-4168200   1542   294285   GCF_016774415.1   NZ_CP068706.1          100.000   3     100.000   1542      100.000   1        1542   3596202   3597743   -      4682670   24
NC_000913.3:4166659-4168200   1542   294285   GCF_016774415.1   NZ_CP068706.1          100.000   4     100.000   1542      100.000   1        1542   4398994   4400535   -      4682670   24
NC_000913.3:4166659-4168200   1542   294285   GCF_016774415.1   NZ_CP068706.1          100.000   5     100.000   1542      100.000   1        1542   4440484   4442025   -      4682670   24
NC_000913.3:4166659-4168200   1542   294285   GCF_016774415.1   NZ_CP068706.1          100.000   6     100.000   1542      100.000   1        1542   4571610   4573151   -      4682670   24
NC_000913.3:4166659-4168200   1542   294285   GCF_016774415.1   NZ_CP068706.1          100.000   7     100.000   1542      100.000   1        1542   4665333   4666874   -      4682670   24
NC_000913.3:4166659-4168200   1542   294285   GCA_003203355.1   BGNE01000097.1         100.000   1     100.000   1542      100.000   1        1542   30        1571      +      1585      25
NC_000913.3:4166659-4168200   1542   294285   GCF_015644465.1   NZ_JADOBZ010000093.1   100.000   1     100.000   1542      100.000   1        1542   15        1556      -      1702      32
NC_000913.3:4166659-4168200   1542   294285   GCF_015644465.1   NZ_JADOBZ010000115.1   100.000   2     7.328     113       100.000   1430     1542   1         113       +      617       1
```


{{< /tab>}}

{{< tab "A plasmid" >}}


```plain
query        qlen    hits    sgenome           sseqid          qcovGnm   hsp   qcovHSP   alenHSP   pident    qstart   qend    sstart    send      sstr   slen      seeds
----------   -----   -----   ---------------   -------------   -------   ---   -------   -------   -------   ------   -----   -------   -------   ----   -------   -----
CP115019.1   52830   58930   GCF_022759845.1   NZ_CP086533.1   97.473    1     97.473    51495     99.998    1        52830   1         51479     +      51479     553
CP115019.1   52830   58930   GCF_022759845.1   NZ_CP086533.1   97.473    2     1.552     820       100.000   9049     9868    23092     23911     +      51479     8
CP115019.1   52830   58930   GCF_022759845.1   NZ_CP086535.1   97.473    3     1.687     891       75.196    51686    52779   26583     27675     -      34058     1
CP115019.1   52830   58930   GCF_022759845.1   NZ_CP086534.1   97.473    4     0.502     265       100.000   19788    20052   29842     30106     +      47185     2
CP115019.1   52830   58930   GCF_022759845.1   NZ_CP086533.1   97.473    5     0.159     84        89.286    8348     8431    19574     19657     +      51479     3
CP115019.1   52830   58930   GCF_022759905.1   NZ_CP086545.1   97.473    1     97.473    51495     99.998    1        52830   1         51479     +      51479     553
CP115019.1   52830   58930   GCF_022759905.1   NZ_CP086545.1   97.473    2     1.552     820       100.000   9049     9868    23092     23911     +      51479     8
CP115019.1   52830   58930   GCF_022759905.1   NZ_CP086547.1   97.473    3     1.687     891       75.196    51686    52779   3843      4935      +      34058     1
CP115019.1   52830   58930   GCF_022759905.1   NZ_CP086546.1   97.473    4     0.502     265       100.000   19788    20052   29842     30106     +      47185     2
CP115019.1   52830   58930   GCF_022759905.1   NZ_CP086545.1   97.473    5     0.159     84        89.286    8348     8431    19574     19657     +      51479     3
CP115019.1   52830   58930   GCF_014826015.1   NZ_CP058621.1   97.473    1     97.473    51495     99.998    1        52830   1         51480     +      51480     504
CP115019.1   52830   58930   GCF_014826015.1   NZ_CP058618.1   97.473    2     2.497     1319      100.000   25153    26471   3019498   3020816   -      4718403   13
CP115019.1   52830   58930   GCF_014826015.1   NZ_CP058620.1   97.473    3     1.554     821       100.000   9049     9869    14895     15715     -      67321     5
CP115019.1   52830   58930   GCF_014826015.1   NZ_CP058620.1   97.473    4     1.554     821       100.000   9049     9869    20250     21070     -      67321     5
CP115019.1   52830   58930   GCF_014826015.1   NZ_CP058620.1   97.473    5     1.554     821       100.000   23721    24541   9689      10509     +      67321     5
CP115019.1   52830   58930   GCF_014826015.1   NZ_CP058620.1   97.473    6     1.554     821       100.000   23721    24541   14896     15716     -      67321     5
CP115019.1   52830   58930   GCF_014826015.1   NZ_CP058620.1   97.473    7     1.554     821       100.000   23721    24541   20251     21071     -      67321     5
CP115019.1   52830   58930   GCF_014826015.1   NZ_CP058620.1   97.473    8     1.552     820       100.000   9049     9868    9690      10509     +      67321     5
CP115019.1   52830   58930   GCF_014826015.1   NZ_CP058621.1   97.473    9     1.552     820       100.000   9049     9868    21166     21985     +      51480     5
CP115019.1   52830   58930   GCF_014826015.1   NZ_CP058622.1   97.473    10    1.552     820       100.000   23722    24541   3424      4243      +      9135      5
```

{{< /tab>}}

{{< tab "Long reads" >}}


Queries are a few Nanopore Q20 reads from a mock metagenomic community.

```plain
    query                qlen   hits   sgenome           sseqid                 qcovGnm   hsp   qcovHSP   alenHSP   pident    qstart   qend   sstart    send      sstr   slen      seeds
    ------------------   ----   ----   ---------------   --------------------   -------   ---   -------   -------   -------   ------   ----   -------   -------   ----   -------   -----
    ERR5396170.1000006   796    4      GCF_013394085.1   NZ_CP040910.1          96.231    1     96.231    766       97.389    31       796    1138938   1139706   +      1887974   24
    ERR5396170.1000006   796    4      GCF_013394085.1   NZ_CP040910.1          96.231    2     95.226    758       96.834    39       796    1280352   1282817   +      1887974   21
    ERR5396170.1000006   796    4      GCF_013394085.1   NZ_CP040910.1          96.231    3     95.226    758       97.493    39       796    32646     33406     +      1887974   24
    ERR5396170.1000006   796    4      GCF_013394085.1   NZ_CP040910.1          96.231    4     95.226    758       97.493    39       796    134468    135228    -      1887974   24
    ERR5396170.1000006   796    4      GCF_013394085.1   NZ_CP040910.1          96.231    5     95.226    758       97.361    39       796    1768935   1769695   +      1887974   24
    ERR5396170.1000006   796    4      GCF_013394085.1   NZ_CP040910.1          96.231    6     95.226    758       97.230    39       796    242012    242772    -      1887974   24
    ERR5396170.1000006   796    4      GCF_013394085.1   NZ_CP040910.1          96.231    7     95.226    758       96.834    39       796    154380    155137    -      1887974   23
    ERR5396170.1000006   796    4      GCF_003344625.1   NZ_QPKJ02000188.1      81.910    1     81.910    652       99.540    66       717    1         653       -      826       21
    ERR5396170.1000006   796    4      GCF_009663775.1   NZ_RDBR01000008.1      91.332    1     91.332    727       86.657    39       796    21391     22151     -      52610     7
    ERR5396170.1000006   796    4      GCF_001591685.1   NZ_BCVJ01000102.1      87.940    1     87.940    700       88.429    66       796    434       1165      +      1933      4
    ERR5396170.1000017   516    1      GCF_013394085.1   NZ_CP040910.1          94.574    1     94.574    488       100.000   27       514    293509    293996    +      1887974   11
    ERR5396170.1000012   848    1      GCF_013394085.1   NZ_CP040910.1          95.165    1     95.165    807       93.804    22       828    190329    191136    -      1887974   21
    ERR5396170.1000052   330    1      GCF_013394085.1   NZ_CP040910.1          90.000    1     90.000    297       100.000   27       323    1161955   1162251   -      1887974   3
    ERR5396170.1000000   698    1      GCF_001457615.1   NZ_LN831024.1          85.673    1     85.673    598       97.157    53       650    4452083   4452685   +      6316979   10
```

{{< /tab>}}

{{< /tabs >}}

Search results (TSV format) above are formatted with [csvtk pretty](https://github.com/shenwei356/csvtk).
