---
title: Step 2. Searching
weight: 10
---

## Table of contents

{{< toc format=html >}}

## TL;DR

1. [Build](https://bioinf.shenwei.me/LexicMap/tutorials/index/) a LexicMap index.

1. Run:

    - For short queries like genes or long reads, returning top N hits.

          lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
               --min-qcov-per-hsp 70 --min-qcov-per-genome 70 --top-n-genomes 1000

    - For longer queries like plasmids, returning all hits.

          lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
              --min-qcov-per-hsp 0   --min-qcov-per-genome 0  --top-n-genomes 0

## Input

{{< hint type=note >}}
**Query length**\
LexicMap is mainly designed for sequence alignment with a small number of queries (gene/plasmid/virus/phage sequences) longer than 200 bp by default.
However, short queries can also be aligned.
{{< /hint >}}

Input should be (gzipped) FASTA or FASTQ records from files or STDIN.


## Hardware requirements

See [benchmark of index building](https://bioinf.shenwei.me/LexicMap/introduction/#searching).

LexicMap is designed to provide fast and low-memory sequence alignment against millions of prokaryotic genomes.

- **CPU:**
    - No specific requirements on CPU type and instruction sets. Both x86 and ARM chips are supported.
    - More is better as LexicMap is a CPU-intensive software. It uses all CPUs by default (`-j/--threads`).
- **RAM**
    - More RAM (> 16 GB) is preferred. The memory usage in searching is mainly related to:
        - The number of matched genomes and sequences.
        - The length of query sequences.
        - Similarities between query and target sequences.
        - The number of threads. It uses all CPUs by default (`-j/--threads`).
- **Disk**
    - SSD disks are preferred to store the index size, while HDD disks are also fast enough.
    - No temporary files are generated during searching.

## Algorithm

<img src="/LexicMap/searching.svg" alt="" width="900"/>

{{< expand "Click to show details." "..." >}}

1. **Masking:**
   Query sequence is masked by the masks of the index. In other words, each mask captures the most similar k-mer which shares the longest prefix with the mask, and stores its position and strand information.
1. **Seeding:**
   For each mask, the captured k-mer is used to search seeds (captured k-mers in reference genomes) sharing **prefixes or suffixes** of at least *p* bases.
    1. Prefix matching
        1. **Setting the search range**: Since the seeded k-mers are stored in lexicographic order, the k-mer matching turns into a range query.
        For example, for a query `CATGCT` requiring matching at least 4-bp prefix is equal to extract k-mers ranging from `CATGAA`, `CATGAC`, `CATGAG`, ...,  to `CATGTT`.
        2. **Retrieving search start point**: The index file of each seed data file stores some k-mers' offsets in the data file, and the index is loaded in RAM.
        3. **Retrieving seed data**: Seed k-mers are read from the file and checked one by one, and k-mers in the search range are returned, along with the k-mer information (genome batch, genome number, location, and strand).
    1. Suffix matching
        1. Reversing the query k-mer and performing prefix matching, returning seeds of reversed k-mers (see indexing algorithm).
1. **Chaining:**
    1. Seeding results, i.e., anchors (matched k-mers from the query and subject sequence), are summarized by genome, and deduplicated.
    2. Performing chaining (see the paper).
1. **Alignment** for each chain.
    1. Extending the anchor region. for extracting sequences from the query and reference genome. For example, extending 1 kb in upstream and downstream of anchor region.
    1. Performing pseudo-alignment with extended query and subject sequences, for find similar regions.
       - For these similar regions that accross more than one reference sequences, splitting them into multiple ones.
    2. Fast alignment of query and subject sequence regions with [our implementation](https://github.com/shenwei356/wfa) of [Wavefront alignment algorithm](https://doi.org/10.1093/bioinformatics/btaa777).
    3. Filtering alignments based on user options.

{{</ expand >}}

## Parameters

**Flags in bold text** are important and frequently used.

{{< tabs "t1" >}}

{{< tab "General" >}}

|Flag                       |Value                      |Function                                                   |Comment                                                                                                                                                                                                                                                                 |
|:--------------------------|:--------------------------|:----------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|**`-j/--threads`**         |Default: all available cpus|Number of CPU cores to use.                                |The value should be >= the number of seed chunk files (“chunks” in info.toml, set by `-c/--chunks` in `lexicmap index`).                                                                                                                                                |
|**`-w/--load-whole-seeds`**|                           |Load the whole seed data into memory for faster search     |Use this if the index is not big and many queries are needed to search.                                                                                                                                                                                                 |
|**`-n/--top-n-genomes`**   |Default 0, 0 for all       |Keep top N genome matches for a query in the chaining phase|Value 1 is not recommended as the best chaining result does not always bring the best alignment, so it better be >= 5. The final number of genome hits might be smaller than this number as some chaining results might fail to pass the criteria in the alignment step.|
|**`-a/--all`**             |                           |Output more columns, e.g., matched sequences.              |Use this if you want to output blast-style format with "lexicmap utils 2blast"                                                                                                                                                                                          |
|`-J/--max-query-conc`      |Default 12, 0 for all      |Maximum number of concurrent queries                       |Bigger values do not improve the batch searching speed and consume much memory.                                                                                                                                                                                         |
|`--max-open-files`         |Default: 1024              |Maximum number of open files                               |It mainly affects candidate subsequence extraction. Increase this value if you have hundreds of genome batches or have multiple queries, and do not forgot to set a bigger `ulimit -n` in shell if the value is > 1024.                                                 |

{{< /tab>}}

{{< tab "Chaining" >}}

|Flag                              |Value       |Function                                                                                        |Comment                                                       |
|:---------------------------------|:-----------|:-----------------------------------------------------------------------------------------------|:-------------------------------------------------------------|
|**`-p, --seed-min-prefix`**       |Default 15  |Minimum (prefix) length of matched seeds.                                                       |Smaller values produce more results at the cost of slow speed.|
|**`-P, --seed-min-single-prefix`**|Default 17  |Minimum (prefix) length of matched seeds if there's only one pair of seeds matched.             |Smaller values produce more results at the cost of slow speed.|
|`--seed-max-dist`                 |Default 1000|Max distance between seeds in seed chaining. It should be <= contig interval length in database.|                                                              |
|`--seed-max-gap`                  |Default 200 |Max gap in seed chaining.                                                                       |                                                              |
                                                 |

{{< /tab>}}

{{< tab "Alignment" >}}

|Flag                             |Value       |Function                                                                                                                                                              |Comment|
|:--------------------------------|:-----------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------|
|**`-Q/--min-qcov-per-genome`**   |Default 0   |Minimum query coverage (percentage) per genome.                                                                                                                       |       |
|**`-q/--min-qcov-per-hsp`**      |Default 0   |Minimum query coverage (percentage) per HSP.                                                                                                                          |       |
|**`-l/--align-min-match-len`**   |Default 50  |Minimum aligned length in a HSP segment.                                                                                                                              |       |
|**`-i/--align-min-match-pident`**|Default 70  |Minimum base identity (percentage) in a HSP segment.                                                                                                                  |       |
|`--align-band`                   |Default 100 |Band size in backtracking the score matrix.                                                                                                                           |       |
|`--align-ext-len`                |Default 1000|Extend length of upstream and downstream of seed regions, for extracting query and target sequences for alignment. It should be <= contig interval length in database.|       |
|`--align-max-gap`                |Default 20  |Maximum gap in a HSP segment.                                                                                                                                         |       |


{{< /tab>}}

{{< /tabs >}}

### Improving searching speed

LexicMap's searching speed is related to many factors:
- **The number of similar sequences in the index/database**. More genome hits cost more time, e.g., 16S rRNA gene.
- **Similarity between query and subject sequences**. Alignment of diverse sequences is slower than that of highly similar sequences.
- **The length of query sequence**. Longer queries run with more time.
- **The I/O performance and load**. LexicMap is I/O bound, because seeds matching and extracting candidate subsequences for alignment require a large number of file readings in parallel.
- **CPU frequency and the number of threads**. Faster CPUs and more threads cost less time.


Here are some tips to improve the search speed.

- **Increasing the concurrency number**
    - Make sure that the value of `-j/--threads` (default: all available CPUs) is ≥ than the number of seed chunk file (default: all available CPUs in the indexing step), which can be found in `info.toml` file, e.g,

          # Seeds (k-mer-value data) files
          chunks = 48

    - Increasing the value of `--max-open-files` (default 512). You might also need to [change the open files limit](https://stackoverflow.com/questions/34588/how-do-i-change-the-number-of-open-files-limit-in-linux).
    - (If you have many queries) Increase the value of `-J/--max-query-conc` (default 12), it will increase the memory.
- (If you have many queries) **Loading the entire seed data into memoy** (It's unnecessary if the index is stored in SSD)
    - Setting `-w/--load-whole-seeds` to load the whole seed data into memory for faster search. For example, for ~85,000 GTDB representative genomes, the memory would be ~260 GB with default parameters.
- **Returning less results**
    - Setting `-n/--top-n-genomes` to keep top N genome matches for a query (0 for all) in chaining phase. For queries with a large number of genome hits, a resonable value such as 1000 would reduce the computation time.

## Steps


- For short queries like genes or long reads, returning top N hits.

      lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
          --min-match-pident 70 \
          --min-qcov-per-hsp 70 \
          --min-qcov-per-genome 70 \
          --top-n-genomes 1000

- For longer queries like plasmids, returning all hits.

      lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
          --min-match-pident 70 \
          --min-qcov-per-hsp 0 \
          --min-qcov-per-genome 0  \
          --top-n-genomes 0


{{< expand "Click to show the log of a demo run." "..." >}}

        $ lexicmap search -d demo.lmi/  q.gene.fasta -o q.gene.fasta.lexicmap.tsv
        10:53:20.200 [INFO] LexicMap v0.6.0 (3e249a2)
        10:53:20.200 [INFO]   https://github.com/shenwei356/LexicMap
        10:53:20.200 [INFO] 
        10:53:20.200 [INFO] checking input files ...
        10:53:20.200 [INFO]   1 input file given: q.gene.fasta
        10:53:20.200 [INFO] 
        10:53:20.200 [INFO] loading index: demo.lmi/
        10:53:20.200 [INFO]   reading masks...
        10:53:20.202 [INFO]   reading indexes of seeds (k-mer-value) data...
        10:53:20.929 [INFO]   creating genome reader pools, each batch with 16 readers...
        10:53:20.930 [INFO] index loaded in 729.488175ms
        10:53:20.930 [INFO] 
        10:53:20.930 [INFO] searching with 16 threads...

        10:53:20.968 [INFO] 
        10:53:20.968 [INFO] processed queries: 1, speed: 1546.205 queries per minute
        10:53:20.968 [INFO] 100.0000% (1/1) queries matched
        10:53:20.968 [INFO] done searching
        10:53:20.968 [INFO] search results saved to: q.gene.fasta.lexicmap.tsv
        10:53:20.969 [INFO] 
        10:53:20.969 [INFO] elapsed time: 768.5836ms
        10:53:20.969 [INFO]

{{< /expand >}}

- Extracting similar sequences for a query gene.

    ```text
    # search matches with query coverage >= 90%
    lexicmap search -d gtdb_complete.lmi/ b.gene_E_faecalis_SecY.fasta --min-qcov-per-hsp 90 --all -o results.tsv

    # extract matched sequences as FASTA format
    sed 1d results.tsv | awk -F'\t' '{print ">"$5":"$14"-"$15":"$16"\n"$20;}' | seqkit seq -g > results.fasta

    seqkit head -n 1 results.fasta | head -n 3
    >NZ_JALSCK010000007.1:39224-40522:-
    TTGTTCAAGCTATTAAAGAACGCCTTTAAAGTCAAAGACATTAGATCAAAAATCTTATTT
    ACAGTTTTAATCTTGTTTGTATTTCGCCTAGGTGCGCACATTACTGTGCCCGGGGTGAAT
    ```

- Exporting blast-like alignment text.

    From file:


    ```text
    lexicmap utils 2blast results.tsv -o results.txt

    ```
    Add genome annotation

    ```
    lexicmap utils 2blast results.tsv -o results.txt --kv-file-genome ass2species.map
    ```

    From stdin:

    ```text
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


## Output

### Alignment result relationship

    Query
    ├── Subject genome                             # A query might have one or more genome hits,
        ├── Subject sequence                       # in different sequences.
            ├── High-Scoring segment Pair (HSP)    # HSP is an alignment segment.

Here, the defination of HSP is similar with that in BLAST. Actually there are small gaps in HSPs.

> A High-scoring Segment Pair (HSP) is a local alignment with no gaps that achieves one of the highest alignment scores in a given search.
> https://www.ncbi.nlm.nih.gov/books/NBK62051/


### Output format

Tab-delimited format with 17+ columns, with 1-based positions.

    1.  query,    Query sequence ID.
    2.  qlen,     Query sequence length.
    3.  hits,     Number of subject genomes.
    4.  sgenome,  Subject genome ID.
    5.  sseqid,   Subject sequence ID.
    6.  qcovGnm,  Query coverage (percentage) per genome: $(aligned bases in the genome)/$qlen.
    7.  hsp,      Nth HSP in the genome. (just for improving readability)
    8.  qcovHSP   Query coverage (percentage) per HSP: $(aligned bases in a HSP)/$qlen.
    9.  alenHSP,  Aligned length in the current HSP.
    10. pident,   Percentage of identical matches in the current HSP.
    11. gaps,     Gaps in the current HSP.
    12. qstart,   Start of alignment in query sequence.
    13. qend,     End of alignment in query sequence.
    14. sstart,   Start of alignment in subject sequence.
    15. send,     End of alignment in subject sequence.
    16. sstr,     Subject strand.
    17. slen,     Subject sequence length.
    18. cigar,    CIGAR string of the alignment.                      (optional with -a/--all)
    19. qseq,     Aligned part of query sequence.                     (optional with -a/--all)
    20. sseq,     Aligned part of subject sequence.                   (optional with -a/--all)
    21. align,    Alignment text ("|" and " ") between qseq and sseq. (optional with -a/--all)

Result ordering:

  1. Within each subject genome, alignments (HSP) are sorted by qcovHSP*pident.
  2. Results of multiple subject genomes are sorted by qcovHSP*pident of the best alignment.


### Examples

{{< tabs "t2" >}}

{{< tab "A single-copy gene (SecY)" >}}

```plain
query                                      qlen   hits   sgenome           sseqid                 qcovGnm   hsp   qcovHSP   alenHSP   pident    gaps   qstart   qend   sstart   send     sstr   slen
----------------------------------------   ----   ----   ---------------   --------------------   -------   ---   -------   -------   -------   ----   ------   ----   ------   ------   ----   -------
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   3580   GCF_000395405.1   NZ_KB947497.1          100.000   1     100.000   1299      100.000   0      1        1299   232279   233577   +      274511
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   3580   GCF_019731615.1   NZ_JAASJA010000010.1   100.000   1     100.000   1299      100.000   0      1        1299   2798     4096     +      42998
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   3580   GCA_004103085.1   RPCL01000012.1         100.000   1     100.000   1299      100.000   0      1        1299   44095    45393    +      84242
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   3580   GCF_023571745.1   NZ_JAMKBS010000014.1   100.000   1     100.000   1299      100.000   0      1        1299   44077    45375    +      84206
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   3580   GCF_013248625.1   NZ_JABTDK010000002.1   100.000   1     100.000   1299      100.000   0      1        1299   9609     10907    +      49787
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   3580   GCF_900092155.1   NZ_FLUS01000006.1      100.000   1     100.000   1299      100.000   0      1        1299   63161    64459    +      77366
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   3580   GCF_902165815.1   NZ_CABHHZ010000005.1   100.000   1     100.000   1299      100.000   0      1        1299   39386    40684    -      200163
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   3580   GCF_014243495.1   NZ_SJAV01000002.1      100.000   1     100.000   1299      100.000   0      1        1299   39085    40383    -      256772
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   3580   GCF_900148695.1   NZ_FRXS01000009.1      100.000   1     100.000   1299      100.000   0      1        1299   39230    40528    -      96692
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   3580   GCF_902164645.1   NZ_LR607334.1          100.000   1     100.000   1299      100.000   0      1        1299   236677   237975   +      3380663
```
{{< /tab>}}

{{< tab "A 16S rRNA gene" >}}


```plain
query                         qlen   hits     sgenome           sseqid              qcovGnm   hsp   qcovHSP   alenHSP   pident    gaps   qstart   qend   sstart    send      sstr   slen
---------------------------   ----   ------   ---------------   -----------------   -------   ---   -------   -------   -------   ----   ------   ----   -------   -------   ----   -------
NC_000913.3:4166659-4168200   1542   293398   GCF_002248685.1   NZ_NQBE01000079.1   100.000   1     100.000   1542      100.000   0      1        1542   40        1581      -      99259
NC_000913.3:4166659-4168200   1542   293398   GCF_017164795.1   NZ_CP062702.1       100.000   1     100.000   1542      100.000   0      1        1542   1270211   1271752   +      5483624
NC_000913.3:4166659-4168200   1542   293398   GCF_017164795.1   NZ_CP062702.1       100.000   2     100.000   1542      100.000   0      1        1542   5466287   5467828   -      5483624
NC_000913.3:4166659-4168200   1542   293398   GCF_017164795.1   NZ_CP062702.1       100.000   3     100.000   1543      99.546    2      1        1542   557008    558549    +      5483624
NC_000913.3:4166659-4168200   1542   293398   GCF_017164795.1   NZ_CP062702.1       100.000   4     100.000   1543      99.482    2      1        1542   4473658   4475199   -      5483624
NC_000913.3:4166659-4168200   1542   293398   GCF_017164795.1   NZ_CP062702.1       100.000   5     100.000   1543      99.482    2      1        1542   5154150   5155691   -      5483624
NC_000913.3:4166659-4168200   1542   293398   GCF_017164795.1   NZ_CP062702.1       100.000   6     100.000   1543      99.482    2      1        1542   5195176   5196717   -      5483624
NC_000913.3:4166659-4168200   1542   293398   GCF_017164795.1   NZ_CP062702.1       100.000   7     100.000   1543      99.482    2      1        1542   5369865   5371406   -      5483624
NC_000913.3:4166659-4168200   1542   293398   GCF_000460355.1   NZ_KE701684.1       100.000   1     100.000   1542      100.000   0      1        1542   1108651   1110192   -      1914390
NC_000913.3:4166659-4168200   1542   293398   GCF_000460355.1   NZ_KE701686.1       100.000   2     100.000   1542      99.741    0      1        1542   100680    102221    +      102235
```


{{< /tab>}}

{{< tab "A plasmid" >}}


```plain
query        qlen    hits    sgenome           sseqid          qcovGnm   hsp   qcovHSP   alenHSP   pident    gaps   qstart   qend    sstart    send      sstr   slen
----------   -----   -----   ---------------   -------------   -------   ---   -------   -------   -------   ----   ------   -----   -------   -------   ----   -------
CP115019.1   52830   58744   GCF_022759845.1   NZ_CP086533.1   97.473    1     75.792    40041     99.995    0      12069    52109   11439     51479     +      51479
CP115019.1   52830   58744   GCF_022759845.1   NZ_CP086533.1   97.473    2     20.316    10733     100.000   0      1        10733   722       11454     +      51479
CP115019.1   52830   58744   GCF_022759845.1   NZ_CP086533.1   97.473    3     1.365     721       100.000   0      52110    52830   1         721       +      51479
CP115019.1   52830   58744   GCF_022759845.1   NZ_CP086535.1   97.473    4     0.916     484       91.116    0      51686    52169   27192     27675     -      34058
CP115019.1   52830   58744   GCF_022759845.1   NZ_CP086535.1   97.473    5     0.829     438       90.868    1      52342    52779   26583     27019     -      34058
CP115019.1   52830   58744   GCF_022759845.1   NZ_CP086533.1   97.473    6     1.552     820       100.000   0      9049     9868    23092     23911     +      51479
CP115019.1   52830   58744   GCF_022759845.1   NZ_CP086534.1   97.473    7     0.502     265       100.000   0      19788    20052   29842     30106     +      47185
CP115019.1   52830   58744   GCF_022759845.1   NZ_CP086533.1   97.473    8     0.159     84        97.619    0      8348     8431    19574     19657     +      51479
CP115019.1   52830   58744   GCF_022759905.1   NZ_CP086545.1   97.473    1     75.792    40041     99.995    0      12069    52109   11439     51479     +      51479
CP115019.1   52830   58744   GCF_022759905.1   NZ_CP086545.1   97.473    2     20.316    10733     100.000   0      1        10733   722       11454     +      51479
CP115019.1   52830   58744   GCF_022759905.1   NZ_CP086545.1   97.473    3     1.365     721       100.000   0      52110    52830   1         721       +      51479
CP115019.1   52830   58744   GCF_022759905.1   NZ_CP086547.1   97.473    4     0.916     484       91.116    0      51686    52169   3843      4326      +      34058
CP115019.1   52830   58744   GCF_022759905.1   NZ_CP086547.1   97.473    5     0.829     438       90.868    1      52342    52779   4499      4935      +      34058
CP115019.1   52830   58744   GCF_022759905.1   NZ_CP086545.1   97.473    6     1.552     820       100.000   0      9049     9868    23092     23911     +      51479
CP115019.1   52830   58744   GCF_022759905.1   NZ_CP086546.1   97.473    7     0.502     265       100.000   0      19788    20052   29842     30106     +      47185
CP115019.1   52830   58744   GCF_022759905.1   NZ_CP086545.1   97.473    8     0.159     84        97.619    0      8348     8431    19574     19657     +      51479
CP115019.1   52830   58744   GCF_014826015.1   NZ_CP058621.1   97.473    1     77.157    40762     99.993    0      12069    52830   9513      50274     +      51480
CP115019.1   52830   58744   GCF_014826015.1   NZ_CP058621.1   97.473    2     18.033    9528      99.990    1      1207     10733   1         9528      +      51480
CP115019.1   52830   58744   GCF_014826015.1   NZ_CP058621.1   97.473    3     2.283     1206      100.000   0      1        1206    50275     51480     +      51480
CP115019.1   52830   58744   GCF_014826015.1   NZ_CP058618.1   97.473    4     2.497     1319      100.000   0      25153    26471   3019498   3020816   -      4718403
```

{{< /tab>}}

{{< tab "Long reads" >}}


Queries are a few Nanopore Q20 reads from a mock metagenomic community.

```plain
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

{{< /tab>}}

{{< /tabs >}}

Search results (TSV format) above are formatted with [csvtk pretty](https://github.com/shenwei356/csvtk).

### Summarizing results

If you would like to summarize alignment results, e.g., the number of species, here's the method.

1. Prepare a two-column tab-delimited file for mapping reference (genome) or sequence IDs to any information (such as species name).
   
        # for GTDB/GenBank/RefSeq genomes downloaded with genome_updater
        cut -f 1,8 assembly_summary.txt > ass2species.tsv

        head -n 3 ass2species.tsv
        GCF_002287175.1 Methanobacterium bryantii
        GCF_000762265.1 Methanobacterium formicicum
        GCF_029601605.1 Methanobacterium formicicum

2. Add information to the alignment result with [csvtk](https://github.com/shenwei356/csvtk) or other tools.

        # add species
        cat b.gene_E_coli_16S.fasta.lexicmap.tsv \
            | csvtk mutate -t --after slen -n species -f sgenome \
            | csvtk replace -t -f species -p "(.+)" -r "{kv}" -k ass2species.tsv \
            > result.with_species.tsv

        # filter result with query coverage >= 80 and count the species
        cat result.with_species.tsv \
            | csvtk uniq -t -f sgenome \
            | csvtk filter2 -t -f "\$qcovHSP >= 80" \
            | csvtk freq -t -f species -nr \
            > result.with_species.tsv.stats.tsv

        csvtk head -t -n 5 result.with_species.tsv.stats.tsv \
            | csvtk pretty -t

        species                    frequency
        ------------------------   ---------
        Salmonella enterica        135065   
        Escherichia coli           128071   
        Streptococcus pneumoniae   51971    
        Staphylococcus aureus      44215    
        Pseudomonas aeruginosa     34254
