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
LexicMap is mainly designed for sequence alignment with a small number of queries (gene/plasmid/virus/phage sequences) longer than 100 bp by default.
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
        - The number and length of query sequences.
        - The number of matched genomes and sequences.
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
- **The I/O performance and load**. LexicMap is I/O bound, because seeds matching (serial access) and extracting candidate subsequences for alignment (random access) require a large number of file readings in parallel.
- **Similarity between query and subject sequences**. Alignment of diverse sequences is slightly slower than that of highly similar sequences.
- **The length of query sequence**. Longer queries run with more time.
- **CPU frequency and the number of threads**. Faster CPUs and more threads cost less time.


Here are some tips to improve the search speed.

- **Returning less results**
    - Bigger `-p/--seed-min-prefix` (default 15) and `-P/--seed-min-single-prefix` (default 17)
      increase the search speed at the cost of decreased sensitivity for distant matches (similarity < 90%).
      Don't worry if you only search highly similar matches.
    - Setting `-n/--top-n-genomes` to keep top N genome matches for a query (0 for all) in chaining phase. 
      For queries with a large number of genome hits, a resonable value such as 1000 would reduce the computation time.
    - **Note that**: alignment result filtering is performed in the final phase, so stricter filtering criteria,
     including `-q/--min-qcov-per-hsp`, `-Q/--min-qcov-per-genome`, and `-i/--align-min-match-pident`,
     do not significantly accelerate the search speed. Hence, you can search with default
     parameters and then filter the result with tools like `awk` or `csvtk`.
- **Increasing the concurrency number**
    - Make sure that the value of `-j/--threads` (default: all available CPUs) is ≥ than the number of seed chunk file (default: all available CPUs in the indexing step), which can be found in `info.toml` file, e.g,

          # Seeds (k-mer-value data) files
          chunks = 48

    - Increasing the value of `--max-open-files` (default 512). You might also need to [change the open files limit](https://stackoverflow.com/questions/34588/how-do-i-change-the-number-of-open-files-limit-in-linux).
    - (If you have many queries) Increase the value of `-J/--max-query-conc` (default 12), it will increase the memory.
- **Loading the entire seed data into memoy** (If you have many queries and the index is not very big. It's unnecessary if the index is stored in SSD)
    - Setting `-w/--load-whole-seeds` to load the whole seed data into memory for faster search. For example, for ~85,000 GTDB representative genomes, the memory would be ~260 GB with default parameters.

## Steps


- For short queries like genes or long reads, returning top N hits.

      lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
          --min-match-pident 70 \
          --min-qcov-per-hsp 70 \
          --min-qcov-per-genome 70 \
          --top-n-genomes 10000

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
    sed 1d results.tsv | awk -F'\t' '{print ">"$5":"$14"-"$15":"$16"\n"$22;}' | seqkit seq -g > results.fasta

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
     Score = 280 bits, Expect = 9.66e-75
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

Tab-delimited format with 19+ columns, with 1-based positions.

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
    18. evalue,   E value.
    19. bitscore, bit score.
    20. cigar,    CIGAR string of the alignment.                      (optional with -a/--all)
    21. qseq,     Aligned part of query sequence.                     (optional with -a/--all)
    22. sseq,     Aligned part of subject sequence.                   (optional with -a/--all)
    23. align,    Alignment text ("|" and " ") between qseq and sseq. (optional with -a/--all)

Result ordering:

  1. Within each subject genome, alignments (HSP) are sorted by qcovHSP*pident.
  2. Results of multiple subject genomes are sorted by qcovHSP*pident of the best alignment.


### Examples

{{< tabs "t2" >}}

{{< tab "A single-copy gene (SecY)" >}}

```plain
query                                      qlen   hits    sgenome           sseqid                 qcovGnm   hsp   qcovHSP   alenHSP   pident    gaps   qstart   qend   sstart   send     sstr   slen     evalue     bitscore
----------------------------------------   ----   -----   ---------------   --------------------   -------   ---   -------   -------   -------   ----   ------   ----   ------   ------   ----   ------   --------   --------
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   41718   GCA_003962555.1   QNGV01000003.1         100.000   1     100.000   1299      100.000   0      1        1299   39229    40527    -      271174   0.00e+00   2343    
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   41718   GCA_900148695.1   FRXS01000009.1         100.000   1     100.000   1299      100.000   0      1        1299   39230    40528    -      96692    0.00e+00   2343    
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   41718   GCA_015919505.1   AAXTGP010000012.1      100.000   1     100.000   1299      100.000   0      1        1299   39064    40362    -      103670   0.00e+00   2343    
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   41718   GCA_028641435.1   DANIKN010000009.1      100.000   1     100.000   1299      100.000   0      1        1299   39290    40588    -      103925   0.00e+00   2343    
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   41718   GCF_902165815.1   NZ_CABHHZ010000005.1   100.000   1     100.000   1299      100.000   0      1        1299   39386    40684    -      200163   0.00e+00   2343    
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   41718   GCA_032868635.1   JAWKAC010000007.1      100.000   1     100.000   1299      100.000   0      1        1299   39280    40578    -      195040   0.00e+00   2343    
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   41718   GCA_031737035.1   JAVGUQ010000004.1      100.000   1     100.000   1299      100.000   0      1        1299   39191    40489    -      292887   0.00e+00   2343    
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   41718   GCA_014353145.1   JACIKK010000004.1      100.000   1     100.000   1299      100.000   0      1        1299   192734   194032   +      233180   0.00e+00   2343    
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   41718   GCA_005237765.1   QPWL01000004.1         100.000   1     100.000   1299      100.000   0      1        1299   10449    11747    -      228528   0.00e+00   2343    
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   41718   GCA_030342345.1   JATAAU010000007.1      100.000   1     100.000   1299      100.000   0      1        1299   134344   135642   +      154082   0.00e+00   2343  
```
{{< /tab>}}

{{< tab "A 16S rRNA gene" >}}


```plain
query                         qlen   hits      sgenome           sseqid              qcovGnm   hsp   qcovHSP   alenHSP   pident    gaps   qstart   qend   sstart    send      sstr   slen      evalue     bitscore
---------------------------   ----   -------   ---------------   -----------------   -------   ---   -------   -------   -------   ----   ------   ----   -------   -------   ----   -------   --------   --------
NC_000913.3:4166659-4168200   1542   1955162   GCA_907175715.1   CAJSIH010000102.1   100.000   1     100.000   1542      100.000   0      1        1542   147       1688      +      1702      0.00e+00   2782    
NC_000913.3:4166659-4168200   1542   1955162   GCF_900174625.1   NZ_LT838196.1       100.000   1     100.000   1542      100.000   0      1        1542   1162946   1164487   +      4809920   0.00e+00   2782    
NC_000913.3:4166659-4168200   1542   1955162   GCF_900174625.1   NZ_LT838196.1       100.000   2     100.000   1542      100.000   0      1        1542   3790963   3792504   -      4809920   0.00e+00   2782    
NC_000913.3:4166659-4168200   1542   1955162   GCF_900174625.1   NZ_LT838196.1       100.000   3     100.000   1542      100.000   0      1        1542   4484052   4485593   -      4809920   0.00e+00   2782    
NC_000913.3:4166659-4168200   1542   1955162   GCF_900174625.1   NZ_LT838196.1       100.000   4     100.000   1542      100.000   0      1        1542   4656667   4658208   -      4809920   0.00e+00   2782    
NC_000913.3:4166659-4168200   1542   1955162   GCF_900174625.1   NZ_LT838196.1       100.000   5     100.000   1542      99.741    0      1        1542   460127    461668    +      4809920   0.00e+00   2764    
NC_000913.3:4166659-4168200   1542   1955162   GCF_900174625.1   NZ_LT838196.1       100.000   6     100.000   1542      99.676    0      1        1542   4525540   4527081   -      4809920   0.00e+00   2758    
NC_000913.3:4166659-4168200   1542   1955162   GCF_900174625.1   NZ_LT838196.1       100.000   7     100.000   1542      99.676    0      1        1542   4750295   4751836   -      4809920   0.00e+00   2758    
NC_000913.3:4166659-4168200   1542   1955162   GCF_013421655.1   NZ_WVHY01000001.1   100.000   1     100.000   1542      100.000   0      1        1542   279590    281131    -      2103828   0.00e+00   2782    
NC_000913.3:4166659-4168200   1542   1955162   GCF_013421655.1   NZ_WVHY01000001.1   100.000   2     100.000   1542      99.741    0      1        1542   1116560   1118101   -      2103828   0.00e+00   2764
```


{{< /tab>}}

{{< tab "A plasmid" >}}


```plain
query        qlen    hits     sgenome           sseqid          qcovGnm   hsp   qcovHSP   alenHSP   pident    gaps   qstart   qend    sstart   send    sstr   slen    evalue     bitscore
----------   -----   ------   ---------------   -------------   -------   ---   -------   -------   -------   ----   ------   -----   ------   -----   ----   -----   --------   --------
CP115019.1   52830   561754   GCF_030863845.1   NZ_CP115019.1   100.000   1     100.000   52830     100.000   0      1        52830   1        52830   +      52830   0.00e+00   95273   
CP115019.1   52830   561754   GCF_030863845.1   NZ_CP115019.1   100.000   2     1.562     827       99.758    2      9044     9868    23715    24541   +      52830   0.00e+00   1496    
CP115019.1   52830   561754   GCF_030863845.1   NZ_CP115019.1   100.000   3     1.565     827       99.758    2      23715    24541   9044     9868    +      52830   0.00e+00   1496    
CP115019.1   52830   561754   GCF_030863845.1   NZ_CP115019.1   100.000   4     0.208     110       98.182    0      29827    29936   35931    36040   -      52830   4.13e-43   190     
CP115019.1   52830   561754   GCF_030863845.1   NZ_CP115018.1   100.000   5     5.397     2851      100.000   0      14179    17029   30028    32878   -      84227   0.00e+00   5142    
CP115019.1   52830   561754   GCF_030863845.1   NZ_CP115018.1   100.000   6     2.996     1583      100.000   0      12069    13651   28264    29846   -      84227   0.00e+00   2856    
CP115019.1   52830   561754   GCF_030863845.1   NZ_CP115018.1   100.000   7     1.556     822       100.000   0      23722    24543   24031    24852   +      84227   0.00e+00   1483    
CP115019.1   52830   561754   GCF_030863845.1   NZ_CP115018.1   100.000   8     1.552     820       100.000   0      9049     9868    24031    24850   +      84227   0.00e+00   1480    
CP115019.1   52830   561754   GCF_030863845.1   NZ_CP115018.1   100.000   9     1.552     820       99.878    0      23722    24541   27427    28246   -      84227   0.00e+00   1474    
CP115019.1   52830   561754   GCF_030863845.1   NZ_CP115018.1   100.000   10    1.554     821       99.756    0      9048     9868    38304    39124   -      84227   0.00e+00   1472
```

{{< /tab>}}

{{< tab "Long reads" >}}


Queries are a few Nanopore Q20 reads from a mock metagenomic community.

```plain
query                qlen   hits   sgenome                   sseqid              qcovGnm   hsp   qcovHSP   alenHSP   pident   gaps   qstart   qend   sstart    send      sstr   slen      evalue      bitscore
------------------   ----   ----   -----------------------   -----------------   -------   ---   -------   -------   ------   ----   ------   ----   -------   -------   ----   -------   ---------   --------
ERR5396170.1000004   190    1      GCF_000227465.1_genomic   NC_016047.1         84.211    1     84.211    165       89.091   5      14       173    4189372   4189536   -      4207222   1.93e-63    253     
ERR5396170.1000006   796    3      GCF_013394085.1_genomic   NZ_CP040910.1       99.623    1     99.623    801       97.628   9      4        796    1138907   1139706   +      1887974   0.00e+00    1431    
ERR5396170.1000006   796    3      GCF_013394085.1_genomic   NZ_CP040910.1       99.623    2     99.623    801       97.628   9      4        796    32607     33406     +      1887974   0.00e+00    1431    
ERR5396170.1000006   796    3      GCF_013394085.1_genomic   NZ_CP040910.1       99.623    3     99.623    801       97.628   9      4        796    134468    135267    -      1887974   0.00e+00    1431    
ERR5396170.1000006   796    3      GCF_013394085.1_genomic   NZ_CP040910.1       99.623    4     99.623    801       97.503   9      4        796    1768896   1769695   +      1887974   0.00e+00    1427    
ERR5396170.1000006   796    3      GCF_013394085.1_genomic   NZ_CP040910.1       99.623    5     99.623    801       97.378   9      4        796    242012    242811    -      1887974   0.00e+00    1422    
ERR5396170.1000006   796    3      GCF_013394085.1_genomic   NZ_CP040910.1       99.623    6     99.623    801       96.879   12     4        796    154380    155176    -      1887974   0.00e+00    1431    
ERR5396170.1000006   796    3      GCF_013394085.1_genomic   NZ_CP040910.1       99.623    7     57.915    469       95.736   9      4        464    1280313   1280780   +      1887974   3.71e-236   829     
ERR5396170.1000006   796    3      GCF_013394085.1_genomic   NZ_CP040910.1       99.623    8     42.839    341       99.120   0      456      796    1282477   1282817   +      1887974   6.91e-168   601     
ERR5396170.1000006   796    3      GCF_009663775.1_genomic   NZ_RDBR01000008.1   99.623    1     99.623    801       93.383   9      4        796    21391     22190     -      52610     0.00e+00    1278    
ERR5396170.1000006   796    3      GCF_003344625.1_genomic   NZ_QPKJ02000188.1   97.362    1     87.437    700       98.143   5      22       717    1         699       -      826       0.00e+00    1249    
ERR5396170.1000006   796    3      GCF_003344625.1_genomic   NZ_QPKJ02000423.1   97.362    2     27.889    222       99.550   0      575      796    1         222       +      510       3.47e-106   396     
ERR5396170.1000000   698    2      GCF_001457615.1_genomic   NZ_LN831024.1       92.264    1     92.264    656       96.341   13     53       696    4452083   4452737   +      6316979   0.00e+00    1169    
ERR5396170.1000000   698    2      GCF_000949385.2_genomic   NZ_JYKO02000001.1   91.977    1     91.977    654       78.135   13     55       696    5638788   5639440   -      5912440   2.68e-176   630     
ERR5396170.1000012   848    2      GCF_013394085.1_genomic   NZ_CP040910.1       98.585    1     98.585    841       96.671   10     13       848    190308    191143    -      1887974   0.00e+00    1472    
ERR5396170.1000012   848    2      GCF_001293735.1_genomic   NZ_BCAH01000003.1   90.212    1     90.212    782       77.110   23     51       815    8230      9005      -      40321     3.19e-214   756
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
