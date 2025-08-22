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
              --align-min-match-pident 80 --min-qcov-per-hsp 70 --min-qcov-per-genome 70 \
              --top-n-genomes 10000

    - For longer queries like plasmids, returning all hits.

          lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
              --align-min-match-pident 70 --min-qcov-per-hsp 0  --min-qcov-per-genome 50 \
              --align-min-match-len 1000 \
              --top-n-genomes 0

## Input

{{< hint type=note >}}
**Query length**\
LexicMap is mainly designed for sequence alignment with a small number of queries (gene/plasmid/virus/phage sequences) longer than 150 bp by default.
{{< /hint >}}

Input should be (gzipped) FASTA or FASTQ records from files or STDIN.


## Hardware requirements

See [benchmark of index building](https://bioinf.shenwei.me/LexicMap/introduction/#searching).

LexicMap is designed to provide fast and low-memory sequence alignment against millions of prokaryotic genomes.

- **CPU:**
    - No specific requirements on CPU type and instruction sets. Both x86 and ARM chips are supported.
    - More is better as LexicMap is a CPU-intensive software. It uses all CPUs by default (`-j/--threads`).
- **RAM**
    - More RAM (>= 16 GB) is preferred. The memory usage in searching is mainly related to:
        - The number and length of query sequences.
        - The number of matched genomes and sequences.
        - Similarities between query and target sequences.
        - The number of threads. It uses all CPUs by default (`-j/--threads`).
        - (Batch searching) The number of concurrent queries (`-J/--max-query-conc`, default 8).
        - (Batch searching) Garbage collection interval (`--gc-interval`, default 64, 0 for disable).
- **Disk**
    - SSD disks are preferred to store the index size, while HDD disks are also fast enough.
    - No temporary files are generated during searching.

## Algorithm

<img src="/LexicMap/searching.svg" alt="" width="900"/>

See the [paper](https://bioinf.shenwei.me/LexicMap/introduction/#citation) for details.

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

|Flag                    |Value                      |Function                                                       |Comment                                                                                                                                                                                                                                                                 |
|:-----------------------|:--------------------------|:--------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|**`-d/--index`**        |                           |Index directory created by "lexicmap index".                   |                                                                                                                                                                                                                                                                        |
|**`-o/--out-file`**     |Default: - (stdout)        |Out file, supports a ".gz" suffix ("-" for stdout).            |                                                                                                                                                                                                                                                                        |
|**`-j/--threads`**      |Default: all available cpus|Number of CPU cores to use.                                    |The value should be >= the number of seed chunk files (“chunks” in info.toml, set by `-c/--chunks` in `lexicmap index`).                                                                                                                                                |
|**`-a/--all`**          |                           |Output more columns, e.g., matched sequences.                  |Use this if you want to output blast-style format with "lexicmap utils 2blast"                                                                                                                                                                                          |
|**`-n/--top-n-genomes`**|Default 0, 0 for all       |Keep top N genome matches for a query in the chaining phase    |Value 1 is not recommended as the best chaining result does not always bring the best alignment, so it better be >= 5. The final number of genome hits might be smaller than this number as some chaining results might fail to pass the criteria in the alignment step.|
|`-J/--max-query-conc`   |Default 8, 0 for all       |Maximum number of concurrent queries                           |Bigger values do not improve the batch searching speed and consume much memory.                                                                                                                                                                                         |
|`--gc-interval`         |Default 64, 0 for disable  |Force garbage collection every N queries.                      |The value can't be too small.                                                                                                                                                                                                                                           |
|`--max-open-files`      |Default: 1024              |Maximum number of open files                                   |It mainly affects candidate subsequence extraction. Increase this value if you have hundreds of genome batches or have multiple queries, and do not forgot to set a bigger `ulimit -n` in shell if the value is > 1024.                                                 |
|`-w/--load-whole-seeds` |                           |Load the whole seed data into memory for faster batch searching|Use this if the index is not big and many queries are needed to search.                                                                                                                                                                                                 |
|`--debug`               |                           |Print debug information, including a progress bar.             |Recommended when searching with one query.                                                                                                                                                                                                                              |

{{< /tab>}}

{{< tab "Chaining" >}}

|Flag                              |Value       |Function                                                                                        |Comment                                                       |
|:---------------------------------|:-----------|:-----------------------------------------------------------------------------------------------|:-------------------------------------------------------------|
|**`-p, --seed-min-prefix`**       |Default 15  |Minimum (prefix) length of matched seeds.                                                       |Smaller values produce more results at the cost of slow speed.|
|**`-P, --seed-min-single-prefix`**|Default 17  |Minimum (prefix) length of matched seeds if there's only one pair of seeds matched.             |Smaller values produce more results at the cost of slow speed.|
|`--seed-max-dist`                 |Default 1000|Max distance between seeds in seed chaining. It should be <= contig interval length in database.|                                                              |
|`--seed-max-gap`                  |Default 50  |Max gap in seed chaining.                                                                       |                                                              |
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
- **The I/O performance and load**. LexicMap is I/O bound, because seeds matching (sequential reading) and extracting candidate subsequences for alignment (**random access**) require a large number of file readings in parallel.
- **Similarity between query and subject sequences**. Alignment of diverse sequences is slightly slower than that of highly similar sequences.
- **The length of query sequence**. Longer queries run with more time.
- **CPU frequency and the number of threads**. Faster CPUs and more threads cost less time.


Here are some tips to improve the search speed.

- **Storing the index on SSD** (It would be very fast!)
- **Returning less results**
    - Bigger `-p/--seed-min-prefix` (default 15) and `-P/--seed-min-single-prefix` (default 17),
      e.g., `-p 17 -P 19`,
      increase the search speed at the cost of decreased sensitivity for distant matches (similarity < 90%) or short queries.
      Don't worry if you only search highly similar matches or long queries.
    - Setting `-n/--top-n-genomes` to keep top N genome matches for a query (0 for all) in chaining phase. 
      For queries with a large number of genome hits, a resonable value such as 1000 would significantly reduce the computation time.
    - **Note that**: alignment result filtering is performed in the final phase, so stricter filtering criteria,
     including `-q/--min-qcov-per-hsp`, `-Q/--min-qcov-per-genome`, and `-i/--align-min-match-pident`,
     do not significantly accelerate the search speed. Hence, you can search with default
     parameters and then filter the result with tools such as [csvtk](https://github.com/shenwei356/csvtk).
- **Increasing the concurrency number**
    - Make sure that the value of `-j/--threads` (default: all available CPUs) is ≥ than the number of seed chunk file (default: all available CPUs in the indexing step), which can be found in `info.toml` file, e.g,
        ```
        # Seeds (k-mer-value data) files
        chunks = 48
        ```
    - Increasing the value of `--max-open-files` (default 1024). You might also need to [change the open files limit](https://stackoverflow.com/questions/34588/how-do-i-change-the-number-of-open-files-limit-in-linux).
    - (If you have many queries) Increase the value of `-J/--max-query-conc` (default 8), which might help. This will increase the memory.
- **Loading the entire seed data into memoy** (*If you have many queries and the index is not very big*. It's unnecessary if the index is stored on SSD)
    - Setting `-w/--load-whole-seeds` to load the whole seed data into memory for faster seed matching. For example, for ~85,000 GTDB representative genomes, the memory would be ~260 GB with default parameters.


### Searching with plasmids or other longer queries

For long queries, such as plasmids, a few parameters can be adjusted for better performance.

- Bigger `-p/--seed-min-prefix` (default 15) and `-P/--seed-min-single-prefix` (default 17),
      e.g., `-p 19 -P 21`. The search sensitivity will not be affected for long queries or high similarity subjects.
- Bigger `-l/--align-min-match-len` (default 50), such as `1000`, because small matches are less informative.

When searching with plasmids, it's recommended to use a strict criterion of `-Q/--min-qcov-per-genome` (`qcovGnm`, default 0), such as 80,
and further filter results with a loose criterion of `-q/--min-qcov-per-hsp` (`qcovHSP`, default 0) after searching, such as 50/60/70.
The reasons are:

- Plasmids are circular, while they are stored linearly. 
  The different starting positions in query and subject sequences would result in two alignment segments (small `qcovHSP`).
- Assemblies can be fragmented, with many contigs, especially these assembled from short reads.
  Therefore, a plasmid might be aligned to multiple contigs with small `qcovHSP`.
    
## Steps


- For short queries like genes or long reads, returning top N hits.

      lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
          --min-match-pident 80 \
          --min-qcov-per-hsp 70 \
          --min-qcov-per-genome 70 \
          --top-n-genomes 10000

- For longer queries like plasmids, returning all hits.

      lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
          --min-match-pident 70 \
          --min-qcov-per-hsp 0 \
          --min-qcov-per-genome 50  \
          --align-min-match-len 1000 \
          --top-n-genomes 0


{{< expand "Click to show the log of a demo run." "..." >}}

        $ lexicmap search -d demo.lmi/  q.gene.fasta -o q.gene.fasta.lexicmap.tsv --debug
        09:56:48.464 [INFO] LexicMap v0.7.0
        09:56:48.464 [INFO]   https://github.com/shenwei356/LexicMap
        09:56:48.464 [INFO] 
        09:56:48.464 [INFO] checking input files ...
        09:56:48.464 [INFO]   1 input file given: q.gene.fasta
        09:56:48.464 [INFO] 
        09:56:48.464 [INFO] loading index: demo.lmi/
        09:56:48.464 [INFO]   reading masks...
        09:56:48.467 [INFO]   reading indexes of seeds (k-mer-value) data...
        09:56:49.434 [INFO]   creating reader pools for 1 genome batches, each with 16 readers...
        09:56:49.434 [INFO] index loaded in 969.583422ms
        09:56:49.434 [INFO] 
        09:56:49.434 [INFO] searching with 16 threads...
        09:56:49.435 [DEBU] NC_000913.3:4166659-4168200 (1542 bp): start to search
        09:56:49.440 [DEBU] NC_000913.3:4166659-4168200 (1542 bp): finished seed-matching (15 genome hits) in 5.354981ms
        09:56:49.441 [DEBU] NC_000913.3:4166659-4168200 (1542 bp): finished chaining (15 genome hits) in 1.045575ms
        checked genomes:  15 / 15 [======================================] ETA: 0s. done
        09:56:49.473 [DEBU] NC_000913.3:4166659-4168200 (1542 bp): finished alignment (15 genome hits) in 32.005224ms
        09:56:49.473 [DEBU] NC_000913.3:4166659-4168200 (1542 bp): finished sorting alignment results (15 genome hits) in 13.758µs

        09:56:49.473 [INFO] 
        09:56:49.473 [INFO] processed queries: 1, speed: 1512.560 queries per minute
        09:56:49.473 [INFO] 100.0000% (1/1) queries matched
        09:56:49.473 [INFO] done searching
        09:56:49.473 [INFO] search results saved to: q.gene.fasta.lexicmap.tsv
        09:56:49.474 [INFO] 
        09:56:49.474 [INFO] elapsed time: 1.009536088s
        09:56:49.474 [INFO]

{{< /expand >}}

- Extracting similar sequences for a query gene.

    ```plain
    # Extracting similar sequences for a query gene.

    # search matches with query coverage >= 90%
    lexicmap search -d demo.lmi/ bench/b.gene_E_faecalis_SecY.fasta -o results.tsv \
        --min-qcov-per-hsp 90

    # extract matched sequences as FASTA format
    lexicmap utils subseq -d demo.lmi -f results.tsv -o results.tsv.aligned.fasta

    seqkit head -n 1 results.tsv.aligned.fasta | head -n 3
    >NZ_KB944588.1:228637-229935:+ query=lcl|NZ_CP064374.1_cds_WP_002359350.1_906 sgenome=GCF_000392875.1 sseqid=NZ_KB944588.1 qcovGnm=100.000 cls=1 hsp=1 qcovHSP=100.000 alenHSP=1299 pident=100.000 gaps=0 qstart=1 qend=1299 sstart=228637 send=229935 sstr=+ slen=274762 evalue=0.00e+00 bitscore=2343
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
    # here, we only align <=200 bp queries and show one medium-similarity result.
    
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

- Merge a query's search results in multiple indexes with `lexicmap utils merge-search-results`
  ([usage](https://bioinf.shenwei.me/LexicMap/usage/utils/merge-search-results/)).
  
    ```text
    lexicmap util merge-search-results -o query.lexicmap.tsv.gz \
        query.lexicmap@genbank.tsv.gz \
        query.lexicmap@refseq.tsv.gz \
        query.lexicmap@allthebacteria.tsv.gz
    ```

## Output

Two output formats are supported:

- The default format is tab-delimited (see details below),
and a query's search results in multiple indexes can be merged with `lexicmap utils merge-search-results`
  ([usage](https://bioinf.shenwei.me/LexicMap/usage/utils/merge-search-results/)).
- Blast-style pair-wise alignment format can be achieved by searching with `-a/--all` and converting with [lexicmap utils 2blast](https://bioinf.shenwei.me/LexicMap/usage/utils/2blast/).

### Alignment result relationship

    Query
    ├── Subject genome
        ├── Subject sequence
            ├── HSP cluster (a cluster of neighboring HSPs)
                ├── High-Scoring segment Pair (HSP)

Here, the defination of HSP is similar with that in BLAST. Actually there are small gaps in HSPs.

> A High-scoring Segment Pair (HSP) is a local alignment with no gaps that achieves one of the highest alignment scores in a given search.
> https://www.ncbi.nlm.nih.gov/books/NBK62051/


### Output format

Tab-delimited format with 20+ columns, with 1-based positions.

    1.  query,    Query sequence ID.
    2.  qlen,     Query sequence length.
    3.  hits,     Number of subject genomes.
    4.  sgenome,  Subject genome ID.
    5.  sseqid,   Subject sequence ID.
    6.  qcovGnm,  Query coverage (percentage) per genome: $(aligned bases in the genome)/$qlen.
    7.  cls,      Nth HSP cluster in the genome. (just for improving readability)
                  It's useful to show if multiple adjacent HSPs are collinear.
    8.  hsp,      Nth HSP in the genome.         (just for improving readability)
    9.  qcovHSP   Query coverage (percentage) per HSP: $(aligned bases in a HSP)/$qlen.
    10. alenHSP,  Aligned length in the current HSP.
    11. pident,   Percentage of identical matches in the current HSP.
    12. gaps,     Gaps in the current HSP.
    13. qstart,   Start of alignment in query sequence.
    14. qend,     End of alignment in query sequence.
    15. sstart,   Start of alignment in subject sequence.
    16. send,     End of alignment in subject sequence.
    17. sstr,     Subject strand.
    18. slen,     Subject sequence length.
    19. evalue,   Expect value.
    20. bitscore, Bit score.
    21. cigar,    CIGAR string of the alignment.                      (optional with -a/--all)
    22. qseq,     Aligned part of query sequence.                     (optional with -a/--all)
    23. sseq,     Aligned part of subject sequence.                   (optional with -a/--all)
    24. align,    Alignment text ("|" and " ") between qseq and sseq. (optional with -a/--all)

**Result ordering:**

  For a HSP cluster, `SimilarityScore = max(bit_score * pident)`.
  1. Within each HSP cluster, HSPs are sorted by `sstart`.
  2. Within each subject genome, HSP clusters are sorted in descending order by `SimilarityScore`.
  3. Results of multiple subject genomes are sorted by the highest `SimilarityScore` of HSP clusters.



### Examples

{{< tabs "t2" >}}

{{< tab "A single-copy gene (SecY)" >}}

```plain
query                                      qlen   hits   sgenome           sseqid                 qcovGnm   cls   hsp   qcovHSP   alenHSP   pident    gaps   qstart   qend   sstart    send      sstr   slen      evalue     bitscore
----------------------------------------   ----   ----   ---------------   --------------------   -------   ---   ---   -------   -------   -------   ----   ------   ----   -------   -------   ----   -------   --------   --------
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   6255   GCF_017641985.1   NZ_SIYB01000008.1      100.000   1     1     100.000   1299      100.000   0      1        1299   38843     40141     -      140750    0.00e+00   2343    
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   6255   GCF_020405635.1   NZ_JAIZEZ010000006.1   100.000   1     1     100.000   1299      100.000   0      1        1299   217200    218498    +      257915    0.00e+00   2343    
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   6255   GCF_902163655.1   NZ_CABHAB010000005.1   100.000   1     1     100.000   1299      100.000   0      1        1299   39259     40557     -      173125    0.00e+00   2343    
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   6255   GCF_020881975.1   NZ_CP086411.1          100.000   1     1     100.000   1299      100.000   0      1        1299   212962    214260    +      2995874   0.00e+00   2343    
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   6255   GCF_900148695.1   NZ_FRXS01000009.1      100.000   1     1     100.000   1299      100.000   0      1        1299   39230     40528     -      96692     0.00e+00   2343    
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   6255   GCF_009735055.1   NZ_WMGL01000019.1      100.000   1     1     100.000   1299      100.000   0      1        1299   38954     40252     -      58628     0.00e+00   2343    
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   6255   GCA_021838685.1   JAKCDA010000012.1      100.000   1     1     100.000   1299      100.000   0      1        1299   44127     45425     +      84305     0.00e+00   2343    
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   6255   GCF_000742975.1   NZ_CP008816.1          100.000   1     1     100.000   1299      100.000   0      1        1299   2311184   2312482   +      2939973   0.00e+00   2343    
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   6255   GCF_925298485.1   NZ_CAKMCI010000007.1   100.000   1     1     100.000   1299      100.000   0      1        1299   56108     57406     +      96359     0.00e+00   2343    
lcl|NZ_CP064374.1_cds_WP_002359350.1_906   1299   6255   GCF_021120765.1   NZ_JADPWI010000009.1   100.000   1     1     100.000   1299      100.000   0      1        1299   38881     40179     -      84306     0.00e+00   2343
```
{{< /tab>}}

{{< tab "A 16S rRNA gene" >}}


```plain
query                         qlen   hits     sgenome           sseqid              qcovGnm   cls   hsp   qcovHSP   alenHSP   pident    gaps   qstart   qend   sstart    send      sstr   slen      evalue     bitscore
---------------------------   ----   ------   ---------------   -----------------   -------   ---   ---   -------   -------   -------   ----   ------   ----   -------   -------   ----   -------   --------   --------
NC_000913.3:4166659-4168200   1542   306060   GCF_000468515.1   NC_022364.1         100.000   1     1     100.000   1542      100.000   0      1        1542   224243    225784    +      4835601   0.00e+00   2782    
NC_000913.3:4166659-4168200   1542   306060   GCF_000468515.1   NC_022364.1         100.000   2     2     100.000   1542      100.000   0      1        1542   2804586   2806127   -      4835601   0.00e+00   2782    
NC_000913.3:4166659-4168200   1542   306060   GCF_000468515.1   NC_022364.1         100.000   3     3     100.000   1542      100.000   0      1        1542   4350745   4352286   +      4835601   0.00e+00   2782    
NC_000913.3:4166659-4168200   1542   306060   GCF_000468515.1   NC_022364.1         100.000   4     4     100.000   1542      99.676    0      1        1542   4391290   4392831   +      4835601   0.00e+00   2758    
NC_000913.3:4166659-4168200   1542   306060   GCF_000468515.1   NC_022364.1         100.000   5     5     100.000   1542      99.611    0      1        1542   4083001   4084542   +      4835601   0.00e+00   2755    
NC_000913.3:4166659-4168200   1542   306060   GCF_000468515.1   NC_022364.1         100.000   6     6     100.000   1542      99.481    0      1        1542   4177970   4179511   +      4835601   0.00e+00   2746    
NC_000913.3:4166659-4168200   1542   306060   GCF_000468515.1   NC_022364.1         100.000   7     7     100.000   1542      99.157    0      1        1542   3548937   3550478   -      4835601   0.00e+00   2722    
NC_000913.3:4166659-4168200   1542   306060   GCA_020564915.1   JAJBZO010000001.1   100.000   1     1     100.000   1542      100.000   0      1        1542   507976    509517    +      4637373   0.00e+00   2782    
NC_000913.3:4166659-4168200   1542   306060   GCA_020564915.1   JAJBZO010000001.1   100.000   2     2     100.000   1542      100.000   0      1        1542   1180104   1181645   +      4637373   0.00e+00   2782    
NC_000913.3:4166659-4168200   1542   306060   GCA_020564915.1   JAJBZO010000001.1   100.000   3     3     100.000   1542      100.000   0      1        1542   3674442   3675983   -      4637373   0.00e+00   2782
```


{{< /tab>}}

{{< tab "A plasmid" >}}


```plain
query        qlen    hits    sgenome           sseqid                 qcovGnm   cls   hsp   qcovHSP   alenHSP   pident    gaps   qstart   qend    sstart    send      sstr   slen      evalue      bitscore
----------   -----   -----   ---------------   --------------------   -------   ---   ---   -------   -------   -------   ----   ------   -----   -------   -------   ----   -------   ---------   --------
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000005.1   100.000   1     1     77.157    40762     99.995    0      12069    52830   5194      45955     +      51479     0.00e+00    73501   
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000005.1   100.000   2     2     10.456    5524      100.000   0      1        5524    45956     51479     +      51479     0.00e+00    9963    
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000005.1   100.000   3     3     9.860     5209      100.000   0      5525     10733   1         5209      +      51479     0.00e+00    9395    
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000005.1   100.000   4     4     1.562     827       99.758    2      9044     9868    16840     17666     +      51479     0.00e+00    1496    
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000005.1   100.000   5     5     1.565     827       99.758    2      23715    24541   3520      4344      +      51479     0.00e+00    1496    
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000005.1   100.000   6     6     0.208     110       98.182    0      29827    29936   29056     29165     -      51479     6.76e-44    190     
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000002.1   100.000   7     7     0.661     349       74.212    0      10092    10440   15029     15377     +      203534    1.04e-53    224     
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000002.1   100.000   7     8     2.531     1337      94.091    0      10733    12069   14193     15529     +      203534    0.00e+00    2055    
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000002.1   100.000   8     9     1.554     821       99.878    0      9049     9869    80484     81304     -      203534    0.00e+00    1476    
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000002.1   100.000   9     10    1.552     820       99.878    0      23722    24541   80485     81304     -      203534    0.00e+00    1474    
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000002.1   100.000   10    11    1.105     584       99.829    0      7767     8350    84382     84965     +      203534    1.41e-301   1049    
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000002.1   100.000   10    12    0.320     169       100.000   0      9133     9301    83857     84025     +      203534    1.87e-78    306     
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000002.1   100.000   11    13    0.439     232       79.310    0      23403    23634   84597     84828     +      203534    2.26e-47    203     
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000002.1   100.000   11    14    0.996     526       100.000   0      23806    24331   83857     84382     +      203534    9.16e-272   949     
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000002.1   100.000   12    15    0.509     269       99.628    0      8082     8350    84697     84965     +      203534    6.55e-131   480     
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000002.1   100.000   12    16    0.996     526       100.000   0      9133     9658    83857     84382     +      203534    9.16e-272   949     
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000002.1   100.000   13    17    0.502     265       100.000   0      19788    20052   88528     88792     +      203534    2.25e-130   479     
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000002.1   100.000   14    18    0.409     224       71.429    8      8827     9042    45434     45657     +      203534    5.28e-36    165     
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000001.1   100.000   15    19    1.554     821       100.000   0      9049     9869    2597473   2598293   +      5332228   0.00e+00    1481    
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000001.1   100.000   16    20    1.554     821       100.000   0      23721    24541   231535    232355    +      5332228   0.00e+00    1481    
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000001.1   100.000   17    21    1.554     821       100.000   0      23721    24541   2597472   2598292   +      5332228   0.00e+00    1481    
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000001.1   100.000   18    22    1.552     820       100.000   0      9049     9868    231536    232355    +      5332228   0.00e+00    1480    
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000001.1   100.000   19    23    0.502     265       100.000   0      19788    20052   4258498   4258762   +      5332228   2.25e-130   479     
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000003.1   100.000   20    24    0.507     268       99.627    0      19785    20052   58030     58297     -      141533    2.28e-130   479     
CP115019.1   52830   65311   GCF_021502915.1   NZ_JAJAAP010000003.1   100.000   21    25    0.424     224       100.000   0      19788    20011   109732    109955    +      141533    3.45e-108   405     
CP115019.1   52830   65311   GCF_016803855.1   NZ_CP048009.1          98.158    1     1     77.157    40762     99.995    0      12069    52830   7768      48529     -      51479     0.00e+00    73501   
CP115019.1   52830   65311   GCF_016803855.1   NZ_CP048009.1          98.158    2     2     14.702    7767      100.000   0      1        7767    1         7767      -      51479     0.00e+00    14008   
CP115019.1   52830   65311   GCF_016803855.1   NZ_CP048009.1          98.158    3     3     5.614     2966      100.000   0      7768     10733   48514     51479     -      51479     0.00e+00    5350    
CP115019.1   52830   65311   GCF_016803855.1   NZ_CP048009.1          98.158    4     4     1.565     827       99.758    2      23715    24541   49379     50203     -      51479     0.00e+00    1496
```

{{< /tab>}}


{{< tab "A prophage" >}}

```plain
query         qlen    hits   sgenome           sseqid          qcovGnm   cls   hsp   qcovHSP   alenHSP   pident   gaps   qstart   qend    sstart    send      sstr   slen      evalue      bitscore   species             
-----------   -----   ----   ---------------   -------------   -------   ---   ---   -------   -------   ------   ----   ------   -----   -------   -------   ----   -------   ---------   --------   --------------------
NC_001895.1   33593   2      GCF_003697165.2   NZ_CP033092.2   77.588    1     1     27.890    9371      97.716   2      1        9369    1864411   1873781   +      4903501   0.00e+00    15953      Escherichia coli    
NC_001895.1   33593   2      GCF_003697165.2   NZ_CP033092.2   77.588    1     2     0.301     101       98.020   0      10308    10408   1873846   1873946   +      4903501   1.72e-43    174        Escherichia coli    
NC_001895.1   33593   2      GCF_003697165.2   NZ_CP033092.2   77.588    2     3     20.665    6942      96.528   4      17441    24382   1882011   1888948   +      4903501   0.00e+00    11459      Escherichia coli    
NC_001895.1   33593   2      GCF_003697165.2   NZ_CP033092.2   77.588    3     4     17.685    5941      97.980   0      24355    30295   1853098   1859038   +      4903501   0.00e+00    10174      Escherichia coli    
NC_001895.1   33593   2      GCF_003697165.2   NZ_CP033092.2   77.588    4     5     8.993     3021      91.526   0      10308    13328   1873846   1876866   +      4903501   0.00e+00    4295       Escherichia coli    
NC_001895.1   33593   2      GCF_003697165.2   NZ_CP033092.2   77.588    5     6     2.438     820       84.390   1      14540    15358   1878798   1879617   +      4903501   1.29e-264   911        Escherichia coli    
NC_001895.1   33593   2      GCF_002949675.1   NZ_CP026774.1   0.976     1     1     0.976     331       85.801   3      13919    14246   3704319   3704649   -      4395762   6.35e-112   403        Shigella dysenteriae
```

{{< /tab>}}

{{< tab "Long reads" >}}


Queries are a few Nanopore Q20 reads from a mock metagenomic community.

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

2. Add information to the alignment results with [csvtk](https://github.com/shenwei356/csvtk) or other tools.

        # add species
        cat b.gene_E_coli_16S.fasta.lexicmap.tsv \
            | csvtk mutate -t --after slen -n species -f sgenome \
            | csvtk replace -t -f species -p "(.+)" -r "{kv}" -k ass2species.tsv \
            > result.with_species.tsv

        # filter results with query coverage >= 80 and count the species
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
