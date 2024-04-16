---
title: Searching
weight: 10
---

## Table of contents

{{< toc format=html >}}

## TL;DR

    # For short queries like genes or long reads, returning top N hits.
    lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
        --min-qcov-per-genome 70 --min-match-pident 70 --min-qcov-per-hsp 70 --top-n-genomes 500

    # For longer queries like plasmids, returning all hits.
    lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
        --min-qcov-per-genome 50 --min-match-pident 70 --min-qcov-per-hsp 0  --top-n-genomes 0

## Input

{{< hint type=note >}}
**Query length**\
LexicMap only supports long (>=500 bp) reads or gene/genome/virus/plasmid/phage sequences.\
However, some short queries can also be aligned.
{{< /hint >}}

Input should be (gzipped) FASTA or FASTQ records from files or STDIN.

## Algorithm

<img src="/LexicMap/searching.svg" alt="" width="900"/>

1. **Masking:**
   Query sequence is masked by the masks of the index. In other way, each mask captures the most similar k-mer and stores its posistion and strand information.
1. **Seeding:**
   For each mask, the captured k-mer is used to search seeds (captured k-mers in reference genomes) sharing prefixes of at least *p* bases.
    1. **Setting the search range**: Since the seeded k-mers are stored in lexicographic order, the k-mer matching turns into a range query.
       For example, for a query `CATGCT` requiring matching at least 4-bp prefix is equal to extract k-mers from `CATGAA` to `CATGTT`.
    2. **Finding the nearest offset**: The index file of each seed data file stores a list (default 512) of k-mers and offsets in the data file, and is read in RAM.
       The nearest k-mer smaller than the range start k-mer (`CATGAA`) is found by binary searching, i.e., `CATCAC` (blue text in the fingure), and the offset is returned.
    3. **Retrieving seed data**: Seed k-mers are read from the file and checked one by one, and k-mers in the search range are returned, along with the k-mer information (genome batch, genome number, location, and strand).
1. **Chaining:**
    1. Seeding results, i.e., anchors (matched k-mers from the query and subject sequence), are summarized by genome, and deduplicated.
    2. Performing chaining.
1. **Alignment** for each chain.
    1. Extending the anchor region. for extracting sequences from the query and reference genome. For example, extending 2 kb in upstream and downstream of anchor region.
    2. Fast alignment of query and subject sequences.
    3. Filtering aligned segments and the whole HSPs (all alignment segments) based on user options.
       - For these HSPs that accross more than one reference sequences, splitting them into multiple HSPs.

## Hardware requirements

LexicMap is designed to provide fast and low-memory sequence alignment against millions of prokaryotic genomes.

- **CPU:**
    - No specific requirements on CPU type and instruction sets. Both x86 and ARM chips are supported.
    - More is better as LexicMap is a CPU-intensive software. It uses all CPUs by default (`-j/--threads`).
- **RAM**
    - More RAM (> 50 GB) is preferred. The memory usage in searching is related to:
    - If the RAM is not sufficient (< 10 GB). Please:
- **Disk**


## Parameters

{{< tabs "t1" >}}

{{< tab "General" >}}


{{< /tab>}}

{{< tab "Chaining" >}}


{{< /tab>}}

{{< tab "Alignment" >}}


{{< /tab>}}

{{< /tabs >}}


## Steps

## Output

### Alignment result relationship

    Query
    ├── Subject genome                             # A query might have one or more genome hits,
        ├── Subject sequence                       # in different sequences.
            ├── High-Scoring segment Pairs (HSP)   # HSP is a cluster of alignment segments.
                ├── HSP segment                    # a local alignment with no gaps.

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

{{< tab "A gene" >}}


{{< /tab>}}

{{< tab "A plasmid" >}}


{{< /tab>}}

{{< tab "Long reads" >}}


{{< /tab>}}

{{< /tabs >}}
