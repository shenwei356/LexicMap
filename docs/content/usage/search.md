---
title: search
weight: 20
---

```plain
$ lexicmap search  -h
Search sequences against an index

Attention:
  1. Input should be (gzipped) FASTA or FASTQ records from files or stdin.
  2. For multiple queries, the order of queries in output might be different from the input.

Tips:
  1. When using -a/--all, the search result would be formatted to Blast-style format
     with 'lexicmap utils 2blast'. And the search speed would be slightly slowed down.
  2. Alignment result filtering is performed in the final phase, so stricter filtering criteria,
     including -q/--min-qcov-per-hsp, -Q/--min-qcov-per-genome, and -i/--align-min-match-pident,
     do not significantly accelerate the search speed. Hence, you can search with default
     parameters and then filter the result with tools like awk or csvtk.

Alignment result relationship:

  Query
  ├── Subject genome
      ├── Subject sequence
          ├── HSP cluster (a cluster of neighboring HSPs)
              ├── High-Scoring segment Pair (HSP)

  Here, the defination of HSP is similar with that in BLAST. Actually there are small gaps in HSPs.

  > A High-scoring Segment Pair (HSP) is a local alignment with no gaps that achieves one of the
  > highest alignment scores in a given search. https://www.ncbi.nlm.nih.gov/books/NBK62051/

Output format:
  Tab-delimited format with 20+ columns, with 1-based positions.

    1.  query,    Query sequence ID.
    2.  qlen,     Query sequence length.
    3.  hits,     Number of subject genomes.
    4.  sgenome,  Subject genome ID.
    5.  sseqid,   Subject sequence ID.
    6.  qcovGnm,  Query coverage (percentage) per genome: $(aligned bases in the genome)/$qlen.
    7.  cls,      Nth HSP cluster in the genome. (just for improving readability)
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

Result ordering:
  For a HSP cluster, SimilarityScore = max(bitscore*pident)
  1. Within each HSP cluster, HSPs are sorted by sstart.
  2. Within each subject genome, HSP clusters are sorted in descending order by SimilarityScore.
  3. Results of multiple subject genomes are sorted by the highest SimilarityScore of HSP clusters.

Usage:
  lexicmap search [flags] -d <index path> [query.fasta.gz ...] [-o query.tsv.gz]

Flags:
      --align-band int                 ► Band size in backtracking the score matrix (pseudo alignment
                                       phase). (default 100)
      --align-ext-len int              ► Extend length of upstream and downstream of seed regions, for
                                       extracting query and target sequences for alignment. It should be
                                       <= contig interval length in database. (default 1000)
      --align-max-gap int              ► Maximum gap in a HSP segment. (default 20)
  -l, --align-min-match-len int        ► Minimum aligned length in a HSP segment. (default 50)
  -i, --align-min-match-pident float   ► Minimum base identity (percentage) in a HSP segment. (default 70)
  -a, --all                            ► Output more columns, e.g., matched sequences. Use this if you
                                       want to output blast-style format with "lexicmap utils 2blast".
      --debug                          ► Print debug information, including a progress bar.
                                       (recommended when searching with one query).
  -h, --help                           help for search
  -d, --index string                   ► Index directory created by "lexicmap index".
  -w, --load-whole-seeds               ► Load the whole seed data into memory for faster seed
                                       matching. It will consume a lot of RAM.
  -e, --max-evalue float               ► Maximum evalue of a HSP segment. (default 10)
      --max-open-files int             ► Maximum opened files. It mainly affects candidate subsequence
                                       extraction. Increase this value if you have hundreds of genome
                                       batches or have multiple queries, and do not forgot to set a
                                       bigger "ulimit -n" in shell if the value is > 1024. (default 1024)
  -J, --max-query-conc int             ► Maximum number of concurrent queries. Bigger values do not
                                       improve the batch searching speed and consume much memory.
                                       (default 12)
  -Q, --min-qcov-per-genome float      ► Minimum query coverage (percentage) per genome.
  -q, --min-qcov-per-hsp float         ► Minimum query coverage (percentage) per HSP.
  -o, --out-file string                ► Out file, supports a ".gz" suffix ("-" for stdout). (default "-")
      --seed-max-dist int              ► Minimum distance between seeds in seed chaining. It should be
                                       <= contig interval length in database. (default 1000)
      --seed-max-gap int               ► Minimum gap in seed chaining. (default 50)
  -p, --seed-min-prefix int            ► Minimum (prefix/suffix) length of matched seeds (anchors).
                                       (default 15)
  -P, --seed-min-single-prefix int     ► Minimum (prefix/suffix) length of matched seeds (anchors) if
                                       there's only one pair of seeds matched. (default 17)
  -n, --top-n-genomes int              ► Keep top N genome matches for a query (0 for all) in chaining
                                       phase. Value 1 is not recommended as the best chaining result
                                       does not always bring the best alignment, so it better be >= 100.
                                       (default 0)

Global Flags:
  -X, --infile-list string   ► File of input file list (one file per line). If given, they are
                             appended to files from CLI arguments.
      --log string           ► Log file.
      --quiet                ► Do not print any verbose information. But you can write them to a file
                             with --log.
  -j, --threads int          ► Number of CPU cores to use. By default, it uses all available cores.
                             (default 16)
```


## Examples

See {{< button size="small" relref="tutorials/search" >}}Searching{{< /button >}}
