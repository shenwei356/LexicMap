---
title: search
weight: 20
---

```plain
$ lexicmap search  -h
Search sequences against an index

Attention:
  1. Input should be (gzipped) FASTA or FASTQ records from files or stdin.
  2. For multiple queries, the order of queries might be different from the input.

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
          ├── High-Scoring segment Pair (HSP)

  Here, the defination of HSP is similar with that in BLAST. Actually there are small gaps in HSPs.

  > A High-scoring Segment Pair (HSP) is a local alignment with no gaps that achieves one of the
  > highest alignment scores in a given search. https://www.ncbi.nlm.nih.gov/books/NBK62051/

Output format:
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

Usage:
  lexicmap search [flags] -d <index path> [query.fasta.gz ...] [-o query.tsv.gz]

Flags:
      --align-band int                 ► Band size in backtracking the score matrix (pseduo alignment
                                       phase). (default 50)
      --align-ext-len int              ► Extend length of upstream and downstream of seed regions, for
                                       extracting query and target sequences for alignment. It should be
                                       <= contig interval length in database. (default 1000)
      --align-max-gap int              ► Maximum gap in a HSP segment. (default 20)
  -l, --align-min-match-len int        ► Minimum aligned length in a HSP segment. (default 50)
  -i, --align-min-match-pident float   ► Minimum base identity (percentage) in a HSP segment. (default 70)
  -a, --all                            ► Output more columns, e.g., matched sequences. Use this if you
                                       want to output blast-style format with "lexicmap utils 2blast".
      --debug                          ► Print debug information.
  -h, --help                           help for search
  -d, --index string                   ► Index directory created by "lexicmap index".
  -w, --load-whole-seeds               ► Load the whole seed data into memory for faster search.
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
      --pseudo-align                   ► Only perform pseudo alignment, alignment metrics, including
                                       qcovGnm, qcovSHP and pident, will be less accurate.
      --seed-max-dist int              ► Max distance between seeds in seed chaining. It should be <=
                                       contig interval length in database. (default 1000)
      --seed-max-gap int               ► Max gap in seed chaining. (default 200)
  -p, --seed-min-prefix int            ► Minimum (prefix) length of matched seeds. (default 15)
  -P, --seed-min-single-prefix int     ► Minimum (prefix) length of matched seeds if there's only one
                                       pair of seeds matched. (default 17)
  -n, --top-n-genomes int              ► Keep top N genome matches for a query (0 for all) in chaining
                                       phase. Value 1 is not recommended as the best chaining result
                                       does not always bring the best alignment, so it better be >= 5.

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
