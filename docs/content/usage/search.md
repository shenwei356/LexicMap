---
title: search
weight: 20
---

```plain
$ lexicmap search  -h
Search sequences against an index

Attention:
  1. Input should be (gzipped) FASTA or FASTQ records from files or stdin.

Alignment result relationship:

  Query
  ├── Subject genome
      ├── Subject sequence
          ├── High-Scoring Pair (HSP)
              ├── HSP fragment

Output format:
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
    11. alen,     Aligned length in current HSP, a HSP might have >=1 HSP fragments.
    12. alenFrag, Aligned length in current HSP fragment.
    13. pident,   Percentage of identical matches in current HSP fragment.
    14. slen,     Subject sequence length.
    15. sstart,   Start of HSP fragment in subject sequence.
    16. send,     End of HSP fragment in subject sequence.
    17. sstr,     Subject strand.
    18. seeds,    Number of seeds in current HSP.

Usage:
  lexicmap search [flags] -d <index path> [query.fasta.gz ...] [-o query.tsv.gz]

Flags:
      --align-band int                 ► Band size in backtracking the score matrix (default 100)
      --align-ext-len int              ► Extend length of upstream and downstream of seed regions, for
                                       extracting query and target sequences for alignment (default 2000)
      --align-max-gap int              ► Maximum gap in a HSP fragment (default 50)
      --align-max-mismatch int         ► Maximum mismatch in a HSP fragment (default 50)
  -l, --align-min-match-len int        ► Minimum aligned length in a HSP fragment (default 50)
  -i, --align-min-match-pident float   ► Minimum base identity (percentage) in a HSP fragment. (default 70)
  -h, --help                           help for search
  -d, --index string                   ► Index directory created by "lexicmap index".
  -w, --load-whole-seeds               ► Load the whole seed data into memory for faster search.
      --max-open-files int             ► Maximum opened files. (default 512)
  -Q, --min-qcov-per-genome float      ► Minimum query coverage (percentage) per genome. (default 50)
  -q, --min-qcov-per-hsp float         ► Minimum query coverage (percentage) per HSP.
  -o, --out-file string                ► Out file, supports and recommends a ".gz" suffix ("-" for
                                       stdout). (default "-")
      --seed-max-dist int              ► Max distance between seeds in seed chaining. (default 10000)
      --seed-max-gap int               ► Max gap in seed chaining. (default 2000)
  -m, --seed-max-mismatch int          ► Minimum mismatch between non-prefix regions of shared
                                       substrings. (default -1)
  -p, --seed-min-prefix int            ► Minimum length of shared substrings (seeds). (default 15)
  -P, --seed-min-single-prefix int     ► Minimum length of shared substrings if there's only one pair.
                                       (default 20)
  -n, --top-n-genomes int              ► Keep top N genome matches for a query (0 for all). (default 500)

Global Flags:
  -X, --infile-list string   ► File of input files list (one file per line). If given, they are
                             appended to files from CLI arguments.
      --log string           ► Log file.
      --quiet                ► Do not print any verbose information. But you can write them to a file
                             with --log.
  -j, --threads int          ► Number of CPUs cores to use. By default, it uses all available cores.
                             (default 16)
```
