---
title: search
weight: 20
---

```plain
$ lexicmap search -h
Search sequences against an index

Attention:
  1. Input should be (gzipped) FASTA or FASTQ records from files or stdin.
  2. One or more input files are accepted, via positional parameters
     and/or a file list via the flag -X/--infile-list.
  3. For multiple queries, the order of queries in output might be different from the input.

Tips:
  1. When using -a/--all, the search result would be formatted to Blast-style format
     with 'lexicmap utils 2blast'. And the search speed would be slightly slowed down.
  2. Alignment result filtering is performed in the final phase, so stricter filtering criteria,
     including -q/--min-qcov-per-hsp, -Q/--min-qcov-per-genome, and -i/--align-min-match-pident,
     do not significantly accelerate the search speed. Hence, you can search with default
     parameters and then filter the result with tools like awk or csvtk.
  3. For searches with -a/--all, LexicMap retains up to 1 GiB of CIGAR strings, aligned
     query/subject sequences, and alignment text in memory by default. When the global budget
     is exceeded, each affected query stores these fields in one file in the system temporary
     directory and removes it after output. Concurrent queries use separate files. On Unix,
     set TMPDIR to choose a temporary directory on a fast disk with sufficient free space.
     Use --max-align-result-memory to change the budget or set it to 0 to disable spilling.
     This limit applies only to these output fields and is not a total process memory limit.
  4. Queries with very many seed matches in large indexes can run out of memory
     during seed collection. Reduce -J/--max-query-conc to lower the memory used
     by concurrent queries. --max-seed-memory enables experimental seed spilling
     (disabled by default) to limit anchor collection buffers.
     The budget is divided among up to -J concurrent collection slots. Decoded seed data
     are delivered in small batches; anchors spill as sorted runs before a buffer growth
     would exceed its query share. Runs are merged by genome for the existing chaining.
     Final Top-N selection, including cutoff ties, still uses chaining scores.
     Files use the system temporary directory (TMPDIR) and are removed after chaining.
     This budget excludes one genome's complete anchor array, chaining scratch,
     fixed I/O buffers, candidate metadata, loaded index data, and alignment memory.
     It is not a total RSS limit.

Taxonomic operations:
  1. Taxonomy data, including NCBI-format taxdump files (-T/--taxdump) and a genome-ID-to-TaxId
     mapping file (-G/--genome2taxid), are needed for filtering genomes by TaxId(s) and/or
     showing taxonomic names.

     Taxdump files can be createted from any taxonomy data with TaxonKit, see
     https://bioinf.shenwei.me/taxonkit/usage/#create-taxdump

     If -T and -G are not provided, it will try to find them in the index directory,
     in the subdirectory "taxdump" and the file "taxdump/taxid.map", respectively.

  2. Users can limit search by TaxId(s) via -t/--taxids or --taxid-file.
     Only genomes with descendant TaxIds of the specific ones or themselves are searched,
     in a similar way with BLAST+ 2.15.0 or later versions.
     Negative values are allowed as a black list.

     For example, searching non-Escherichia (561) genera of Enterobacteriaceae (543) family with
     -t 543,-561.

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

Result ordering:
  For a HSP cluster, SimilarityScore = max(bitscore*pident)
  1. Within each HSP cluster, HSPs are sorted by sstart.
  2. Within each subject genome, HSP clusters are sorted in descending order by SimilarityScore.
  3. Results of multiple subject genomes are sorted by the highest SimilarityScore of HSP clusters.

Usage:
  lexicmap search [flags] -d <index path> [query.fasta[.gz] ...] [-o result.tsv[.gz]]

Flags:
      --align-band int                   ► Band size in backtracking the score matrix (pseudo
                                         alignment phase). (default 100)
      --align-ext-len int                ► Extend length of upstream and downstream of seed regions,
                                         for extracting query and target sequences for alignment. It
                                         should be <= contig interval length in database. (default 1000)
      --align-max-gap int                ► Maximum gap in a HSP segment. (default 20)
  -l, --align-min-match-len int          ► Minimum aligned length in a HSP segment. (default 50)
  -i, --align-min-match-pident float     ► Minimum base identity (percentage) in a HSP segment.
                                         (default 70)
  -a, --all                              ► Output more columns, e.g., matched sequences. Use this if
                                         you want to output blast-style format with "lexicmap utils 2blast".
      --debug                            ► Print debug information, including a progress bar.
                                         (recommended when searching with one query).
      --gc-interval int                  ► Force garbage collection every N queries (0 for disable).
                                         The value can't be too small. (default 64)
  -G, --genome2taxid string              ► Two-column tabular file for mapping genome ID to TaxId,
                                         needed for filtering results with TaxIds. Genome IDs in the
                                         index can be exported via "lexicmap utils genomes -d db.lmi/ |
                                         csvtk cut -t -f 1 | csvtk uniq -Ut". (default: <taxdump
                                         path>/taxid.map)
  -h, --help                             help for search
  -d, --index string                     ► Index directory created by "lexicmap index".
  -k, --keep-genomes-without-taxid       ► Keep genome hits without TaxId, i.e., those without TaxId
                                         in the --genome2taxid file.
  -w, --load-whole-seeds                 ► Load the whole seed data into memory for faster seed
                                         matching. It will consume a lot of RAM.
      --max-align-result-memory string   ► Maximum memory for retaining CIGAR, query sequence, subject
                                         sequence, and alignment text across concurrent queries. Values
                                         support K/M/G/T suffixes. When the global budget is exceeded,
                                         the affected query spills these fields to a temporary file.
                                         This is not a total RSS limit (0 disables spilling). (default "1G")
  -e, --max-evalue float                 ► Maximum evalue of a HSP segment. (default 10)
      --max-open-files int               ► Maximum opened files. It mainly affects candidate
                                         subsequence extraction. Increase this value if you have
                                         hundreds of genome batches or have multiple queries, and do not
                                         forgot to set a bigger "ulimit -n" in shell if the value is >
                                         1024. (default 1024)
  -J, --max-query-conc int               ► Maximum number of concurrent queries. Bigger values do not
                                         improve the batch searching speed and consume much memory.
                                         Reduce this value when memory is limited. (default 8)
      --max-seed-memory string           ► Experimental anchor collection buffer budget shared across
                                         up to -J query slots (K/M/G/T suffixes; 0 disables spilling).
                                         Uses TMPDIR for temporary files. Not a total memory limit. See
                                         Tips in --help for details. (default "0")
  -Q, --min-qcov-per-genome float        ► Minimum query coverage (percentage) per genome.
  -q, --min-qcov-per-hsp float           ► Minimum query coverage (percentage) per HSP.
  -o, --out-file string                  ► Out file, supports a ".gz" suffix ("-" for stdout).
                                         (default "-")
      --seed-max-dist int                ► Minimum distance between seeds in seed chaining. It should
                                         be <= contig interval length in database. (default 1000)
      --seed-max-gap int                 ► Minimum gap in seed chaining. (default 50)
  -p, --seed-min-prefix int              ► Minimum (prefix/suffix) length of matched seeds (anchors).
                                         (default 15)
  -P, --seed-min-single-prefix int       ► Minimum (prefix/suffix) length of matched seeds (anchors)
                                         if there's only one pair of seeds matched. (default 17)
      --show-avg-qual                    ► Add average quality of the aligned region as a suffix to
                                         alenHSP field, e.g., 128:21.6, where 21.6 is the average
                                         quality of the aligned region.
  -g, --show-genome-name                 ► Add the taxonomic name as a prefix to sgenome fied. Flags
                                         -T/--taxdump and -G/--genome2taxid are needed.
  -s, --show-species-name                ► Add the species name as a prefix to sgenome fied. Flags
                                         -T/--taxdump and -G/--genome2taxid are needed.
      --show-sseq-idx                    ► Add 1-based genome chunk and subject-sequence index
                                         prefixes to sseqid values, e.g., c2/3:s1/10:contig00001, where
                                         c2/3 denotes chunk 2 of 3 and s1/10 denotes sequence 1 of 10.
  -T, --taxdump string                   ► Directory containing taxdump files (nodes.dmp, names.dmp,
                                         etc.), needed for filtering results with TaxIds. For other
                                         non-NCBI taxonomy data, please use 'taxonkit create-taxdump' to
                                         create taxdump files. (default: <index path>/taxdump)
      --taxid-file string                ► TaxIds from a file for filtering results, where the taxids
                                         are equal to or are the children of the given taxids. Negative
                                         values are allowed as a black list.
  -t, --taxids strings                   ► TaxIds(s) for filtering results, where the taxids are equal
                                         to or are the children of the given taxids. Negative values are
                                         allowed as a black list.
  -N, --top-n-chains int                 ► Keep the top N chains in a genome for the query (0 for all)
                                         in the chaining phase. Value 1 is not recommended as the best
                                         chaining result does not always bring the best alignment, so
                                         it's better be >= 10. (default 0)
  -n, --top-n-genomes int                ► Keep the top N genome matches for a query (0 for all) in
                                         the chaining phase, including all matches tied at the cutoff
                                         score. Value 1 is not recommended as the best chaining result
                                         does not always bring the best alignment, so it's better be >=
                                         100. (default 0)

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
