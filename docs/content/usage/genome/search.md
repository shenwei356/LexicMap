---
title: lexicmap genome search
linkTitle: search
weight: 30
---

## Usage

```plain
$ lexicmap genome search -h
Search genomes against an index, with ANI and AF computed

Algorithm:
  1. Genome screening: candidate genomes are screened by the total length of shared seeds
     (k-mers longer than the value of -p/--seed-min-prefix) between the query genome and
     the genomes in the index, with lexichash masking.
  2. Alignment:
     a) (default): cut the query genome into non-overlapping fragments and align them to
        the candidate genomes with an approach similar to that used in 'lexicmap search',
        with steps of seed matching, chaining, pseudoalignment, and base-level alignment.
     b) (OrthoANI mode): cut both the query genome and candidate genomes into non-overlapping
        fragments and only orthologous fragment pairs are used for calculating the ANI and
        AF values, which is similar to the algorithm of OrthoANI.

Attention:
  1. Input should be (gzipped) FASTA records from files or stdin, with one genome per file.
     Experimental feature: Long-read FASTQ files are also accepted. For these inputs,
     --adjust-ani-by-quality can be used to estimate ANI corrected for sequencing errors:
         e = mean_i(10^(-Q_i/10)) and ANIadj = (ANI - e/3) / (1 - 4e/3),
     where Q_i are the Phred scores in aligned query regions and ANI is expressed as a
     fraction. The correction assumes independent, symmetric substitution errors.
  2. One or more input files are accepted, via positional parameters
     and/or a file list via the flag -X/--infile-list.
  3. For multiple queries, the order of queries in output might be different from the input.

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

  3. Flag --show-genome-name and --show-species-name can be used to show the taxonomic/species
     name of the subject genome in the output.

Output format:
  Tab-delimited format with 10 columns. --adjust-ani-by-quality appends
  ANIadj as the last column.

    1.  query,    Query genome ID.
    2.  subject,  Subject genome ID.
    3.  ANI,      Average nucleotide identity.
    4.  qAF,      Align fraction of the query genome (sum_aligned_length / sum_fragment_length).
    5.  sAF,      Align fraction of the subject genome.
                  Default mode : sum_aligned_length / sum_fragment_length)
                  OrthoANI mode: sum_aligned_length / genome size)
    6.  qctgs,    Number of contigs in the query genome.
    7.  qsize,    Size of the query genome.
    8.  sctgs,    Number of contigs in the subject genome.
    9.  ssize,    Size of the subject genome.
    10. sname,    Taxonomic name of the subject name.
    11. ANIadj,   Quality-adjusted ANI (optional, experimental feature).

Usage:
  lexicmap genome search [flags] -d <index path> [query.fasta[.gz] ...] [-o result.tsv[.gz]]

Flags:
  -O, --OrthoANI                       ► Compute OrthoANI using reciprocal best hit of fragment pairs.
                                       Type 'lexicmap genome search --help' for details.
      --adjust-ani-by-quality          ► (Experimental feature) Estimate ANI corrected for Phred+33
                                       sequencing errors in the aligned query regions and add the
                                       ANIAdjusted column. FASTQ input is required. The raw ANI is still
                                       used for filtering and sorting.
      --align-band int                 ► Band size in backtracking the score matrix (pseudo alignment
                                       phase). (default 100)
      --align-max-gap int              ► Maximum gap in a HSP segment. (default 100)
  -l, --align-min-match-len int        ► Minimum aligned length in a HSP segment. (default 30)
  -i, --align-min-match-pident float   ► Minimum base identity (percentage) in a HSP segment. (default 70)
      --debug                          ► Print debug information, including a progress bar.
                                       (recommended when searching with one query).
      --extra                          ► Show extra columns for -S/--only-genome-screening.
      --frag-size int                  ► The size of non-overlap fragments cut for ANI computation.
                                       (default 1020)
      --gc-interval int                ► Force garbage collection every N queries (0 for disable). The
                                       value can't be too small. (default 4)
  -G, --genome2taxid string            ► Two-column tabular file for mapping genome ID to TaxId,
                                       needed for filtering results with TaxIds. Genome IDs in the index
                                       can be exported via "lexicmap utils genomes -d db.lmi/ | csvtk
                                       cut -t -f 1 | csvtk uniq -Ut". (default: <taxdump path>/taxid.map)
  -h, --help                           help for search
  -d, --index string                   ► Index directory created by "lexicmap index".
  -k, --keep-genomes-without-taxid     ► Keep genome hits without TaxId, i.e., those without TaxId in
                                       the --genome2taxid file.
      --kmer-scale int                 ► Using 1/scale of k-mers for seeding (default mode) or
                                       fragment comparison (OrthoANI mode). Available values: 2, 4, 8.
                                       (default 4)
  -w, --load-whole-seeds               ► Load the whole seed data into memory for faster seed
                                       matching. It will consume a lot of RAM.
  -m, --masks int                      ► Only use seed data of N masks for genome
                                       prefiltering/screening. It should be 0 (for all masks in the
                                       index) or power of 4 (needs to be >= 256, e.g., 256, 1024, 4096,
                                       16384). (default 1024)
  -e, --max-evalue float               ► Maximum evalue of a HSP segment. (default 1e-15)
      --max-open-files int             ► Maximum opened files. It mainly affects candidate genome
                                       extraction. Increase this value if you have hundreds of genome
                                       batches or have multiple queries, and do not forgot to set a
                                       bigger "ulimit -n" in shell if the value is > 1024. (default 1024)
  -J, --max-query-conc int             ► Maximum number of concurrent queries. (default 8)
      --max-subject-genome-size int    ► Maximum size of subject genomes to be considered (in MB).
                                       (default 20)
  -F, --min-af float                   ► Only output results where one genome has aligned fraction >
                                       than this value (percentage). (default 15)
  -I, --min-ani float                  ► Only output results where one genome has ANI > than this
                                       value (percentage). (default 70)
      --min-frag-size int              ► The minimum length of fragments in the end of a sequence
                                       during cutting fragments. (default 100)
  -q, --min-qcov-per-hsp float         ► Minimum query coverage (percentage) per HSP. (default 30)
  -S, --only-genome-screening          ► Only perform genome screening, no ANI computation.
  -o, --out-file string                ► Out file, supports a ".gz" suffix ("-" for stdout). (default "-")
      --ref-name-regexp string         ► Regular expression (must contains "(" and ")") for extracting
                                       the reference name from the input filename. Attention: use double
                                       quotation marks for patterns containing commas, e.g., -p
                                       '"A{2,}"'. (default
                                       "(?i)(.+)\\.(f[aq](st[aq])?|fna)(\\.gz|\\.xz|\\.zst|\\.bz2)?$")
  -p, --seed-min-prefix int            ► Minimum prefix length of matched seeds in the genome
                                       filtering phase. (default 21)
  -g, --show-genome-name               ► Show the taxonomic name of subject genome in field 'sname'.
                                       Flags -T/--taxdump and -G/--genome2taxid are needed.
  -s, --show-species-name              ► Show the species name of subject genome in field 'sname'.
                                       Flags -T/--taxdump and -G/--genome2taxid are needed.
  -T, --taxdump string                 ► Directory containing taxdump files (nodes.dmp, names.dmp,
                                       etc.), needed for filtering results with TaxIds. For other
                                       non-NCBI taxonomy data, please use 'taxonkit create-taxdump' to
                                       create taxdump files. (default: <index path>/taxdump)
      --taxid-file string              ► TaxIds from a file for filtering results, where the taxids
                                       are equal to or are the children of the given taxids. Negative
                                       values are allowed as a black list.
  -t, --taxids strings                 ► TaxIds(s) for filtering results, where the taxids are equal
                                       to or are the children of the given taxids. Negative values are
                                       allowed as a black list.
  -N, --top-n-chains int               ► Keep the top N chains in a genome for the query (0 for all)
                                       in the chaining phase. Value 1 is not recommended as the best
                                       chaining result does not always bring the best alignment. (default 5)
  -n, --top-n-genomes int              ► Keep the top N genome matches for a query (0 for all).
                                       Screening includes all candidates tied with the Nth score; ANI
                                       ranking then returns at most N matches. (default 10)
      --windows int                    ► The number of windows in lexichash masking, for genome
                                       screening (1 is sufficient). (default 1)

Global Flags:
  -X, --infile-list string   ► File of input file list (one file per line). If given, they are
                             appended to files from CLI arguments.
      --log string           ► Log file.
      --quiet                ► Do not print any verbose information. But you can write them to a file
                             with --log.
  -j, --threads int          ► Number of CPU cores to use. By default, it uses all available cores.
                             (default 16)
```

Examples:

1. Simple search

        $ lexicmap genome search -d demo.lmi/ refs/GCF_003697165.2.fa.gz \
            | csvtk pretty -t
    
        query             subject           ANI       qAF       sAF       qcontigs   qsize     scontigs   ssize     sname
        ---------------   ---------------   -------   -------   -------   --------   -------   --------   -------   -----
        GCF_003697165.2   GCF_003697165.2   100.000   100.000   100.000   2          5034834   2          5034834        
        GCF_003697165.2   GCF_002949675.1   96.821    67.230    73.932    2          5034834   2          4578459        
        GCF_003697165.2   GCF_002950215.1   96.735    72.831    74.255    2          5034834   3          4938295        
        GCF_003697165.2   GCF_000006945.2   81.614    54.733    55.655    2          5034834   2          4951383        
        GCF_003697165.2   GCF_000742135.1   80.588    45.868    41.642    2          5034834   5          5545784        
        
2. Show genome names
    
        $ lexicmap genome search -d demo.lmi/ refs/GCF_003697165.2.fa.gz \
              --taxdump ~/.taxonkit --genome2taxid taxid.map --show-genome-name  \
            | csvtk pretty -t

        query             subject           ANI       qAF       sAF       qcontigs   qsize     scontigs   ssize     sname                
        ---------------   ---------------   -------   -------   -------   --------   -------   --------   -------   ---------------------
        GCF_003697165.2   GCF_003697165.2   100.000   100.000   100.000   2          5034834   2          5034834   Escherichia coli     
        GCF_003697165.2   GCF_002949675.1   96.821    67.230    73.932    2          5034834   2          4578459   Shigella dysenteriae 
        GCF_003697165.2   GCF_002950215.1   96.735    72.831    74.255    2          5034834   3          4938295   Shigella flexneri    
        GCF_003697165.2   GCF_000006945.2   81.614    54.733    55.655    2          5034834   2          4951383   Salmonella enterica  
        GCF_003697165.2   GCF_000742135.1   80.588    45.868    41.642    2          5034834   5          5545784   Klebsiella pneumoniae
