---
title: lexicmap genome compare
linkTitle: compare
weight: 10
---

## Usage

```plain
$ lexicmap genome compare -h
Compare genome pairs and compute ANI and AF

Input:
  - Option 1:
    Two or more FASTA files. All combinations of 2 genomes will be compared.
  - Option 2:
    Two or more genome IDs, with the LexicMap index given by '-d/--index'.
    All combinations of 2 genomes will be compared.
  - Option 3:
    One or more tab-delimited files given by '--pair-file', with genome IDs in
    the first two columns. The files can be the output of 'lexicmap genome pair'.
    The LexicMap index must be given by '-d/--index', and '-H' skips header lines.

Attention:
  1. Fine alignment prioritizes the highest-scoring pseudoalignment chain.
     Repeated sequences can produce many k-mer matches, so alignment time can differ
     substantially between genomes of similar sizes.
  2. Prepared genomes and fragment indexes are reused across pairs within
     --max-genome-cache-memory. This budget does not limit total process memory.
  3. Both genomes must meet --min-af. In the default mode, at least one direction
     must meet --min-ani; OrthoANI uses its single reciprocal-best-hit ANI.
     Both thresholds are inclusive.
  4. FASTA input is case-insensitive. Soft masking does not remove lowercase
     A/C/G/T from stored genomes. Two-bit encoding replaces ambiguous bases
     (including N) with A/C/G/T in indexed genomes.
     ANI and AF can therefore differ between file and index modes in ambiguous regions.

Output format:
  Tab-delimited format with 11 columns. ANI values are printed with four decimal places
  and AF values with three.

    1.  genome1,  Genome 1.
    2.  genome2,  Genome 2.
    3.  tANI,     Total Average nucleotide identity, calculated by dividing the total matched
                  bases in pairwise alignments (genome 1 vs genome 2 and genome 2 vs genome 1)
                  by the total genome sizes of two genomes.
    4.  ANI1,     Average nucleotide identity when aligning genome 1 to genome 2.
    5.  ANI2,     Average nucleotide identity when aligning genome 2 to genome 1.
    6.  AF1,      Align fraction of genome 1 (sum_aligned_length / genome1_size).
    7.  AF2,      Align fraction of genome 2 (sum_aligned_length / genome2_size).
                  Terminal pieces shorter than half a fragment join the preceding fragment
                  of the same contig. Contigs shorter than the alignment seed length are
                  not aligned but remain included in genome sizes for AF.
    8.  ctgs1,    Number of contigs in genome 1.
    9.  size1,    Size of the genome 1.
    10. ctgs2,    Number of contigs in genome 2.
    11. size2,    Size of the genome 1.

Usage:
  lexicmap genome compare [flags] 

Flags:
  -O, --OrthoANI                         ► Compute OrthoANI using reciprocal best hit of fragment pairs.
      --align-band int                   ► Band size in backtracking the score matrix (pseudo
                                         alignment phase). (default 100)
      --align-max-gap int                ► Maximum gap in a HSP segment. (default 100)
  -l, --align-min-match-len int          ► Minimum aligned length in a HSP segment. (default 30)
  -i, --align-min-match-pident float     ► Minimum base identity (percentage) in a HSP segment.
                                         (default 70)
      --debug                            ► Print debug information, including a progress bar.
                                         (recommended when searching with one query).
      --frag-size int                    ► The size of non-overlap fragments cut for ANI computation.
                                         (default 1020)
  -f, --full-input-path                  ► Use the full path of input file as the genome name.
      --gc-interval int                  ► Force garbage collection every N queries (0 for disable).
                                         The value can't be too small. (default 128)
  -h, --help                             help for compare
  -d, --index string                     ► Index directory created by "lexicmap index". When given,
                                         positional arguments are genome IDs.
      --kmer-scale int                   ► Using 1/scale of k-mers for seeding (default mode) or
                                         fragment comparison (OrthoANI mode). Available values: 2, 4, 8.
                                         (default 4)
  -e, --max-evalue float                 ► Maximum evalue of a HSP segment. (default 1e-15)
      --max-genome-cache-memory string   ► Memory budget for reusing prepared genomes and
                                         fine-alignment indexes across pairs (supports K/M/G; 0 disables
                                         reuse). Uses conservative cache accounting; active uncached
                                         comparisons, alignment scratch, readers, and runtime memory are
                                         additional. (default "4G")
      --max-genome-size int              ► Maximum size of genomes to be considered (in MB) when
                                         reading genome from an index. (default 20)
      --max-open-files int               ► Maximum opened files. It mainly affects candidate genome
                                         extraction. Increase this value if you have hundreds of genome
                                         batches or have multiple queries, and do not forgot to set a
                                         bigger "ulimit -n" in shell if the value is > 1024. (default 1024)
  -F, --min-af float                     ► Minimum aligned fraction for both genomes (percentage,
                                         inclusive).
  -I, --min-ani float                    ► Minimum ANI in at least one direction (percentage,
                                         inclusive). (default 70)
  -q, --min-qcov-per-hsp float           ► Minimum query coverage (percentage) per HSP. (default 30)
  -o, --out-file string                  ► Out file, supports a ".gz" suffix ("-" for stdout).
                                         (default "-")
  -P, --pair-file strings                ► Tab-delimited file(s) containing genome-ID pairs in the
                                         first two columns. Requires -d/--index; can be repeated.
      --ref-name-regexp string           ► Regular expression (must contains "(" and ")") for
                                         extracting the reference name from the input filename.
                                         Attention: use double quotation marks for patterns containing
                                         commas, e.g., -p '"A{2,}"'. (default
                                         "(?i)(.+)\\.(f[aq](st[aq])?|fna)(\\.gz|\\.xz|\\.zst|\\.bz2)?$")
  -H, --skip-header-line                 ► Skip the header line in every --pair-file.

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

1. From two or more sequence files

        $ lexicmap genome compare refs/GCF_003697165.2.fa.gz refs/GCF_002949675.1.fa.gz refs/GCF_002950215.1.fa.gz \
            | csvtk pretty -t
    
        genome1           genome2           tANI      ANI1      ANI2      AF1      AF2      ctgs1   size1     ctgs2   size2  
        ---------------   ---------------   -------   -------   -------   ------   ------   -----   -------   -----   -------
        GCF_003697165.2   GCF_002949675.1   71.8008   96.8012   96.5448   67.613   81.014   2       5034834   2       4578459
        GCF_003697165.2   GCF_002950215.1   73.4098   96.7220   96.8348   73.252   77.978   2       5034834   3       4938295
        GCF_002949675.1   GCF_002950215.1   77.0516   97.0923   97.6772   82.972   74.934   2       4578459   3       4938295

2. From an index and some genome IDs

        $ lexicmap genome compare -d demo.lmi/ GCF_003697165.2 GCF_002949675.1 GCF_002950215.1 \
            | csvtk pretty -t
  
        genome1           genome2           tANI      ANI1      ANI2      AF1      AF2      ctgs1   size1     ctgs2   size2  
        ---------------   ---------------   -------   -------   -------   ------   ------   -----   -------   -----   -------
        GCF_003697165.2   GCF_002949675.1   71.8008   96.8012   96.5448   67.613   81.014   2       5034834   2       4578459
        GCF_003697165.2   GCF_002950215.1   73.4098   96.7220   96.8348   73.252   77.978   2       5034834   3       4938295
        GCF_002949675.1   GCF_002950215.1   77.0516   97.0923   97.6772   82.972   74.934   2       4578459   3       4938295
