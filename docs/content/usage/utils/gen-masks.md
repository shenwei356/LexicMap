---
title: gen-masks
weight: 40
---

## Usage

```plain
$ lexicmap utils gen-masks  -h
Generate masks from the top N largest genomes

How:

    |ATTATAACGCCACGGGGAGCCGCGGGGTTTC One k-bp mask
    |--------========_______________
    |
    |-------- Prefixes for covering all possible P-bp DNA.
    |         The length is the largest number for 4^P <= #masks
    |
    |--------======== Extend prefixes, chosen from the most frequent extended prefixes
    |                 of which the prefix-derived k-mers do not overlap masked k-mers.
    |
    |                _______________ Randomly generated bases

Usage:
  lexicmap utils gen-masks [flags] [-k <k>] [-n <masks>] [-n <top-n>] [-D <seeds.tsv.gz>] [-o masks.txt] { -I <seqs dir> | -X <file list>}

Flags:
  -G, --big-genomes string        ► Out file of skipped files with genomes >= -G/--max-genome
  -r, --file-regexp string        ► Regular expression for matching sequence files in -I/--in-dir,
                                  case ignored. (default "\\.(f[aq](st[aq])?|fna)(.gz)?$")
  -h, --help                      help for gen-masks
  -I, --in-dir string             ► Directory containing FASTA/Q files. Directory symlinks are followed.
  -k, --kmer int                  ► Maximum k-mer size. K needs to be <= 32. (default 31)
  -m, --masks int                 ► Number of masks. (default 40000)
  -g, --max-genome int            ► Maximum genome size. Extremely large genomes (non-isolate
                                  assemblies) will be skipped. (default 15000000)
  -o, --out-file string           ► Out file of generated masks. The ".gz" suffix is not recommended.
                                  ("-" for stdout). (default "-")
  -P, --prefix-ext int            ► Extension length of prefixes, higher values -> smaller maximum
                                  seed distances. (default 8)
  -s, --rand-seed int             ► Rand seed for generating random masks. (default 1)
  -N, --ref-name-regexp string    ► Regular expression (must contains "(" and ")") for extracting the
                                  reference name from the filename. (default
                                  "(?i)(.+)\\.(f[aq](st[aq])?|fna)(.gz)?$")
  -D, --seed-pos string           ► Out file of seed postions and distances, supports and recommends a
                                  ".gz" suffix.
  -B, --seq-name-filter strings   ► List of regular expressions for filtering out sequences by
                                  header/name, case ignored.
  -S, --skip-file-check           ► Skip input file checking when given files or a file list.
  -n, --top-n int                 ► Select the top N largest genomes for generating masks. (default 20)

Global Flags:
  -X, --infile-list string   ► File of input files list (one file per line). If given, they are
                             appended to files from CLI arguments.
      --log string           ► Log file.
      --quiet                ► Do not print any verbose information. But you can write them to a file
                             with --log.
  -j, --threads int          ► Number of CPUs cores to use. By default, it uses all available cores.
                             (default 16)

```

## Examples

This command only generate masks, the flags are the same as these in `lexicmap index`.

    $ lexicmap utils gen-masks -I refs/ -k 31 --top-n 3 --masks 40000 -o masks.txt

    $ head -n 10 masks.txt
    AAAAAAACAGGGGCACGACGGTTTATGAGCT
    AAAAAACCGCCGGCCAGGGAAAGCTCATGCC
    AAAAAAGCCCCGCGAGACAAAGCGAAAATAG
    AAAAAATACTGTCTGGCTAGCGTTTCGGGGA
    AAAAACAACCTTGCTTCGTCGTTCCAGTGCC
    AAAAACCCATCACTTCTGTGGCCGGTGTATC
    AAAAACGCAACTGAAGACGCCGACAAACTCT
    AAAAACTCCATGAAACCCACAATATAAATAG
    AAAAAGAAACGCTCAATAATTTGTATGGTTT
    AAAAAGCCGCGCGGTTCAGATGGACGGCATT

The generated masks could be used in `lexicmap index` via the flag `-M/--mask-file`.
