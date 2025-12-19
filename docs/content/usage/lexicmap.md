---
title: lexicmap
weight: 0
---

```plain
$ lexicmap -h

   LexicMap: efficient sequence alignment against millions of prokaryotic genomes

    Version: v0.8.1
  Documents: https://bioinf.shenwei.me/LexicMap
Source code: https://github.com/shenwei356/LexicMap
Please cite: https://doi.org/10.1038/s41587-025-02812-8 Nature Biotechnology (2025)

Usage:
  lexicmap [command] 

Available Commands:
  autocompletion Generate shell autocompletion scripts
  index          Generate an index from FASTA/Q sequences
  search         Search sequences against an index
  utils          Some utilities
  version        Print version information and check for update

Flags:
  -h, --help                 help for lexicmap
  -X, --infile-list string   ► File of input file list (one file per line). If given, they are
                             appended to files from CLI arguments.
      --log string           ► Log file.
      --quiet                ► Do not print any verbose information. But you can write them to a file
                             with --log.
  -j, --threads int          ► Number of CPU cores to use. By default, it uses all available cores.
                             (default 16)

Use "lexicmap [command] --help" for more information about a command.

```
