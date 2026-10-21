---
title: lexicmap
linkTitle: lexicmap
weight: 0
---

```plain
$ lexicmap -h

   LexicMap: efficient sequence/genome alignment against millions of prokaryotic genomes

    Version: v0.10.0
  Documents: https://bioinf.shenwei.me/LexicMap
Source code: https://github.com/shenwei356/LexicMap
Please cite: https://doi.org/10.1038/s41587-025-02812-8

Usage:
  lexicmap [command]

Available Commands:
  autocompletion Generate shell autocompletion scripts
  genome         Commands for genome searching, prefilter, and comparison
  index          Generate an index from FASTA/Q sequences
  search         Search sequences against an index
  utils          Some utilities
  version        Print version information and check for update

Flags:
      --force-gzip-output    ► Gzip-compress result output regardless of the filename, including stdout.
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

Output and log files must differ from input files, including symbolic links and
hard links. This also applies to file lists, auxiliary inputs, and index data.
Log files must be outside output directories and temporary merge directories
that the command may replace, including when `--force` is used.
Commands exit with an error if writing, flushing, or closing their result output
fails.

Files from `-X/--infile-list` are appended to positional inputs, including an
explicit `-` for stdin. With no positional inputs, a non-empty file list supplies
the inputs without adding stdin.

Stdin (`-`) may be used only once across positional files, file-list entries, and
auxiliary inputs. For example, `-X -` cannot be combined with a sequence input
from stdin. In indexed `genome compare`, a genome ID named `-` remains an ID.

Result output is gzip-compressed when the output filename ends in `.gz`.
Use the global flag `--force-gzip-output` to enable gzip compression for other
filenames or stdout, for example:

```bash
lexicmap search -d db.lmi query.fasta --force-gzip-output > results.tsv.gz
```
