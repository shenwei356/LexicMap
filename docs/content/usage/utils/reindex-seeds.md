---
title: reindex-seeds
weight: 50
---

## Usage

```plain
$ lexicmap utils reindex-seeds -h
Recreate indexes of k-mer-value (seeds) data

Usage:
  lexicmap utils reindex-seeds [flags]

Flags:
  -h, --help             help for reindex-seeds
  -d, --index string     ► Index directory created by "lexicmap index".
  -p, --partitions int   ► Number of partitions for re-indexing seeds (k-mer-value data) files.
                         (default 512)

Global Flags:
  -X, --infile-list string   ► File of input files list (one file per line). If given, they are
                             appended to files from CLI arguments.
      --log string           ► Log file.
      --quiet                ► Do not print any verbose information. But you can write them to a file
                             with --log.
  -j, --threads int          ► Number of CPUs cores to use. By default, it uses all available cores.
                             (default 16
```

## Examples
