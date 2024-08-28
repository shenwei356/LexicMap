---
title: genomes
weight: 20
---

## Usage

```plain
$ lexicmap utils genomes -h
View genome IDs in the index

Usage:
  lexicmap utils genomes [flags]

Flags:
  -h, --help              help for genomes
  -d, --index string      ► Index directory created by "lexicmap index".
  -o, --out-file string   ► Out file, supports the ".gz" suffix ("-" for stdout). (default "-")

Global Flags:
  -X, --infile-list string   ► File of input file list (one file per line). If given, they are
                             appended to files from CLI arguments.
      --log string           ► Log file.
      --quiet                ► Do not print any verbose information. But you can write them to a file
                             with --log.
  -j, --threads int          ► Number of CPU cores to use. By default, it uses all available cores.
                             (default 8)
```

## Examples


```
$ lexicmap utils genomes -d demo.lmi/
GCF_000148585.2
GCF_001457655.1
GCF_900638025.1
GCF_001096185.1
GCF_006742205.1
GCF_001544255.1
GCF_000392875.1
GCF_001027105.1
GCF_009759685.1
GCF_002949675.1
GCF_002950215.1
GCF_000006945.2
GCF_003697165.2
GCF_000742135.1
GCF_000017205.1
```
