---
title: lexicmap genome
linkTitle: genome
weight: 50
geekdocCollapseSection: true
---

Subcommands:

- [compare](compare/)
- [pair](pair/)
- [search](search/)

Usage

```plain
$ lexicmap genome
Commands for genome searching, prefilter, and comparison

Usage:
  lexicmap genome [command] 

Available Commands:
  compare     Compare genome pairs and compute ANI and AF
  pair        Find similar genome pairs in the index
  search      Search genomes against an index, with ANI and AF computed

Flags:
  -h, --help   help for genome

Global Flags:
  -X, --infile-list string   ► File of input file list (one file per line). If given, they are
                             appended to files from CLI arguments.
      --log string           ► Log file.
      --quiet                ► Do not print any verbose information. But you can write them to a file
                             with --log.
  -j, --threads int          ► Number of CPU cores to use. By default, it uses all available cores.
                             (default 16)

Use "lexicmap genome [command] --help" for more information about a command.
```
