---
title: utils
weight: 40
geekdocCollapseSection: true
---

```plain
$ lexicmap utils
Some utilities

Usage:
  lexicmap utils [command]

Available Commands:
  2blast          Convert the default search output to blast-style format
  edit-genome-ids Edit genome IDs in the index via a regular expression
  genomes         View genome IDs in the index
  kmers           View k-mers captured by the masks
  masks           View masks of the index or generate new masks randomly
  reindex-seeds   Recreate indexes of k-mer-value (seeds) data
  remerge         Rerun the merging step for an unfinished index
  seed-pos        Extract and plot seed positions via reference name(s)
  subseq          Extract subsequence via reference name, sequence ID, position and strand

Flags:
  -h, --help   help for utils

Global Flags:
  -X, --infile-list string   ► File of input file list (one file per line). If given, they are
                             appended to files from CLI arguments.
      --log string           ► Log file.
      --quiet                ► Do not print any verbose information. But you can write them to a file
                             with --log.
  -j, --threads int          ► Number of CPU cores to use. By default, it uses all available cores.
                             (default 16)
```


Subcommands:

- [2blast](2blast/)
- [masks](masks/)
- [kmers](kmers/)
- [genomes](genomes/)
- [subseq](subseq/)
- [seed-pos](seed-pos/)
- [reindex-seeds](reindex-seeds/)
- [remerge](remerge/)
- [edit-genome-ids](edit-genome-ids/)
