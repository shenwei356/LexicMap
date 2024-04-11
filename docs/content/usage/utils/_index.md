---
title: utils
weight: 40
---

```plain
$ lexicmap utils
Some utilities

Usage:
  lexicmap utils [command]

Available Commands:
  gen-masks     Generate masks from the top N largest genomes
  kmers         View k-mers captured by the masks
  masks         View masks of the index or generate new masks randomly
  reindex-seeds Recreate indexes of k-mer-value (seeds) data
  seed-pos      Extract seed positions via reference names
  subseq        Extract subsequence via reference name, sequence ID, position and strand

Flags:
  -h, --help   help for utils

Global Flags:
  -X, --infile-list string   ► File of input files list (one file per line). If given, they are
                             appended to files from CLI arguments.
      --log string           ► Log file.
      --quiet                ► Do not print any verbose information. But you can write them to a file
                             with --log.
  -j, --threads int          ► Number of CPUs cores to use. By default, it uses all available cores.
                             (default 16)
```
