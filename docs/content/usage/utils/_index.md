---
title: lexicmap utils
linkTitle: utils
weight: 60
geekdocCollapseSection: true
---

Subcommands:

- [2blast](2blast/)
- [2sam](2sam/)
- [merge-search-results](merge-search-results/)
- [masks](masks/)
- [kmers](kmers/)
- [genomes](genomes/)
- [genome-details](genomes-details/)
- [genome-seqs](genome-seqs)
- [subseq](subseq/)
- [seed-pos](seed-pos/)
- [reindex-seeds](reindex-seeds/)
- [reindex-seeds2](reindex-seeds2/)
- [remerge](remerge/)
- [edit-genome-ids](edit-genome-ids/)

Usage

```plain
$ lexicmap utils
Some utilities

Usage:
  lexicmap utils [command] 

Available Commands:
  2blast               Convert the default search output to blast-style format
  2sam                 Convert the default search output to SAM format
  edit-genome-ids      Edit genome IDs in the index via a regular expression
  genome-details       Extract or view genome details in the index
  genome-seqs          Extract all sequences of a given genome
  genomes              View genome IDs in the index
  kmers                View k-mers captured by the masks
  masks                View masks of the index or generate new masks randomly
  merge-search-results Merge a query's search results from multiple indexes
  reindex-seeds        Recreate primary indexes of seeds data
  reindex-seeds2       Recreate adaptive two-level indexes of seeds data
  remerge              Rerun the merging step for an unfinished index
  seed-pos             Extract and plot seed positions via reference name(s)
  subseq               Extract subsequence via 1) reference name, sequence ID, position and strand, or 2) search result

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
