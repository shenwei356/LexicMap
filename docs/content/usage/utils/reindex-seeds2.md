---
title: lexicmap utils reindex-seeds
linkTitle: reindex-seeds
weight: 55
---

## Usage

```plain
$ lexicmap utils reindex-seeds -h
Recreate adaptive two-level indexes of seeds data.

The seeds files are not changed. The primary prefix length is the mask prefix
plus log4(partitions), and the secondary prefix adds 2 bases. With the standard
7-bp mask prefix and 4096 partitions, these lengths are 13 and 15 bp. Large
primary-prefix blocks receive a secondary index, while smaller blocks retain
the original linear-scan path. The command rewrites each primary .idx file and
creates its matching .idx15 file. Different seeds files are processed in
parallel using -j/--threads.

This index primarily improves seed-matching performance for batch queries with
"lexicmap search". The end-to-end speedup is most noticeable when
-n/--top-n-genomes limits the number of candidates passed to downstream
alignment; searches dominated by alignment, especially single-query searches,
may see only a small overall gain.

Usage:
  lexicmap utils reindex-seeds2 [flags] 

Flags:
  -h, --help               help for reindex-seeds2
  -d, --index string       ► Index directory created by "lexicmap index".
      --partitions int     ► Number of anchor partitions. The value needs to be a power of 4 and at
                           least 4. The primary prefix is mask-prefix + log4(partitions), and the
                           secondary prefix adds 2 bases; the 5-byte checkpoint suffix also requires K -
                           primary-prefix <= 20. The default 4096 produces 13/15-bp prefixes with the
                           standard 7-bp mask prefix. (default 4096)
      --threshold string   ► Minimum size of a primary-prefix seeds block for creating its secondary
                           index. The minimum is 4K because seed readers use 4 KiB buffers; the default
                           is 8K. (default "8K")

Global Flags:
  -X, --infile-list string   ► File of input file list (one file per line). If given, they are
                             appended to files from CLI arguments.
      --log string           ► Log file.
      --quiet                ► Do not print any verbose information. But you can write them to a file
                             with --log.
  -j, --threads int          ► Number of CPU cores to use. By default, it uses all available cores.
                             (default 16)
```

## Examples


    $ lexicmap utils reindex-seeds2 -d demo.lmi/
    16:28:30.100 [INFO] LexicMap v0.10.0
    16:28:30.100 [INFO]   https://github.com/shenwei356/LexicMap
    16:28:30.100 [INFO] 
    16:28:30.100 [INFO]  CWD: /home/shenwei/go/src/github.com/shenwei356/LexicMap/demo
    16:28:30.100 [INFO]  CMD: lexicmap utils reindex-seeds2 -d demo.lmi/
    16:28:30.100 [INFO] DATE: 2026-10-09
    16:28:30.100 [INFO] 
    16:28:30.100 [INFO] recreating adaptive seed indexes with 4096 partitions and a 8.0 KiB block threshold for: demo.lmi/
    processed files:  16 / 16 [======================================] ETA: 0s. done
    16:28:30.286 [INFO] primary-prefix blocks: 81920000 total, 2107085 non-empty, 0 indexed (0.00% of non-empty)
    16:28:30.286 [INFO] idx15 size: 768 B
    16:28:30.286 [INFO] 
    16:28:30.286 [INFO] elapsed time: 185.895153ms
    16:28:30.286 [INFO] 
