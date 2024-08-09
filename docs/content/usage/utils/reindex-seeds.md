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
      --partitions int   ► Number of partitions for re-indexing seeds (k-mer-value data) files. The
                         value needs to be the power of 4. (default 1024)

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


    $ lexicmap utils reindex-seeds -d demo.lmi/ --partitions 1024
    10:20:29.150 [INFO] recreating seed indexes with 1024 partitions for: demo.lmi/
    processed files:  16 / 16 [======================================] ETA: 0s. done
    10:20:29.166 [INFO] update index information file: demo.lmi/info.toml
    10:20:29.166 [INFO]   finished updating the index information file: demo.lmi/info.toml
    10:20:29.166 [INFO]
    10:20:29.166 [INFO] elapsed time: 15.981266ms
    10:20:29.166 [INFO]
