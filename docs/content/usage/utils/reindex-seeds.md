---
title: reindex-seeds
weight: 50
---

## Usage

```plain
$ lexicmap utils reindex-seeds -h
Recreate indexes of k-mer-value (seeds) data

Experimental feature:

  The flag --plain-format can save indexes of seed data in plain format,
  so marker/anchor k-mers and their offsets in the seed file can be accessed with mmap.
  This reduces the startup time (1-6 seconds).
  
  This flag is usually used along with a bigger value of --partition, such as 65536 (4^8),
  to reduce the seed matching time, by omitting the reading of some unwanted seed data.
  However, larger values of --partition would result in bigger .idx files. 
  E.g., the default 4096 requires < 1 GB, while 655536 needs 20 GB.

  Attention:
    This feature only benefits searching a small number of queries against big databases.
  For a lot of queries, the speed would be slower, and the memory would be too high,
  as more and more seed index data will be mapped into memory.

Usage:
  lexicmap utils reindex-seeds [flags] 

Flags:
  -h, --help             help for reindex-seeds
  -d, --index string     ► Index directory created by "lexicmap index".
      --partitions int   ► Number of partitions for re-indexing seeds (k-mer-value data) files. The
                         value needs to be the power of 4 (default 4096)
      --plain-format     ► Save indexes of seed data in plain format for faster quering with mmap, at
                         the cost of bigger index size.

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
