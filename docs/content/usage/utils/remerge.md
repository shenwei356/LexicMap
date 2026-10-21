---
title: lexicmap utils remerge
linkTitle: remerge
weight: 60
---

```plain
$ lexicmap utils remerge -h

Rerun the merging step for an unfinished index

When to use this command?

- Only one thread is used for merging indexes, which happens when there are
  a lot (>200 batches) of batches ($inpu_files / --batch-size) and the value
  of --max-open-files is not big enough. E.g.,

  22:54:24.420 [INFO] merging 297 indexes...
  22:54:24.455 [INFO]   [round 1]
  22:54:24.455 [INFO]     batch 1/1, merging 297 indexes to xxx.lmi.tmp/r1_b1 with 1 threads...

  ► Then you can run this command with a bigger --max-open-files (e.g., 4096) and
  -J/--seed-data-threads (e.g., 12. 12 needs to be <= (4096-8)/(297+2)=13.7).
  And you need to set a bigger 'ulimit -n' if the value of --max-open-files is bigger than 1024.

- The Slurm/PBS job time limit is almost reached and the merging step won't be finished before that.

- Disk quota is reached in the merging step.

Large batch counts are merged in multiple rounds within --max-open-files,
including 8 files reserved for standard streams, logging, and runtime use.
Original batches are retained until the final index is published, so failed
merges can be retried. Genome data use hard links where supported; otherwise,
retaining the batches requires extra disk space for copies.
Log paths must not alias files in the original temporary batches.
Neither the output index directory nor its temporary directory (<index>.tmp)
may be the current directory or an alias of it.

Use -2/--seed-index2 to write adaptive two-level seed indexes during the final
merge, as in "lexicmap index". Set --seed-index2-threshold to the threshold used
for the original build. Prefix lengths are read from the existing batch indexes.

Usage:
  lexicmap utils remerge [flags] [flags] -d <index path>

Flags:
  -h, --help                           help for remerge
  -d, --index string                   ► Index directory created by "lexicmap index".
      --max-open-files int             ► Maximum open files for merging indexes, including 8 reserved
                                       files (minimum: 12 for multiple batches, or 13 with
                                       --seed-index2). Large batch counts are merged in multiple rounds.
                                       Set "ulimit -n" at least this high. (default 1024)
  -J, --seed-data-threads int          ► Number of threads for writing seed data and merging seed
                                       chunks from all batches, in range [1, -c/--chunks]. Merging
                                       threads are limited by --max-open-files. (default 8)
  -2, --seed-index2                    ► Write adaptive two-level seed indexes (.idx and .idx15)
                                       alongside the final seeds, as in "lexicmap index --seed-index2".
                                       Prefix lengths are read from the existing batch indexes.
      --seed-index2-threshold string   ► Minimum seed block size for creating a secondary index with
                                       --seed-index2. Use the threshold from the original build.
                                       Minimum: 4K. The default 8K is recommended. (default "8K")

Global Flags:
      --force-gzip-output    ► Gzip-compress result output regardless of the filename, including stdout.
  -X, --infile-list string   ► File of input file list (one file per line). If given, they are
                             appended to files from CLI arguments.
      --log string           ► Log file.
      --quiet                ► Do not print any verbose information. But you can write them to a file
                             with --log.
  -j, --threads int          ► Number of CPU cores to use. By default, it uses all available cores.
                             (default 16)
```
