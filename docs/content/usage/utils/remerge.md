---
title: remerge
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
  -J/--seed-data-threads (e.g., 12. 12 needs be <= 4096/(297+2)=13.7).
  And you need to set a bigger 'ulimit -n' if the value of --max-open-files is bigger than 1024.

- The Slurm/PBS job time limit is almost reached and the merging step won't be finished before that.

- Disk quota is reached in the merging step.

Usage:
  lexicmap utils remerge [flags] [flags] -d <index path>

Flags:
  -h, --help                    help for remerge
  -d, --index string            ► Index directory created by "lexicmap index".
      --max-open-files int      ► Maximum opened files, used in merging indexes. If there are >100
                                batches, please increase this value and set a bigger "ulimit -n" in
                                shell. (default 768)
  -J, --seed-data-threads int   ► Number of threads for writing seed data and merging seed chunks from
                                all batches, the value should be in range of [1, -c/--chunks]. If there
                                are >100 batches, please also increase the value of --max-open-files and
                                set a bigger "ulimit -n" in shell. (default 8)

Global Flags:
  -X, --infile-list string   ► File of input file list (one file per line). If given, they are
                             appended to files from CLI arguments.
      --log string           ► Log file.
      --quiet                ► Do not print any verbose information. But you can write them to a file
                             with --log.
  -j, --threads int          ► Number of CPU cores to use. By default, it uses all available cores.
                             (default 16)
```
