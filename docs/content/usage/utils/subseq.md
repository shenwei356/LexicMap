---
title: subseq
weight: 20
---

## Usage

```plain
$ lexicmap utils subseq -h
Exextract subsequence via reference name, sequence ID, position and strand

Attention:
  1. The option -s/--seq-id is optional.
     1) If given, the positions are that in the original sequence.
         2) If not given, the positions are that in the concatenated sequence.

Usage:
  lexicmap utils subseq [flags]

Flags:
  -h, --help              help for subseq
  -d, --index string      ► Index directory created by "lexicmap index".
  -w, --line-width int    ► Line width of sequence (0 for no wrap). (default 60)
  -o, --out-file string   ► Out file, supports and recommends a ".gz" suffix ("-" for stdout).
                          (default "-")
  -n, --ref-name string   ► Reference name.
  -r, --region string     ► Region of the subsequence (1-based).
  -R, --revcom            ► Extract subsequence on the negative strand.
  -s, --seq-id string     ► Sequence ID. If the value is empty, the positions in the region are
                          treated as that in the concatenated sequence.

Global Flags:
  -X, --infile-list string   ► File of input files list (one file per line). If given, they are
                             appended to files from CLI arguments.
      --log string           ► Log file.
      --quiet                ► Do not print any verbose information. But you can write them to a file
                             with --log.
  -j, --threads int          ► Number of CPUs cores to use. By default, it uses all available cores.
                             (default 16)
```

## Examples

