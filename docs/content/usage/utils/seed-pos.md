---
title: seed-pos
weight: 30
---

## Usage

```plain
$ lexicmap utils seed-pos -h
Extract seed positions via reference names

Attentions:
  0. This command requires the index to be created with the flag --save-seed-pos in lexicmap index.
  1. Seed/K-mer positions (column pos) are 1-based.
     For reference genomes with multiple sequences, the sequences were
     concatenated to a single sequence with intervals of N's.

Usage:
  lexicmap utils seed-pos [flags]

Flags:
  -a, --all-refs           ► Output for all reference genomes. This would take a long time for an
                           index with a lot of genomes.
  -b, --bins int           ► Number of bins in histograms. (default 100)
      --color-index int    ► Color index (1-7). (default 1)
      --force              ► Overwrite existing output directory.
      --height float       ► Histogram height (unit: inch). (default 4)
  -h, --help               help for seed-pos
  -d, --index string       ► Index directory created by "lexicmap index".
  -o, --out-file string    ► Out file, supports and recommends a ".gz" suffix ("-" for stdout).
                           (default "-")
  -O, --plot-dir string    ► Output directory for histograms of seed distances.
      --plot-ext string    ► Histogram plot file extention. (default ".png")
  -t, --plot-title         ► Plot genome ID as the title.
  -n, --ref-name strings   ► Reference name(s).
      --width float        ► Histogram width (unit: inch). (default 6)

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


