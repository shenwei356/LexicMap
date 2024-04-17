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

1. Adding the flag `--save-seed-pos` in index building.

        $ lexicmap index -I refs/ -O demo.lmi --top-n 3 --save-seed-pos --force

2. Listing seed position of one genome.

        $ lexicmap utils seed-pos -d demo.lmi/ -n GCF_000006945.2 -o seed_distance.tsv

        $ head -n 10 seed_distance.tsv | csvtk pretty -t
        ref               pos   strand   distance
        ---------------   ---   ------   --------
        GCF_000006945.2   113   +        112
        GCF_000006945.2   291   -        178
        GCF_000006945.2   297   -        6
        GCF_000006945.2   299   -        2
        GCF_000006945.2   322   +        23
        GCF_000006945.2   340   +        18
        GCF_000006945.2   389   +        49
        GCF_000006945.2   533   -        144
        GCF_000006945.2   598   -        65

    Check the biggest seed distances.

        $ csvtk freq -t -f distance seed_distance.tsv \
            | csvtk sort -t -k distance:nr \
            | head -n 20 \
            | csvtk pretty -t

        distance   frequency
        --------   ---------
        1586       1
        1434       1
        1418       1
        1398       1
        1282       2
        1268       1
        1261       1
        1233       1
        1175       1
        1166       1
        1158       1
        1129       1
        1106       1
        1079       2
        1074       1
        1066       1
        1060       1
        1041       1
        1034       2

    Plot the histogram of distances between seeds.

        $ lexicmap utils seed-pos -d demo.lmi/ -n GCF_000006945.2 -o seed_distance.tsv  --plot-dir seed_distance

    <img src="/LexicMap/GCF_000006945.2.png" alt="" width="600"/>

3. Listing seed position of all genomes.

        $ lexicmap utils seed-pos -d demo.lmi/ --all-refs -o seed-pos.tsv.gz

    Show the number of seed positions in each genome.

        $ csvtk freq -t -f ref -nr seed-pos.tsv.gz | csvtk pretty -t
        ref               frequency
        ---------------   ---------
        GCF_002950215.1   43149
        GCF_002949675.1   42888
        GCF_001457655.1   42444
        GCF_006742205.1   42050
        GCF_900638025.1   41939
        GCF_001027105.1   41925
        GCF_000392875.1   41487
        GCF_009759685.1   41058
        GCF_000148585.2   41029
        GCF_000017205.1   41012
        GCF_000006945.2   41002
        GCF_003697165.2   40785
        GCF_000742135.1   40632
        GCF_001096185.1   40234
        GCF_001544255.1   40188

    Plot the histograms of distances between seeds for all genomes.

        $ lexicmap utils seed-pos -d demo.lmi/ --all-refs -o seed-pos.tsv.gz --plot-dir seed_distance


The output (TSV format) is formatted with [csvtk pretty](https://github.com/shenwei356/csvtk).
