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
     The positions can be used to extract subsequence with 'lexicmap utils subseq'.
  2. A distance between seeds (column distance) with a value of "-1" means it's the first seed
     in that sequence, and the distance can't be computed currently.
  3. All degenerate bases in reference genomes were converted to the lexicographic first bases.
     E.g., N was converted to A. Therefore, consecutive A's in output might be N's in the genomes.

Extra columns:
  Using -v/--verbose will output more columns:
     pre_pos,  the position of the previous seed.
     len_aaa,  length of consecutive A's.
     seq,      sequence between the previous and current seed.

Usage:
  lexicmap utils seed-pos [flags]

Flags:
  -a, --all-refs             ► Output for all reference genomes. This would take a long time for an
                             index with a lot of genomes.
  -b, --bins int             ► Number of bins in histograms. (default 100)
      --color-index int      ► Color index (1-7). (default 1)
      --force                ► Overwrite existing output directory.
      --height float         ► Histogram height (unit: inch). (default 4)
  -h, --help                 help for seed-pos
  -d, --index string         ► Index directory created by "lexicmap index".
      --max-open-files int   ► Maximum opened files, used for extracting sequences. (default 512)
  -D, --min-dist int         ► Only output records with seed distance >= this value.
  -o, --out-file string      ► Out file, supports and recommends a ".gz" suffix ("-" for stdout).
                             (default "-")
  -O, --plot-dir string      ► Output directory for histograms of seed distances.
      --plot-ext string      ► Histogram plot file extention. (default ".png")
  -t, --plot-title           ► Plot genome ID as the title.
  -n, --ref-name strings     ► Reference name(s).
  -v, --verbose              ► Show more columns including position of the previous seed and sequence
                             between the two seeds. Warning: it's slow to extract the sequences,
                             recommend set -D 1000 or higher values to filter results
      --width float          ► Histogram width (unit: inch). (default 6)

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

1. Adding the flag `--save-seed-pos` in index building.

        $ lexicmap index -I refs/ -O demo.lmi --save-seed-pos --force

2. Listing seed position of one genome.

        $ lexicmap utils seed-pos -d demo.lmi/ -n GCF_000017205.1 -o seed_distance.tsv

        $ head -n 10 seed_distance.tsv | csvtk pretty -t
        ref               pos   strand   distance
        ---------------   ---   ------   --------
        GCF_000017205.1   2     +        1
        GCF_000017205.1   41    -        39
        GCF_000017205.1   45    +        4
        GCF_000017205.1   74    -        29
        GCF_000017205.1   85    -        11
        GCF_000017205.1   119   -        34
        GCF_000017205.1   130   -        11
        GCF_000017205.1   185   +        55
        GCF_000017205.1   269   -        84

    Or only list records with seed distance longer than a threshold.

        $ lexicmap utils seed-pos -d demo.lmi/ -n GCF_000017205.1 -D 850 | csvtk pretty -t | head -n 3
        15:34:58.669 [INFO] seed positions of 1 genomes(s) saved to -
        ref               pos       strand   distance
        ---------------   -------   ------   --------
        GCF_000017205.1   30713     -        850

    Check the biggest seed distances.

        $ csvtk freq -t -f distance seed_distance.tsv \
            | csvtk sort -t -k distance:nr \
            | head -n 10 \
            | csvtk pretty -t

        distance   frequency
        --------   ---------
        899        2
        898        4
        897        1
        896        5
        895        4
        894        3
        893        3
        892        4
        891        2


    Plot the histogram of distances between seeds.

        $ lexicmap utils seed-pos -d demo.lmi/ -n GCF_000017205.1 -o seed_distance.tsv  --plot-dir seed_distance

    In the plot below, there's a peak at 200 bp, because LexicMap fills sketching deserts with extra k-mers (seeds) of which their distance is 200 bp by default.

    <img src="/LexicMap/GCF_000017205.1.png" alt="" width="600"/>

2. More columns including sequences between two seeds.

        $ lexicmap utils seed-pos -d demo.lmi/  -n GCF_000017205.1 -v \
            | head -n4 | csvtk pretty -t -W 50 --clip
        ref               pos   strand   distance   pre_pos   len_aaa   seq
        ---------------   ---   ------   --------   -------   -------   ---------------------------------------
        GCF_000017205.1   2     +        1          0         0         T
        GCF_000017205.1   41    -        39         1         5         TAAAGAGACCGGCGATTCTAGTGAAATCGAACGGGCAGG
        GCF_000017205.1   45    +        4          40        1         TCAA

    Or only list records with seed distance longer than a threshold.

        $ lexicmap utils seed-pos -d demo.lmi/ -n GCF_000017205.1 -v -D 890 \
            | head -n 2 \
            | csvtk pretty -t -W 50 --clip
        ref               pos      strand   distance   pre_pos   len_aaa   seq
        ---------------   ------   ------   --------   -------   -------   --------------------------------------------------
        GCF_000017205.1   152018   -        892        151125    21        CGCGGCCCAGCCATGCCTACTGGGACCTCTCGCCGGGGATCGATTTC...


3. Listing seed position of all genomes.

        $ lexicmap utils seed-pos -d demo.lmi/ --all-refs -o seed-pos.tsv.gz

    Show the number of seed positions in each genome.
    Frequencies larger than 40000 (the number of masks) means some k-mers can be foud in more than one positions in a genome.

        $ csvtk freq -t -f ref -nr seed-pos.tsv.gz | csvtk pretty -t
        ref               frequency
        ---------------   ---------
        GCF_000017205.1   45737
        GCF_002950215.1   43617
        GCF_002949675.1   43469
        GCF_001457655.1   42112
        GCF_006742205.1   42102
        GCF_900638025.1   42008
        GCF_001027105.1   41855
        GCF_000742135.1   41419
        GCF_000392875.1   41391
        GCF_003697165.2   41194
        GCF_009759685.1   41137
        GCF_000006945.2   41114
        GCF_000148585.2   41075
        GCF_001096185.1   40233
        GCF_001544255.1   40165

    Plot the histograms of distances between seeds for all genomes.

        $ lexicmap utils seed-pos -d demo.lmi/ --all-refs -o seed-pos.tsv.gz --plot-dir seed_distance --force
        processed files:  15 / 15 [======================================] ETA: 0s. done
        11:48:31.346 [INFO] seed positions of 15 genomes(s) saved to seed-pos.tsv.gz
        11:48:31.346 [INFO] histograms of 15 genomes(s) saved to seed_distance

        $ ls seed_distance/
        GCF_000006945.2.png  GCF_000392875.1.png  GCF_001096185.1.png  GCF_002949675.1.png  GCF_006742205.1.png
        GCF_000017205.1.png  GCF_000742135.1.png  GCF_001457655.1.png  GCF_002950215.1.png  GCF_009759685.1.png
        GCF_000148585.2.png  GCF_001027105.1.png  GCF_001544255.1.png  GCF_003697165.2.png  GCF_900638025.1.png

    Some genomes, e.g., GCF_000392875.1, might have a few big seed distances around gaps (N's). In LexicMap, the N's are converted to A's.


    ```text
    $ lexicmap utils seed-pos -d demo.lmi/ -n GCF_000392875.1 -v -D 1000 | csvtk pretty -t -W 80
    ref               pos       strand   distance   pre_pos   len_aaa   seq
    ---------------   -------   ------   --------   -------   -------   --------------------------------------------------------------------------------
    GCF_000392875.1   503031    +        1161       501869    1116      ATGAGCCAACAGTAGAAGGTGAAAAAGTAGAAATCGGTGGTAAAGTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGTGTCAT
    GCF_000392875.1   2640078   +        1344       2638733   1150      CAACTCCTGTACTAGTATTTAAGTGTCCATTATTCCCCCCATTTTTTTGCTCCTTTTTATTTTCCCCACTATTTTTCAAT
                                                                        GTTAATTGCTTCACTGCCGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                                                                        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGCTTGTTTCAGTGCTTCGCTGTAGGCTTTCCAGCTGCT
                                                                        TGCGGTGTAATCTTTTTCTTGGTGTTCTTTTTGTTCCTGAATTAATTTTTCTAACGCTTCTTTC
```

The output (TSV format) is formatted with [csvtk pretty](https://github.com/shenwei356/csvtk).
