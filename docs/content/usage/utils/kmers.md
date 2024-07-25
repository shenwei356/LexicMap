---
title: kmers
weight: 10
---

```plain
$ lexicmap utils kmers -h
View k-mers captured by the masks

Attention:
  1. Mask index (column mask) is 1-based.
  2. Prefix means the length of shared prefix between a k-mer and the mask.
  3. K-mer positions (column pos) are 1-based.
     For reference genomes with multiple sequences, the sequences were
     concatenated to a single sequence with intervals of N's.
  4. Reversed means if the k-mer is reversed for suffix matching.

Usage:
  lexicmap utils kmers [flags] -d <index path> [-m <mask index>] [-o out.tsv.gz]

Flags:
  -h, --help              help for kmers
  -d, --index string      ► Index directory created by "lexicmap index".
  -m, --mask int          ► View k-mers captured by Xth mask. (0 for all) (default 1)
  -o, --out-file string   ► Out file, supports and recommends a ".gz" suffix ("-" for stdout).
                          (default "-")

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

1. The default output is captured k-mers of the first mask.

        $ lexicmap utils kmers --quiet -d demo.lmi/ | head -n 20 | csvtk pretty -t
        mask   kmer                              prefix   number   ref               pos       strand   reversed
        ----   -------------------------------   ------   ------   ---------------   -------   ------   --------
        1      AAAAAAAAAATACTAAGAGGTGACAAAAGAG   4        1        GCF_001544255.1   2142679   +        yes
        1      AAAAAAAAACAAAGCGGACTTGGACATTGTC   4        1        GCF_000006945.2   3170307   +        yes
        1      AAAAAAAAACCAGTAAAAAAAGGGGAGTAGA   4        1        GCF_000392875.1   771896    +        yes
        1      AAAAAAAAACGACTTACCATTAACGTTCAAG   4        1        GCF_003697165.2   803728    +        yes
        1      AAAAAAAAACTAGGGTTAAATGCCTTATGTT   4        1        GCF_009759685.1   442423    +        yes
        1      AAAAAAAAAGAGATGAAAAAGGGTGTATTCG   4        1        GCF_001544255.1   1493451   -        yes
        1      AAAAAAAAATAAAATATCTAACGAGCAAATT   4        1        GCF_001096185.1   2065540   +        yes
        1      AAAAAAAAATACCATAGACTATGCTCTTAGT   4        1        GCF_000392875.1   134079    -        yes
        1      AAAAAAAAATAGAGTTTTTTTTCTGGATAAG   4        1        GCF_000392875.1   795189    +        yes
        1      AAAAAAAAATGTTAACAGAAGGTCCCTACCT   4        1        GCF_002950215.1   2765957   +        yes
        1      AAAAAAAACAAAAGCTATACTGGTCATGTTC   4        1        GCF_000006945.2   3635995   +        yes
        1      AAAAAAAACAAAGATACATTTAGGACGGTTA   4        1        GCF_000006945.2   616481    -        yes
        1      AAAAAAAACAGCCCACCGCCGATTGCGGAAT   4        1        GCF_000742135.1   1208620   +        yes
        1      AAAAAAAACAGGGTGTCGTGCCCTTGTCAGT   4        1        GCF_003697165.2   627153    -        yes
        1      AAAAAAAACAGGGTGTTCTTAGATAAAAGGG   4        1        GCF_000742135.1   1723387   -        yes
        1      AAAAAAAACATATAGTTGTGAAGGCATTGGA   4        1        GCF_001027105.1   2508079   -        yes
        1      AAAAAAAACCAGTAAAAAAAGGGGAGTAGAA   4        1        GCF_000392875.1   771895    +        yes
        1      AAAAAAAACCATATTATGTCCGATCCTCACA   4        1        GCF_000392875.1   1060650   +        yes
        1      AAAAAAAACCCTTCGTCAAGCATTATGGAAT   4        1        GCF_000392875.1   1139573   -        yes


1. Specify the mask.

        $ lexicmap utils kmers --quiet -d demo.lmi/ --mask 12345 | csvtk pretty -t
        mask    kmer                              prefix   number   ref               pos       strand   reversed
        -----   -------------------------------   ------   ------   ---------------   -------   ------   --------
        12345   CATGTTACAAAAGGTGGGTCAGGCAACGTAT   7        1        GCF_001457655.1   335112    -        yes
        12345   CATGTTACCAAGGTTAGTCGTATGGCGCTAC   7        1        GCF_001457655.1   23755     -        yes
        12345   CATGTTACGCGTATTTTAGCGGCTCGCGGAC   7        1        GCF_000006945.2   702224    +        yes
        12345   CATGTTATAACGGCCTATGAATCGGCATTAC   9        1        GCF_009759685.1   2591866   +        no
        12345   CATGTTATACGTTGAAACTGTCTTGTTAATA   9        1        GCF_001096185.1   1142460   +        yes
        12345   CATGTTATACTTTAGATACTTATTTTTAGGA   9        1        GCF_000392875.1   1524553   +        no
        12345   CATGTTATAGAAGGACGTCGACATCTTGTGG   10       1        GCF_000017205.1   3140677   +        no
        12345   CATGTTATAGAATTACATACATTGTAACATG   10       1        GCF_006742205.1   704431    -        no
        12345   CATGTTATAGCACGCTTAATCGCTTGATCCC   13       1        GCF_001027105.1   2655846   +        no
        12345   CATGTTATAGCATCCTTTTACGTGAAAAGGT   12       1        GCF_000742135.1   4136093   +        no
        12345   CATGTTATAGCCAGCAAATGGAAGCATCGCG   11       1        GCF_009759685.1   492828    -        no
        12345   CATGTTATAGCCATTGATGGTAACTTTGATG   11       1        GCF_001096185.1   536843    +        no
        12345   CATGTTATAGCCTGAAAGGTGCTAAACAACT   11       1        GCF_000006945.2   4876155   +        no
        12345   CATGTTATAGCCTTCTCCAAGACCAATCAAA   11       1        GCF_000148585.2   1667015   +        no
        12345   CATGTTATAGCGTAAATCAGCACCGCGCGCC   11       3        GCF_002949675.1   1871326   +        no
        12345   CATGTTATAGCGTAAATCAGCACCGCGCGCC   11       3        GCF_002950215.1   2326544   +        no
        12345   CATGTTATAGCGTAAATCAGCACCGCGCGCC   11       3        GCF_003697165.2   3996124   +        no
        12345   CATGTTATAGCTAACTGCGACTTGTGGCACA   11       1        GCF_900638025.1   991007    -        no
        12345   CATGTTATAGTCGTGAGGTTCTAAAAAAACT   10       1        GCF_001544255.1   1091256   -        no
        12345   CATGTTATAGTTTGTCTTACCGCTACTGAAA   10       1        GCF_002950215.1   1457055   +        yes
        12345   CATGTTATATCCTTCTTGAATACGAGCAATA   9        1        GCF_000392875.1   1963573   +        no
        12345   CATGTTATATGAACCTTCAACCTTATTTGAC   9        1        GCF_001457655.1   1510084   +        no
        12345   CATGTTATCCAGGTATTTCACCAGCGCACGC   8        1        GCF_000006945.2   836525    +        no
        12345   CATGTTATCGAATATTATAACATCGGCTCCC   8        1        GCF_000148585.2   1372855   +        yes
        12345   CATGTTATCGATAAGGCTATATATGACCTTA   8        1        GCF_002950215.1   878140    -        no
        12345   CATGTTATCGCTCAGGGTCTGCGGGTATATC   8        1        GCF_002950215.1   1880029   +        yes
        12345   CATGTTATGCGTATAAAGACGAGTAAAGGTT   8        1        GCF_009759685.1   3827118   +        no
        12345   CATGTTATGCTGGGACATTTAGCACCGCTAC   8        1        GCF_000006945.2   1988134   +        yes

    "reversed" means means if the k-mer is reversed for suffix matching.
    E.g., `CATGTTACAAAAGGTGGGTCAGGCAACGTAT` is reversed, so you need to reverse it before searching in the genome.


        $ seqkit locate -p $(echo CATGTTACAAAAGGTGGGTCAGGCAACGTAT | rev) refs/GCF_001457655.1.fa.gz -M | csvtk pretty -t
        seqID           patternName                       pattern                           strand   start    end
        -------------   -------------------------------   -------------------------------   ------   ------   ------
        NZ_LN831035.1   TATGCAACGGACTGGGTGGAAAACATTGTAC   TATGCAACGGACTGGGTGGAAAACATTGTAC   -        335112   335142


1. For all masks. The result might be very big, therefore, writing to gzip format is recommended.


        $ lexicmap utils kmers -d demo.lmi/ --mask 0 -o kmers.tsv.gz

        $ zcat kmers.tsv.gz | csvtk freq -t -f mask -nr | head -n 10
        mask    frequency
        1       610
        40000   568
        31      435
        20      432
        39997   423
        28      419
        30018   415
        30027   403
        79      396

1. K-mers of a specific mask

        $ lexicmap utils kmers -d demo.lmi/ -m 12345 | head -n 20 | csvtk pretty -t
        mask    kmer                              prefix   number   ref               pos       strand   reversed
        -----   -------------------------------   ------   ------   ---------------   -------   ------   --------
        12345   CATGTTACAAAAGGTGGGTCAGGCAACGTAT   7        1        GCF_001457655.1   335112    -        yes
        12345   CATGTTACCAAGGTTAGTCGTATGGCGCTAC   7        1        GCF_001457655.1   23755     -        yes
        12345   CATGTTACGCGTATTTTAGCGGCTCGCGGAC   7        1        GCF_000006945.2   702224    +        yes
        12345   CATGTTATAACGGCCTATGAATCGGCATTAC   9        1        GCF_009759685.1   2591866   +        no
        12345   CATGTTATACGTTGAAACTGTCTTGTTAATA   9        1        GCF_001096185.1   1142460   +        yes
        12345   CATGTTATACTTTAGATACTTATTTTTAGGA   9        1        GCF_000392875.1   1524553   +        no
        12345   CATGTTATAGAAGGACGTCGACATCTTGTGG   10       1        GCF_000017205.1   3140677   +        no
        12345   CATGTTATAGAATTACATACATTGTAACATG   10       1        GCF_006742205.1   704431    -        no
        12345   CATGTTATAGCACGCTTAATCGCTTGATCCC   13       1        GCF_001027105.1   2655846   +        no
        12345   CATGTTATAGCATCCTTTTACGTGAAAAGGT   12       1        GCF_000742135.1   4136093   +        no
        12345   CATGTTATAGCCAGCAAATGGAAGCATCGCG   11       1        GCF_009759685.1   492828    -        no
        12345   CATGTTATAGCCATTGATGGTAACTTTGATG   11       1        GCF_001096185.1   536843    +        no
        12345   CATGTTATAGCCTGAAAGGTGCTAAACAACT   11       1        GCF_000006945.2   4876155   +        no
        12345   CATGTTATAGCCTTCTCCAAGACCAATCAAA   11       1        GCF_000148585.2   1667015   +        no
        12345   CATGTTATAGCGTAAATCAGCACCGCGCGCC   11       3        GCF_002949675.1   1871326   +        no
        12345   CATGTTATAGCGTAAATCAGCACCGCGCGCC   11       3        GCF_002950215.1   2326544   +        no
        12345   CATGTTATAGCGTAAATCAGCACCGCGCGCC   11       3        GCF_003697165.2   3996124   +        no
        12345   CATGTTATAGCTAACTGCGACTTGTGGCACA   11       1        GCF_900638025.1   991007    -        no
        12345   CATGTTATAGTCGTGAGGTTCTAAAAAAACT   10       1        GCF_001544255.1   1091256   -        no

    Lengths of shared prefixes between probes and captured k-mers.

        zcat kmers.tsv.gz \
          | csvtk filter2 -t -f '$reversed == "no"'\
          | csvtk plot hist -t -f prefix -o prefix.hist.png \
              --xlab "length of common prefixes between captured k-mers and masks"


    <img src="/LexicMap/prefix.hist.png" alt="" width="400"/>

The output (TSV format) is formatted with [csvtk pretty](https://github.com/shenwei356/csvtk).
