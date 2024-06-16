---
title: kmers
weight: 10
---

```plain
$ lexicmap utils kmers -h
View k-mers captured by the masks

Attentions:
  1. Mask index (column mask) is 1-based.
  2. K-mer positions (column pos) are 1-based.
     For reference genomes with multiple sequences, the sequences were
     concatenated to a single sequence with intervals of N's.

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
                             (default 16
```

## Examples

1. The default output is captured k-mers of the first mask.

        $ lexicmap utils kmers --quiet -d demo.lmi/ | csvtk pretty -t
        mask   kmer                              number   ref               pos       strand
        ----   -------------------------------   ------   ---------------   -------   ------
        1      AAAACACATGCTTTCACTGACTTGGAATGCA   1        GCF_001457655.1   389653    +
        1      AAAACACATGGATTGTTAAAAGGTAGTTGGC   1        GCF_900638025.1   2061446   -
        1      AAAACACATGTAAGCCCCAACCAGGCGGCTT   1        GCF_000742135.1   2569538   -
        1      AAAACACATGTCTAAAATTATCGGTATTGAC   2        GCF_000148585.2   326139    +
        1      AAAACACATGTCTAAAATTATCGGTATTGAC   2        GCF_001096185.1   34675     -
        1      AAAACACATGTGAGGCAGGCGCTCGCCTGTC   1        GCF_001544255.1   938768    -
        1      AAAACACATGTGCAAATCCATATGTGTTTAG   1        GCF_002950215.1   2793719   +
        1      AAAACACATGTGTTGTTTAAATCAAATTATG   1        GCF_001027105.1   1413381   +
        1      AAAACACATGTGTTTAATCACCTTAATTCAA   1        GCF_006742205.1   729899    +
        1      AAAACACATGTTCACGGCGGCAGGCTGCAAT   1        GCF_003697165.2   1581455   +
        1      AAAACACATGTTCTCAATACTCGCCTGACGC   1        GCF_000006945.2   1274137   -
        1      AAAACACATGTTGATCATCATAAATACAGCG   1        GCF_002949675.1   3925773   -
        1      AAAACACATGTTGATCTATTCTTATAGCTCA   1        GCF_009759685.1   3295037   -
        1      AAAACACATGTTTCAAACATTTTAGCAAAAC   1        GCF_000392875.1   2491283   -
        1      AAAACACATGTTTCACACAACTTCACCCAAT   1        GCF_000017205.1   4394137   +


1. Specify the mask.

        $ lexicmap utils kmers --quiet -d demo.lmi/ --mask 12345 | csvtk pretty -t
        mask    kmer                              number   ref               pos       strand
        -----   -------------------------------   ------   ---------------   -------   ------
        12345   CATGTTATAGAAGGACGTCGACATCTTGTGG   1        GCF_000017205.1   3140677   +
        12345   CATGTTATAGAATTACATACATTGTAACATG   1        GCF_006742205.1   704431    -
        12345   CATGTTATAGCACGCTTAATCGCTTGATCCC   1        GCF_001027105.1   2655846   +
        12345   CATGTTATAGCATCCTTTTACGTGAAAAGGT   1        GCF_000742135.1   4136093   +
        12345   CATGTTATAGCCAGCAAATGGAAGCATCGCG   1        GCF_009759685.1   492828    -
        12345   CATGTTATAGCCATTGATGGTAACTTTGATG   1        GCF_001096185.1   536843    +
        12345   CATGTTATAGCCTGAAAGGTGCTAAACAACT   1        GCF_000006945.2   4876155   +
        12345   CATGTTATAGCCTTCTCCAAGACCAATCAAA   1        GCF_000148585.2   1667015   +
        12345   CATGTTATAGCGTAAATCAGCACCGCGCGCC   3        GCF_003697165.2   3996124   +
        12345   CATGTTATAGCGTAAATCAGCACCGCGCGCC   3        GCF_002949675.1   1871326   +
        12345   CATGTTATAGCGTAAATCAGCACCGCGCGCC   3        GCF_002950215.1   2326544   +
        12345   CATGTTATAGCTAACTGCGACTTGTGGCACA   1        GCF_900638025.1   991007    -
        12345   CATGTTATAGTAAACAAAAGTAGTGACGAAT   1        GCF_000392875.1   1539455   -
        12345   CATGTTATAGTCGTGAGGTTCTAAAAAAACT   1        GCF_001544255.1   1091256   -
        12345   CATGTTATATGAACCTTCAACCTTATTTGAC   1        GCF_001457655.1   1510084   +


1. For all masks. The result might be very big, therefore, writing to gzip format is recommended.


        $ lexicmap utils kmers -d demo.lmi/ --mask 0 -o kmers.tsv.gz

        $ zcat kmers.tsv.gz | csvtk freq -t -f 1 -nr | head -n 10
        mask    frequency
        38388   142
        38030   138
        31609   130
        15643   129
        35924   129
        35160   89
        8960    88
        36439   86
        19132   84

        $ lexicmap utils kmers -d demo.lmi/ -m 38388 | head -n 20 | csvtk pretty -t
        mask    kmer                              number   ref               pos       strand
        -----   -------------------------------   ------   ---------------   -------   ------
        38388   TTCCATCAGATATTGCAGTTGCCGCGCCAGC   1        GCF_000017205.1   2157565   -
        38388   TTCCATCAGATCCTCCTACTTCAAAGACCAG   1        GCF_001096185.1   1159072   +
        38388   TTCCATCAGATCCTTTGTCACTACCTGAAGC   1        GCF_001027105.1   1108238   -
        38388   TTCCATCAGATGATGACCGGTGGACGAACAC   1        GCF_001544255.1   1315664   -
        38388   TTCCATCAGATGCTTCTGGTTCTTTATTTAA   1        GCF_000392875.1   503146    -
        38388   TTCCATCAGATGGAAGGTCTGATTGTTGATA   1        GCF_000006945.2   1416284   +
        38388   TTCCATCAGATGTCCCACTTGTTCAGCTACC   9        GCF_000742135.1   38864     -
        38388   TTCCATCAGATGTCCCACTTGTTCAGCTACC   9        GCF_000742135.1   818406    +
        38388   TTCCATCAGATGTCCCACTTGTTCAGCTACC   9        GCF_000742135.1   1622053   -
        38388   TTCCATCAGATGTCCCACTTGTTCAGCTACC   9        GCF_000742135.1   2064518   +
        38388   TTCCATCAGATGTCCCACTTGTTCAGCTACC   9        GCF_000742135.1   2227218   -
        38388   TTCCATCAGATGTCCCACTTGTTCAGCTACC   9        GCF_000742135.1   3381540   -
        38388   TTCCATCAGATGTCCCACTTGTTCAGCTACC   9        GCF_000742135.1   5301234   +
        38388   TTCCATCAGATGTCCCACTTGTTCAGCTACC   9        GCF_000742135.1   5431183   -
        38388   TTCCATCAGATGTCCCACTTGTTCAGCTACC   9        GCF_000742135.1   5433483   -
        38388   TTCCATCAGATGTCCTTCCTGCTCCGCTACT   122      GCF_002950215.1   17655     +
        38388   TTCCATCAGATGTCCTTCCTGCTCCGCTACT   122      GCF_002950215.1   70503     -
        38388   TTCCATCAGATGTCCTTCCTGCTCCGCTACT   122      GCF_002950215.1   81837     +
        38388   TTCCATCAGATGTCCTTCCTGCTCCGCTACT   122      GCF_002950215.1   83307     +

The output (TSV format) is formatted with [csvtk pretty](https://github.com/shenwei356/csvtk).
