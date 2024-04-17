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
  -X, --infile-list string   ► File of input files list (one file per line). If given, they are
                             appended to files from CLI arguments.
      --log string           ► Log file.
      --quiet                ► Do not print any verbose information. But you can write them to a file
                             with --log.
  -j, --threads int          ► Number of CPUs cores to use. By default, it uses all available cores.
                             (default 16)
```

## Examples

1. The default output is captured k-mers of the first mask.

        $ lexicmap utils kmers --quiet -d demo.lmi/ | csvtk pretty -t
        mask   kmer                              number   ref               pos       strand
        ----   -------------------------------   ------   ---------------   -------   ------
        1      AAAAAAACGTAATTGTGTGGGTACAGATTGG   1        GCF_000148585.2   1175190   +
        1      AAAAAAACGTCAATCAATTTGACCGCACTTT   1        GCF_900638025.1   774645    +
        1      AAAAAAACGTCACAACAACAAGACAAAGTGA   1        GCF_000392875.1   2565649   -
        1      AAAAAAACGTCACTTGTCCTTTCCTACCCTT   1        GCF_001457655.1   1711468   +
        1      AAAAAAACGTCATAAACCGCAACAGCCGCGC   1        GCF_000742135.1   144039    -
        1      AAAAAAACGTCATAACGAATTTCCGGGCGCG   1        GCF_000006945.2   3120524   -
        1      AAAAAAACGTCATATCCCAATGGCGATTACT   1        GCF_006742205.1   640624    +
        1      AAAAAAACGTCATCGCTTGCATTAGAAAGGT   3        GCF_003697165.2   339475    -
        1      AAAAAAACGTCATCGCTTGCATTAGAAAGGT   3        GCF_002949675.1   2873885   +
        1      AAAAAAACGTCATCGCTTGCATTAGAAAGGT   3        GCF_002950215.1   2692249   +
        1      AAAAAAACGTCATGAAAAACGGCGATTGACA   1        GCF_001544255.1   1283395   -
        1      AAAAAAACGTCCATTCTATCAAATCGTTGTG   1        GCF_009759685.1   970548    -
        1      AAAAAAACGTCCCAATATTATTGGGACGTTT   1        GCF_001027105.1   2044210   -
        1      AAAAAAACGTCGGAGCCACATCTGGCTCTTC   1        GCF_000017205.1   2674197   +
        1      AAAAAAACGTCTAGCCAATGACCCATTCCTG   1        GCF_001096185.1   2127010   +


1. Specify the maks.

        $ lexicmap utils kmers --quiet -d demo.lmi/ --mask 12345 | csvtk pretty -t
        mask    kmer                              number   ref               pos       strand
        -----   -------------------------------   ------   ---------------   -------   ------
        12345   CATGGTGCTTAAAAAAGTAGGGGCAAAGTCC   1        GCF_009759685.1   1468773   +
        12345   CATGGTGCTTAACCTCGACAAATGTATCGGT   1        GCF_002949675.1   357987    -
        12345   CATGGTGCTTAACCTCGATAAATGTATCGGC   1        GCF_002950215.1   532721    -
        12345   CATGGTGCTTAATCGCCGCCTGCGCCGCCAC   1        GCF_000742135.1   3801879   +
        12345   CATGGTGCTTACGCAGGTCAGCGGGTATGAC   1        GCF_000006945.2   3792773   -
        12345   CATGGTGCTTACGGTGGATCGCGGTCCATTC   1        GCF_000017205.1   383787    -
        12345   CATGGTGCTTACTAGTTATGAAATCATCAGT   1        GCF_001457655.1   1887005   -
        12345   CATGGTGCTTACTTTACGGTTTTGCTGGCCT   1        GCF_003697165.2   1579125   +
        12345   CATGGTGCTTATGCTACGGGTGATGCATTCT   1        GCF_001096185.1   1060887   -
        12345   CATGGTGCTTCAACCGAACCAAAATTTTGAA   1        GCF_000392875.1   2836858   +
        12345   CATGGTGCTTCATAACCAGGTACTAAACGTT   1        GCF_001027105.1   1305877   -
        12345   CATGGTGCTTCATCAAGTGAAAACTCTTATA   1        GCF_000148585.2   456813    -
        12345   CATGGTGCTTCGTAACCTGGTACAAGACGTT   1        GCF_006742205.1   1534889   +
        12345   CATGGTGCTTGAGGATTGATCAACGGCACTG   1        GCF_001544255.1   1428518   +
        12345   CATGGTGCTTGATCATGGAAAATCACTTGTG   1        GCF_900638025.1   468185    +


1. For all masks. The result might be very big, therefore, writing to gzip format is recommended.


        $ lexicmap utils kmers --quiet -d demo.lmi/ --mask 0 -o kmers.tsv.gz

The output (TSV format) is formatted with [csvtk pretty](https://github.com/shenwei356/csvtk).
