---
title: masks
weight: 0
---

```plain
$ lexicmap utils masks -h
View masks of the index or generate new masks randomly

Usage:
  lexicmap utils masks [flags] { -d <index path> | [-k <k>] [-n <masks>] [-s <seed>] } [-o out.tsv.gz]

Flags:
  -h, --help              help for masks
  -d, --index string      ► Index directory created by "lexicmap index".
  -k, --kmer int          ► Maximum k-mer size. K needs to be <= 32. (default 31)
  -m, --masks int         ► Number of masks. (default 1000)
  -o, --out-file string   ► Out file, supports and recommends a ".gz" suffix ("-" for stdout).
                          (default "-")
  -p, --prefix int        ► Length of mask k-mer prefix for checking low-complexity (0 for no
                          checking). (default 15)
  -s, --seed int          ► The seed for generating random masks. (default 1)

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

```plain
$ lexicmap utils masks --quiet -d demo.lmi/ | head -n 10
1       AAAAAAACGTCATCGAACCGTGACTTTCTCC
2       AAAAAACCATTCAACAAGCTGTTCGACCCGG
3       AAAAAAGCCATCGTTGCAAGTACCATTGATC
4       AAAAAATAAGACTCTCGAGCATTTCTCACAA
5       AAAAACACCATCCCCGACTCGTATATGCCTC
6       AAAAACCACGCTGCCTGTCCAGCCATGGCCC
7       AAAAACGCCCCGCCAGGACTAATTCCGATTT
8       AAAAACTAGTGTATGCACAAGACGTACCTTC
9       AAAAAGAATAGTGCGCCGTCCTTTACACAGC
10      AAAAAGCCGGCGAAGAAATGCTTGCATGTTC

$ lexicmap utils masks --quiet -d demo.lmi/ | tail -n 10
39991   TTTTTCGCCAGAGGATAGGTCTAAGAATAGC
39992   TTTTTCTAACATCACCGCGGATATTCTATAC
39993   TTTTTGAAGGCGCGGGCTCAGTGGGCTTTCT
39994   TTTTTGCATAATAGTTGTCGGGAGAAGTTTA
39995   TTTTTGGAAGGGGCACGCGGAGGAGCTGCAT
39996   TTTTTGTCCATGAGTCTAACAACTAGCCACC
39997   TTTTTTAAAATGCAAGATCTTATGTACAATT
39998   TTTTTTCCAGTGACCACGATGGGCACTTATA
39999   TTTTTTGCCGGAACCACGCTAGCGCCAGACG
40000   TTTTTTTATTGAAACGACAACCTAACGCGGG
```
