---
title: masks
weight: 5
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
  -m, --masks int         ► Number of masks. (default 40000)
  -o, --out-file string   ► Out file, supports and recommends a ".gz" suffix ("-" for stdout).
                          (default "-")
  -p, --prefix int        ► Length of mask k-mer prefix for checking low-complexity (0 for no
                          checking). (default 15)
  -s, --seed int          ► The seed for generating random masks. (default 1)

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

```plain
$ lexicmap utils masks --quiet -d demo.lmi/ | head -n 10
1       AAAAAAAAGTCACTTGACAATCCACACGGTG
2       AAAAAAACTGCTTGCACCTTTCTCGCCTCTC
3       AAAAAAATTCTCGGCGGTGTTTCCAGGCGCA
4       AAAAAACCCAAGCGCGAAAGCCTGAACAACC
5       AAAAAACGTGGCGTCCCCTGTATAACGGCTA
6       AAAAAAGAGGGGAAGCAAGCTGAAGGATATG
7       AAAAAAGCTTAGTGTGAATGAATGGCTTCCG
8       AAAAAATCCAGGGTTCCGTTAAGGATCTGTC
9       AAAAAATGCCTCGCAGAGCAGGCTATGCTGA
10      AAAAAATTGATTCTTAGAGCGTTCCCGCCCA

$ lexicmap utils masks --quiet -d demo.lmi/ | tail -n 10
39991   TTTTTTACACGCTGTGACTGCATTACAAAAA
39992   TTTTTTAGCCAGGGTTCACAGCGCCAAAACA
39993   TTTTTTATCGGACGCCAAGTTTGTAATCGTC
39994   TTTTTTCACTCGCATCTAGGAAGGAAGCATA
39995   TTTTTTCTTGCATCGTATTCAGCACGTTCCT
39996   TTTTTTGCCGAGTGACCCCGAAAAGCTCACA
39997   TTTTTTGGCGTGAGGCATTGTTTACTGCCTT
39998   TTTTTTTAAGTGGTCGTGGTAGGAGCCTCAC
39999   TTTTTTTCCGTAACTAGGTTCTGGCGATTCC
40000   TTTTTTTGAGGGTATAAGATAGAGAAAAGCT

# check a specific mask

$ lexicmap utils masks --quiet -d demo.lmi/ -m 12345
12345   CATTAGTAGAAGAAGGCACAATGTATCGTCG
```

Freqency of prefixes.

```
$ lexicmap utils masks --quiet -d demo.lmi/ \
  | csvtk mutate -Ht -f 2 -p '^(.{7})' \
  | csvtk freq -Ht -f 3 -nr \
  | head -n 10
AAAAAAA 3
AAAAAAT 3
AAAAACA 3
AAAAACC 3
AAAAACG 3
AAAAACT 3
AAAAAGC 3
AAAAAGG 3
AAAAAGT 3
AAAAATT 3

$ lexicmap utils masks --quiet -d demo.lmi/ \
  | csvtk mutate -Ht -f 2 -p '^(.{7})' \
  | csvtk freq -Ht -f 3 -n \
  | head -n 10
AAAAAAC 2
AAAAAAG 2
AAAAAGA 2
AAAAATA 2
AAAAATC 2
AAAAATG 2
AAAACAC 2
AAAACAT 2
AAAACCG 2
AAAACGC 2
```

Frequency of frequencies. i.e., for 40,000 masks, 4<sup>*7*</sup> = 16384.
All 16,384 masks are duplicated twice, and 7,232 of them are duplicated 3 times.

```
$ lexicmap utils masks --quiet -d demo.lmi/ | csvtk mutate -Ht -f 2 -p '^(.{7})' | csvtk freq -Ht -f 3 -n | csvtk freq -Ht -f 2 -k
2       9152
3       7232
```
