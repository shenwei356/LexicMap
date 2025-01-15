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
1       AAAAAAATTCTCGGCGGTGTTTCCAGGCGCA
2       AAAAAACGTGGCGTCCCCTGTATAACGGCTA
3       AAAAAAGAGGGGAAGCAAGCTGAAGGATATG
4       AAAAAATACAGGCTGGCATCTTTAACCCACC
5       AAAAAATCCAGGGTTCCGTTAAGGATCTGTC
6       AAAAACATTCATGCTAGCATACCTTGGCAAC
7       AAAAACCACAATGTGGAAGCACGAGAGGATT
8       AAAAACCTGTACCCACCCGACGTGGATCCTC
9       AAAAACGTAGGCGTACCTCTCATAGCTTGTA
10      AAAAACTATGGATACTTGCCGTAAATCACCT

$ lexicmap utils masks --quiet -d demo.lmi/ | tail -n 10
19991   TTTTTGAACTTGTGAAAAAGGCAGATGTGTG
19992   TTTTTGCGTTTATGCTGCCCTCAAACCATCT
19993   TTTTTGGATCCACTGTACGAGCACACTACCC
19994   TTTTTGTGGCTCATCGGGATCGGGAGCAGTC
19995   TTTTTTACATGTTGGGCTAGGGGCGGTTCAC
19996   TTTTTTATCGGACGCCAAGTTTGTAATCGTC
19997   TTTTTTCTTGCATCGTATTCAGCACGTTCCT
19998   TTTTTTGCCGAGTGACCCCGAAAAGCTCACA
19999   TTTTTTTATCGAGGCATGGTTGAAGACGGGT
20000   TTTTTTTCCGTAACTAGGTTCTGGCGATTCC

# check a specific mask

$ lexicmap utils masks --quiet -d demo.lmi/ -m 12345
12345   GCTGCACACGCAAAGACTCACGTCTTCAACG
```

Freqency of prefixes.

```
$ lexicmap utils masks --quiet -d demo.lmi/ \
  | csvtk mutate -Ht -f 2 -p '^(.{7})' \
  | csvtk freq -Ht -f 3 -nr \
  | head -n 10
AAAAAAT 2
AAAAACC 2
AAAAACT 2
AAAAAGG 2
AAAAAGT 2
AAAAATT 2
AAAACCA 2
AAAACCC 2
AAAACGA 2
AAAACTA 2

$ lexicmap utils masks --quiet -d demo.lmi/ \
  | csvtk mutate -Ht -f 2 -p '^(.{7})' \
  | csvtk freq -Ht -f 3 -n \
  | head -n 10
AAAAAAA 1
AAAAAAC 1
AAAAAAG 1
AAAAACA 1
AAAAACG 1
AAAAAGA 1
AAAAAGC 1
AAAAATA 1
AAAAATC 1
AAAAATG 1
```

Frequency of frequencies. i.e., for 20,000 masks, 4<sup>*7*</sup> = 16384.
In them, 3,616 of them are duplicated 2 times. 12768 + 2 * 3616 = 20000.

```
$ lexicmap utils masks --quiet -d demo.lmi/ \
  | csvtk mutate -Ht -f 2 -p '^(.{7})' \
  | csvtk freq -Ht -f 3 -n \
  | csvtk freq -Ht -f 2 -k
1       12768
2       3616
```
