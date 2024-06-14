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
  -m, --masks int         ► Number of masks. (default 1000)
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
1       AAAACACATGTGGCCGTTCCGATAATGATCC
2       AAAACACGCGATCAGGTCGTCCGTCCCAGTG
3       AAAACACTATGGCCTGATTACCCCATCCCGA
4       AAAACAGGACCGTCCTAGGGTCAATGGTTCG
5       AAAACAGTCTTGTATTATGTACTTCACATTC
6       AAAACATTCTTTCAGATCCATGCACGAGCAG
7       AAAACATTGGTCCTATTGGCGTCACTCGATA
8       AAAACCAGCTCTGTAAGCACTAACAACGCTA
9       AAAACCAGTTACTACGGTTTTCCGCAATTGG
10      AAAACCATGTGCGGCCCATAACTGGTGGTAA

$ lexicmap utils masks --quiet -d demo.lmi/ | tail -n 10
39991   TTTTGGTCTACAGAAAGTGCGTTATAGATTT
39992   TTTTGGTGCTCATACAGTACCGCTGTCTCCA
39993   TTTTGTAGCGCCATACCGGAGTTCACGGGTC
39994   TTTTGTATGGAATACTTTACAGTCATCAGTT
39995   TTTTGTCATACAGGAAGCACGGTCGAACGTT
39996   TTTTGTCATCAGTCGGCACTTAGGGGAACCG
39997   TTTTGTCCACCCTATTGAAGATGGAGTCCCA
39998   TTTTGTGAATAAGAGATCCTGTCGCAGGAAA
39999   TTTTGTGGACGGTAGAACCCCGTGGTCACCA
40000   TTTTGTGTAGCATAAAGGGCTGGCGACCTGT
```
