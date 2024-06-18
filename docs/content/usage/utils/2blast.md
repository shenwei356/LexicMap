---
title: 2blast
weight: 0
---

## Usage

```plain
$ lexicmap utils 2blast -h
Convert the default search output to blast-style format

Input:
   - Output of 'lexicmap search' with the flag -a/--all.

Usage:
  lexicmap utils 2blast [flags]

Flags:
  -b, --buffer-size string   ► Size of buffer, supported unit: K, M, G. You need increase the value
                             when "bufio.Scanner: token too long" error reported (default "20M")
  -h, --help                 help for 2blast
  -o, --out-file string      ► Out file, supports and recommends a ".gz" suffix ("-" for stdout).
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


From stdin.

```text
$ seqkit seq -M 500 q.long-reads.fasta.gz \
    | seqkit head -n 2 \
    | lexicmap search -d demo.lmi/ -a \
    | lexicmap utils 2blast

Query = GCF_000017205.1_r160
Length = 478

[Subject genome #1/1] = GCF_000017205.1
Query coverage per genome = 95.188%

>NC_009656.1
Length = 6588339

 HSP #1
 Query coverage per seq = 95.188%, Aligned length = 463, Identities = 95.680%, Gaps = 12
 Query range = 13-467, Subject range = 4866862-4867320, Strand = Plus/Plus

Query  13       CCTCAAACGAGTCC-AACAGGCCAACGCCTAGCAATCCCTCCCCTGTGGGGCAGGGAAAA  71
                |||||||||||||| |||||||| ||||||  | ||||||||||||| ||||||||||||
Sbjct  4866862  CCTCAAACGAGTCCGAACAGGCCCACGCCTCACGATCCCTCCCCTGTCGGGCAGGGAAAA  4866921

Query  72       TCGTCCTTTATGGTCCGTTCCGGGCACGCACCGGAACGGCGGTCATCTTCCACGGTGCCC  131
                |||||||||||||||||||||||||||||||||||||||||||||| |||||||||||||
Sbjct  4866922  TCGTCCTTTATGGTCCGTTCCGGGCACGCACCGGAACGGCGGTCAT-TTCCACGGTGCCC  4866980

Query  132      GCCCACGGCGGACCCGCGGAAACCGACCCGGGCGCCAAGGCGCCCGGGAACGGAGTA-CA  190
                ||| ||||||||||| ||||||||||||||||||||||||||||||||||||||||| ||
Sbjct  4866981  GCC-ACGGCGGACCC-CGGAAACCGACCCGGGCGCCAAGGCGCCCGGGAACGGAGTATCA  4867038

Query  191      CTCGGCGTTCGGCCAGCGACAGC---GACGCGTTGCCGCCCACCGCGGTGGTGTTCACCG  247
                |||||||| ||||||||||||||   ||||||||||||||||||||||||||||||||||
Sbjct  4867039  CTCGGCGT-CGGCCAGCGACAGCAGCGACGCGTTGCCGCCCACCGCGGTGGTGTTCACCG  4867097

Query  248      AGGTGGTGCGCTCGCTGAC-AAACGCAGCAGGTAGTTCGGCCCGCCGGCCTTGGGACCG-  305
                ||||||||||||||||||| |||||||||||||||||||||||||||||||||||||||
Sbjct  4867098  AGGTGGTGCGCTCGCTGACGAAACGCAGCAGGTAGTTCGGCCCGCCGGCCTTGGGACCGG  4867157

Query  306      TGCCGGACAGCCCGTGGCCGCCGAACAGTTGCACGCCCACCACCGCGCCGAT-TGGTTTC  364
                |||||||||||||||||||||||||| ||||||||||||||||||||||||| ||||| |
Sbjct  4867158  TGCCGGACAGCCCGTGGCCGCCGAACGGTTGCACGCCCACCACCGCGCCGATCTGGTTGC  4867217

Query  365      GGTTGACGTAGAGGTTGCCGACCCGCGCCAGCTCTTGGATGCGGCGGGCGGTTTCCTCGT  424
                |||||||||||||||||||||||||||||||||||| |||||||||||||||||||||||
Sbjct  4867218  GGTTGACGTAGAGGTTGCCGACCCGCGCCAGCTCTTCGATGCGGCGGGCGGTTTCCTCGT  4867277

Query  425      TGCGGCTGTGGACCCCCATGGTCAGGCCGAAACCGGTGGCGTT  467
                |||||||||||||||||||||||||||||||||||||||||||
Sbjct  4867278  TGCGGCTGTGGACCCCCATGGTCAGGCCGAAACCGGTGGCGTT  4867320


Query = GCF_006742205.1_r100
Length = 431

[Subject genome #1/1] = GCF_006742205.1
Query coverage per genome = 92.575%

>NZ_AP019721.1
Length = 2422602

 HSP #1
 Query coverage per seq = 92.575%, Aligned length = 402, Identities = 98.507%, Gaps = 4
 Query range = 33-431, Subject range = 1321677-1322077, Strand = Plus/Minus

Query  33       TAAAACGATTGCTAATGAGTCACGTATTTCATCTGGTTCGGTAACTATACCGTCTACTAT  92
                ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1322077  TAAAACGATTGCTAATGAGTCACGTATTTCATCTGGTTCGGTAACTATACCGTCTACTAT  1322018

Query  93       GGACTCAGTGTAACCCTGTAATAAAGAGATTGGCGTACGTAATTCATGTG-TACATTTGC  151
                |||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||
Sbjct  1322017  GGACTCAGTGTAACCCTGTAATAAAGAGATTGGCGTACGTAATTCATGTGATACATTTGC  1321958

Query  152      TATAAAATCTTTTTTCATTTGATCAAGATTATGTTCATTTGTCATATCACAGGATGACCA  211
                |||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||
Sbjct  1321957  TATAAAATCTTTTTTCATTTGATCAAGATTATGTTCATTTGTCATATCAC-GGATGACCA  1321899

Query  212      TGACAATACCACTTCTACCATTTGTTTGAATTCTATCTATATAACTGGAGATAAATACAT  271
                ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1321898  TGACAATACCACTTCTACCATTTGTTTGAATTCTATCTATATAACTGGAGATAAATACAT  1321839

Query  272      AGTACCTTGTATTAATTTCTAATTCTAA-TACTCATTCTGTTGTGATTCAAATGGTGCTT  330
                |||||||||||||||||||||||||||| ||||||||||||||||||||||||| |||||
Sbjct  1321838  AGTACCTTGTATTAATTTCTAATTCTAAATACTCATTCTGTTGTGATTCAAATGTTGCTT  1321779

Query  331      CAATTTGCTGTTCAATAGATTCTTTTGAAAAATCATCAATGTGACGCATAATATAATCAG  390
                |||||||||||||||||||||||||||||||||||||||||||||||||||||| |||||
Sbjct  1321778  CAATTTGCTGTTCAATAGATTCTTTTGAAAAATCATCAATGTGACGCATAATATCATCAG  1321719

Query  391      CCATCTTGTT-GACAATATGATTTCACGTTGATTATTAATGC  431
                |||||||||| |||||||||||||||||||||||||||||||
Sbjct  1321718  CCATCTTGTTTGACAATATGATTTCACGTTGATTATTAATGC  1321677


```


From file.

    $ lexicmap utils 2blast r.lexicmap.tsv -o r.lexicmap.txt
