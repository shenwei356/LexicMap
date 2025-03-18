---
title: 2blast
weight: 0
---

## Usage

```plain
$ lexicmap utils 2blast -h
Convert the default search output to blast-style format

LexicMap only stores genome IDs and sequence IDs, without description information.
But the option -g/--kv-file-genome enables adding description data after the genome ID
with a tabular key-value mapping file.

Input:
   - Output of 'lexicmap search' with the flag -a/--all.

Usage:
  lexicmap utils 2blast [flags]

Flags:
  -b, --buffer-size string      ► Size of buffer, supported unit: K, M, G. You need increase the value
                                when "bufio.Scanner: token too long" error reported (default "20M")
  -h, --help                    help for 2blast
  -i, --ignore-case             ► Ignore cases of sgenome and sseqid
  -g, --kv-file-genome string   ► Two-column tabular file for mapping the target genome ID (sgenome)
                                to the corresponding value
  -s, --kv-file-seq string      ► Two-column tabular file for mapping the target sequence ID (sseqid)
                                to the corresponding value
  -o, --out-file string         ► Out file, supports and recommends a ".gz" suffix ("-" for stdout).
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
    | lexicmap utils 2blast --kv-file-genome ass2species.map

Query = GCF_000017205.1_r160
Length = 478

[Subject genome #1/1] = GCF_000017205.1 Pseudomonas aeruginosa
Query coverage per genome = 98.536%

>NC_009656.1 
Length = 6588339

 HSP cluster #1, HSP #1
 Score = 883 bits, Expect = 3.60e-256
 Query coverage per seq = 98.536%, Aligned length = 479, Identities = 94.990%, Gaps = 15
 Query range = 7-477, Subject range = 4866857-4867328, Strand = Plus/Plus

Query  7        GGTGGCCCTCAAACGAGTCC-AACAGGCCAACGCCTAGCAATCCCTCCCCTGTGGGGCAG  65
                ||||||| |||||||||||| |||||||| ||||||  | ||||||||||||| ||||||
Sbjct  4866857  GGTGGCC-TCAAACGAGTCCGAACAGGCCCACGCCTCACGATCCCTCCCCTGTCGGGCAG  4866915

Query  66       GGAAAATCGTCCTTTATGGTCCGTTCCGGGCACGCACCGGAACGGCGGTCATCTTCCACG  125
                |||||||||||||||||||||||||||||||||||||||||||||||||||| |||||||
Sbjct  4866916  GGAAAATCGTCCTTTATGGTCCGTTCCGGGCACGCACCGGAACGGCGGTCAT-TTCCACG  4866974

Query  126      GTGCCCGCCCACGGCGGACCCGCGGAAACCGACCCGGGCGCCAAGGCGCCCGGGAACGGA  185
                ||||||||| ||||||||||| ||||||||||||||||||||||||||||||||||||||
Sbjct  4866975  GTGCCCGCC-ACGGCGGACCC-CGGAAACCGACCCGGGCGCCAAGGCGCCCGGGAACGGA  4867032

Query  186      GTA-CACTCGGCGTTCGGCCAGCGACAGC---GACGCGTTGCCGCCCACCGCGGTGGTGT  241
                ||| |||||||||| ||||||||||||||   ||||||||||||||||||||||||||||
Sbjct  4867033  GTATCACTCGGCGT-CGGCCAGCGACAGCAGCGACGCGTTGCCGCCCACCGCGGTGGTGT  4867091

Query  242      TCACCGAGGTGGTGCGCTCGCTGAC-AAACGCAGCAGGTAGTTCGGCCCGCCGGCCTTGG  300
                ||||||||||||||||||||||||| ||||||||||||||||||||||||||||||||||
Sbjct  4867092  TCACCGAGGTGGTGCGCTCGCTGACGAAACGCAGCAGGTAGTTCGGCCCGCCGGCCTTGG  4867151

Query  301      GACCG-TGCCGGACAGCCCGTGGCCGCCGAACAGTTGCACGCCCACCACCGCGCCGAT-T  358
                ||||| |||||||||||||||||||||||||| ||||||||||||||||||||||||| |
Sbjct  4867152  GACCGGTGCCGGACAGCCCGTGGCCGCCGAACGGTTGCACGCCCACCACCGCGCCGATCT  4867211

Query  359      GGTTTCGGTTGACGTAGAGGTTGCCGACCCGCGCCAGCTCTTGGATGCGGCGGGCGGTTT  418
                |||| ||||||||||||||||||||||||||||||||||||| |||||||||||||||||
Sbjct  4867212  GGTTGCGGTTGACGTAGAGGTTGCCGACCCGCGCCAGCTCTTCGATGCGGCGGGCGGTTT  4867271

Query  419      CCTCGTTGCGGCTGTGGACCCCCATGGTCAGGCCGAAACCGGTGGCGTTTGATGGCCCT  477
                ||||||||||||||||||||||||||||||||||||||||||||||||| ||| ||| |
Sbjct  4867272  CCTCGTTGCGGCTGTGGACCCCCATGGTCAGGCCGAAACCGGTGGCGTT-GATCGCC-T  4867328


Query = GCF_006742205.1_r100
Length = 431

[Subject genome #1/1] = GCF_006742205.1 Staphylococcus epidermidis
Query coverage per genome = 93.968%

>NZ_AP019721.1 
Length = 2422602

 HSP cluster #1, HSP #1
 Score = 740 bits, Expect = 2.39e-213
 Query coverage per seq = 93.968%, Aligned length = 408, Identities = 98.284%, Gaps = 4
 Query range = 27-431, Subject range = 1321677-1322083, Strand = Plus/Minus

Query  27       TTCATTTAAAACGATTGCTAATGAGTCACGTATTTCATCTGGTTCGGTAACTATACCGTC  86
                ||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1322083  TTCATCTAAAACGATTGCTAATGAGTCACGTATTTCATCTGGTTCGGTAACTATACCGTC  1322024

Query  87       TACTATGGACTCAGTGTAACCCTGTAATAAAGAGATTGGCGTACGTAATTCATGTG-TAC  145
                |||||||||||||||||||||||||||||||||||||||||||||||||||||||| |||
Sbjct  1322023  TACTATGGACTCAGTGTAACCCTGTAATAAAGAGATTGGCGTACGTAATTCATGTGATAC  1321964

Query  146      ATTTGCTATAAAATCTTTTTTCATTTGATCAAGATTATGTTCATTTGTCATATCACAGGA  205
                |||||||||||||||||||||||||||||||||||||||||||||||||||||||| |||
Sbjct  1321963  ATTTGCTATAAAATCTTTTTTCATTTGATCAAGATTATGTTCATTTGTCATATCAC-GGA  1321905

Query  206      TGACCATGACAATACCACTTCTACCATTTGTTTGAATTCTATCTATATAACTGGAGATAA  265
                ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1321904  TGACCATGACAATACCACTTCTACCATTTGTTTGAATTCTATCTATATAACTGGAGATAA  1321845

Query  266      ATACATAGTACCTTGTATTAATTTCTAATTCTAA-TACTCATTCTGTTGTGATTCAAATG  324
                |||||||||||||||||||||||||||||||||| |||||||||||||||||||||||||
Sbjct  1321844  ATACATAGTACCTTGTATTAATTTCTAATTCTAAATACTCATTCTGTTGTGATTCAAATG  1321785

Query  325      GTGCTTCAATTTGCTGTTCAATAGATTCTTTTGAAAAATCATCAATGTGACGCATAATAT  384
                 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1321784  TTGCTTCAATTTGCTGTTCAATAGATTCTTTTGAAAAATCATCAATGTGACGCATAATAT  1321725

Query  385      AATCAGCCATCTTGTT-GACAATATGATTTCACGTTGATTATTAATGC  431
                 ||||||||||||||| |||||||||||||||||||||||||||||||
Sbjct  1321724  CATCAGCCATCTTGTTTGACAATATGATTTCACGTTGATTATTAATGC  1321677

```


From file.

    $ lexicmap utils 2blast r.lexicmap.tsv -o r.lexicmap.txt
