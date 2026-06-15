---
title: genome-seqs
weight: 22
---

## Usage

```plain
Extract all sequences of a given genome

Attention:
  1. All degenerate bases in reference genomes were converted to the lexicographic first bases.
     E.g., N was converted to A. Therefore, consecutive A's in output might be N's in the genomes.
  2. Large genomes fragmented into multiple chunks during indexing (total size > 15 Mb by default,
     configurable with -g/--max-genome in 'lexicmap index'), such as many fungal genomes,
     may have their sequence order rearranged relative to the original input files.

Usage:
  lexicmap utils genome-seqs [flags] 

Flags:
  -h, --help              help for genome-seqs
  -d, --index string      ► Index directory created by "lexicmap index".
  -w, --line-width int    ► Line width of sequence (0 for no wrap). (default 60)
  -o, --out-file string   ► Out file, supports the ".gz" suffix ("-" for stdout). (default "-")
  -n, --ref-name string   ► Reference name.

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

Find a genome without N's in the demo, such as `demo/refs/GCF_001544255.1.fa.gz`.

```text
$ seqkit stats demo/refs/*.fa.gz -a -G Nn -T | csvtk cut -t -f 1-12 | csvtk pretty -t
processed files:  15 / 15 [======================================] ETA: 0s. done
file                              format   type   num_seqs   sum_len   min_len   avg_len     max_len   Q1        Q2        Q3        sum_gap
-------------------------------   ------   ----   --------   -------   -------   ---------   -------   -------   -------   -------   -------
demo/refs/GCF_000006945.2.fa.gz   FASTA    DNA    2          4951383   93933     2475691.5   4857450   93933     2475692   4857450   0      
demo/refs/GCF_000017205.1.fa.gz   FASTA    DNA    1          6588339   6588339   6588339.0   6588339   6588339   6588339   6588339   0      
demo/refs/GCF_000148585.2.fa.gz   FASTA    DNA    1          1868883   1868883   1868883.0   1868883   1868883   1868883   1868883   0      
demo/refs/GCF_000392875.1.fa.gz   FASTA    DNA    3          2881400   274762    960466.7    1924212   478594    682426    1303319   4965   
demo/refs/GCF_000742135.1.fa.gz   FASTA    DNA    5          5545784   16331     1109156.8   5284261   42420     95930     106842    1100   
demo/refs/GCF_001027105.1.fa.gz   FASTA    DNA    2          2782562   27490     1391281.0   2755072   27490     1391281   2755072   0      
demo/refs/GCF_001096185.1.fa.gz   FASTA    DNA    24         2117177   366       88215.7     543880    9252      48418     133692    4      
demo/refs/GCF_001457655.1.fa.gz   FASTA    DNA    1          1890645   1890645   1890645.0   1890645   1890645   1890645   1890645   0      
demo/refs/GCF_001544255.1.fa.gz   FASTA    DNA    38         2484851   563       65390.8     386577    1364      20536     115935    0      
demo/refs/GCF_002949675.1.fa.gz   FASTA    DNA    2          4578459   182697    2289229.5   4395762   182697    2289230   4395762   0      
demo/refs/GCF_002950215.1.fa.gz   FASTA    DNA    3          4938295   113130    1646098.3   4659463   139416    165702    2412582   0      
demo/refs/GCF_003697165.2.fa.gz   FASTA    DNA    2          5034834   131333    2517417.0   4903501   131333    2517417   4903501   0      
demo/refs/GCF_006742205.1.fa.gz   FASTA    DNA    2          2427041   4439      1213520.5   2422602   4439      1213520   2422602   0      
demo/refs/GCF_009759685.1.fa.gz   FASTA    DNA    2          3990388   9540      1995194.0   3980848   9540      1995194   3980848   0      
demo/refs/GCF_900638025.1.fa.gz   FASTA    DNA    1          2062405   2062405   2062405.0   2062405   2062405   2062405   2062405   0      
```

Extract contigs and check the content.

```text
$ lexicmap utils genome-seqs -d demo/demo.lmi/ -n GCF_001544255.1 | seqkit stats 
file  format  type  num_seqs    sum_len  min_len   avg_len  max_len
-     FASTA   DNA         38  2,484,851      563  65,390.8  386,577

$ lexicmap utils genome-seqs -d demo/demo.lmi/ -n GCF_001544255.1 | seqkit sum
seqkit.v0.1_DLS_k0_4369d22ee7050db2833ae66e78aa28a6     -

# original genome

$ seqkit stats demo/refs/GCF_001544255.1.fa.gz
file                             format  type  num_seqs    sum_len  min_len   avg_len  max_len
demo/refs/GCF_001544255.1.fa.gz  FASTA   DNA         38  2,484,851      563  65,390.8  386,577

$ seqkit sum demo/refs/GCF_001544255.1.fa.gz
seqkit.v0.1_DLS_k0_4369d22ee7050db2833ae66e78aa28a6     demo/refs/GCF_001544255.1.fa.gz
```

But note that, Lexicmap does not store description in the header line.

```text
$ lexicmap utils genome-seqs -d demo/demo.lmi/ -n GCF_001544255.1 | seqkit seq -n | head -n 3
NZ_BCQD01000001.1
NZ_BCQD01000002.1
NZ_BCQD01000003.1

$ seqkit seq -n demo/refs/GCF_001544255.1.fa.gz | head -n 3
NZ_BCQD01000001.1 Enterococcus faecium NBRC 100486, whole genome shotgun sequence
NZ_BCQD01000002.1 Enterococcus faecium NBRC 100486, whole genome shotgun sequence
NZ_BCQD01000003.1 Enterococcus faecium NBRC 100486, whole genome shotgun sequence
```
