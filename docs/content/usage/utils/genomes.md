---
title: genomes
weight: 20
---

## Usage

```plain
$ lexicmap utils genomes -h
View genome IDs in the index

Usage:
  lexicmap utils genomes [flags]

Flags:
  -e, --extra             ► Show extra columns, such as where the genome is stored (genome_batch,
                          genome_index)
  -h, --help              help for genomes
  -d, --index string      ► Index directory created by "lexicmap index".
  -o, --out-file string   ► Out file, supports the ".gz" suffix ("-" for stdout). (default "-")

Global Flags:
  -X, --infile-list string   ► File of input file list (one file per line). If given, they are
                             appended to files from CLI arguments.
      --log string           ► Log file.
      --quiet                ► Do not print any verbose information. But you can write them to a file
                             with --log.
  -j, --threads int          ► Number of CPU cores to use. By default, it uses all available cores.
                             (default 8)
```

## Examples


```
$ lexicmap utils genomes -d demo.lmi/ \
    | csvtk pretty -t
ref               chunked
---------------   -------
GCF_000148585.2          
GCF_001457655.1          
GCF_900638025.1          
GCF_001096185.1          
GCF_006742205.1          
GCF_001544255.1          
GCF_001027105.1          
GCF_000392875.1          
GCF_009759685.1          
GCF_002949675.1          
GCF_000006945.2          
GCF_000017205.1          
GCF_003697165.2          
GCF_002950215.1          
GCF_000742135.1
```

Extra columns

```
$ lexicmap utils genomes -d demo.lmi/ --extra \
    | csvtk pretty -t
ref               chunked   genome_batch   genome_index
---------------   -------   ------------   ------------
GCF_001457655.1             0              0           
GCF_900638025.1             0              1           
GCF_000148585.2             0              2           
GCF_001096185.1             0              3           
GCF_006742205.1             0              4           
GCF_001027105.1             0              5           
GCF_001544255.1             0              6           
GCF_000392875.1             0              7           
GCF_009759685.1             0              8           
GCF_000006945.2             0              9           
GCF_002949675.1             0              10          
GCF_002950215.1             0              11          
GCF_003697165.2             0              12          
GCF_000742135.1             0              13          
GCF_000017205.1             0              14          
```

