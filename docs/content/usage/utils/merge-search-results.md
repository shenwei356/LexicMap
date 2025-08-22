---
title: merge-search-results
weight: 1
---

## Usage

```plain
$ lexicmap utils merge-search-results -h
Merge a query's search results from multiple indexes

Attention:
  1. These search results should come from the same ONE query.
     If not, please specify one query with the flag -q/--query
  2. We assume that genome IDs are distinct across all indexes.
  3. One or more input files are accepted, via positional parameters
     and/or a file list via the flag -X/--infile-list.
  4. Both the default 20- and 24-column formats are supported,
     and formats better be consistent across all input files.
     If not, the output format would be the one with a valid record.
  5. The column "hits" in the output will be set to 0.

Usage:
  lexicmap utils merge-search-results [flags] 

Flags:
  -b, --buffer-size string   ► Size of buffer, supported unit: K, M, G. You need increase the value
                             when "bufio.Scanner: token too long" error reported (default "20M")
  -h, --help                 help for merge-search-results
  -o, --out-file string      ► Out file, supports the ".gz" suffix ("-" for stdout). (default "-")
  -q, --query string         ► Query ID to merge

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

All search results belong to one query.

```text
lexicmap util merge-search-results -o query.lexicmap.tsv.gz \
    query.lexicmap@genbank.tsv.gz \
    query.lexicmap@refseq.tsv.gz \
    query.lexicmap@allthebacteria.tsv.gz
```

Some files contain search results of multiple queries, then specify one query to merge.

```text
lexicmap util merge-search-results -o query.lexicmap.tsv.gz \
    query.lexicmap@genbank.tsv.gz \
    query.lexicmap@refseq.tsv.gz \
    query.lexicmap@allthebacteria.tsv.gz \
    --query NC_000913.3:4166659-4168200
```


