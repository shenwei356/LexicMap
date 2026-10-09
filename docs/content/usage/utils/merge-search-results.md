---
title: lexicmap utils merge-search-results
linkTitle: merge-search-results
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
  5. Flag -n/--top-n-genomes ranks genomes by their highest alignment
     bitscore * pident, retaining all ties at the cutoff score.
     Unlike 'lexicmap search -n', this does not use chaining scores.
     Selected genome results are buffered in memory to count retained genome hits.

Usage:
  lexicmap utils merge-search-results [flags] 

Flags:
  -b, --buffer-size string   ► Size of buffer, supported unit: K, M, G. You need increase the value
                             when "bufio.Scanner: token too long" error reported (default "20M")
  -h, --help                 help for merge-search-results
  -o, --out-file string      ► Out file, supports the ".gz" suffix ("-" for stdout). (default "-")
  -q, --query string         ► Query ID to merge
  -n, --top-n-genomes int    ► Keep the top N genome matches by the highest alignment bitscore *
                             pident per genome, including all matches tied at the cutoff score (0 for
                             all). This filters alignment results, not chaining scores.

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

Building 3 indexes with demo data.

```text
lexicmap index -O demo_part1.lmi -X <(ls refs/*.fa.gz | awk 'NR<=5')
lexicmap index -O demo_part2.lmi -X <(ls refs/*.fa.gz | awk 'NR>=6 && NR<=10')
lexicmap index -O demo_part3.lmi -X <(ls refs/*.fa.gz | awk 'NR>10')
```

Searching a query in these 3 indexes.

```text
query=bench/b.gene_E_coli_16S.fasta

# in practice, one can search in multiple cluster nodes
for index in demo_*.lmi; do
    lexicmap search $query -d $index -o t.lexicmap@$index.tsv.gz --debug
done
```

Merge search results

```text
lexicmap utils merge-search-results -o t.lexicmap.tsv.gz t.lexicmap@*.tsv.gz

22:41:03.963 [INFO] 15 genome hits merged from 3 files for query: NC_000913.3:4166659-4168200
```

To retain the top N genomes after merging, use `-n/--top-n-genomes`:

```text
lexicmap utils merge-search-results -n 100 -o t.lexicmap.top100.tsv.gz t.lexicmap@*.tsv.gz
```

This ranks genomes by the highest alignment `bitscore * pident` among their HSPs,
using the values in the input TSV files. All genomes tied at the Nth score are retained,
along with all of their HSPs, so more than N genomes may be returned. The `hits` column
reports the number of retained genomes. The default, `-n 0`, keeps all genomes.

`lexicmap search -n` filters candidates by chaining score before alignment.
The merge command filters completed alignment results and cannot recover candidates
discarded during the original searches.

With `-n` enabled, selected genome results are buffered in memory until their count
is known, then written with the correct `hits` value. No temporary files are used.
This also supports standard input. The buffer includes all ties at the cutoff score,
so its size depends on the number and size of retained alignments, not only N.

If some files contain search results of multiple queries, then specify one query to merge.

```text
# search with multiple queries
for index in demo_*.lmi; do
    lexicmap search bench/*.fasta -d $index -o t.lexicmap2@$index.tsv.gz --debug
done

# this would failed
lexicmap utils merge-search-results -o t.lexicmap2.tsv.gz t.lexicmap2@*.tsv.gz

22:46:07.807 [ERRO] inconsistent queries: 'blaLAQ__NG_076677.1' in file 't.lexicmap2@demo_part2.lmi.tsv.gz' and 'blaSHV__NG_242606.1' in file 't.lexicmap2@demo_part1.lmi.tsv.gz. Please specify one query with flag -q/--query

# specify one query

lexicmap utils merge-search-results -o t.lexicmap2.tsv.gz t.lexicmap2@*.tsv.gz \
    -q NC_000913.3:4166659-4168200
    
22:49:01.076 [INFO] 15 genome hits merged from 3 files for query: NC_000913.3:4166659-4168200
```
