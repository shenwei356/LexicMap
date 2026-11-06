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
     If not, please specify one query with the flag -q/--query.
     Query selection also applies when only one input file is given.
     Merging without -q rejects mixed query IDs, including within one subject group.
     All records for the selected query must have the same positive query length.
     Hit counts must be positive integers, identities finite values in [0,100],
     and bit scores integers (negative scores are supported).
  2. We assume that genome IDs are distinct across all indexes.
  3. One or more input files are accepted, via positional parameters
     and/or a file list via the flag -X/--infile-list.
  4. Both the default 20- and 24-column formats are supported,
     and the first retained result selects the output format.
     In 24-column output, missing alignment fields are left empty.
     In 20-column output, alignment fields are omitted.
  5. Flag -n/--top-n-genomes ranks genomes by their highest alignment
     bitscore * pident, retaining all ties at the cutoff score.
     Unlike 'lexicmap search -n', this does not use chaining scores.
     Ranking uses the rounded TSV values, which can change the original search order.
     Only the top N genomes and cutoff ties are buffered for final output.
  6. Reader buffers start at up to 64 KiB and grow for long alignment rows up to
     -b/--buffer-size (default 20M). Increase -b if a row exceeds this limit.
  7. Empty inputs are skipped with a warning naming each file (unless --quiet).
     If all inputs are empty, the command succeeds without writing results.
     A completed search writes a header even when no matches are found.
     Header-only inputs are accepted.
  8. Inputs are read directly when their count fits --max-open-files minus 9
     reserved descriptors. Otherwise, decompressed inputs are first copied to one
     file in the system temporary directory (TMPDIR on Unix), then merged.
     The temporary file is removed on completion or error.
     Allow temporary disk space for the decompressed inputs, including with -n.

Usage:
  lexicmap utils merge-search-results [flags]

Flags:
  -b, --buffer-size string   ► Size of buffer, supported unit: K, M, G. You need increase the value
                             when "bufio.Scanner: token too long" error reported (default "20M")
  -h, --help                 help for merge-search-results
      --max-open-files int   ► File budget including 9 reserved descriptors (minimum: 10). Read inputs
                             directly when they fit; otherwise use one temporary file. Set "ulimit -n"
                             at least this high. (default 1024)
  -o, --out-file string      ► Out file, supports the ".gz" suffix ("-" for stdout). (default "-")
  -q, --query string         ► Query ID to merge
  -n, --top-n-genomes int    ► Keep the top N genome matches by the highest alignment bitscore *
                             pident per genome, including all matches tied at the cutoff score (0 for
                             all). This filters alignment results, not chaining scores.

Global Flags:
      --force-gzip-output    ► Gzip-compress result output regardless of the filename, including stdout.
  -X, --infile-list string   ► File of input file list (one file per line). If given, they are
                             appended to files from CLI arguments.
      --log string           ► Log file.
      --quiet                ► Do not print any verbose information. But you can write them to a file
                             with --log.
  -j, --threads int          ► Number of CPU cores to use. By default, it uses all available cores.
                             (default 16)
```

## Examples

Without `-q`, all records must have the same query ID, even when only one file
is supplied and `-n` is left at its default of `0`. Use `-q` to select one query
from a file containing results for multiple queries.

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
is known, then written with the correct `hits` value.
This also supports standard input. The buffer includes all ties at the cutoff score,
so its size depends on the number and size of retained alignments, not only N.

`--max-open-files` defaults to `1024`, including 9 reserved descriptors for standard
streams, runtime files, logging, output, and temporary storage. Inputs are read
directly when their count is at most `--max-open-files - 9`. For example, 200 inputs
are read directly with the default budget, while `--max-open-files 64` allows 55
direct inputs. Set `ulimit -n` at least as high as the configured budget.

When the input count exceeds the available budget, all inputs are copied to one
file in the system temporary directory (`TMPDIR` on Unix) before merging. Allow
disk space for all decompressed input data, including with `-n`. The file is removed
on completion or error.

If some files contain search results for multiple queries, specify one query to merge.

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
