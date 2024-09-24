---
title: FAQs
weight: 60
---
## Table of contents

{{< toc format=html >}}

## Does LexicMap support short reads?

LexicMap is mainly designed for sequence alignment with a small number of queries (gene/plasmid/virus/phage sequences) longer than 200 bp by default.
However, short queries can also be aligned.

If you just want to search long (>1kb) queries for highy similar (>95%) targets, you can build an index with a bigger `-D/--seed-max-desert` (200 by default), e.g.,

    --seed-max-desert 450 --seed-in-desert-dist 150

Bigger values decrease the search sensitivity for distant targets, speed up the indexing
speed, decrease the indexing memory occupation and decrease the index size. While the
alignment speed is almost not affected.


## Does LexicMap support fungi genomes?

Yes. LexicMap mainly supports small genomes including prokaryotic, viral, and plasmid genomes.
**Fungi can also be supported, just remember to increase the value of `-g/--max-genome` when running `lexicmap index`,
which is used to skip genomes larger than 15Mb by default**.

```
  -g, --max-genome int            ► Maximum genome size. Extremely large genomes (e.g., non-isolate
                                  assemblies from Genbank) will be skipped. (default 15000000)
```

Maximum genome size is about 268 Mb (268,435,456). More precisely:

    $total_bases + ($num_contigs - 1) * 1000 <= 268,435,456

as we concatenate contigs with 1000-bp intervals of N’s to reduce the sequence scale to index.

For big and complex genomes, like the human genome (chr1 is ~248 Mb) which has many repetitive sequences, LexicMap would be slow to align.


## How's the hardware requirement?

- For index building. See details [hardware requirement](https://bioinf.shenwei.me/LexicMap/tutorials/index/#hardware-requirements).
- For seaching. See details [hardware requirement](https://bioinf.shenwei.me/LexicMap/tutorials/search/#hardware-requirements).


## Can I extract the matched sequences?

Yes, `lexicmap search` has a flag

```
  -a, --all                            ► Output more columns, e.g., matched sequences. Use this if you
                                       want to output blast-style format with "lexicmap utils 2blast".
```

to output CIGAR string, aligned query and subject sequences.

```
18. cigar,    CIGAR string of the alignment                       (optional with -a/--all)
19. qseq,     Aligned part of query sequence.                     (optional with -a/--all)
20. sseq,     Aligned part of subject sequence.                   (optional with -a/--all)
21. align,    Alignment text ("|" and " ") between qseq and sseq. (optional with -a/--all)
```

An example:

    # Extracting similar sequences for a query gene.

    # search matches with query coverage >= 90%
    lexicmap search -d gtdb_complete.lmi/ b.gene_E_faecalis_SecY.fasta -o results.tsv \
        --min-qcov-per-hsp 90 --all

    # extract matched sequences as FASTA format
    sed 1d results.tsv | awk -F'\t' '{print ">"$5":"$14"-"$15":"$16"\n"$20;}' \
        | seqkit seq -g > results.fasta

    seqkit head -n 1 results.fasta | head -n 3
    >NZ_JALSCK010000007.1:39224-40522:-
    TTGTTCAAGCTATTAAAGAACGCCTTTAAAGTCAAAGACATTAGATCAAAAATCTTATTT
    ACAGTTTTAATCTTGTTTGTATTTCGCCTAGGTGCGCACATTACTGTGCCCGGGGTGAAT


And `lexicmap util 2blast` can help to convert the tabular format to Blast-style format,
see [examples](https://bioinf.shenwei.me/LexicMap/usage/utils/2blast/#examples).

## How can I extract the upstream and downstream flanking sequences of matched regions?

[lexicmap utils subseq](https://bioinf.shenwei.me/LexicMap/usage/utils/subseq/)
can extract subsequencess via genome ID, sequence ID and positions.
So you can use these information from the search result and expand the region positions to extract flanking sequences.



## Why isn't the pident 100% when aligning with a sequence from the reference genomes?

It happens if there are some degenerate bases (e.g., `N`) in the query sequence.
In the indexing step, all degenerate bases are converted to their lexicographic first bases. E.g., `N` is converted to `A`.
While for the query sequences, we don't convert them.


## Why is LexicMap slow for batch searching?

LexicMap is mainly designed for sequence alignment with a small number of queries against a database with a huge number (up to 17 million) of genomes.
There are some ways to improve the search speed of `lexicmap search`.

- **Increasing the concurrency number**
    - Make sure that the value of `-j/--threads` (default: all available CPUs) is ≥ than the number of seed chunk file (default: all available CPUs in the indexing step), which can be found in `info.toml` file, e.g,

          # Seeds (k-mer-value data) files
          chunks = 48

    - Increasing the value of `--max-open-files` (default 512). You might also need to [change the open files limit](https://stackoverflow.com/questions/34588/how-do-i-change-the-number-of-open-files-limit-in-linux).
    - (If you have many queries) Increase the value of `-J/--max-query-conc` (default 12), it will increase the memory.
- **Loading the entire seed data into memoy** (It's unnecessary if the index is stored in SSD)
    - Setting `-w/--load-whole-seeds` to load the whole seed data into memory for faster search. For example, for ~85,000 GTDB representative genomes, the memory would be ~260 GB with default parameters.
- **Returning less results**
    - Setting `-n/--top-n-genomes` to keep top N genome matches for a query (0 for all) in chaining phase. For queries with a large number of genome hits, a resonable value such as 1000 would reduce the computation time.
- Sacrificing accuracy
    - Setting `--pseudo-align` to only perform pseudo alignment, which is slightly faster and uses less memory.
    It can be used in searching with long and divergent query sequences like nanopore long-reads.

{{< button relref="/usage/search"  >}}Click{{< /button >}}  to read more detail of the usage.

