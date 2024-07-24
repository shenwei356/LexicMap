---
title: FAQs
weight: 60
---

### Does LexicMap support short reads?

LexicMap only supports long (>=500 bp, by default) reads or gene/genome/virus/plasmid/phage sequences.
However, some short queries can also be aligned.

For shorter queries like 200 bp, just build an index with smaller `-D/--seed-max-desert`, e.g.,

    --seed-max-desert 200 --seed-in-desert-dist 50

which would generate densers seeds and provide more sensitive results for less similar targets.
The costs are slower (~3X) indexing speed, higher (~2X) indexing memory and bigger (~50%) index size.


### Does LexicMap support fungi genomes?

Yes. LexicMap mainly supports small genomes including prokaryotic, viral, and plasmid genomes.
Fungi can also be supported, just remember to increase the value of `-g/--max-genome` when running `lexicmap index`,
which is used to skip genomes larger than 15Mb by default.

```
  -g, --max-genome int            ► Maximum genome size. Extremely large genomes (e.g., non-isolate
                                  assemblies from Genbank) will be skipped. (default 15000000)
```

Maximum genome size is about 268 Mb (268,435,456). More precisely:

    $total_bases + ($num_contigs - 1) * 1000 <= 268,435,456

as we insert 1000-bp intervals of N’s between contigs to reduce the sequence scale to index.

For big and complex genomes, like the human genome (chr1 is ~248 Mb) which has many repetitive sequences, LexicMap would be slow to align.


### How's the hardware requirement?

- For index building. See details [hardware requirement](https://bioinf.shenwei.me/LexicMap/tutorials/index/#hardware-requirements).
- For seaching. See details [hardware requirement](https://bioinf.shenwei.me/LexicMap/tutorials/search/#hardware-requirements).


### Can I extract the matched sequences?

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

And `lexicmap util 2blast` can help to convert the tabular format to Blast-style format,
see [examples](https://bioinf.shenwei.me/LexicMap/usage/utils/2blast/#examples).

## How can I extract the upstream and downstream flanking sequences of matched regions?

[lexicmap utils subseq](https://bioinf.shenwei.me/LexicMap/usage/utils/subseq/)
can extract subsequencess via genome ID, sequence ID and positions.
So you can use these information from the search result and expand the region positions to extract flanking sequences.



### Why isn't the pident 100% when aligning with a sequence from the reference genomes?

It happens if there are some degenerate bases (e.g., `N`) in the query sequence.
In the indexing step, all degenerate bases are converted to their lexicographic first bases. E.g., `N` is converted to `A`.
While for the query sequences, we don't convert them.


### Why is LexicMap slow for batch searching?

- LexicMap is mainly designed for sequence alignment with a small number of queries against a database with a huge number (up to 16 million) of genomes.

- `lexicmap search` has a flag `-w/--load-whole-seeds` to load the whole seed data into memory for
faster search.
    - For example, for ~85,000 GTDB representative genomes, the memory would be ~170 GB with default parameters.
- `lexicmap search` also has a flag `--pseudo-align` to only perform pseudo alignment, which is slightly faster and uses less memory.
It can be used in searching with long and divergent query sequences like nanopore long-reads.

{{< button relref="/usage/search"  >}}Click{{< /button >}}  to read more detail of the usage.

