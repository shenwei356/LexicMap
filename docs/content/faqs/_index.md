---
title: FAQs
weight: 60
---

### Does LexicMap support fungi genomes?

Yes. LexicMap mainly supports small genomes including prokaryotic, viral, and plasmid genomes.
Fungi can also be supported, just remember to increase the value of `-g/--max-genome` when running `lexicmap index`,
which is used to skip genomes larger than 15Mb by default.

```
  -g, --max-genome int            ► Maximum genome size. Extremely large genomes (e.g., non-isolate
                                  assemblies from Genbank) will be skipped. (default 15000000)
```

For big and complex genomes, like the human genome which has many repetitive sequences, LexicMap would be slow to align.

### Does LexicMap support short reads?

No. LexicMap only supports long (>=500 bp) reads or gene/genome/virus/plasmid/phage sequences.
However, some short queries can also be aligned.

### How's the hardware requirement?

- For index building. See details [hardware requirement](https://bioinf.shenwei.me/LexicMap/tutorials/index/#hardware-requirements).
- For seaching. See details [hardware requirement](https://bioinf.shenwei.me/LexicMap/tutorials/search/#hardware-requirements).

### Can I extract the matched sequences?

Yes, `lexicmap search` has a flag

```
  -a, --all                            ► Output more columns, e.g., matched sequences.
```

to output CIGAR string, aligned query and subject sequences.

```
18. cigar,    CIGAR string of the alignment                       (optional with -a/--all)
19. qseq,     Aligned part of query sequence.                     (optional with -a/--all)
20. sseq,     Aligned part of subject sequence.                   (optional with -a/--all)
21. align,    Alignment text ("|" and " ") between qseq and sseq. (optional with -a/--all)
```

### Why isn't the pident 100% when aligning with a sequence from the reference genomes?

It happens if there are some degenerate bases (e.g., `N`) in the query sequence.
In the indexing step, all degenerate bases are converted to their lexicographic first bases. E.g., `N` is converted to `A`.
While for the query sequences, we don't convert them.


### Why is LexicMap slow for batch searching?

- LexicMap is mainly designed for sequence alignment with a small number of queries against a database with a huge number (up to 16 million) of genomes.

- `lexicmap search` has a flag `-w/--load-whole-seeds` to load the whole seed data into memory for
faster search.
    - For example, for ~85,000 GTDB representative genomes, searching on an index built with
20,000 masks, the memory would be ~100 GB with default parameters.
- `lexicmap search` also has a flag `--pseudo-align` to only perform pseudo alignment, which is faster and uses less memory.
It can be used in searching with long and divergent query sequences like nanopore long-reads.

{{< button relref="/usage/search"  >}}Click{{< /button >}}  to read more detail of the usage.

