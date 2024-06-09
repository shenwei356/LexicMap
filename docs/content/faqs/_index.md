---
title: FAQs
weight: 60
---

### Does LexicMap support fungi genomes?

Yes. LexicMap mainly supports small genomes including prokaryotic, viral, and plasmid genomes.
Fungi can also be supported.

For big and complex genomes, like the human genome which has many repetitive sequences, LexicMap would be slow to align.

### Does LexicMap support short reads?

No. LexicMap only supports long (>=500 bp, >=1kb is prefered) reads or gene/genome/virus/plasmid/phage sequences.
However, some short queries can also be aligned.

### How's the hardware requirement?

For index building. See details [hardware requirement](https://bioinf.shenwei.me/LexicMap/tutorials/index/#hardware-requirements).
- More CPUs would accelerate indexing.
- The memory occupation is linear with:
    - The size of the genome batch.
    - The number of LexicHash masks.

For seaching. See details [hardware requirement](https://bioinf.shenwei.me/LexicMap/tutorials/search/#hardware-requirements).
- More CPUs would accelerate searching.
- The memory occupation mainly depends on the length of queries, the number of hits, and the distance between query and target sequences.

### Can I extract the matched sequences?

Yes, `lexicmap search` has a flag

```
  -a, --all                            ► Output more columns, e.g., matched sequences.
```

to ouput aligned query and subject sequences.

```
18. qseq,     Aligned part of query sequence.   (optional with -a/--all)
19. sseq,     Aligned part of subject sequence. (optional with -a/--all)
```

### What is not the pident 100% when aligning with a sequence from the reference genomes?

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
It can be used insearching with long and divergent query sequences like nanopore long-reads.

{{< button relref="/usage/search"  >}}Click{{< /button >}}  to read more detail of the usage.

