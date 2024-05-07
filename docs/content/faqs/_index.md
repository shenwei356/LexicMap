---
title: FAQs
weight: 60
---

### Does LexicMap support fungi genomes?

Yes. LexicMap mainly supports small genomes including prokaryotic, viral, and plasmid genomes.
Fungi can be supported with more masks.

### Does LexicMap support short reads?

No. LexicMap only supports long (>=500 bp, >=1kb is prefered) reads or gene/genome/virus/plasmid/phage sequences.
However, some short queries can also be aligned.

### How's the hardware requirement?

For index building. See details [hardware requirement](https://bioinf.shenwei.me/LexicMap/tutorials/index/#hardware-requirements).
- More CPUs would accelerate indexing.
- The memory occupation is linear with
    - The size of the genome batch.
    - The number of LexicHash masks.

For seaching. See details [hardware requirement](https://bioinf.shenwei.me/LexicMap/tutorials/search/#hardware-requirements).
- More CPUs would accelerate searching.
- The memory occupation mainly depends on the length of queries, the number of hits, and band size in backtracking.

### Why is LexicMap slow for batch searching?

- LexicMap is mainly designed for sequence alignment with a small number of queries against a database with a huge number (up to 16 million) of genomes.

- `lexicmap search` has a flag `-w/--load-whole-seeds` to load the whole seed data into memory for
faster search.
    - For example, for ~85,000 GTDB representative genomes, searching on an index built with
20,000 masks, the memory would be ~100 GB with default parameters.
    -  {{< button relref="/usage/search"  >}}Click{{< /button >}}  to read more detail of the usage.
