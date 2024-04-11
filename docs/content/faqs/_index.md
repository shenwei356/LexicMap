---
title: FAQs
weight: 60
---

### Does LexicMap support fungi genomes?

No. LexicMap only supports short genomes including prokaryotic, viral, and plasmid genomes.

### Does LexicMap support short reads?

No. LexicMap only supports long (>=500 bp) reads or gene/genome/plasmid sequences.
However, some short queries can also be aligned.

### How's the hardware requirement?

For index building.
- More CPUs would accelerate indexing.
- The memory occupation is linear with
    - The size of the genome batch
    - The number of LexicHash masks.

For seaching.
- More CPUs would accelerate searching.
- The memory occupation depends on the number of hits.

### Why is LexicMap slow for batch searching?

- LexicMap is designed to scale to alignment against a large number of genomes.

- `lexicmap search` has a flag `-w/--load-whole-seeds` to load the whole seed data into memory for
faster search.
    - For example, for ~85,000 GTDB representative genomes, searching on an index built with
20,000 masks, the memory would be ~75 GB.
