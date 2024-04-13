---
title: Releases
weight: 30
---



## Latest Version

### v0.3.0

[v0.3.0](https://github.com/shenwei356/LexicMap/releases/tag/v0.3.0) - 2024-04-xx [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/LexicMap/v0.3.0/total.svg)](https://github.com/shenwei356/LexicMap/releases/tag/v0.3.0)

- `lexicmap index`:
    - **Generate masks from the top N biggest genomes instead of randomly generation**.
    - **Use longer (1000bp N's, previous: k-1) intervals between contigs**.
    - Fix a concurrency bug between genome data writing and k-mer-value data collecting.
    - Change the format of k-mer-value index file, and fix the computation of index partitions.
- `lexicmap search`:
    - **Better support of long queries**.
    - **Add a new flag `-w/--load-whole-seeds` for loading the whole seed data into memory for faster search**.
    - **Parallelize alignment in each query**, so it's faster for a single query.
    - Change output format.
    - Add output of query start and end positions.
    - Fix a seed-chaining bug.
    - Fix a target sequence extracting bug.
    - Keep indexes of genome data in memory.
- `lexicmap utils kmers`:
    - Fix a little bug, wrong number of k-mers for the second k-mer in each k-mer pair.
- New commands:
    - `lexicmap utils gen-masks` for generating masks from the top N largest genomes.
    - `lexicmap utils seed-pos` for extracting seed positions via reference names.
    - `lexicmap utils reindex-seeds` for recreating indexes of k-mer-value (seeds) data.


{{< hint type=note >}}

- Please run `lexicmap version` to check update !!!
- Please run `lexicmap autocompletion` to update shell autocompletion script !!!
{{< /hint >}}

## Previous Versions


### v0.2.0

[v0.2.0](https://github.com/shenwei356/LexicMap/releases/tag/v0.2.0) - 2024-02-02 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/LexicMap/v0.2.0/total.svg)](https://github.com/shenwei356/LexicMap/releases/tag/v0.2.0)

- Software architecture and index formats are redesigned to reduce searching memory occupation.
- Indexing: genomes are processed in batches to reduce RAM usage, then indexes of all batches are merged.
- Searching: seeds matching is performed on disk yet it's ultra-fast.

### v0.1.0

[v0.1.0](https://github.com/shenwei356/LexicMap/releases/tag/v0.1.0) - 2024-01-15 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/LexicMap/v0.1.0/total.svg)](https://github.com/shenwei356/LexicMap/releases/tag/v0.1.0)

- The first release.
- Seed indexing and querying are performed in RAM.
- GTDB r214 with 10k masks: index size 75GB, RAM: 130GB.
