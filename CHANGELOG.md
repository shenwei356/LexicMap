# Changelog

### v0.3.0 - 2024-04-xx

- `lexicmap index`:
    - Generating masks from the top N biggest genomes instead of randomly generating.
    - Fix a concurrency bug between genome data writing and k-mer-value data collecting.
    - Change the format of k-mer-value index file.
- `lexicmap search`:
    - Fix a seed-chaining bug.
    - Fix a target sequence extracting bug.
    - Add output of query start and end positions.
    - Add a new flag `-w/--load-whole-seeds` for loading the whole seed data into memory for faster search.
    - Keep index of genome data in memory.
- `lexicmap utils kmers`:
    - Fix a little bug, wrong number of k-mers for the second k-mer in each k-mer pair.
- New commands:
    - `lexicmap utils gen-masks` for generating masks from the top N largest genomes.
    - `lexicmap utils seed-pos` for extracting seed positions via reference names.
    - `lexicmap utils reindex-seeds` for recreating indexes of k-mer-value (seeds) data.

### v0.2.0 - 2024-02-02

- Software architecture and index formats are redesigned to reduce searching memory occupation.
- Indexing: genomes are processed in batches to reduce RAM usage, then indexes of all batches are merged.
- Searching: seeds matching is performed on disk yet it's ultra-fast.

### v0.1.0 - 2024-01-15

- The first release.
- Seed indexing and querying are performed in RAM.
- GTDB r214 with 10k masks: index size 75GB, RAM: 130GB.
