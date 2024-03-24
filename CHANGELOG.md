# Changelog

### v0.3.0 - 2024-03-24

- `lexicmap index`:
    - Generating masks from the top N biggest genomes instead of randomly generating.
- `lexicmap search`:
    - Fix a seed-chaining bug.
    - Fix a target sequence extracting bug.
- New command: `lexicmap utils gen-masks`.

### v0.2.0 - 2024-02-02

- Software architecture and index formats are redesigned to reduce searching memory occupation.
- Indexing: genomes are processed in batches to reduce RAM usage, then indexes of all batches are merged.
- Searching: seeds matching is performed on disk yet it's ultra-fast.

### v0.1.0 - 2024-01-15

- The first release.
- Seed indexing and querying are performed in RAM.
- GTDB r214 with 10k masks: index size 75GB, RAM: 130GB.
