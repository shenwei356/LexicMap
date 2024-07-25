---
title: Releases
weight: 30
---



## Latest version

### v0.4.0

[v0.4.0](https://github.com/shenwei356/LexicMap/releases/tag/v0.4.0) - 2024-07-xx [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/LexicMap/v0.4.0/total.svg)](https://github.com/shenwei356/LexicMap/releases/tag/v0.4.0)

- New commands:
    - **`lexicmap utils 2blast`: Convert the default search output to blast-style format**.
- `lexicmap index`:
    - **Support suffix matching of seeds, now seeds are immune to any single SNPs!!!**, at the cost of doubled seed data.
    - **Better sketching desert filling for highly-repetitive regions**.
    - **Change the default value of `--seed-max-desert` from 900 to 200 to increase alignment sensitivity**.
    - **Mask gap regions (N's)**.
    - Fix skipping interval regions by further including the last k-1 bases of contigs.
    - Fix a bug in indexing small genomes.
    - Change the default value of `-b, --batch-size` from 10,000 to 5,000.
    - Improve lexichash data structure.
    - Write and merge seed data in parallel, new flag `-J/--seed-data-threads`.
    - Improve the log.
- `lexicmap search`:
    - **Fix chaining for highly-repetitive regions**.
    - **Perform more accurate alignment with [WFA](https://github.com/shenwei356/wfa)**.
    - Fix object recycling and reduce memory usage.
    - Fix alignment against genomes with many short contigs.
    - Fix early quit when meeting a sequence shorter than k.
    - Add a new option `-J/--max-query-conc` to limit the miximum number of concurrent queries,
      with a default valule of 12 instead of the number of CPUs, which reduces the memory usage
      in batch searching.
    - Result format:
        - Cluster alignments of each target sequence.
        - Remove the column `seeds`.
        - Add columns `gaps`, `cigar`, `align`, which can be reformated with `lexicmap utils 2blast`.
- `lexicmap utils kmers`:
    - Fix the progress bar.
    - Fix a bug where some masks do not have any k-mer.
    - Add a new column `prefix` to show the length of common prefix between the seed and the probe.
    - Add a new column `reversed` to indicate if the k-mer is reversed for suffix matching.
- `lexicmap utils masks`:
    - Add the support of only outputting a specific mask.
- `lexicmap utils seed-pos`:
    - New columns: `sseqid` and `pos_seq`.
    - More accurate seed distance.
    - Add histograms of numbers of seed in sliding windows.
- `lexicmap utils subseq`:
    - Fix a bug when the given end position is larger than the sequence length.
    - Add the strand ("+" or "-") in the sequence header.

{{< hint type=note >}}

- Please run `lexicmap version` to check update !!!
- Please run `lexicmap autocompletion` to update shell autocompletion script !!!
{{< /hint >}}

## Previous versions

### v0.3.0

[v0.3.0](https://github.com/shenwei356/LexicMap/releases/tag/v0.3.0) - 2024-05-14 [![Github Releases (by Release)](https://img.shields.io/github/downloads/shenwei356/LexicMap/v0.3.0/total.svg)](https://github.com/shenwei356/LexicMap/releases/tag/v0.3.0)

- `lexicmap index`:
    - **Better seed coverage by filling sketching deserts**.
    - **Use longer (1000bp N's, previous: k-1) intervals between contigs**.
    - Fix a concurrency bug between genome data writing and k-mer-value data collecting.
    - Change the format of k-mer-value index file, and fix the computation of index partitions.
    - Optionally save seed positions which can be outputted by `lexicmap utils seed-pos`.
- `lexicmap search`:
    - **Improved seed-chaining algorithm**.
    - **Better support of long queries**.
    - **Add a new flag `-w/--load-whole-seeds` for loading the whole seed data into memory for faster search**.
    - **Parallelize alignment in each query**, so it's faster for a single query.
    - **Optional outputing matched query and subject sequences**.
    - 2-5X searching speed with a faster masking method.
    - Change output format.
    - Add output of query start and end positions.
    - Fix a target sequence extracting bug.
    - Keep indexes of genome data in memory.
- `lexicmap utils kmers`:
    - Fix a little bug, wrong number of k-mers for the second k-mer in each k-mer pair.
- New commands:
    - `lexicmap utils gen-masks` for generating masks from the top N largest genomes.
    - `lexicmap utils seed-pos` for extracting seed positions via reference names.
    - `lexicmap utils reindex-seeds` for recreating indexes of k-mer-value (seeds) data.
    - `lexicmap utils genomes` for list genomes IDs in the index.

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

