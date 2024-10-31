# Changelog

### v0.5.0 - 2024-xx-xx

This version generates the same alignment results as v0.4.0.

- New commands:
    - **`lexicmap utils remerge`: Rerun the merging step for an unfinished index**.
- `lexicmap index`:
    - **Big genomes with thousands of contigs (big yet fragmented assemblies) are automatically split into multiple chunks, and alignments from these chunks will be merged.**
    - **Change the default value of `--partitions` from 1024 to 4096, which increases the seed-matching speed at the cost of 2 GiB more memory occupation**.
      For existing lexicmap indexes, just run `lexicmap utils reindex-seeds --partitions 4096` to re-create seed indexes.
    - Change the default value of `-c/--chunks` from all available CPUs to the value of `-j/--threads`.
    - Change the default value of `--max-open-files` from 512 to 1024.
    - Add a new flag `--debug`.
- `lexicmap search`:
    - Fix positions and alignment texts for queries with highly repetitive sequences in flanking regions. [#9](https://github.com/shenwei356/LexicMap/issues/9)
    - Automatically adjust arguments `--seed-max-dist` and `--align-ext-len` for indexes with a smaller contig interval size.
    - More accurate `-n/--top-n-genomes`, and add new help message.
    - Change the default value of `--max-open-files` from 512 to 1024.
    - Improve the speed of anchor deduplication, genome information extraction, and result ordering.
    - Improve the speed of seed matching when using `-w/--load-whole-seeds`.
    - Improve the speed of alignment, and reduce the memory usage.
    - Remain compatible after the change of `lexicmap index`.
    - Add a new flag `--debug`.
- `lexicmap utils genomes`:
    - Do not sort genome ids.
    - Add a header line and add another column to show if the reference genome is chunked.
- `lexicmap utils subseq`:
    - Remain compatible after the change of `lexicmap index`.
- `lexicmap utils seed-pos`:
    - Remain compatible after the change of `lexicmap index`, while histograms are plotted separately for multiple genome chunks.
- `lexicmap utils reindex-seeds`:
    - Change the default value of `--partitions` from 1024 to 4096.

### v0.4.0 - 2024-08-15

- New commands:
    - **`lexicmap utils 2blast`: Convert the default search output to blast-style format**.
- `lexicmap index`:
    - **Support suffix matching of seeds, now seeds are immune to any single SNP!!!**, at the cost of doubled seed data.
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
    - Use buffered reader for seeds file reading.
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

### v0.3.0 - 2024-05-14

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

### v0.2.0 - 2024-02-02

- Software architecture and index formats are redesigned to reduce searching memory occupation.
- Indexing: genomes are processed in batches to reduce RAM usage, then indexes of all batches are merged.
- Searching: seeds matching is performed on disk yet it's ultra-fast.

### v0.1.0 - 2024-01-15

- The first release.
- Seed indexing and querying are performed in RAM.
- GTDB r214 with 10k masks: index size 75GB, RAM: 130GB.
