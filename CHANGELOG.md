# Changelog

### v0.7.1 - 2025-xx-xx

- `lexicmap index`:
    - Reduce memory usage for huge genomes like Logan contigs.
- `lexicmap search`:
    - **Reduce memory usage, especially for batch searching with long queries**.
    - Slightly improve the search speed.
    - **Support limiting search by TaxId(s)** via `-t/--taxids` or `--taxid-file`.
      Only genomes with descendant TaxIds of the specific ones or themselves are searched,
      in a similar way with BLAST+ 2.15.0 or later versions.
      Negative values are allowed as a black list.
      For example, searching non-Escherichia (561) genera of Enterobacteriaceae (543) family with `-t 543,-561`.
      Users only need to provide NCBI-format taxdump files (`-T/--taxdump`, can also create from
      any taxonomy data with [TaxonKit](https://bioinf.shenwei.me/taxonkit/usage/#create-taxdump))
      and a genome-ID-to-TaxId mapping file (`-G/--genome2taxid`).
      There's no need to rebuild the index.
    - Check if the output file and the log file are the same.

### v0.7.0 - 2025-04-11

Please rebuild the index, as some seeds in the genome end regions were missed during computation.

- `lexicmap index`:
    - **Fix a little bug in seed desert filling** -- forgot to fill the region (a few hundred bases) behind the last seed.
- `lexicmap search`:
    - **Improve seed chaining** -- more accurate for complex anchors.
    - **Improve pseudoalignment in repetitive regions**.
    - Change the default value of `--seed-max-gap` from 200 to 50.

### v0.6.1 - 2025-03-31

- `lexicmap search`:
    - Fix the program hang in the debug mode when no chaining result is returned.
- `lexicmap version`:
    - Do not show commit hash by default.

### v0.6.0 - 2025-03-25

This version is compatible with indexes created by previous versions (requires a one-time, automatic preprocessing),
but rebuilding the index is recommended for more accurate results on short queries (<500bp).
However, indexes created by this version are not compatible with previous versions when the number of batches is <= 512.

- `lexicmap index`:
    - **Change default option values to bring a higher sensitivity for short (<=500, especially <=250) queries,
      faster indexing speed, and faster seed-matching speed<s>, at a cost of slightly larger index</s>**.
        - `-m/--masks`: 40,000 -> 20,000. 
           40k is unnecessary especially for small genomes, where seeds would be very crowded,
           with a big proportion of seed distance being between 0-50 bp.
        - `-D/--seed-max-desert`: 200 -> 100. This provides a smaller seed window guarantee.
    - **Reduce index size by using 3 bytes rather than 4 for saving seed data when the number of batches is <= 512**,
      which requires only 9 (17 minus 8) bits to store the batch index. 
      We also [recommend controlling the number of batches for better performance](https://bioinf.shenwei.me/LexicMap/tutorials/index/#notes-for-indexing-with-large-datasets).
    - **Fix seed desert filling near gap regions**.
- `lexicmap search`:
    - **Improve pseudoalignment to produce longer alignment regions**.
    - **Add 3 extra columns: `cls`, `evalue` and `bitscore`**, and a new option `-e/--max-evalue`.
    - Reduce memory usage.
    - Remove flag `--pseudo-align`.
    - Add a progress bar for `--debug`.
- `lexicmap utils seed-pos`:
    - Change default option values of sliding window.

### v0.5.0 - 2024-12-18

This version is compatible with indexes created by LexicMap v0.4.0, but rebuilding the index is recommended for more accurate results.

- New commands:
    - **`lexicmap utils remerge`: Rerun the merging step for an unfinished index**.
- `lexicmap index`:
    - **Big genomes with thousands of contigs (big yet fragmented assemblies) are automatically split into multiple chunks, and alignments from these chunks will be merged.**
    - **Change the default value of `--partitions` from 1024 to 4096, which increases the seed-matching speed at the cost of 2 GiB more memory occupation**.
      For existing lexicmap indexes, just run `lexicmap utils reindex-seeds --partitions 4096` to re-create seed indexes.
    - **Do not save seeds of low-complexity**.
    - Fix high memory usage in writing seed data.
    - Change the default value of `-c/--chunks` from all available CPUs to the value of `-j/--threads`.
    - Change the default value of `--max-open-files` from 512 to 1024.
    - Add a new flag `--debug`.
- `lexicmap search`:
    - **Improving chaining, pseudoalignment, and alignment for highly repetitive sequences**.
    - **More accurate chaining score with better chaining of overlapped anchors, this produces more accurate results with `-n/--top-n-genomes`**: 
         - Merging two overlapped non-gapped anchors into a longer one.
         - For these with gaps, only the non-overlapped part of the second anchor is used to compute the weight.
         - Using the score of the best chain (rather than the sum) for sorting genomes when using `-n`.
    - Fix positions and alignment texts for queries with highly repetitive sequences in end regions. [#9](https://github.com/shenwei356/LexicMap/issues/9)
    - Skip seeds of low-complexity.
    - Change the default value of `--max-open-files` from 512 to 1024.
    - Change the default value of `--align-band` from 50 to 100.
    - Improve the speed of anchor deduplication, genome information extraction, and result ordering.
    - Improve the speed of chaining for long queries.
    - Improve the speed of seed matching when using `-w/--load-whole-seeds`.
    - **Improve the speed of alignment, and reduce the memory usage**.
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
