# Changelog

### v0.10.0 - 2026-xx-xx

There is a small change in the seed computation, but re-indexing is unnecessary.

- `lexicmap`:
    - All commands log the invocation, working directory, software version, build commit (when available), and project URL at startup, and total elapsed time on completion.
    - All commands check that the log file differs from the output file before opening it.
    - Progress bars in `search`, `genome search`, and `genome compare` respect `--quiet`, including when `--debug` is enabled.
- New commands:
    - **`lexicmap utils reindex-seeds2`: Recreate adaptive two-level indexes of seed data to speed up prefix matching in large indexes.**
    - **`lexicmap genome search`: Search genomes against an index, with ANI and AF computed**.
    - **`lexicmap genome pair`: Find similar genome pairs in the index**.
    - **`lexicmap genome compare`: Compare genome pairs and compute ANI and AF**.
    - `lexicmap utils genome-details`: Extract or view genome details in the index.
    - `lexicmap utils genome-seqs`: Extract all sequences of a given genome.
- `lexicmap index`:
    - **Faster indexing (30-45% less time) and ~30% lower memory usage by optimizing seed computation and merging**.
    - **Fixed a strand bias in seed computation that skipped some negative-strand k-mers during the first round of probe capture (k-mer masking)**.
      This caused more k-mers to be captured on the positive strand, but had a negligible effect
      on alignment sensitivity after seed deserts were filled. Only a small fraction of seeds change
      when rebuilding an index.
      To preserve compatibility, the old algorithm is used for index formats v3.0 to v3.4.
    - **Add `-2/--seed-index2` to write adaptive two-level seed indexes alongside the final seed files**,
      with `--seed-index2-threshold` (default `8K`, minimum `4K`) controlling the minimum block size.
      These indexes speed up prefix matching in large indexes by using longer prefixes to locate seeds within large data blocks.
      Overall gains are most noticeable for batch queries with `lexicmap search` when `-n/--top-n-genomes` limits downstream alignments,
      making seed matching a larger share of runtime. Without this limit, sequence alignment often dominates runtime.
    - Changed the default value of `-g/--max-genome` from 15Mb to 20Mb,
      as a few genomes in RefSeq are larger than 15Mb (e.g., GCA_051525975.1).
    - Fixed data races in parallel seed computation and progress reporting.
    - Fixed a panic when a custom mask file contained a different number of masks from `-m/--masks`.
      The mask file now determines the number of masks.
    - Keep batch merging within `--max-open-files`, including output files and reserved descriptors, by limiting input groups and merge threads. Applies to `utils remerge` as well.
- `lexicmap search`:
    - **Faster searching speed, mainly for batch querying with `-n/--top-n-genomes`**.
        - **Parallelize anchor generation and collection from seed-matching results to reduce collector bottlenecks for high-hit queries**.
        - Optimize chaining to reduce memory use and garbage collection overhead.
        - Release seed anchors after chaining, and retain only output fields after alignment. Inspired by @d-callan's proposal in [#38](https://github.com/shenwei356/LexicMap/pull/38).
        - Faster pseudoalignment for long queries.
    - **Add two flags to limit the memory usage when searching in huge indexes** such as Logan Project,
      where the large number of seed matches can cause memory exhaustion during seed collection.
      See how to [trade off speed for memory usage](https://bioinf.shenwei.me/LexicMap/tutorials/search/#trade-off-speed-for-memory-usage).
        - **Add `--max-seed-memory` (default `0`, disabled) to stream seed data and spill sorted anchors before their collection buffers exceed the budget divided among query slots**. Use this and set a large value (such as 1/2 to 3/4 of the available free RAM). See the help message for more details. Inspired by @d-callan's proposal in [#37](https://github.com/shenwei356/LexicMap/pull/37).
        - **Add `--max-align-result-memory` (default `1G`) to spill large `-a/--all` output fields to temporary files once their global in-memory budget is exhausted**. Useful when long queries or many alignment hits consume substantial memory for retained CIGAR strings, aligned sequences, and alignment text, especially with concurrent queries. Inspired by @d-callan's proposal in [#38](https://github.com/shenwei356/LexicMap/pull/38).
    - **Updated the WFA implementation to follow standard end-to-end global alignment semantics and WFA2-compatible tie-breaking**.
      This may slightly change CIGAR strings and derived statistics for some low-similarity hits.
      In limited tests, the resulting alignments tended to be slightly shorter and contain fewer gaps.
    - Keep all genome matches tied at the Nth chaining score when using `-n/--top-n-genomes`; the number of retained candidates may exceed N.
    - Flag `-T/--taxdump`: set a default value `<index path>/taxdump`.
    - Flag `-G/--genome2taxid`: set a default value `<taxdump path>/taxid.map`.
    - Added new flags `-g/--show-genome-name` and `-s/--show-species-name` to add the taxonomic/species name as a prefix to the `sgenome` field.
    - Added a new flag `--show-sseq-idx` to add 1-based genome chunk and subject-sequence index prefixes to `sseqid` values.
    - Added a new flag `--show-avg-qual` to add the average quality of the aligned region as a suffix to the `alenHSP` field.
    - Fixed TaxId filtering with only negative TaxIds, which discarded the first seed hit from each allowed genome.
    - Fixed a data race when extending the pseudoalignment region.
- `lexicmap index, lexicmap utils edit-genome-ids/genome-details`:
    - Truncate genome/sequence IDs longer than 65,535 characters.
- `lexicmap utils subseq`:
    - Fixed a hang with `-f/--search-result` when `--max-open-files` was smaller than the number of genome batches.
- `lexicmap utils merge-search-results`:
    - Add `-n/--top-n-genomes` to retain genomes by their highest alignment `bitscore * pident`, including all ties at the cutoff score. Unlike `search -n`, this filters alignment results rather than chaining scores.
- `lexicmap utils seed-pos`:
    - Skip empty seed results and finish reading result fields before recycling them.
    - Fixed data races in seed-position reader creation and error handling.
- `lexicmap util kmers`:
    - Faster speed for printing all seed data (`--mask 0`).
- `lexicmap utils reindex-seeds`:
    - Remove obsolete `.idx15` files after successfully rebuilding their single-level seed indexes.
    - Added progress bars to all seed chunks.

### v0.9.0 - 2026-03-13

- New commands:
    - **`lexicmap utils 2sam`: Convert the default search output to SAM format** ([#26](https://github.com/shenwei356/LexicMap/issues/26)).
       Attention: This command requires search results generated by the current LexicMap version.
- `lexicmap index`:
    - **Added a new flag `--soft-masking` to support soft-masked genomes (usually eukaryotic genomes)**. 
      Lowercase bases in soft-masked low-complexity regions will be treated as A's and won't be seeded,
      while they will be saved for base-level alignment.
      The value of this flag is appended to the index information file `info.toml`.
    - Added a new option `--max-kmer-freq`. If a mask captures the same k-mer at more than *N* positions in a genome,
      only the first N positions will be retained. This option may reduce search sensitivity, but it can be useful
      when simply checking whether a query matches any position in a genome that contains many tandem repeat sequences.
      The value of this flag is appended to the index information file `info.toml`.
    - Added sequence validity checking. [#30](https://github.com/shenwei356/LexicMap/issues/30)
- `lexicmap search`:
    - **Fixed the computation of bitscore and evalue**. Previous values were slightly overestimated.
    - **Fixed CIGAR to follow the SAM spec**. Previously, 'D' and 'I' were inverted.
    - Fixed merging search results from genome chunks. Some were not merged.
    - **New flag `-N/--top-n-chains` for keeping the top N chains in a genome for the query** (0 for all) in the chaining phase.
      It reduces search time when one needs only the most similar matches in a genome.
    - **Improved seed chaining speed and reduced the memory usage, especially for genomes with lots of repeat sequences**.
- `lexicmap utils subseq`:
    - Fixed a concurrency bug when using search result as the input.
- `lexicmap utils genomes`:
    - Added a new flag `-e/--extra` to show more information.
- `lexicmap utils masks`:
    - Changed the default values of `-m/--masks` and `-p/--prefix` to match those in `lexicmap index`.

### v0.8.1 - 2025-12-19

- `lexicmap search`:
    - Fix a runtime panic. [#23](https://github.com/shenwei356/LexicMap/issues/23).

### v0.8.0 - 2025-09-10

No changes to the index format (see [Index format changelog](https://bioinf.shenwei.me/LexicMap/tutorials/index/#index-format-changelog)).

- New commands:
    - **`lexicmap utils merge-search-results`: Merge a query's search results from multiple indexes**.
    - **`lexicmap utils edit-genome-ids`: Edit genome IDs in the index via a regular expression**.
      It's helpful when users forgot to use the flag `-N/--ref-name-regexp`
      to extract the genome ID from the sequence file during indexing.
      This command helps fix this without rebuilding the index.
- `lexicmap index`:
    - **Significantly reduce the memory usage (by up to 25%)** in the merge step.
      Also reduce memory usage for huge datasets, such as long reads or contigs in the [Logan project](https://github.com/IndexThePlanet/Logan).
- `lexicmap search`:
    - **Reduce memory usage, particularly for batch searching (by up to 50%)**.
    - **Improve search speed, mainly for batch searching**.
    - **Support limiting search by TaxId(s)** via `-t/--taxids` or `--taxid-file`.
      Only genomes with descendant TaxIds of the specific ones or themselves are searched,
      in a similar way to BLAST+ 2.15.0 or later versions.
      Negative values are allowed as a blacklist.
      For example, searching non-Escherichia (561) genera of the Enterobacteriaceae (543) family with `-t 543,-561`.
      Users only need to provide NCBI-format taxdump files (`-T/--taxdump`, can also create from
      any taxonomy data with [TaxonKit](https://bioinf.shenwei.me/taxonkit/usage/#create-taxdump))
      and a genome-ID-to-TaxId mapping file (`-G/--genome2taxid`).
      There's no need to rebuild the index.
    - Check if the output file and the log file are the same.
    - Reduce the time of seed matching when using `-w`.
    - Change the default value of `--max-query-conc` from 12 to 8.
    - New flag `--gc-interval` (default 64, 0 to disable) for forcing garbage collection every N queries. This decreases memory usage a lot.
- `lexicmap utils subseq`:
    - **Accept the output file of `lexicmap search` as the input**.
      So one can extract matched sequences (including flanking regions) from the index, after alignment with `lexicmap search` with or without using the flag `-a/--all`.
    - Support extending aligned regions with `-U/--upstream` and/or `-D/--downstream`.

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
    - **Change default option values to provide higher sensitivity for short (<=500, especially <=250) queries,
      faster indexing speed, and faster seed-matching speed<s>, at the cost of a slightly larger index</s>**.
        - `-m/--masks`: 40,000 -> 20,000. 
           40k is unnecessary especially for small genomes, where seeds would be very crowded,
           with a large proportion of seed distances between 0 and 50 bp.
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
    - **Change the default value of `--partitions` from 1024 to 4096, which increases the seed-matching speed at the cost of 2 GiB of additional memory**.
      For existing lexicmap indexes, just run `lexicmap utils reindex-seeds --partitions 4096` to re-create seed indexes.
    - **Do not save low-complexity seeds**.
    - Fix high memory usage in writing seed data.
    - Change the default value of `-c/--chunks` from all available CPUs to the value of `-j/--threads`.
    - Change the default value of `--max-open-files` from 512 to 1024.
    - Add a new flag `--debug`.
- `lexicmap search`:
    - **Improving chaining, pseudoalignment, and alignment for highly repetitive sequences**.
    - **More accurate chaining score with better chaining of overlapped anchors, this produces more accurate results with `-n/--top-n-genomes`**: 
         - Merging two overlapped non-gapped anchors into a longer one.
         - For those with gaps, only the non-overlapping part of the second anchor is used to compute the weight.
         - Using the score of the best chain (rather than the sum) for sorting genomes when using `-n`.
    - Fix positions and alignment texts for queries with highly repetitive sequences in end regions. [#9](https://github.com/shenwei356/LexicMap/issues/9)
    - Skip low-complexity seeds.
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
    - Use a buffered reader for reading seed files.
    - Fix object recycling and reduce memory usage.
    - Fix alignment against genomes with many short contigs.
    - Fix early quit when meeting a sequence shorter than k.
    - Add a new option `-J/--max-query-conc` to limit the maximum number of concurrent queries,
      with a default value of 12 instead of the number of CPUs, which reduces the memory usage
      in batch searching.
    - Result format:
        - Cluster alignments of each target sequence.
        - Remove the column `seeds`.
        - Add columns `gaps`, `cigar`, `align`, which can be reformatted with `lexicmap utils 2blast`.
- `lexicmap utils kmers`:
    - Fix the progress bar.
    - Fix a bug where some masks do not have any k-mer.
    - Add a new column `prefix` to show the length of the common prefix between the seed and the probe.
    - Add a new column `reversed` to indicate if the k-mer is reversed for suffix matching.
- `lexicmap utils masks`:
    - Add support for outputting only a specific mask.
- `lexicmap utils seed-pos`:
    - New columns: `sseqid` and `pos_seq`.
    - More accurate seed distance.
    - Add histograms of numbers of seeds in sliding windows.
- `lexicmap utils subseq`:
    - Fix a bug when the given end position is larger than the sequence length.
    - Add the strand ("+" or "-") in the sequence header.

### v0.3.0 - 2024-05-14

- `lexicmap index`:
    - **Better seed coverage by filling sketching deserts**.
    - **Use longer (1000bp N's, previous: k-1) intervals between contigs**.
    - Fix a concurrency bug between genome data writing and k-mer-value data collecting.
    - Change the format of k-mer-value index file, and fix the computation of index partitions.
    - Optionally save seed positions which can be output by `lexicmap utils seed-pos`.
- `lexicmap search`:
    - **Improved seed-chaining algorithm**.
    - **Better support for long queries**.
    - **Add a new flag `-w/--load-whole-seeds` for loading the whole seed data into memory for faster search**.
    - **Parallelize alignment in each query**, so it's faster for a single query.
    - **Optional output of matched query and subject sequences**.
    - 2-5X searching speed with a faster masking method.
    - Change output format.
    - Add output of query start and end positions.
    - Fix a target sequence extraction bug.
    - Keep indexes of genome data in memory.
- `lexicmap utils kmers`:
    - Fix a little bug, wrong number of k-mers for the second k-mer in each k-mer pair.
- New commands:
    - `lexicmap utils gen-masks` for generating masks from the top N largest genomes.
    - `lexicmap utils seed-pos` for extracting seed positions via reference names.
    - `lexicmap utils reindex-seeds` for recreating indexes of k-mer-value (seeds) data.
    - `lexicmap utils genomes` for listing genome IDs in the index.

### v0.2.0 - 2024-02-02

- Software architecture and index formats are redesigned to reduce searching memory usage.
- Indexing: genomes are processed in batches to reduce RAM usage, then indexes of all batches are merged.
- Searching: seed matching is performed on disk, yet it's ultra-fast.

### v0.1.0 - 2024-01-15

- The first release.
- Seed indexing and querying are performed in RAM.
- GTDB r214 with 10k masks: index size 75GB, RAM: 130GB.
