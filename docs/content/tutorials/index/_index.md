---
title: Building an index
weight: 0
---

## Input

LexicMap is designed for small genomes like Archaea, Bacteria, Viruses and plasmids.

**Sequences of each reference genome should be saved in a separate FASTA/Q file, with the identifier in the file name**.

- **File type**: FASTA/Q files, in plain text or gzip-compressed format.
- **File name**: "Genome ID" + "File extention". E.g., `GCF_000006945.2.fa.gz`.
    - **Genome ID**: they should be distinctive for accurate result interpretation, which will be shown in the search result.
    - File extention: a regular expression set by the flag `-N/--ref-name-regexp` is used to extract genome IDs from the file name.
      The default value supports common sequence file extentions, including `.fa`, `.fasta`, `.fna`, `.fa.gz`, `.fasta.gz`, and `.fna.gz`.
- **Sequences**:
    - **Only DNA or RNA sequences are supported**.
    - **Sequence IDs** should be distinctive for accurate result interpretation, which will be shown in the search result.
    - One or more sequences are allowed.
        - Unwanted sequences can be filtered out by regular expressions from the flag `-B/--seq-name-filter`.
    - **Genome size limit**. Some none-isolate assemblies might have extremely large genomes, e.g., [GCA_000765055.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000765055.1/) has >150 Mb.
     The flag `-g/--max-genome` (default 15 Mb) is used to skip these input files, and the file list would be written to a file
     via the flag `-G/--big-genomes`..
- At most 17,179,869,184 genomes are supported.

Input files can be given via one of the below ways:

- **Positional arguments**. For a few input files.
- A **file list** via the flag `-X/--infile-list`  with one file per line.
  It can be STDIN (`-`), e.g., you can filter a file list and pass it to `lexicmap`.
- A **directory** containing input files via the flag `-I/--in-dir`.
    - Multiple-level directories are supported.
    - Directory and file symlinks are followed.

## How

1. Generating [LexicHash masks](https://doi.org/10.1093/bioinformatics/btad652).

        |ATTATAACGCCACGGGGAGCCGCGGGGTTTC One k-bp mask
        |--------========_______________
        |
        |-------- Prefixes for covering all possible P-bp DNA.
        |         The length is the largest number for 4^P <= #masks
        |
        |--------======== Extend prefixes, chosen from the most frequent extended prefixes
        |                 of which the prefix-derived k-mers do not overlap masked k-mers.
        |
        |                _______________ Randomly generated bases


    1. Reading all input genomes and getting the genome size information.
    2. Recording k-mers on the positive and negative strands of the top *N* biggest genomes.
    3. Generating all permutations of *p*-bp prefixes that can cover all possible k-mers, *p* is the biggest value for 4<sup>*p*</sup> <= *M* (desired number of masks), e.g., *p*=7 for 40,000 masks.
    4. For each prefix.
        1. Counting the k-mer frequencies of extend prefixes (`-P/--prefix-ext`, default 8) for each top *N* genome,
           and randomly choose one from the most frequent one.
        2. Randomly generating left *k*-*p*-*P* bases.
        3. Capturing the most similar k-mer in each genome, and record k-mer locations in a interval tree.
        4. Go to 1.4.1, and only counting k-mers which do not overlap existing regions.
    5. Generate left *M* - 4<sup>*p*</sup> masks similary, with a few differences.
        - Only generating masks with non-low-complexity prefixes.
        - If there's no available k-mers for a extend prefix, just randomly generate one.

2. For each genome batch (`-b/--batch-size`, default 10000, max 131072).
    1. For each genome file in a genome batch.
        1. Optionally discarding sequences via regular expression (`-B/--seq-name-filter`).
        2. Skipping genomes bigger than the value of `-g/--max-genome`.
        3. Concatenating all sequences, with intervals of 1000-bp N's.
        4. Capureing the most similar k-mer for each mask and saving k-mer and its location(s) and strand information.
        5. Saving the concatenated genome sequence (bit-packed, 2 bits for one base) and genome information (genome ID, size, and lengths of all sequences) into the genome data file, and creating an index file for the genome data file for fast random subsequence extraction.
    2. Compressing k-mer and the corresponding data (k-mer-data, or seeds data) into chunks of files, and creating an index file for each k-mer-data file for fast seeding.
    3. Writing summary information into `info.toml` file.
3. Merging indexes of multiple batches.



## Parameters

{{< hint type=note >}}
**Genome size**\
LexicMap is only suitable for small genomes like Archaea, Bacteria, Viruses and plasmids.
{{< /hint >}}

## Steps

## Output
