---
title: Building an index
weight: 0
---

## Input

LexicMap is designed for small genomes like Archaea, Bacteria, Viruses and plasmids.

**Sequences of each reference genome should be saved in a separate FASTA/Q file, with the identifier in the file name**.

- **File type**: FASTA format, supported compression formats: gzip (`.gz`), xz (`.xz`), zstd (`.zst`), and bzip2 (`.bz2`).
- **File name**: "Genome ID" + "File extention". E.g., `GCF_000006945.2.fa.gz`.
    - **Genome ID**: they should be distinctive for result interpretation , which will be shown in the search result.
    - File extention: a regular expression given by the flag `-N/--ref-name-regexp` is used to extract genome IDs.
      The default value supports common file extentions, including `.fa`, `.fasta`, `.fna`, `.fa.gz`, `.fasta.gz`, and `.fna.gz`.
- **Sequences**:
    - **Only DNA or RNA sequences are supported**.
    - **Sequence IDs** should be distinctive for result interpretation, which will be shown in the search result.
    - One or more sequences are allowed.
        - Unwanted sequences can be filtered out by regular expressions from the flag `-B/--seq-name-filter`.
    - **Genome size limit**. Some none-isolate assemblies might have extremely large genomes, e.g., [GCA_000765055.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000765055.1/) has >150 Mb.
     The flag `-g/--max-genome` (default 15 Mb) is used to skip these input files, and the file list would be written to a file
     via the flag `-G/--big-genomes`..

Input files can be given via one of the ways below:

- **Positional arguments**. For a few input files.
- A **file list** via the flag `-X/--infile-list`  with one file per line.
  It can be STDIN (`-`), e.g., you can filter a file list and pass it to `lexicmap`.
- A **directory** containing input files via the flag `-I/--in-dir`.
    - Multiple-level directories are supported.
    - Directory and file symlinks are followed.

## How


## Parameters

{{< hint type=note >}}
**Genome size**\
LexicMap is only suitable for small genomes like Archaea, Bacteria, Viruses and plasmids.
{{< /hint >}}


## Output
