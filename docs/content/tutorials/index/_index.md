---
title: Step 1. Building a database
weight: 0
---

{{< hint type=note >}}
Terminology differences:

- On this page and in the LexicMap command line options, the term **"mask"** is used, following the terminology in the LexicHash paper.
- In the LexicMap manuscript, however, we use **"probe"** as it is easier to understand.
  These masks consist of thousands of k-mers and capture k-mers from sequences through prefix matching, functioning similarly to DNA probes in molecular biology.
{{< /hint >}}

## Table of contents

{{< toc format=html >}}

## TL;DR

1. Prepare input files:
    - **Sequences of each reference genome should be saved in a separate FASTA file, with identifiers (no tab symbols) in the file names**.
      E.g., GCF_000006945.2.fna.gz
        - A regular expression is also available to extract the reference ID from the file name.
          E.g., `--ref-name-regexp '^(\w{3}_\d{9}\.\d+)'` extracts `GCF_000006945.2` from GenBank assembly file `GCF_000006945.2_ASM694v2_genomic.fna.gz`.
        - Even if you forgot to use `-N/--ref-name-regexp`,
          [lexicmap utils edit-genome-ids](https://bioinf.shenwei.me/LexicMap/usage/utils/edit-genome-ids/) (available since v0.8.0) can fix this without rebuilding the index.
    - If you save *a few* **small** (viral) **complete** genomes (one sequence per genome) in each file, it's feasible as sequence IDs in search results can help distinguish target genomes.
2. Run:
    - From a directory with multiple genome files:

          lexicmap index -I genomes/ -O db.lmi

    - From a file list with one file per line:

          lexicmap index --skip-file-check -X files.txt -O db.lmi


## Input

{{< hint type=note >}}
**Genome size**\
LexicMap is mainly suitable for small genomes like Archaea, Bacteria, Viruses, fungi, and plasmids.

Maximum genome size: 268 Mb (268,435,456).
More precisely:

    $total_bases + ($num_contigs - 1) * 1000 <= 268,435,456

as we concatenate contigs with 1000-bp intervals of N’s to reduce the number of sequences to index.

{{< /hint >}}


<font color="red">**Sequences of each reference genome should be saved in a separate FASTA file, with identifiers in the file names**.</font>
If you save *a few* **small** (viral) **complete** genomes (one sequence per genome) in each file, it's feasible as sequence IDs in search results can help distinguish target genomes.

- **File type**: FASTA/Q files, in **plain text or gzip/xz/zstd/bzip2/lz4 compressed** formats.
- **File name**: "Genome ID" + "File extension". E.g., `GCF_000006945.2.fna.gz`.
    - **Genome ID**: **IDs must not contain tab ("\t") symbols and should be distinct for accurate interpretation of results**. They will be shown in the search results.
        - A regular expression is also available to extract the reference ID from the file name.
          E.g., `--ref-name-regexp '^(\w{3}_\d{9}\.\d+)'` extracts `GCF_000006945.2` from GenBank assembly file `GCF_000006945.2_ASM694v2_genomic.fna.gz`.
        - **If you forgot to use `-N/--ref-name-regexp`,
        [lexicmap utils edit-genome-ids](https://bioinf.shenwei.me/LexicMap/usage/utils/edit-genome-ids/) (available since v0.8.0) can fix this without rebuilding the index**.
- **File extension**: a regular expression set by the flag `-r/--file-regexp` is used to match input files.
      The default value supports common sequence file extensions, e.g., `.fa`, `.fasta`, `.fna`, `.fa.gz`, `.fasta.gz`, `.fna.gz`, `fasta.xz`, `fasta.zst`, and `fasta.bz2`.
- **Sequences**:
    - **Only DNA or RNA sequences are supported**.
       - **Soft-masked sequences** have been supported since v0.9.0 with the flag `--soft-masking`.
         Lowercase bases in soft-masked low-complexity regions will be treated as A's and won't be seeded,
         while they will be saved for base-level alignment.
    - **Sequence IDs** should be distinct for accurate interpretation of results. They will be shown in the search results.
    - Sequence descriptions (text after the sequence ID) are not saved. If you do need them, you can create a mapping file
      (`cat files.txt | seqkit seq -n -X - | sed -E 's/\s+/\t/' > id2desc.tsv`) and use it to [add descriptions to search results](https://bioinf.shenwei.me/LexicMap/tutorials/search/#summarizing-results).
    - **One or more sequences (contigs) in each file are allowed**.
        - **Unwanted sequences** (such as plasmids) can be filtered out by regular expressions from the flag `-B/--seq-name-filter`.
    - **Genome size limit**. Some non-isolate assemblies might have extremely large genomes, e.g., [GCA_000765055.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000765055.1/) has >150 Mb.
     The flag `-g/--max-genome` (default 15 Mb) is used to skip these input files, and the file list would be written to a file
     via the flag `-G/--big-genomes`.
        - **Changes since v0.5.0**:
            - Genomes with any single contig larger than the threshold will be skipped as before.
            - However, **fragmented (with many contigs) genomes with a total number of bases larger than the threshold will
              be split into chunks** and alignments from these chunks will be merged in `lexicmap search`.
        - **For fungal genomes, please increase the value of `-g/--max-genome`**.
    - **Minimum sequence length**. A flag `-l/--min-seq-len` can filter out sequences shorter than the threshold (default is the `k` value).
- **At most 17,179,869,184 (2<sup>34</sup>) genomes are supported**. For more genomes, please create a file list and split it into multiple parts, and build an index for each part.

**Input files can be given in one of the following ways:**

- **Positional arguments**. For a few input files.
- A **file list** via the flag `-X/--infile-list`  with one file per line.
  **It can be STDIN (`-`)**, e.g., you can filter a file list and pass it to `lexicmap index`.
    - **The flag `-S/--skip-file-check` is optional and skips input file checks if you are sure these files exist**.
    By default, LexicMap checks the existence of all input files, which would take tens of minutes for >1M files.
- A **directory** containing input files via the flag `-I/--in-dir`.
    - **Multiple-level directories are supported**. So you don't need to save hundreds of thousands of files in one directory.
    - **Directory and file symlinks are followed**.

## Hardware requirements

See [benchmark of index building](https://bioinf.shenwei.me/LexicMap/introduction/#indexing).

LexicMap is designed to provide fast and low-memory sequence alignment against millions of prokaryotic genomes.

- **CPU:**
    - No specific requirements on CPU type and instruction sets. Both x86 and ARM chips are supported.
    - More is better as LexicMap is CPU-intensive software. **It uses all CPUs by default (`-j/--threads`)**.
- **RAM**
    - More RAM (> 200 GB for >2 million genomes) is preferred. The memory usage in index building is mainly related to:
        - **The number of masks** (`-m/--masks`, default 20,000). Bigger values improve the search sensitivity slightly, increase the index size, and slow down the search speed. For smaller genomes like phages/viruses, m=5,000 is high enough.
        - **The number of genomes**. Generally, more genomes consume more memory, but we can reduce it by using smaller genome batch sizes (see below).
        - **The genome batch size**  (`-b/--batch-size`, default 5,000). <font color="red">This is the main parameter to adjust **memory usage**</font>. Bigger values increase indexing memory usage.
        - **The divergence between genome sequences in each batch**. Diverse genomes consume more memory.
        - **The maximum seed distance** or **the maximum sketching desert size** (`-D/--seed-max-desert`, default 100),
          and the distance of k-mers to fill deserts (`-d/--seed-in-desert-dist`, default 50).
          <font color="#ff5733">These are the main parameters to adjust **search sensitivity**.</font>
          Bigger `-D/--seed-max-desert` values decrease the search sensitivity, speed up indexing,
          decrease indexing memory usage, and decrease index size. Alignment speed is largely unaffected.
    - **If the RAM is not sufficient**. Please:
        - **Use a smaller genome batch size**. It decreases indexing memory usage and has little effect on search performance.
- **Disk**
    - More is better. LexicMap index size is related to the number of input genomes, the divergence between genome sequences, the number of masks, and the maximum seed distance. See [some examples](#index-size).
        - **Note that index size grows sublinearly with the number of genomes**. Because seed data are compressed with the VARINT-GB algorithm, more genomes lead to smaller compression ratios.
    - SSD disks are preferred, while HDD disks are also fast enough.

## Algorithm

<img src="/LexicMap/indexing.svg" alt="" width="900"/>

See the [paper](https://bioinf.shenwei.me/LexicMap/introduction/#citation) for details.

{{< expand "Click to show details." "..." >}}

1. **Generating *m* [LexicHash masks](https://doi.org/10.1093/bioinformatics/btad652)**.

    1. Generate *m* prefixes.
        1. Generating all permutations of *p*-bp prefixes that can cover all possible k-mers; *p* is the largest value for 4<sup>*p*</sup> <= *m* (desired number of masks), e.g., *p*=7 for 20,000 masks. (4<sup>*7*</sup> = 16384)
        3. Duplicating these prefixes to *m* prefixes.
    2. For each prefix,
        1. Randomly generating the remaining *k*-*p* bases.
        3. If the mask is duplicated, re-generating.

2. **Building an index for each genome batch** (`-b/--batch-size`, default 5,000, max 131,072).

    1. For each genome file in a genome batch.
        1. Optionally discarding sequences via regular expression (`-B/--seq-name-filter`).
        2. Skipping genomes bigger than the value of `-g/--max-genome`.
        3. Concatenating all sequences, with intervals of 1000-bp N's.
        4. Capturing the most similar k-mer (in non-gap and non-interval regions) for each mask and recording the k-mer and its location(s) and strand information. Base N is treated as A.
        5. Filling sketching deserts (genome regions longer than `--seed-max-desert` [default 100] without any captured k-mers/seeds).
           In a sketching desert, not a single k-mer is captured because there's another k-mer in another place which shares a longer prefix with the mask.
           As a result, for a query similar to sequences in this region, none of the captured k-mers can match the correct seeds.
            1. For a desert region (`start`, `end`), masking the extended region (`start-1000`, `end+1000`) with the masks.
            2. Starting from `start`, approximately every `--seed-in-desert-dist` (default 50) bp, finding a k-mer which is captured by some mask, and adding the k-mer and its position information into the index of that mask.
        6. Saving the concatenated genome sequence (bit-packed, 2 bits for one base, N is treated as A) and genome information (genome ID, size, and lengths of all sequences) into the genome data file, and creating an index file for the genome data file for fast random subsequence extraction.
    2. Duplicate and reverse all k-mers, and save each reversed k-mer along with the duplicated position information in the seed data of the closest (sharing the longest prefix) mask. This is for suffix matching of seeds.
    2. Compressing k-mers and the corresponding data (k-mer-data, or seed data, including genome batch, genome number, location, and strand) into chunks of files, and creating an index file for each k-mer-data file for fast seeding.
    3. Writing summary information into the `info.toml` file.

3. **Merging indexes of multiple batches**.
    1. For each k-mer-data chunk file (belonging to a list of masks), serially reading data of each mask from all batches,
      merging them and writing them to a new file.
    2. For genome data files, just move them.
    3. Concatenating `genomes.map.bin`, which maps each genome ID to its batch ID and index in the batch.
    4. Update the index summary file.

{{</ expand >}}

## Parameters

{{< hint type=note >}}
**Query length**

**LexicMap is mainly designed for sequence alignment with a small number of queries (gene/plasmid/virus/phage sequences) longer than 150 bp by default**.

**If you want to search some short reads, you need to build the index with small values of `-D/--seed-max-desert` (default 100) and `-d/--seed-in-desert-dist` (default 50), e.g., `-D 60 -d 30` for 125bp reads, or `-D 50 -D 25` for 100bp reads**. This will increase indexing time and index size. This is less of a concern with a small number of genomes, such as < 10,000.

If you just want to search long (>1kb) queries for highly similar (>95%) targets, you can build an index with a bigger `-D/--seed-max-desert` (default 100) and `-d/--seed-in-desert-dist` (default 50), e.g., `-D 300 -d 150`. Bigger values decrease search sensitivity for distant targets, speed up
indexing, decrease indexing memory usage, and decrease index size. Alignment speed is largely unaffected.

Note that **LexicMap is slow for ultra-long (>1Mb) queries, and the alignment might be fragmented**.
{{< /hint >}}



**Flags in bold text** are important and frequently used.

{{< tabs "t1" >}}

{{< tab "General" >}}

|Flag              |Value                      |Function                     |Comment                                                                                                                                                 |
|:-----------------|:--------------------------|:----------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------|
|**`-j/--threads`**|Default: all available CPUs|Number of CPU cores to use.  |► If the value is smaller than the number of available CPUs, make sure to use the same value for `-c/--chunks`.                                            |
|`--soft-masking`  |Default: false             |Support soft-masked sequences|► Lowercase bases in soft-masked low-complexity regions will be treated as A's and won't be seeded, while they will be saved for base-level alignment.  |

{{< /tab>}}

{{< tab "Genome batches" >}}

| Flag                  | Value                      | Function                                | Comment                                                                                                                                                                                                                                                                                                                                                                       |
| :-------------------- | :------------------------- | :-------------------------------------- | :---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **`-b/--batch-size`** | Max: 131072, default: 5000 | Maximum number of genomes in each batch | If the number of input files exceeds this number, input files are split into multiple batches and indexes are built for all batches. In the end, seed files are merged, while genome data files are kept unchanged and collected. ■ Bigger values increase indexing memory usage and increase batch searching speed, while single-query search speed is not affected.    |

{{< /tab>}}

{{< tab "LexicHash mask generation" >}}

| Flag             | Value                | Function               | Comment                                                                                                                                                                                    |
| :--------------- | :------------------- | :--------------------- | :----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-M/--mask-file` | A file               | File with custom masks | File with custom masks, which can be exported from an existing index or newly generated by "lexicmap utils masks". This flag overrides `-k/--kmer`, `-m/--masks`, `-s/--rand-seed`, etc.   |
| **`-k/--kmer`**  | Max: 32, default: 31 | K-mer size             | ■ Bigger values improve the search specificity and do not increase the index size.                                                                                                         |
| **`-m/--masks`** | Default: 20,000      | Number of masks        | ■ Bigger values improve the search sensitivity slightly, increase the index size, and slow down the search speed. For smaller genomes like phages/viruses, m=5,000 is high enough.                 |


{{< /tab>}}


{{< tab "Seeds (k-mer-value) data" >}}

|Flag                        |Value                                       |Function                                                  |Comment                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
|:---------------------------|:-------------------------------------------|:---------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|**`--seed-max-desert`**     |Default: 100                                |Maximum length of distances between seeds                 |The default value of 100 guarantees that queries >=200 bp will match at least two seeds. ► Large regions with no seeds are called sketching deserts. Deserts with seed distance larger than this value will be filled by choosing k-mers roughly every --seed-in-desert-dist (50 by default) bases. ■ Bigger values decrease the search sensitivity for distant targets, speed up indexing, decrease indexing memory usage, and decrease index size. Alignment speed is largely unaffected.    |
|**`-c/--chunks`**           |Maximum: 128, default: value of -j/--threads|Number of seed file chunks                                |Bigger values accelerate the search speed at the cost of a higher disk read load. ► The value should not exceed the maximum number of open files set by the operating system. ► Make sure the value of `-j/--threads` in `lexicmap search` is >= this value.                                                                                                                                                                                                                                   |
|**`-J/--seed-data-threads`**|Maximum: -c/--chunks, default: 8            |Number of threads for merging seed chunks from all batches|Merging reserves 8 open files. The actual number of merging threads is min(--seed-data-threads, (--max-open-files - 8)/($inputs + 2)), where $inputs is the number of batches in the current group (at most --max-open-files - 10). ■ Bigger values increase I/O load in HDDs for large batch counts.                                                                                                                                                                                          |
|`-p/--partitions`           |Default: 4096                               |Number of partitions for indexing each seed file          |Bigger values slightly increase memory usage. ► After indexing, `lexicmap utils reindex-seeds` can be used to reindex the seed data with another value of this flag.                                                                                                                                                                                                                                                                                                                           |
|**`--max-open-files`**      |Default: 1024                               |Maximum number of open files                              |Used in merging indexes of multiple genome batches (minimum: 12), including 8 reserved files. Large batch counts are merged in multiple rounds. Increasing the budget can allow more merging threads; set `ulimit -n` at least this high.                                                                                                                                                                                                                                                      |

{{< /tab>}}

{{< /tabs >}}

Also see the [usage](https://bioinf.shenwei.me/LexicMap/usage/#index) of `lexicmap index`.

### Notes for indexing with large datasets

<font color="#ff5733">If you have hundreds of thousands of input genomes or more, it's **better to control the number of genome batches**</font>, which can be calculated via

    $num_input_files / --batch-size

E.g., for GenBank prokaryotic genomes: 2,340,672 / 5000 (default)  = 468.
The number is too big, and **it would slow down the seed-data merging step in `lexicmap index`** and **candidate sequence extraction in `lexicmap search`**.

Therefore, if you have enough memory, you can set a bigger `--batch-size` (e.g., 2,340,672 / 25000 = 93.6).

If the batch number is still big (e.g. 300), you can set bigger `--max-open-files` (e.g., `4096`) and `-J/--seed-data-threads` (e.g., `12`. 12 <= 4096/300 = 13.6)
to accelerate the merging step. Meanwhile, don't forget to increase the maximum open files per process via `ulimit -n 4096`.

If you forgot these settings, you can rerun the merging step for an unfinished index via [lexicmap utils remerge](https://bioinf.shenwei.me/LexicMap/usage/utils/remerge/)
(available since v0.5.0, also see [FAQ: how to resume the indexing](https://bioinf.shenwei.me/LexicMap/faqs/#how-to-resume-the-indexing-as-slurm-job-limit-is-almost-reached-while-lexicmap-index-is-still-in-the-merging-step)). Other cases to use this command:
- Only one thread is used for merging indexes, which happens when there are
many batches (>200) (`$inpu_files / --batch-size`) and the value
of `--max-open-files` is not large enough.
- The Slurm/PBS job time limit is almost reached and the merging step won't finish in time.
- The disk quota is reached during the merging step.

## Steps

We use a small dataset for demonstration.

1. Preparing the test genomes (15 bacterial genomes) in the `refs` directory.

   Note that the genome files contain the assembly accessions (ID) in the file names.

        git clone https://github.com/shenwei356/LexicMap
        cd LexicMap/demo/

        ls refs/
        GCF_000006945.2.fa.gz  GCF_000392875.1.fa.gz  GCF_001096185.1.fa.gz  GCF_002949675.1.fa.gz  GCF_006742205.1.fa.gz
        GCF_000017205.1.fa.gz  GCF_000742135.1.fa.gz  GCF_001457655.1.fa.gz  GCF_002950215.1.fa.gz  GCF_009759685.1.fa.gz
        GCF_000148585.2.fa.gz  GCF_001027105.1.fa.gz  GCF_001544255.1.fa.gz  GCF_003697165.2.fa.gz  GCF_900638025.1.fa.gz

1. Building an index with genomes from **a directory**.

        lexicmap index -I refs/ -O demo.lmi

    It would take about 2 seconds and 1.5 GB RAM on a 16-CPU PC.

    Optionally, we can also use **a file list** as the input.

        $ head -n 3 files.txt
        refs/GCF_000006945.2.fa.gz
        refs/GCF_000017205.1.fa.gz
        refs/GCF_000148585.2.fa.gz

        lexicmap index -S -X files.txt -O demo.lmi

    {{< expand "Click to show the log of a demo run." "..." >}}

    # here we set a small --batch-size 5
    # memusg is used to measure time and memory usage.
    #    https://github.com/shenwei356/memusg/

    $ memusg -t -s "lexicmap index -I refs/ -O demo.lmi --batch-size 5"
    16:22:52.012 [INFO] LexicMap v0.10.0
    16:22:52.012 [INFO]   https://github.com/shenwei356/LexicMap
    16:22:52.012 [INFO] 
    16:22:52.012 [INFO]  CWD: /home/shenwei/go/src/github.com/shenwei356/LexicMap/demo
    16:22:52.012 [INFO]  CMD: lexicmap index -I refs/ -O demo.lmi --batch-size 5
    16:22:52.012 [INFO] DATE: 2026-10-09
    16:22:52.012 [INFO] 
    16:22:52.013 [INFO] checking input files ...
    16:22:52.013 [INFO]   scanning files from directory: refs/
    16:22:52.013 [INFO]   15 input file(s) given
    16:22:52.013 [INFO] 
    16:22:52.013 [INFO] --------------------- [ main parameters ] ---------------------
    16:22:52.013 [INFO] 
    16:22:52.013 [INFO] input and output:
    16:22:52.013 [INFO]   input directory: refs/
    16:22:52.013 [INFO]     regular expression of input files: (?i)\.(f[aq](st[aq])?|fna)(\.gz|\.xz|\.zst|\.bz2)?$
    16:22:52.013 [INFO]     *regular expression for extracting reference name from file name: (?i)(.+)\.(f[aq](st[aq])?|fna)(\.gz|\.xz|\.zst|\.bz2)?$
    16:22:52.013 [INFO]     *regular expressions for filtering out sequences: []
    16:22:52.013 [INFO]   min sequence length: 31
    16:22:52.013 [INFO]   max genome size: 20000000
    16:22:52.013 [INFO]   output directory: demo.lmi
    16:22:52.013 [INFO] 
    16:22:52.013 [INFO] mask generation:
    16:22:52.013 [INFO]   k-mer size: 31
    16:22:52.013 [INFO]   number of masks: 20000
    16:22:52.013 [INFO]   rand seed: 1
    16:22:52.013 [INFO] 
    16:22:52.013 [INFO] seed data:
    16:22:52.013 [INFO]   maximum sketching desert length: 100
    16:22:52.013 [INFO]   distance of k-mers to fill deserts: 50
    16:22:52.013 [INFO]   seeds data chunks: 16
    16:22:52.013 [INFO]   seeds data indexing partitions: 4096
    16:22:52.013 [INFO] 
    16:22:52.013 [INFO] general:
    16:22:52.013 [INFO]   genome batch size: 5
    16:22:52.013 [INFO]   threads: 16
    16:22:52.013 [INFO]   batch merge threads: 8
    16:22:52.013 [INFO] 
    16:22:52.013 [INFO] 
    16:22:52.013 [INFO] --------------------- [ generating masks ] ---------------------
    16:22:52.017 [INFO] 
    16:22:52.017 [INFO] --------------------- [ building index ] ---------------------
    16:22:52.313 [INFO] 
    16:22:52.313 [INFO]   ------------------------[ batch 1/3 ]------------------------
    16:22:52.313 [INFO]   building index for batch 1 with 5 files...
    processed files:  5 / 5 [======================================] ETA: 0s. done
    16:22:53.162 [INFO]   writing seeds...
    16:22:53.235 [INFO]   finished writing seeds in 72.926598ms
    16:22:53.235 [INFO]   finished building index for batch 1 in: 921.70147ms
    16:22:53.261 [INFO] 
    16:22:53.262 [INFO]   ------------------------[ batch 2/3 ]------------------------
    16:22:53.262 [INFO]   building index for batch 2 with 5 files...
    processed files:  5 / 5 [======================================] ETA: 0s. done
    16:22:54.174 [INFO]   writing seeds...
    16:22:54.263 [INFO]   finished writing seeds in 88.10286ms
    16:22:54.263 [INFO]   finished building index for batch 2 in: 1.001148132s
    16:22:54.290 [INFO] 
    16:22:54.290 [INFO]   ------------------------[ batch 3/3 ]------------------------
    16:22:54.290 [INFO]   building index for batch 3 with 5 files...
    processed files:  5 / 5 [======================================] ETA: 0s. done
    16:22:55.158 [INFO]   writing seeds...
    16:22:55.231 [INFO]   finished writing seeds in 72.762597ms
    16:22:55.231 [INFO]   finished building index for batch 3 in: 940.67608ms
    16:22:55.283 [INFO] 
    16:22:55.283 [INFO] merging 3 indexes...
    16:22:55.283 [INFO]   [round 1]
    16:22:55.283 [INFO]     batch 1/1, merging 3 indexes to demo.lmi.tmp/r1_b1 with 8 threads...
    16:22:55.367 [INFO]   [round 1] finished in 84.472031ms
    16:22:55.367 [INFO] rename demo.lmi.tmp/r1_b1 to demo.lmi
    16:22:55.381 [INFO] 
    16:22:55.381 [INFO] finished building LexicMap index from 15 files with 20000 masks in 3.368742159s
    16:22:55.381 [INFO] LexicMap index saved: demo.lmi
    16:22:55.381 [INFO] 
    16:22:55.381 [INFO] elapsed time: 3.368781691s
    16:22:55.381 [INFO] 

    elapsed time: 3.541s
    peak rss: 1.25 GB
    {{< /expand >}}

From v0.10.0, it is optional to create adaptive two-level indexes of seed data.
It can improve seed-matching performance for batch queries with `lexicmap search`.
See the [usage](https://bioinf.shenwei.me/LexicMap/usage/utils/reindex-seeds2) for more details.

## Output

The LexicMap index is a directory with multiple files.

### File structure

    $ tree demo.lmi/
    demo.lmi/                    # the index directory
    ├── genomes                  # directory of genome data
    │   ├── batch_0000           # genome data of one batch
    │   │   ├── genomes.bin      # genome data file, containing genome ID, size, sequence lengths, bit-packed sequences
    │   │   └── genomes.bin.idx  # index of genome data file, for fast subsequence extraction
    │   └── batch_0001
    │   │   ├── genomes.bin
    │   │   └── genomes.bin.idx
    │   ... ...
    ├── seeds                    # seed data: pairs of k-mer and its location information (genome batch, genome number, location, strand)
    │   ├── chunk_000.bin        # seed data file
    │   ├── chunk_000.bin.idx    # index of seed data file, for fast seed searching and data extraction
    │   ... ...
    │   ├── chunk_015.bin        # the number of chunks is set by flag `-c/--chunks`, default: #cpus
    │   └── chunk_015.bin.idx
    ├── genomes.chunks.bin       # lists of genome chunks which belong to the same genome
    ├── genomes.map.bin          # mapping genome ID to batch number and genome number in the batch
    ├── info.toml                # summary of the index
    └── masks.bin                # mask data

### Index size

LexicMap index size is related to the number of input genomes, the divergence between genome sequences, the number of masks, and the maximum seed distance.

**Note that index size grows sublinearly with the number of genomes**. Because seed data are compressed with the VARINT-GB algorithm, more genomes lead to smaller compression ratios (smaller is better).

{{< tabs "t2" >}}

{{< tab "Demo data" >}}

    # 15 genomes, 54 Mb (54,142,446 bp)
    
    demo.lmi/: 78.36 MiB (82,165,269)
     65.26 MiB      seeds
     12.94 MiB      genomes
    156.28 KiB      masks.bin
         600 B      info.toml
         375 B      genomes.map.bin
           0 B      genomes.chunks.bin

{{< /tab>}}

{{< tab "GTDB repr" >}}

    # 85,205 genomes, 274 Gbp (273,848,490,566 bp)
    
    gtdb_repr.lmi: 213.27 GiB (228,999,914,466)
    146.49 GiB      seeds
     66.78 GiB      genomes
      2.03 MiB      genomes.map.bin
    156.28 KiB      masks.bin
         613 B      info.toml
          48 B      genomes.chunks.bin

{{< /tab>}}

{{< tab "GTDB complete" >}}

    # 402,538 genomes, 1.5 Tbp (1,501,864,915,560 bp)
    
    gtdb_complete.lmi: 905.34 GiB (972,098,200,328)
    542.34 GiB      seeds
    362.99 GiB      genomes
      9.60 MiB      genomes.map.bin
    156.28 KiB      masks.bin
         616 B      info.toml
         168 B      genomes.chunks.bin

{{< /tab>}}


{{< tab "GenBank+RefSeq" >}}

    # 2,340,672 genomes, 9.2 Tbp (9,192,651,650,196 bp)
         
    genbank_refseq.lmi: 4.96 TiB (5,454,659,703,138)
      2.79 TiB      seeds
      2.17 TiB      genomes
     55.81 MiB      genomes.map.bin
    156.28 KiB      masks.bin
      3.59 KiB      genomes.chunks.bin
         619 B      info.toml

{{< /tab>}}


{{< tab "AllTheBacteria HQ" >}}

    # 1,858,610 genomes, 7.5 Tbp (7,493,622,021,123 bp)
    
    atb_hq.lmi: 3.91 TiB (4,304,515,140,156)
      2.15 TiB      seeds
      1.77 TiB      genomes
     39.22 MiB      genomes.map.bin
    156.28 KiB      masks.bin
         619 B      info.toml
          24 B      genomes.chunks.bin


{{< /tab>}}

{{< /tabs >}}

- Directory/file sizes are counted with https://github.com/shenwei356/dirsize v1.2.1 (`dirsize $file`, **base: 1024**).
- Index building parameters: `-k 31 -m 20000 -D 100 -d 50`. Genome batch size: `-b 25000` for GenBank+RefSeq and AllTheBacteria datasets, `-b 5000` (default) for others.


## Explore the index

We provide several commands to explore the index data and extract indexed subsequences:

- Extracting sequences
    1. `lexicmap utils genome-seqs` can extract all sequences of a given genome,
        see the [usage and example](https://bioinf.shenwei.me/LexicMap/usage/utils/genome-seqs/).
    1. `lexicmap utils subseq` can extract subsequences via 1) reference name, sequence ID, position and strand, or 2) search result,
        see the [usage and example](https://bioinf.shenwei.me/LexicMap/usage/utils/subseq/).
- Genome data
    1. `lexicmap utils genomes` can list genome IDs of indexed genomes,
        see the [usage and example](https://bioinf.shenwei.me/LexicMap/usage/utils/genomes/).
    1. `lexicmap utils genome-details` can extract or view genome details in the index,
        see the [usage and example](https://bioinf.shenwei.me/LexicMap/usage/utils/genome-details/).
- Seed data
    1. `lexicmap utils kmers` can list details of all seeds (k-mers), including reference, location(s), the strand, and the k-mer direction.
        see the [usage and example](https://bioinf.shenwei.me/LexicMap/usage/utils/kmers/).
    1. `lexicmap utils seed-pos` can help to explore the seed positions,
        see the [usage and example](https://bioinf.shenwei.me/LexicMap/usage/utils/seed-pos/).
        Before that, the flag `--save-seed-pos` needs to be added to `lexicmap index`.
    1. `lexicmap utils masks` can list masks of the index,
        see the [usage and example](https://bioinf.shenwei.me/LexicMap/usage/utils/masks/).

## Index format changelog

Index version information is available in the `info.toml` file of each LexicMap index.
LexicMap search and other utility commands check compatibility via the main version.

|Index version|LexicMap version|Supported LexicMap versions|Date      |Changes                                                                                                                |
|:------------|:---------------|:--------------------------|:---------|:----------------------------------------------------------------------------------------------------------------------|
|3.5          |0.10.0          |0.6.0 +                    |2026-xx-xx|A small fraction of seeds change after fixing the lexichash computation.                                              |
|3.4          |0.7.0           |0.6.0 +                    |2025-04-11|Fix filling the seed desert region behind the last seed of a genome.                                                   |
|3.3          |0.6.0           |0.6.0 +                    |2025-03-25|Reduce index size for batches <= 512. Add the total bases of index to info.toml for computing the Evalue. Denser seeds.|
|3.1          |0.5.0           |0.4.0 +                    |2024-12-18|Change the default partitions of seed data index.                                                                      |
|3.0          |0.4.0           |0.4.0 +                    |2024-08-15|Support suffix matching of seeds. Better seed desert filling for highly-repetitive regions. Denser seeds.              |
|1.1          |0.3.0           |0.3.0                      |2024-05-14|Change the format of seed data index. Use longer contig intervals.                                                     |
|0.1          |0.1.0           |0.1.0 - 0.2.0              |2024-01-25|First version.                                                                                                         |


**What's next:** {{< button size="small" relref="tutorials/search" >}}Searching{{< /button >}}
