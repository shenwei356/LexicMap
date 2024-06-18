---
title: Building an index
weight: 0
---

## Table of contents

{{< toc format=html >}}

## TL;DR

1. Prepare input files:
    - **Sequences of each reference genome should be saved in separate FASTA/Q files, with identifiers in the file names**.
      E.g., GCF_000006945.2.fna.gz
2. Run:
    - From a directory with multiple genome files:

          lexicmap index -I genomes/ -O db.lmi

    - From a file list with one file per line:

          lexicmap index -X files.txt -O db.lmi


## Input

{{< hint type=note >}}
**Genome size**\
LexicMap is only suitable for small genomes like Archaea, Bacteria, Viruses and plasmids.
{{< /hint >}}

**Sequences of each reference genome should be saved in separate FASTA/Q files, with identifiers in the file names**.

- **File type**: FASTA/Q files, in plain text or gzip/xz/zstd/bzip2 compressed formats.
- **File name**: "Genome ID" + "File extention". E.g., `GCF_000006945.2.fna.gz`.
    - **Genome ID**: they should be distinct for accurate result interpretation, which will be shown in the search result.
    - File extention: a regular expression set by the flag `-N/--ref-name-regexp` is used to extract genome IDs from the file name.
      The default value supports common sequence file extentions, e.g., `.fa`, `.fasta`, `.fna`, `.fa.gz`, `.fasta.gz`, `.fna.gz`, `fasta.xz`, `fasta.zst`, and `fasta.bz2`.
    - [brename](https://github.com/shenwei356/brename) can help to batch rename files safely.
    - If you don't want to change the original file names, you can
        1. Create and change to a new directory.
        2. Create symbolic links (`ln -s`) for all genome files.
        3. Batch rename all the symbolic links with brename.
        4. Use this directory as input via the flag `-I/--in-dir`.
- **Sequences**:
    - **Only DNA or RNA sequences are supported**.
    - **Sequence IDs** should be distinct for accurate result interpretation, which will be shown in the search result.
    - One or more sequences in each file are allowed.
        - Unwanted sequences can be filtered out by regular expressions from the flag `-B/--seq-name-filter`.
    - **Genome size limit**. Some none-isolate assemblies might have extremely large genomes, e.g., [GCA_000765055.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000765055.1/) has >150 Mb.
     The flag `-g/--max-genome` (default 15 Mb) is used to skip these input files, and the file list would be written to a file
     via the flag `-G/--big-genomes`.
- **At most 17,179,869,184 (2<sup>34</sup>) genomes are supported**. For more genomes, just build multiple indexes.

**Input files can be given via one of the following ways:**

- **Positional arguments**. For a few input files.
- A **file list** via the flag `-X/--infile-list`  with one file per line.
  It can be STDIN (`-`), e.g., you can filter a file list and pass it to `lexicmap index`.
  *The flag `-S/--skip-file-check` is optional for skiping input file checking if you believe these files do exist*.
- A **directory** containing input files via the flag `-I/--in-dir`.
    - Multiple-level directories are supported.
    - Directory and file symlinks are followed.

## Hardware requirements

See [benchmark of index building](https://bioinf.shenwei.me/LexicMap/introduction/#indexing).

LexicMap is designed to provide fast and low-memory sequence alignment against millions of prokaryotic genomes.

- **CPU:**
    - No specific requirements on CPU type and instruction sets. Both x86 and ARM chips are supported.
    - More is better as LexicMap is a CPU-intensive software. It uses all CPUs by default (`-j/--threads`).
- **RAM**
    - More RAM (> 50 GB) is preferred. The memory usage in index building is mainly related to:
        - The number of masks (`-m/--masks`, default 40,000).
        - The number of genomes.
        - The divergence between genome sequences.
        - **The genome batch size**  (`-b/--batch-size`, default 10,000). This is the main parameter to adjust memory usage.
    - **If the RAM is not sufficient (< 50 GB)**. Please:
        1. **Use a smaller genome batch size**. It decreases indexing memory occupation and has little affection on searching performance.
        2. Use a smaller number of masks, e.g., 20,000 performs well for small genomes (<=5 Mb). And if the queries are long (>= 2kb), there's little affection for the alignment results.
- **Disk**
    - More (>2 TB) is better. The index size is related to the input genomes and the number of masks. See [some examples](#index-size).
    - SSD disks are preferred, while HDD disks are also fast enough.

## Algorithm

<img src="/LexicMap/indexing.svg" alt="" width="900"/>

1. **Generating *m* [LexicHash masks](https://doi.org/10.1093/bioinformatics/btad652)**.

    1. Generate *m* prefixes.
        1. Generating all permutations of *p*-bp prefixes that can cover all possible k-mers, *p* is the biggest value for 4<sup>*p*</sup> <= *m* (desired number of masks), e.g., *p*=7 for 40,000 masks.
        2. Removing low-complexity prefixes. E.g., 16176 out of 16384 (4^7) prefixes are left.
        3. Duplicating these prefixes to *m* prefixes.
    2. For each prefix,
        1. Randomly generating left *k*-*p* bases.
        2. If the *P*-prefix (`-p/--seed-min-prefix`) is of low-complexity, re-generating. *P* is the minimum length of substring matches, default 15.
        3. If the mask is duplicated, re-generating.

2. **Building an index for each genome batch** (`-b/--batch-size`, default 10,000, max 131,072).

    1. For each genome file in a genome batch.
        1. Optionally discarding sequences via regular expression (`-B/--seq-name-filter`).
        2. Skipping genomes bigger than the value of `-g/--max-genome`.
        3. Concatenating all sequences, with intervals of 1000-bp N's.
        4. Capturing the most similar k-mer for each mask and recording the k-mer and its location(s) and strand information. Base N is treated as A.
        5. Filling sketching deserts (genome regions longer than `--seed-max-desert` without any captured k-mers/seeds).
           In a sketching desert, not a single k-mer is captured because there's another k-mer in another place which shares a longer prefix with the mask.
           As a result, for a query similar to seqs in this region, all captured k-mers can’t match the correct seeds.
            1. For a desert region (`start`, `end`), masking the extended region (`start-1000`, `end+1000`) with the masks.
            2. Starting from `start`, every around `--seed-in-desert-dist` (default 200) bp, finding a k-mer which is captured by some mask, and add the k-mer and its position information into the index of that mask.
        6. Saving the concatenated genome sequence (bit-packed, 2 bits for one base, N is treated as A) and genome information (genome ID, size, and lengths of all sequences) into the genome data file, and creating an index file for the genome data file for fast random subsequence extraction.
    2. Compressing k-mers and the corresponding data (k-mer-data, or seeds data, including genome batch, genome number, location, and strand) into chunks of files, and creating an index file for each k-mer-data file for fast seeding.
    3. Writing summary information into `info.toml` file.

3. **Merging indexes of multiple batches**.
    1. For each k-mer-data chunk file (belonging to a list of masks), serially reading data of each mask from all batches,
      merging them and writting to a new file.
    2. For genome data files, just moving them.
    3. Concatenating `genomes.map.bin`, which maps each genome ID to its batch ID and index in the batch.
    4. Update the index summary file.

## Parameters

**Flags in bold text** are important and frequently used.

{{< tabs "t1" >}}

{{< tab "Genome batches" >}}

|Flag                 |Value                      |Function                               |Comment                                                                                                                                                                                                                                                                                |
|:--------------------|:--------------------------|:--------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|**`-b/--batch-size`**|Max: 131072, default: 10000|Maximum number of genomes in each batch|If the number of input files exceeds this number, input files are split into multiple batches and indexes are built for all batches. In the end, seed files are merged, while genome data files are kept unchanged and collected. ► Bigger values increase indexing memory occupation  |

{{< /tab>}}

{{< tab "LexicHash mask generation" >}}

|Flag                  |Value               |Function                                                  |Comment                                                                                                                                                                                     |
|:---------------------|:-------------------|:---------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|`-M/--mask-file`      |A file              |File with custom masks                                    |File with custom masks, which could be exported from an existing index or newly generated by "lexicmap utils masks". This flag oversides `-k/--kmer`, `-m/--masks`, `-s/--rand-seed`, et al.|
|**`-k/--kmer`**       |Max: 32, default: 31|K-mer size                                                |► Bigger values improve the search specificity and do not increase the index size.                                                                                                          |
|**`-m/--masks`**      |Default: 40000      |Number of masks                                           |► Bigger values improve the search sensitivity, increase the index size, and slow down the search speed.                                                                                    |
|`-p/--seed-min-prefix`|Max: 32, Default: 15|Minimum length of shared substrings (anchors) in searching|This value is used to remove masks with a prefix of low-complexity.                                                                                                                         |


{{< /tab>}}


{{< tab "Seeds (k-mer-value) data" >}}

|Flag                   |Value                       |Function                                        |Comment                                                                                                                                                                                                                                                                                                                                                                                  |
|:----------------------|:---------------------------|:-----------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|**`--seed-max-desert`**|Default: 450                |Maximum length of distances between seeds       |The default value of 450 guarantees queries >450 bp would match at least one seed. ► Smaller values improve the search sensitivity and slightly increase the index size. ►Large regions with no seeds are called sketching deserts. Deserts with seed distance larger than this value will be filled by choosing k-mers roughly every `--seed-in-desert-dist` (150 by default) bases.    |
|`-c/--chunks`          |Maximum: 128, default: #CPUs|Number of seed file chunks                      |Bigger values accelerate the search speed at the cost of a high disk reading load. The maximum number should not exceed the maximum number of open files set by the operating systems.                                                                                                                                                                                                   |
|`-p/--partitions`      |Default: 512                |Number of partitions for indexing each seed file|Bigger values bring a little higher memory occupation. 512 is a good value with high searching speed, larger or smaller values would decrease the speed in `lexicmap search`. ► After indexing, `lexicmap utils reindex-seeds` can be used to reindex the seeds data with  another value of this flag.                                                                                   |
|`--max-open-files`     |Default: 512                |Maximum number of open files                    |It's only used in merging indexes of multiple genome batches.                                                                                                                                                                                                                                                                                                                            |

{{< /tab>}}

{{< /tabs >}}

Also see the [usage](https://bioinf.shenwei.me/LexicMap/usage/#index) of `lexicmap index`.

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

    It would take about 2 seconds and 2 GB RAM in a 16-CPU PC.

    Optionally, we can also use **a file list** as the input.

        $ head -n 3 files.txt
        refs/GCF_000006945.2.fa.gz
        refs/GCF_000017205.1.fa.gz
        refs/GCF_000148585.2.fa.gz

        lexicmap index -X files.txt -O demo.lmi

{{< expand "Click to show the log of a demo run." "..." >}}

    # here we set a small --batch-size 5
    $ lexicmap index -I refs/ -O demo.lmi --batch-size 5
    20:18:51.410 [INFO] LexicMap v0.4.0
    20:18:51.410 [INFO]   https://github.com/shenwei356/LexicMap
    20:18:51.410 [INFO]
    20:18:51.410 [INFO] checking input files ...
    20:18:51.410 [INFO]   15 input file(s) given
    20:18:51.410 [INFO]
    20:18:51.410 [INFO] --------------------- [ main parameters ] ---------------------
    20:18:51.410 [INFO]
    20:18:51.410 [INFO] input and output:
    20:18:51.410 [INFO]   input directory: refs/
    20:18:51.410 [INFO]     regular expression of input files: (?i)\.(f[aq](st[aq])?|fna)(\.gz|\.xz|\.zst|\.bz2)?$
    20:18:51.410 [INFO]     *regular expression for extracting reference name from file name: (?i)(.+)\.(f[aq](st[aq])?|fna)(\.gz|\.xz|\.zst|\.bz2)?$
    20:18:51.410 [INFO]     *regular expressions for filtering out sequences: []
    20:18:51.410 [INFO]   max genome size: 15000000
    20:18:51.410 [INFO]   output directory: demo.lmi
    20:18:51.410 [INFO]
    20:18:51.410 [INFO] mask generation:
    20:18:51.410 [INFO]   k-mer size: 31
    20:18:51.410 [INFO]   number of masks: 40000
    20:18:51.410 [INFO]   rand seed: 1
    20:18:51.410 [INFO]   prefix length for checking low-complexity in mask generation: 15
    20:18:51.410 [INFO]
    20:18:51.410 [INFO] seed data:
    20:18:51.410 [INFO]   maximum sketching desert length: 900
    20:18:51.410 [INFO]   distance of k-mers to fill deserts: 200
    20:18:51.410 [INFO]   seeds data chunks: 16
    20:18:51.410 [INFO]   seeds data indexing partitions: 512
    20:18:51.410 [INFO]
    20:18:51.410 [INFO] general:
    20:18:51.410 [INFO]   genome batch size: 5
    20:18:51.410 [INFO]
    20:18:51.410 [INFO]
    20:18:51.410 [INFO] --------------------- [ generating masks ] ---------------------
    20:18:51.728 [INFO]
    20:18:51.728 [INFO] --------------------- [ building index ] ---------------------
    20:18:51.884 [INFO]
    20:18:51.884 [INFO]   ------------------------[ batch 0 ]------------------------
    20:18:51.884 [INFO]   building index for batch 0 with 5 files...
    processed files:  5 / 5 [======================================] ETA: 0s. done
    20:18:52.463 [INFO]   writing seeds...
    20:18:52.606 [INFO]   finished writing seeds in 143.292864ms
    20:18:52.606 [INFO]   finished building index for batch 0 in: 722.316862ms
    20:18:52.606 [INFO]
    20:18:52.606 [INFO]   ------------------------[ batch 1 ]------------------------
    20:18:52.606 [INFO]   building index for batch 1 with 5 files...
    processed files:  5 / 5 [======================================] ETA: 0s. done
    20:18:54.073 [INFO]   writing seeds...
    20:18:54.251 [INFO]   finished writing seeds in 178.055055ms
    20:18:54.251 [INFO]   finished building index for batch 1 in: 1.645034631s
    20:18:54.251 [INFO]
    20:18:54.251 [INFO]   ------------------------[ batch 2 ]------------------------
    20:18:54.251 [INFO]   building index for batch 2 with 5 files...
    processed files:  5 / 5 [======================================] ETA: 0s. done
    20:18:54.993 [INFO]   writing seeds...
    20:18:55.182 [INFO]   finished writing seeds in 188.702811ms
    20:18:55.182 [INFO]   finished building index for batch 2 in: 931.156662ms
    20:18:55.183 [INFO]
    20:18:55.183 [INFO] merging 3 indexes...
    20:18:55.183 [INFO]   [round 1]
    20:18:55.183 [INFO]     batch 1/1, merging 3 indexes to demo.lmi.tmp/r1_b1
    20:18:55.741 [INFO]   [round 1] finished in 558.606877ms
    20:18:55.741 [INFO] rename demo.lmi.tmp/r1_b1 to demo.lmi
    20:18:55.746 [INFO]
    20:18:55.746 [INFO] finished building LexicMap index from 15 files with 40000 masks in 4.336354688s
    20:18:55.746 [INFO] LexicMap index saved: demo.lmi
    20:18:55.746 [INFO]
    20:18:55.746 [INFO] elapsed time: 4.336386606s
    20:18:55.746 [INFO]

{{< /expand >}}

## Output

The LexicMap index is a directory with multiple files.

### File structure

    $ tree demo.lmi/
    demo.lmi/                    # the index directory
    ├── genomes                  # directory of genome data
    │   └── batch_0000           # genome data of one batch
    │       ├── genomes.bin      # genome data file, containing genome ID, size, sequence lengths, bit-packed sequences
    │       └── genomes.bin.idx  # index of genome data file, for fast subsequence extraction
    ├── seeds                    # seed data: pairs of k-mer and its location information (genome batch, genome number, location, strand)
    │   ├── chunk_000.bin        # seed data file
    │   ├── chunk_000.bin.idx    # index of seed data file, for fast seed searching and data extraction
    ... ... ...
    │   ├── chunk_015.bin        # the number of chunks is set by flag `-c/--chunks`, default: #cpus
    │   └── chunk_015.bin.idx
    ├── genomes.map.bin          # mapping genome ID to batch number of genome number in the batch
    ├── info.toml                # summary of the index
    └── masks.bin                # mask data

### Index size


{{< tabs "t2" >}}

{{< tab "Demo data" >}}

    # 15 genomes
    demo.lmi/: 26.87 MB
      13.64 MB      seeds
      12.93 MB      genomes
     312.53 KB      masks.bin
      375.00 B      genomes.map.bin
      260.00 B      info.toml

{{< /tab>}}

{{< tab "GTDB repr" >}}

    # 85,205 genomes
    gtdb_repr.lmi: 110.37 GB
      66.78 GB      genomes
      43.59 GB      seeds
       2.03 MB      genomes.map.bin
     312.53 KB      masks.bin
      266.00 B      info.toml

{{< /tab>}}

{{< tab "GTDB complete" >}}

    # 402,538 genomes
    gtdb_complete.lmi: 509.99 GB
     362.98 GB      genomes
     147.00 GB      seeds
       9.60 MB      genomes.map.bin
     312.53 KB      masks.bin
      269.00 B      info.toml

{{< /tab>}}


{{< tab "Genbank+RefSeq" >}}

    # 2,340,672 genomes
    genbank_refseq.lmi: 2.91 TB
       2.17 TB      genomes
     764.22 GB      seeds
      55.81 MB      genomes.map.bin
     312.53 KB      masks.bin
      270.00 B      info.toml

{{< /tab>}}


{{< tab "AllTheBacteria HQ" >}}

    # 1,858,610 genomes
    atb_hq.lmi: 2.32 TB
       1.77 TB      genomes
     569.44 GB      seeds
      39.22 MB      genomes.map.bin
     312.53 KB      masks.bin
      270.00 B      info.toml

{{< /tab>}}

{{< /tabs >}}

- Directory/file sizes are counted with https://github.com/shenwei356/dirsize.
- Index building parameters: `-k 31 -m 40000`. Genome batch size: `-b 10000` for GTDB datasets, `-b 50000` for others.

**What's next:** {{< button size="small" relref="tutorials/search" >}}Searching{{< /button >}}
