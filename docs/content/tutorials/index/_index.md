---
title: Building an index
weight: 0
---

## Table of contents

{{< toc format=html >}}

## TL;DR

1. Prepare input files:
    - **Sequences of each reference genome should be saved in separate FASTA/Q files, with identifiers in the file names**.
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
- **At most 17,179,869,184 genomes are supported**. For more genomes, just build multiple indexes.

Input files can be given via one of the following ways:

- **Positional arguments**. For a few input files.
- A **file list** via the flag `-X/--infile-list`  with one file per line.
  It can be STDIN (`-`), e.g., you can filter a file list and pass it to `lexicmap`.
- A **directory** containing input files via the flag `-I/--in-dir`.
    - Multiple-level directories are supported.
    - Directory and file symlinks are followed.

## Hardware requirements

LexicMap is designed to provide fast and low-memory sequence alignment against millions of prokaryotic genomes.

- **CPU:**
    - No specific requirements on CPU type and instruction sets. Both x86 and ARM chips are supported.
    - More is better as LexicMap is a CPU-intensive software. It uses all CPUs by default (`-j/--threads`).
- **RAM**
    - More RAM (> 50 GB) is preferred. The memory usage in index building is mainly related to:
        - The number of masks (`-m, --masks`, default 40,000).
        - The number of genome.
        - The divergence between genome sequences.
        - The genome batch size  (`-b/--batch-size`, default 10,000).
    - If the RAM is not sufficient (< 50 GB). Please:
        1. Use a smaller genome batch size. It decreases indexing memory occupation and has little effect on searching performance.
        2. Use a smaller number of masks, e.g., 20,000 performs well for small genomes (<=5 Mb). And if the queries are long (>= 2kb), there's little affection for the alignment results.
- **Disk**
    - More (>2 TB) is better. The index size is related to the input genomes and the number of masks. See [some examples](#index-size).
    - SSD disks are preferred, while HDD disks are also fast enough.

## Algorithm

<img src="/LexicMap/indexing.svg" alt="" width="900"/>

1. **Generating [LexicHash masks](https://doi.org/10.1093/bioinformatics/btad652)**.

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

2. **Building an index for each genome batch** (`-b/--batch-size`, default 10000, max 131072).

    1. For each genome file in a genome batch.
        1. Optionally discarding sequences via regular expression (`-B/--seq-name-filter`).
        2. Skipping genomes bigger than the value of `-g/--max-genome`.
        3. Concatenating all sequences, with intervals of 1000-bp N's.
        4. Capureing the most similar k-mer for each mask and saving k-mer and its location(s) and strand information.
        5. Saving the concatenated genome sequence (bit-packed, 2 bits for one base) and genome information (genome ID, size, and lengths of all sequences) into the genome data file, and creating an index file for the genome data file for fast random subsequence extraction.
    2. Compressing k-mer and the corresponding data (k-mer-data, or seeds data, including genome batch, genome number, location, and strand) into chunks of files, and creating an index file for each k-mer-data file for fast seeding.
    3. Writing summary information into `info.toml` file.

3. **Merging indexes of multiple batches**.

## Parameters

**Flags in bold text** are important and frequently used.

{{< tabs "t1" >}}

{{< tab "Genome batches" >}}

|Flag                 |Value                          |Function                               |Comment                                                                                                                                                                                                                                                                                    |
|:--------------------|:------------------------------|:--------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|**`-b/--batch-size`**|Maximum: 131072, default: 10000|Maximum number of genomes in each batch|If the number of input files exceeds this number, input files are into multiple batches and indexes are built for all batches. Next, seed files are merged into a big one, while genome data files are kept unchanged and collected. ► Bigger values increase indexing memory occupation.  |

{{< /tab>}}

{{< tab "LexicHash mask generation" >}}

|Flag             |Value                   |Function                                       |Comment                                                                                                                                                   |
|:----------------|:-----------------------|:----------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------|
|`-M/--mask-file` |A file                  |File with custom masks                         |It could be generated by `lexicmap utils gen-masks`. This flag oversides `-k/--kmer`, `-m/--masks`, `-s/--rand-seed`, `-n/--top-n`, and `-P/--prefix-ext`.|
|**`-k/--kmer`**  |Maximum: 32, default: 31|K-mer size                                     |Bigger values improve the search specificity and do not increase the index size.                                                                          |
|**`-m/--masks`** |Default: 40000          |Number of masks                                |Bigger values improve the search sensitivity and increase the index size.                                                                                 |
|**`-n/--top-n`** |Default: 20             |The top N largest genomes for generating masks.|Bigger values increase the indexing memory occupation (1.0~1.5 GB per genome for 20k masks).                                                              |
|`-P/--prefix-ext`|Default: 8              |Extension length of prefixes                   |Bigger values improve the search sensitivity by decreasing the maximum seed distances.                                                                    |                                                               |


{{< /tab>}}


{{< tab "Seeds (k-mer-value) data" >}}

|Flag              |Value                       |Function                                        |Comment                                                                                                                                                                                                                                                                                        |
|:-----------------|:---------------------------|:-----------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|`-c/--chunks`     |Maximum: 128, default: #CPUs|Number of seed file chunks                      |Bigger values accelerate the search speed at the cost of a high disk reading load. The maximum number should not exceed the maximum number of open files set by the operating systems.                                                                                                         |
|`-p/--partitions` |Default: 512                |Number of partitions for indexing each seed file|Bigger values bring higher memory occupation. 512 is a good value with high searching speed, larger or smaller values would decrease the speed in `lexicmap search`. ► After indexing, `lexicmap utils reindex-seeds` can be used to reindex the seeds data with  another value of this flag.  |
|`--max-open-files`|Default: 512                |Maximum number of open files                    |It's only used in merging indexes of multiple genome batches.                                                                                                                                                                                                                                  |

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

    It would take about 3 minutes and 10 GB RAM in a 16-CPU PC. You can use a smaller value of `-n/--top-n` to decrease memory usage.

    Optionally, we can also use **a file list** as the input.

        $ head -n 3 files.txt
        refs/GCF_000006945.2.fa.gz
        refs/GCF_000017205.1.fa.gz
        refs/GCF_000148585.2.fa.gz

        lexicmap index -X files.txt -O demo.lmi

{{< expand "Click to show the log of a demo run." "..." >}}

    $ lexicmap index -I refs/ -O demo.lmi --top-n 4 --batch-size 5
    21:27:50.928 [INFO] LexicMap v0.3.0
    21:27:50.928 [INFO]   https://github.com/shenwei356/LexicMap
    21:27:50.928 [INFO]
    21:27:50.928 [INFO] checking input files ...
    21:27:50.928 [INFO]   15 input file(s) given
    21:27:50.928 [INFO]
    21:27:50.928 [INFO] --------------------- [ main parameters ] ---------------------
    21:27:50.928 [INFO]
    21:27:50.928 [INFO] input and output:
    21:27:50.928 [INFO]   input directory: refs/
    21:27:50.928 [INFO]     regular expression of input files: (?i)\.(f[aq](st[aq])?|fna)(.gz)?$
    21:27:50.928 [INFO]     *regular expression for extracting reference name from file name: (?i)(.+)\.(f[aq](st[aq])?|fna)(.gz)?$
    21:27:50.928 [INFO]     *regular expressions for filtering out sequences: []
    21:27:50.928 [INFO]   max genome size: 15000000
    21:27:50.928 [INFO]   output directory: demo.lmi
    21:27:50.928 [INFO]
    21:27:50.928 [INFO] k-mer size: 31
    21:27:50.928 [INFO] number of masks: 40000
    21:27:50.928 [INFO] rand seed: 1
    21:27:50.928 [INFO] top N genomes for generating mask: 4
    21:27:50.928 [INFO] prefix extension length: 8
    21:27:50.928 [INFO]
    21:27:50.928 [INFO] seeds data chunks: 16
    21:27:50.928 [INFO] seeds data indexing partitions: 512
    21:27:50.928 [INFO]
    21:27:50.928 [INFO] genome batch size: 5
    21:27:50.928 [INFO]
    21:27:50.928 [INFO]
    21:27:50.928 [INFO] --------------------- [ generating masks ] ---------------------
    21:27:50.928 [INFO]
    21:27:50.928 [INFO] generating masks from the top 4 out of 15 genomes...
    21:27:50.928 [INFO]
    21:27:50.928 [INFO]   checking genomes sizes of 15 files...
    processed files:  15 / 15 [======================================] ETA: 0s. done
    21:27:50.979 [INFO]     0 genomes longer than 15000000 are filtered out
    21:27:50.979 [INFO]     genome size range in the top 4 files: [4951383, 6588339]
    21:27:50.979 [INFO]
    21:27:50.979 [INFO]   collecting k-mers from 4 files...
    processed files: 4/4
    21:27:57.633 [INFO]
    21:27:57.633 [INFO]   generating masks...
    21:27:57.638 [INFO]     generating 16384 masks covering all 7-bp prefixes...
    processed prefixes: 16384/16384
    21:28:19.639 [INFO]     generating left 23616 masks...
    processed prefixes: 23616/23616
    21:28:45.343 [INFO]
    21:28:45.343 [INFO]   maximum distance between seeds:
    21:28:45.345 [INFO]     GCF_000006945.2: 1092
    21:28:45.346 [INFO]     GCF_003697165.2: 1151
    21:28:45.347 [INFO]     GCF_000742135.1: 1406
    21:28:45.348 [INFO]     GCF_000017205.1: 1050
    21:28:45.358 [INFO]
    21:28:45.358 [INFO]   finished generating masks in: 54.430001475s
    21:28:45.359 [INFO]
    21:28:45.359 [INFO] --------------------- [ building index ] ---------------------
    21:28:45.518 [INFO]
    21:28:45.518 [INFO]   ------------------------[ batch 0 ]------------------------
    21:28:45.518 [INFO]   building index for batch 0 with 5 files...
    processed files:  5 / 5 [======================================] ETA: 0s. done
    21:28:46.720 [INFO]   writing seeds...
    21:28:46.864 [INFO]   finished writing seeds in 144.463288ms
    21:28:46.864 [INFO]   finished building index for batch 0 in: 1.345988413s
    21:28:46.864 [INFO]
    21:28:46.864 [INFO]   ------------------------[ batch 1 ]------------------------
    21:28:46.864 [INFO]   building index for batch 1 with 5 files...
    processed files:  5 / 5 [======================================] ETA: 0s. done
    21:28:48.637 [INFO]   writing seeds...
    21:28:48.816 [INFO]   finished writing seeds in 178.371487ms
    21:28:48.816 [INFO]   finished building index for batch 1 in: 1.951457718s
    21:28:48.816 [INFO]
    21:28:48.816 [INFO]   ------------------------[ batch 2 ]------------------------
    21:28:48.816 [INFO]   building index for batch 2 with 5 files...
    processed files:  5 / 5 [======================================] ETA: 0s. done
    21:28:50.171 [INFO]   writing seeds...
    21:28:50.347 [INFO]   finished writing seeds in 175.806993ms
    21:28:50.347 [INFO]   finished building index for batch 2 in: 1.53094663s
    21:28:50.347 [INFO]
    21:28:50.347 [INFO] merging 3 indexes...
    21:28:50.347 [INFO]   [round 1]
    21:28:50.347 [INFO]     batch 1/1, merging 3 indexes to demo.lmi.tmp/r1_b1
    21:28:50.919 [INFO]   [round 1] finished in 572.049922ms
    21:28:50.919 [INFO] rename demo.lmi.tmp/r1_b1 to demo.lmi
    21:28:50.924 [INFO]
    21:28:50.924 [INFO] finished building LexicMap index from 15 files with 40000 masks in 59.996543239s
    21:28:50.924 [INFO] LexicMap index saved: demo.lmi
    21:28:50.924 [INFO]
    21:28:50.924 [INFO] elapsed time: 59.996570756s
    21:28:50.924 [INFO]

{{< /expand >}}

## Output

The LexicMap index is a directory with multiple files.

### File structure

    $ tree demo.lmi/
    demo.lmi/                    # the index directory
    ├── genomes                  # directory of genome data
    │   └── batch_0000           # genome data of a batch
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
    ├── masks.bin                # mask data

### Index size


{{< tabs "t2" >}}

{{< tab "Demo data" >}}

    # 15 genomes
    $ dirsize  demo.lmi/
    demo.lmi/: 26.63 MB
      13.39 MB      seeds
      12.93 MB      genomes
     312.53 KB      masks.bin
      375.00 B      genomes.map.bin
      261.00 B      info.toml

{{< /tab>}}

{{< tab "GTDB repr" >}}

    # 85,205 genomes
    gtdb_repr.lmi: 109.16 GB
      66.79 GB      genomes
      42.37 GB      seeds
       2.03 MB      genomes.map.bin
     312.53 KB      masks.bin
      267.00 B      info.toml

{{< /tab>}}

{{< tab "GTDB complete" >}}

    # 402,538 genomes
    gtdb_complete.lmi: 506.57 GB
     362.99 GB      genomes
     143.57 GB      seeds
       9.60 MB      genomes.map.bin
     312.53 KB      masks.bin
      269.00 B      info.toml

{{< /tab>}}


{{< tab "Genbank+RefSeq" >}}

    # 2,340,672 genomes
    genbank_refseq.lmi: 2.90 TB
       2.17 TB      genomes
     754.04 GB      seeds
      55.81 MB      genomes.map.bin
     312.53 KB      masks.bin
      271.00 B      info.toml

{{< /tab>}}


{{< tab "AllTheBacteria HQ" >}}

    # 1,858,610 genomes
    2kk-HQ.lmi: 2.32 TB
       1.77 TB      genomes
     563.16 GB      seeds
      39.22 MB      genomes.map.bin
     312.53 KB      masks.bin
      271.00 B      info.toml

{{< /tab>}}

{{< /tabs >}}

Index building parameters: `-k 31 -m 40000`. Genome batch size: `-b 10000` for GTDB datasets, `-b 131072` for others.

**What's next:** {{< button size="small" relref="tutorials/search" >}}Searching{{< /button >}}
