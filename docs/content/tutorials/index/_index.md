---
title: Step 1. Building a database
weight: 0
---

Terminology differences:

- On this page and in the LexicMap command line options, the term **"mask"** is used, following the terminology in the LexicHash paper.
- In the LexicMap manuscript, however, we use **"probe"** as it is easier to understand.
  Because these masks, which consist of thousands of k-mers and capture k-mers from sequences through prefix matching, function similarly to DNA probes in molecular biology.

## Table of contents

{{< toc format=html >}}

## TL;DR

1. Prepare input files:
    - **Sequences of each reference genome should be saved in separate FASTA/Q files, with identifiers (no tab symbols) in the file names**.
      E.g., GCF_000006945.2.fna.gz
        - A regular expression is also available to extract reference id from the file name.
          E.g., `--ref-name-regexp '^(\w{3}_\d{9}\.\d+)'` extracts `GCF_000006945.2` from GenBank assembly file `GCF_000006945.2_ASM694v2_genomic.fna.gz`
    - While if you save *a few* **small** (viral) **complete** genomes (one sequence per genome) in each file, it's feasible as sequence IDs in search result can help to distinguish targe genomes.
2. Run:
    - From a directory with multiple genome files:

          lexicmap index -I genomes/ -O db.lmi

    - From a file list with one file per line:

          lexicmap index -X files.txt -O db.lmi


## Input

{{< hint type=note >}}
**Genome size**\
LexicMap is mainly suitable for small genomes like Archaea, Bacteria, Viruses and plasmids.

Maximum genome size: 268 Mb (268,435,456).
More precisely:

    $total_bases + ($num_contigs - 1) * 1000 <= 268,435,456

as we concatenate contigs with 1000-bp intervals of N’s to reduce the sequence scale to index.

{{< /hint >}}


**Sequences of each reference genome should be saved in separate FASTA/Q files, with identifiers in the file names**. While if you save *a few* **small** (viral) **complete** genomes (one sequence per genome) in each file, it's feasible as sequence IDs in search result can help to distinguish targe genomes.

- **File type**: FASTA/Q files, in plain text or gzip/xz/zstd/bzip2 compressed formats.
- **File name**: "Genome ID" + "File extention". E.g., `GCF_000006945.2.fna.gz`.
    - **Genome ID**: **they must not contain tab ("\t") symbols, and should be distinct for accurate result interpretation**, which will be shown in the search result.
        - A regular expression is also available to extract reference id from the file name.
          E.g., `--ref-name-regexp '^(\w{3}_\d{9}\.\d+)'` extracts `GCF_000006945.2` from GenBank assembly file `GCF_000006945.2_ASM694v2_genomic.fna.gz`
- **File extention**: a regular expression set by the flag `-r/--file-regexp` is used to match input files.
      The default value supports common sequence file extentions, e.g., `.fa`, `.fasta`, `.fna`, `.fa.gz`, `.fasta.gz`, `.fna.gz`, `fasta.xz`, `fasta.zst`, and `fasta.bz2`.
- **Sequences**:
    - **Only DNA or RNA sequences are supported**.
    - **Sequence IDs** should be distinct for accurate result interpretation, which will be shown in the search result.
    - Sequence description (text behind sequence ID) is not saved. If you do need it, you can create a mapping file
      (`seqkit seq -n ref.fa.gz | sed -E 's/\s+/\t/' > id2desc.tsv`) and use it to [add description in search result](https://bioinf.shenwei.me/LexicMap/tutorials/search/#summarizing-results).
    - **One or more sequences (contigs) in each file are allowed**.
        - Unwanted sequences can be filtered out by regular expressions from the flag `-B/--seq-name-filter`.
    - **Genome size limit**. Some none-isolate assemblies might have extremely large genomes, e.g., [GCA_000765055.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000765055.1/) has >150 Mb.
     The flag `-g/--max-genome` (default 15 Mb) is used to skip these input files, and the file list would be written to a file
     via the flag `-G/--big-genomes`.
        - **Changes since v0.5.0**:
            - Genomes with any single contig larger than the threshold will be skipped as before.
            - However, **fragmented (with many contigs) genomes with the total bases larger than the threshold will
              be split into chunks** and alignments from these chunks will be merged in "lexicmap search".
        - **For fungi genomes, please increase the value**.
    - **Minimum sequence length**. A flag `-l/--min-seq-len` can filter out sequences shorter than the threshold (default is the `k` value).
- **At most 17,179,869,184 (2<sup>34</sup>) genomes are supported**. For more genomes, please create a file list and split it into multiple parts, and build an index for each part.

**Input files can be given via one of the following ways:**

- **Positional arguments**. For a few input files.
- A **file list** via the flag `-X/--infile-list`  with one file per line.
  **It can be STDIN (`-`)**, e.g., you can filter a file list and pass it to `lexicmap index`.
  *The flag `-S/--skip-file-check` is optional for skiping input file checking if you believe these files do exist*.
- A **directory** containing input files via the flag `-I/--in-dir`.
    - **Multiple-level directories are supported**. So you don't need to saved hundreds of thousand files into one directoy.
    - **Directory and file symlinks are followed**.

## Hardware requirements

See [benchmark of index building](https://bioinf.shenwei.me/LexicMap/introduction/#indexing).

LexicMap is designed to provide fast and low-memory sequence alignment against millions of prokaryotic genomes.

- **CPU:**
    - No specific requirements on CPU type and instruction sets. Both x86 and ARM chips are supported.
    - More is better as LexicMap is a CPU-intensive software. **It uses all CPUs by default (`-j/--threads`)**.
- **RAM**
    - More RAM (> 100 GB) is preferred. The memory usage in index building is mainly related to:
        - **The number of masks** (`-m/--masks`, default 40,000). Bigger values improve the search sensitivity, increase the index size, and slow down the search speed. For smaller genomes like phages/viruses, m=10,000 is high enough.
        - **The number of genomes**. More genomes consume more memory.
        - **The divergence between genome sequences in each batch**. Diverse genomes consume more memory.
        - **The genome batch size**  (`-b/--batch-size`, default 5,000). **This is the main parameter to adjust memory usage**. Bigger values increase indexing memory occupation.
        - **The maximum seed distance** or **the maximum sketching desert size** (`-D/--seed-max-desert`, default 200),
          and the distance of k-mers to fill deserts (`-d/--seed-in-desert-dist`, default 50).
          Bigger `-D/--seed-max-desert` values decrease the search sensitivity for distant targets, speed up the indexing speed,
          decrease the indexing memory occupation and decrease the index size. While the alignment speed is almost not affected.
    - **If the RAM is not sufficient**. Please:
        - **Use a smaller genome batch size**. It decreases indexing memory occupation and has little affection on searching performance.
        - Use a smaller number of masks, e.g., 20,000 performs well for small genomes (<=5 Mb). And if the queries are long (>= 2kb), there's little affection for the alignment results.
- **Disk**
    - More is better. LexicMap index size is related to the number of input genomes, the divergence between genome sequences, the number of masks, and the maximum seed distance. See [some examples](#index-size).
        - **Note that the index size is not linear with the number of genomes, it's sublinear**. Because the seed data are compressed with VARINT-GB algorithm, more genomes bring higher compression rates.
    - SSD disks are preferred, while HDD disks are also fast enough.

## Algorithm

<img src="/LexicMap/indexing.svg" alt="" width="900"/>

{{< expand "Click to show details." "..." >}}

1. **Generating *m* [LexicHash masks](https://doi.org/10.1093/bioinformatics/btad652)**.

    1. Generate *m* prefixes.
        1. Generating all permutations of *p*-bp prefixes that can cover all possible k-mers, *p* is the biggest value for 4<sup>*p*</sup> <= *m* (desired number of masks), e.g., *p*=7 for 40,000 masks. (4<sup>*7*</sup> = 16384)
        3. Duplicating these prefixes to *m* prefixes.
    2. For each prefix,
        1. Randomly generating left *k*-*p* bases.
        3. If the mask is duplicated, re-generating.

2. **Building an index for each genome batch** (`-b/--batch-size`, default 5,000, max 131,072).

    1. For each genome file in a genome batch.
        1. Optionally discarding sequences via regular expression (`-B/--seq-name-filter`).
        2. Skipping genomes bigger than the value of `-g/--max-genome`.
        3. Concatenating all sequences, with intervals of 1000-bp N's.
        4. Capturing the most similar k-mer (in non-gap and non-interval regions) for each mask and recording the k-mer and its location(s) and strand information. Base N is treated as A.
        5. Filling sketching deserts (genome regions longer than `--seed-max-desert` [default 200] without any captured k-mers/seeds).
           In a sketching desert, not a single k-mer is captured because there's another k-mer in another place which shares a longer prefix with the mask.
           As a result, for a query similar to seqs in this region, all captured k-mers can’t match the correct seeds.
            1. For a desert region (`start`, `end`), masking the extended region (`start-1000`, `end+1000`) with the masks.
            2. Starting from `start`, every around `--seed-in-desert-dist` (default 50) bp, finding a k-mer which is captured by some mask, and adding the k-mer and its position information into the index of that mask.
        6. Saving the concatenated genome sequence (bit-packed, 2 bits for one base, N is treated as A) and genome information (genome ID, size, and lengths of all sequences) into the genome data file, and creating an index file for the genome data file for fast random subsequence extraction.
    2. Duplicate and reverse all k-mers, and save each reversed k-mer along with the duplicated position information in the seed data of the closest (sharing the longgest prefix) mask. This is for suffix matching of seeds.
    2. Compressing k-mers and the corresponding data (k-mer-data, or seeds data, including genome batch, genome number, location, and strand) into chunks of files, and creating an index file for each k-mer-data file for fast seeding.
    3. Writing summary information into `info.toml` file.

3. **Merging indexes of multiple batches**.
    1. For each k-mer-data chunk file (belonging to a list of masks), serially reading data of each mask from all batches,
      merging them and writting to a new file.
    2. For genome data files, just moving them.
    3. Concatenating `genomes.map.bin`, which maps each genome ID to its batch ID and index in the batch.
    4. Update the index summary file.

{{</ expand >}}

## Parameters

{{< hint type=note >}}
**Query length**\
LexicMap is mainly designed for sequence alignment with a small number of queries (gene/plasmid/virus/phage sequences) longer than 200 bp by default.
However, short queries can also be aligned.

If you just want to search long (>1kb) queries for highy similar (>95%) targets, you can build an index with a bigger `-D/--seed-max-desert` (200 by default), e.g.,

    --seed-max-desert 450 --seed-in-desert-dist 150

Bigger values decrease the search sensitivity for distant targets, speed up the indexing
speed, decrease the indexing memory occupation and decrease the index size. While the
alignment speed is almost not affected.
{{< /hint >}}

**Flags in bold text** are important and frequently used.

{{< tabs "t1" >}}

{{< tab "General" >}}

|Flag              |Value                      |Function                   |Comment                                                                                                      |
|:-----------------|:--------------------------|:--------------------------|:------------------------------------------------------------------------------------------------------------|
|**`-j/--threads`**|Default: all available CPUs|Number of CPU cores to use.|► If the value is smaller than the number of available CPUs, make sure set the same value to `-c/--chunks`.  |

{{< /tab>}}

{{< tab "Genome batches" >}}

|Flag                 |Value                     |Function                               |Comment                                                                                                                                                                                                                                                                                                                                                                        |
|:--------------------|:-------------------------|:--------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|**`-b/--batch-size`**|Max: 131072, default: 5000|Maximum number of genomes in each batch|If the number of input files exceeds this number, input files are split into multiple batches and indexes are built for all batches. In the end, seed files are merged, while genome data files are kept unchanged and collected. ■ Bigger values increase indexing memory occupation and increase batch searching speed, while single query searching speed is not affected.  |

{{< /tab>}}

{{< tab "LexicHash mask generation" >}}

|Flag            |Value               |Function              |Comment                                                                                                                                                                                   |
|:---------------|:-------------------|:---------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|`-M/--mask-file`|A file              |File with custom masks|File with custom masks, which could be exported from an existing index or newly generated by "lexicmap utils masks". This flag oversides `-k/--kmer`, `-m/--masks`, `-s/--rand-seed`, etc.|
|**`-k/--kmer`** |Max: 32, default: 31|K-mer size            |■ Bigger values improve the search specificity and do not increase the index size.                                                                                                        |
|**`-m/--masks`**|Default: 40,000     |Number of masks       |■  Bigger values improve the search sensitivity, increase the index size, and slow down the search speed. For smaller genomes like phages/viruses, m=10,000 is high enough.               |


{{< /tab>}}


{{< tab "Seeds (k-mer-value) data" >}}

|Flag                        |Value                                       |Function                                                                        |Comment                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
|:---------------------------|:-------------------------------------------|:-------------------------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|**`--seed-max-desert`**     |Default: 200                                |Maximum length of distances between seeds                                       |The default value of 200 guarantees queries >200 bp would match at least one seed. ► Large regions with no seeds are called sketching deserts. Deserts with seed distance larger than this value will be filled by choosing k-mers roughly every --seed-in-desert-dist (50 by default) bases. ■ Bigger values decrease the search sensitivity for distant targets, speed up the indexing speed, decrease the indexing memory occupation and decrease the index size. While the alignment speed is almost not affected.    |
|**`-c/--chunks`**           |Maximum: 128, default: value of -j/--threads|Number of seed file chunks                                                      |Bigger values accelerate the search speed at the cost of a high disk reading load. ► The value should not exceed the maximum number of open files set by the operating systems. ► Make sure the value of `-j/--threads` in `lexicmap search` is >= this value.                                                                                                                                                                                                                                                            |
|**`-J/--seed-data-threads`**|Maximum: -c/--chunks, default: 8            |Number of threads for writing seed data and merging seed chunks from all batches|The actual value is min(--seed-data-threads, max(1, --max-open-files/($batches_1_round + 2))), where $batches_1_round = min(int($input_files / --batch-size), --max-open-files). ■ Bigger values increase indexing speed at the cost of slightly higher memory occupation.                                                                                                                                                                                                                                                |
|`-p/--partitions`           |Default: 1024                               |Number of partitions for indexing each seed file                                |Bigger values bring a little higher memory occupation. ► After indexing, `lexicmap utils reindex-seeds` can be used to reindex the seeds data with  another value of this flag.                                                                                                                                                                                                                                                                                                                                           |
|**`--max-open-files`**      |Default: 768                                |Maximum number of open files                                                    |It's only used in merging indexes of multiple genome batches. If there are >100 batches, i.e., ($input_files / --batch-size), please increase this value and set a bigger `ulimit -n` in shell.                                                                                                                                                                                                                                                                                                                           |

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

    It would take about 6 seconds and 3 GB RAM in a 16-CPU PC.

    Optionally, we can also use **a file list** as the input.

        $ head -n 3 files.txt
        refs/GCF_000006945.2.fa.gz
        refs/GCF_000017205.1.fa.gz
        refs/GCF_000148585.2.fa.gz

        lexicmap index -X files.txt -O demo.lmi

    {{< expand "Click to show the log of a demo run." "..." >}}

        # here we set a small --batch-size 5

        $ lexicmap index -I refs/ -O demo.lmi --batch-size 5
        16:22:49.745 [INFO] LexicMap v0.4.0 (14c2606)
        16:22:49.745 [INFO]   https://github.com/shenwei356/LexicMap
        16:22:49.745 [INFO]
        16:22:49.745 [INFO] checking input files ...
        16:22:49.745 [INFO]   15 input file(s) given
        16:22:49.745 [INFO]
        16:22:49.745 [INFO] --------------------- [ main parameters ] ---------------------
        16:22:49.745 [INFO]
        16:22:49.745 [INFO] input and output:
        16:22:49.745 [INFO]   input directory: refs/
        16:22:49.745 [INFO]     regular expression of input files: (?i)\.(f[aq](st[aq])?|fna)(\.gz|\.xz|\.zst|\.bz2)?$
        16:22:49.745 [INFO]     *regular expression for extracting reference name from file name: (?i)(.+)\.(f[aq](st[aq])?|fna)(\.gz|\.xz|\.zst|\.bz2)?$
        16:22:49.745 [INFO]     *regular expressions for filtering out sequences: []
        16:22:49.745 [INFO]   max genome size: 15000000
        16:22:49.745 [INFO]   output directory: demo.lmi
        16:22:49.745 [INFO]
        16:22:49.745 [INFO] mask generation:
        16:22:49.745 [INFO]   k-mer size: 31
        16:22:49.745 [INFO]   number of masks: 40000
        16:22:49.745 [INFO]   rand seed: 1
        16:22:49.745 [INFO]   prefix length for checking low-complexity in mask generation: 15
        16:22:49.745 [INFO]
        16:22:49.745 [INFO] seed data:
        16:22:49.745 [INFO]   maximum sketching desert length: 450
        16:22:49.745 [INFO]   distance of k-mers to fill deserts: 150
        16:22:49.745 [INFO]   seeds data chunks: 16
        16:22:49.745 [INFO]   seeds data indexing partitions: 1024
        16:22:49.745 [INFO]
        16:22:49.745 [INFO] general:
        16:22:49.745 [INFO]   genome batch size: 5
        16:22:49.745 [INFO]   batch merge threads: 8
        16:22:49.745 [INFO]
        16:22:49.745 [INFO]
        16:22:49.745 [INFO] --------------------- [ generating masks ] ---------------------
        16:22:50.180 [INFO]
        16:22:50.180 [INFO] --------------------- [ building index ] ---------------------
        16:22:50.328 [INFO]
        16:22:50.328 [INFO]   ------------------------[ batch 1/3 ]------------------------
        16:22:50.328 [INFO]   building index for batch 1 with 5 files...
        processed files:  5 / 5 [======================================] ETA: 0s. done
        16:22:51.192 [INFO]   writing seeds...
        16:22:51.264 [INFO]   finished writing seeds in 71.756662ms
        16:22:51.264 [INFO]   finished building index for batch 1 in: 935.464336ms
        16:22:51.264 [INFO]
        16:22:51.264 [INFO]   ------------------------[ batch 2/3 ]------------------------
        16:22:51.264 [INFO]   building index for batch 2 with 5 files...
        processed files:  5 / 5 [======================================] ETA: 0s. done
        16:22:53.126 [INFO]   writing seeds...
        16:22:53.212 [INFO]   finished writing seeds in 86.823785ms
        16:22:53.212 [INFO]   finished building index for batch 2 in: 1.948770015s
        16:22:53.212 [INFO]
        16:22:53.212 [INFO]   ------------------------[ batch 3/3 ]------------------------
        16:22:53.212 [INFO]   building index for batch 3 with 5 files...
        processed files:  5 / 5 [======================================] ETA: 0s. done
        16:22:54.350 [INFO]   writing seeds...
        16:22:54.437 [INFO]   finished writing seeds in 87.058101ms
        16:22:54.437 [INFO]   finished building index for batch 3 in: 1.224414126s
        16:22:54.437 [INFO]
        16:22:54.437 [INFO] merging 3 indexes...
        16:22:54.437 [INFO]   [round 1]
        16:22:54.437 [INFO]     batch 1/1, merging 3 indexes to demo.lmi.tmp/r1_b1 with 8 threads...
        16:22:54.613 [INFO]   [round 1] finished in 175.640164ms
        16:22:54.613 [INFO] rename demo.lmi.tmp/r1_b1 to demo.lmi
        16:22:54.620 [INFO]
        16:22:54.620 [INFO] finished building LexicMap index from 15 files with 40000 masks in 4.875616203s
        16:22:54.620 [INFO] LexicMap index saved: demo.lmi
        16:22:54.620 [INFO]
        16:22:54.620 [INFO] elapsed time: 4.875654824s
        16:22:54.620 [INFO]

    {{< /expand >}}

### Notes for indexing with large datasets

If you have hundreds of thousands of input genomes or more, it's **better to control the number of genome batches**, which can be calculated via

    $num_input_files / --batch-size

E.g, for GenBank prokaryotic genomes: 2,340,672 / 5000 (default)  = 468.
The number is too big, and **it would slow down the seed-data merging step in `lexicmap index`** and **candidate sequence extraction in `lexicmap search`**.

Therefore, if you have enough memory, you can set a bigger `--batch-size` (e.g., 2,340,672 / 25000 = 93.6).

If the batch number is still big (e.g. 300), you can set bigger `--max-open-files` (e.g., `4096`) and `-J/--seed-data-threads` (e.g., `12`. 12 <= 4096/300 = 13.6)
to accelerate the merging step. Meanwhile, don't forget to increase the maximum open files per process via `ulimit -n 4096`.

If you forgot these setting, you can rerun the merging step for an unfinished index via [lexicmap utils remerge](https://bioinf.shenwei.me/LexicMap/usage/utils/remerge/)
(available since v0.5.0, also see [FAQ: how to resume the indexing](https://bioinf.shenwei.me/LexicMap/faqs/#how-to-resume-the-indexing-as-slurm-job-limit-is-almost-reached-while-lexicmap-index-is-still-in-the-merging-step)). Other cases to use this command:
- Only one thread is used for merging indexes, which happens when there are
a lot (>200 batches) of batches (`$inpu_files / --batch-size`) and the value
of `--max-open-files` is not big enough.
- The Slurm/PBS job time limit is almost reached and the merging step won't be finished before that.
- Disk quota is reached in the merging step.


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

LexicMap index size is related to the number of input genomes, the divergence between genome sequences, the number of masks, and the maximum seed distance.

**Note that the index size is not linear with the number of genomes, it's sublinear**. Because the seed data are compressed with VARINT-GB algorithm, more genome bring higher compression rates.

{{< tabs "t2" >}}

{{< tab "Demo data" >}}

    # 15 genomes
    demo.lmi: 73.30 MB (73,297,328)
      59.41 MB      seeds
      13.57 MB      genomes
     320.03 kB      masks.bin
         375 B      genomes.map.bin
         323 B      info.toml

{{< /tab>}}

{{< tab "GTDB repr" >}}

    # 85,205 genomes
    gtdb_repr.lmi: 228.15 GB (228,149,871,198)
     156.44 GB      seeds
      71.71 GB      genomes
       2.13 MB      genomes.map.bin
     320.03 kB      masks.bin
         329 B      info.toml

{{< /tab>}}

{{< tab "GTDB complete" >}}

    # 402,538 genomes
    gtdb_complete.lmi: 972.85 GB (972,854,821,322)
     583.10 GB      seeds
     389.74 GB      genomes
      10.06 MB      genomes.map.bin
     320.03 kB      masks.bin
         330 B      info.toml

{{< /tab>}}


{{< tab "Genbank+RefSeq" >}}

    # 2,340,672 genomes
    genbank_refseq.lmi: 5.43 TB (5,428,003,631,182)
       3.04 TB      seeds
       2.38 TB      genomes
      58.52 MB      genomes.map.bin
     320.03 kB      masks.bin
         332 B      info.toml

{{< /tab>}}


{{< tab "AllTheBacteria HQ" >}}

    # 1,858,610 genomes
    atb_hq.lmi: 4.26 TB (4,261,437,129,065)
       2.32 TB      seeds
       1.94 TB      genomes
      41.12 MB      genomes.map.bin
     320.03 kB      masks.bin
         332 B      info.toml


{{< /tab>}}

{{< /tabs >}}

- Directory/file sizes are counted with https://github.com/shenwei356/dirsize v1.2.1 (`dirsize -k`, **base: 1000**).
- Index building parameters: `-k 31 -m 40000`. Genome batch size: `-b 5000` for GTDB datasets, `-b 25000` for others.


## Explore the index

We provide several commands to explore the index data and extract indexed subsequences:

1. `lexicmap utils genomes` can list genome IDs of indexed genomes,
    see the [usage and example](https://bioinf.shenwei.me/LexicMap/usage/utils/genomes/).
1. `lexicmap utils masks` can list masks of the index,
    see the [usage and example](https://bioinf.shenwei.me/LexicMap/usage/utils/masks/).
1. `lexicmap utils kmers` can list details of all seeds (k-mers), including reference, location(s), the strand, and the k-mer direction.
    see the [usage and example](https://bioinf.shenwei.me/LexicMap/usage/utils/kmers/).
1. `lexicmap utils seed-pos` can help to explore the seed positions,
    see the [usage and example](https://bioinf.shenwei.me/LexicMap/usage/utils/seed-pos/).
    Before that, the flag `--save-seed-pos` needs to be added to `lexicmap index`.
1. `lexicmap utils subseq` can extract subsequences via genome ID, sequence ID and positions,
    see the [usage and example](https://bioinf.shenwei.me/LexicMap/usage/utils/subseq/).


**What's next:** {{< button size="small" relref="tutorials/search" >}}Searching{{< /button >}}
