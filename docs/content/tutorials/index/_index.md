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
LexicMap is mainly suitable for small genomes like Archaea, Bacteria, Viruses and plasmids.

Maximum genome size: 268 Mb (268,435,456).
More precisely:

    $total_bases + ($num_contigs - 1) * 1000 <= 268,435,456

as we insert 1000-bp intervals of N's between contigs to reduce the sequence scale to index.
{{< /hint >}}


**Sequences of each reference genome should be saved in separate FASTA/Q files, with identifiers in the file names**.Click to show

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
    - More RAM (> 100 GB) is preferred. The memory usage in index building is mainly related to:
        - The number of masks (`-m/--masks`, default 40,000).
        - The number of genomes.
        - The divergence between genome sequences. Diverse genomes consume more memory.
        - **The genome batch size**  (`-b/--batch-size`, default 5,000). **This is the main parameter to adjust memory usage**.
        - **The maximum seed distance** or **the maximum sketching desert size** (`-D/--seed-max-desert`, default 200),
          and the distance of k-mers to fill deserts (`-d/--seed-in-desert-dist`, default 50).
          Bigger `-D/--seed-max-desert` values decrease the search sensitivity for distant targets, speed up the indexing speed,
          decrease the indexing memory occupation and decrease the index size. While the alignment speed is almost not affected.
    - **If the RAM is not sufficient**. Please:
        1. **Use a smaller genome batch size**. It decreases indexing memory occupation and has little affection on searching performance.
        2. Use a smaller number of masks, e.g., 20,000 performs well for small genomes (<=5 Mb). And if the queries are long (>= 2kb), there's little affection for the alignment results.
- **Disk**
    - More (>2 TB) is better. LexicMap index size is related to the number of input genomes, the divergence between genome sequences, the number of masks, and the maximum seed distance. See [some examples](#index-size).
        - **Note that the index size is not linear with the number of genomes, it's sublinear**. Because the seed data are compressed with VARINT-GB algorithm, more genome bring higher compression rates.
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
        4. Capturing the most similar k-mer (in non-gap and non-interval regions) for each mask and recording the k-mer and its location(s) and strand information. Base N is treated as A.
        5. Filling sketching deserts (genome regions longer than `--seed-max-desert` without any captured k-mers/seeds).
           In a sketching desert, not a single k-mer is captured because there's another k-mer in another place which shares a longer prefix with the mask.
           As a result, for a query similar to seqs in this region, all captured k-mers can’t match the correct seeds.
            1. For a desert region (`start`, `end`), masking the extended region (`start-1000`, `end+1000`) with the masks.
            2. Starting from `start`, every around `--seed-in-desert-dist` (default 150) bp, finding a k-mer which is captured by some mask, and add the k-mer and its position information into the index of that mask.
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

{{< tab "Genome batches" >}}

|Flag                 |Value                     |Function                               |Comment                                                                                                                                                                                                                                                                                                                                                                        |
|:--------------------|:-------------------------|:--------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|**`-b/--batch-size`**|Max: 131072, default: 5000|Maximum number of genomes in each batch|If the number of input files exceeds this number, input files are split into multiple batches and indexes are built for all batches. In the end, seed files are merged, while genome data files are kept unchanged and collected. ■ Bigger values increase indexing memory occupation and increase batch searching speed, while single query searching speed is not affected.  |

{{< /tab>}}

{{< tab "LexicHash mask generation" >}}

|Flag                  |Value               |Function                                                  |Comment                                                                                                                                                                                   |
|:---------------------|:-------------------|:---------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|`-M/--mask-file`      |A file              |File with custom masks                                    |File with custom masks, which could be exported from an existing index or newly generated by "lexicmap utils masks". This flag oversides `-k/--kmer`, `-m/--masks`, `-s/--rand-seed`, etc.|
|**`-k/--kmer`**       |Max: 32, default: 31|K-mer size                                                |■ Bigger values improve the search specificity and do not increase the index size.                                                                                                        |
|**`-m/--masks`**      |Default: 40,000     |Number of masks                                           |■  Bigger values improve the search sensitivity, increase the index size, and slow down the search speed. For smaller genomes like phages/viruses, m=10,000 is high enough.               |
|`-p/--seed-min-prefix`|Max: 32, Default: 15|Minimum length of shared substrings (anchors) in searching|This value is used to remove masks with a prefix of low-complexity.                                                                                                                       |


{{< /tab>}}


{{< tab "Seeds (k-mer-value) data" >}}

|Flag                        |Value                           |Function                                                                        |Comment                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
|:---------------------------|:-------------------------------|:-------------------------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|**`--seed-max-desert`**     |Default: 200                    |Maximum length of distances between seeds                                       |The default value of 200 guarantees queries >200 bp would match at least one seed. ► Large regions with no seeds are called sketching deserts. Deserts with seed distance larger than this value will be filled by choosing k-mers roughly every --seed-in-desert-dist (50 by default) bases. ■ Bigger values decrease the search sensitivity for distant targets, speed up the indexing speed, decrease the indexing memory occupation and decrease the index size. While the alignment speed is almost not affected.    |
|`-c/--chunks`               |Maximum: 128, default: #CPUs    |Number of seed file chunks                                                      |Bigger values accelerate the search speed at the cost of a high disk reading load. The maximum number should not exceed the maximum number of open files set by the operating systems.                                                                                                                                                                                                                                                                                                                                    |
|**`-J/--seed-data-threads`**|Maximum: -c/--chunks, default: 8|Number of threads for writing seed data and merging seed chunks from all batches|■ Bigger values increase indexing speed at the cost of slightly higher memory occupation.                                                                                                                                                                                                                                                                                                                                                                                                                                 |
|`-p/--partitions`           |Default: 512                    |Number of partitions for indexing each seed file                                |Bigger values bring a little higher memory occupation. 512 is a good value with high searching speed, larger or smaller values would decrease the speed in `lexicmap search`. ► After indexing, `lexicmap utils reindex-seeds` can be used to reindex the seeds data with  another value of this flag.                                                                                                                                                                                                                    |
|`--max-open-files`          |Default: 512                    |Maximum number of open files                                                    |It's only used in merging indexes of multiple genome batches.                                                                                                                                                                                                                                                                                                                                                                                                                                                             |

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

    It would take about 3 seconds and 2 GB RAM in a 16-CPU PC.

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
        16:22:49.745 [INFO]   seeds data indexing partitions: 512
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
    demo.lmi/: 59.55 MB
      46.31 MB      seeds
      12.93 MB      genomes
     312.53 KB      masks.bin
      375.00 B      genomes.map.bin
      322.00 B      info.toml

{{< /tab>}}

{{< tab "GTDB repr" >}}

    # 85,205 genomes/
    gtdb_repr.lmi: 212.58 GB
     145.79 GB      seeds
      66.78 GB      genomes
       2.03 MB      genomes.map.bin
     312.53 KB      masks.bin
      328.00 B      info.toml

{{< /tab>}}

{{< tab "GTDB complete" >}}

    # 402,538 genomes
    gtdb_complete.lmi: 905.95 GB
     542.97 GB      seeds
     362.98 GB      genomes
       9.60 MB      genomes.map.bin
     312.53 KB      masks.bin
      329.00 B      info.toml

{{< /tab>}}


{{< tab "Genbank+RefSeq" >}}

    # 2,340,672 genomes
    genbank_refseq.lmi: 4.94 TB
       2.77 TB      seeds
       2.17 TB      genomes
      55.81 MB      genomes.map.bin
     312.53 KB      masks.bin
      331.00 B      info.toml

{{< /tab>}}


{{< tab "AllTheBacteria HQ" >}}

    # 1,858,610 genomes
    atb_hq.lmi: 3.88 TB
       2.11 TB      seeds
       1.77 TB      genomes
      39.22 MB      genomes.map.bin
     312.53 KB      masks.bin
      331.00 B      info.toml

{{< /tab>}}

{{< /tabs >}}

- Directory/file sizes are counted with https://github.com/shenwei356/dirsize.
- Index building parameters: `-k 31 -m 40000`. Genome batch size: `-b 5000` for GTDB datasets, `-b 25000` for others.


## Explore the index

1. `lexicmap utils genomes` can list genome IDs of indexed genomes,
    see the [usage and example](https://bioinf.shenwei.me/LexicMap/usage/utils/genomes/).
1. `lexicmap utils masks` can list masks of the index,
    see the [usage and example](https://bioinf.shenwei.me/LexicMap/usage/utils/masks/).
1. `lexicmap utils kmers` can list details of all seeds (k-mers), including reference, location(s) and the strand.
    see the [usage and example](https://bioinf.shenwei.me/LexicMap/usage/utils/kmers/).
1. `lexicmap utils seed-pos` can help to explore the seed positions,
    see the [usage and example](https://bioinf.shenwei.me/LexicMap/usage/utils/seed-pos/).
    Before that, the flag `--save-seed-pos` needs to be added to `lexicmap index`.
1. `lexicmap utils subseq` can extract subsequences via genome ID, sequence ID and positions,
    see the [usage and example](https://bioinf.shenwei.me/LexicMap/usage/utils/subseq/).


**What's next:** {{< button size="small" relref="tutorials/search" >}}Searching{{< /button >}}
