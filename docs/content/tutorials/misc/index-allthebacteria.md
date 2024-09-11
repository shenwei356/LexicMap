---
title: Indexing AllTheBacteria
weight: 15
---

**Make sure you have enough disk space, at least 8 TB, >10 TB is preferred.**

Tools:

- https://github.com/shenwei356/rush, for running jobs

Info:

- [AllTheBacteria](https://github.com/AllTheBacteria/AllTheBacteria), All WGS isolate bacterial INSDC data to June 2023 uniformly assembled, QC-ed, annotated, searchable.
- Preprint: [AllTheBacteria - all bacterial genomes assembled, available and searchable](https://www.biorxiv.org/content/10.1101/2024.03.08.584059v1)

## Steps for v0.2

1. Downloading assemblies tarballs here (except these starting with `unknown__`) to a directory (like `atb`):
https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/0.2/assembly/

        mkdir -p atb;
        cd atb;

        # assembly file list, 650 files in total
        wget https://bioinf.shenwei.me/LexicMap/AllTheBacteria-v0.2.url.txt

        # download
        #   rush is used: https://github.com/shenwei356/rush
        #   The download.rush file stores finished jobs, which will be skipped in a second run for resuming jobs.
        cat AllTheBacteria-v0.2.url.txt | rush --eta -j 2 -c -C download.rush 'wget {}'


        # list of high-quality samples
        wget https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/0.2/metadata/hq_set.sample_list.txt.gz

1. Decompressing all tarballs. The decompressed genomes are stored in plain text,
   so we use `gzip` (can be replaced with faster `pigz` ) to compress them to save disk space.

        # {^asm.tar.xz} is for removing the suffix "asm.tar.xz"
        ls *.tar.xz | rush --eta -c -C decompress.rush 'tar -Jxf {}; gzip -f {^asm.tar.xz}/*.fa'
        cd ..

    After that, the assemblies directory would have multiple subdirectories.
    When you give the directory to `lexicmap index -I`, it can recursively scan (plain or gz/xz/zstd-compressed) genome files.
    You can also give a file list with selected assemblies.

        $ tree atb | more
        atb
        ├── achromobacter_xylosoxidans__01
        │   ├── SAMD00013333.fa.gz
        │   ├── SAMD00049594.fa.gz
        │   ├── SAMD00195911.fa.gz
        │   ├── SAMD00195914.fa.gz


        # disk usage

        $ du -sh atb
        2.9T    atb

        $ du -sh atb --apparent-size
        2.1T    atb

2. Creating a LexicMap index. (more details: https://bioinf.shenwei.me/LexicMap/tutorials/index/)

        # file paths of all samples
        find atb/ -name "*.fa.gz" > atb_all.txt

        # wc -l atb_all.txt
        # 1876015 atb_all.txt

        # file paths of high-quality samples
        grep -w -f <(zcat atb/hq_set.sample_list.txt.gz) atb_all.txt > atb_hq.txt

        # wc -l atb_hq.txt
        # 1858610 atb_hq.txt



       # index
       lexicmap index -S -X atb_hq.txt -O atb_hq.lmi -b 25000 --log atb_hq.lmi.log

   For 1,858,610 HQ genomes, on a 48-CPU machine, time: 48 h, ram: 85 GB, index size: 3.88 TB.
   If you don't have enough memory, please decrease the value of `-b`.

        # disk usage

        $ du -sh atb_hq.lmi
        4.6T    atb_hq.lmi

        $ du -sh atb_hq.lmi --apparent-size
        3.9T    atb_hq.lmi

        $ dirsize atb_hq.lmi

        atb_hq.lmi: 3.88 TiB (4,261,437,129,065)
          2.11 TiB      seeds
          1.77 TiB      genomes
         39.22 MiB      genomes.map.bin
        312.53 KiB      masks.bin
             332 B      info.toml

    Note that, there's a tmp directory `atb_hq.lmi` being created during indexing.
    In the tmp directory, the seed data would be bigger than the final size of `seeds` directory,
    however, the genome files are simply moved to the final index.
