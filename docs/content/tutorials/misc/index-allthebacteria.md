---
title: Indexing AllTheBacteria
weight: 15
---

Info:

- [AllTheBacteria](https://github.com/AllTheBacteria/AllTheBacteria), All WGS isolate bacterial INSDC data to June 2023uniformly assembled, QC-ed, annotated, searchable.
- Preprint: [AllTheBacteria - all bacterial genomes assembled, available and searchable](https://www.biorxiv.org/content/10.1101/2024.03.08.584059v1)

Steps:

1. Downloading assemblies tarballs here (except these starting with `unknown__`) to a directory (like assemblies):
https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/0.2/assembly/

1. Decompressing all tarballs.

        cd assemblies;
        ls *.tar.xz | parallel --eta 'tar -Jxf {}; gzip {}/*.fa'
        cd ..

    After that, the assemblies directory would have multiple subdirectories.
    When you give the directory to `lexicmap index -I`, it can recursively scan (plain or gz/xz/zstd-compressed) genome files.

2. Creating a LexicMap index. (more details: https://bioinf.shenwei.me/LexicMap/tutorials/index/)

       lexicmap index -I assemblies/ -O atb.lmi -b 25000 --log atb.lmi.log

   For 1,858,610 HQ genomes, on a 48-CPU machine, time: 48 h, ram: 85 GB, index size: 3.88 TB.
   If you don't have enough memory, please decrease the value of `-b`.

