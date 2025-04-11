---
title: Indexing GenBank+RefSeq
weight: 10
---

**Make sure you have enough disk space, >10 TB is preferred.**

Tools:

- https://github.com/pirovc/genome_updater, for downloading genomes
- https://github.com/shenwei356/seqkit, for checking sequence files
- https://github.com/shenwei356/rush, for running jobs

Data:

    time genome_updater.sh -d "refseq,genbank" -g "archaea,bacteria" \
        -f "genomic.fna.gz" -o "genbank" -M "ncbi" -t 12 -m -L curl

    cd genbank/2024-02-15_11-00-51/


    # ----------------- check the file integrity -----------------

    genomes=files

    # corrupted files
    # find $genomes -name "*.gz" \
    fd ".gz$" $genomes \
        | rush --eta 'seqkit seq -w 0 {} > /dev/null; if [ $? -ne 0 ]; then echo {}; fi' \
        > failed.txt

    # empty files
    find $genomes -name "*.gz" -size 0 >> failed.txt

    # delete these files
    cat failed.txt | rush '/bin/rm {}'

    # redownload them:
    # run the genome_updater command again, with the flag -i

Indexing. On a 48-CPU machine, time: 56 h, ram: 181 GB, index size: 4.96 TiB.
If you don't have enough memory, please decrease the value of `-b`.

    lexicmap index \
        -I files/ \
        --ref-name-regexp '^(\w{3}_\d{9}\.\d+)' \
        -O genbank_refseq.lmi --log genbank_refseq.lmi.log \
        -b 25000

    # dirsize genbank_refseq.lmi
    genbank_refseq.lmi: 4.96 TiB (5,454,659,703,138)
      2.79 TiB      seeds
      2.17 TiB      genomes
     55.81 MiB      genomes.map.bin
    156.28 KiB      masks.bin
      3.59 KiB      genomes.chunks.bin
         619 B      info.toml
