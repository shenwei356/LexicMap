---
title: Indexing GTDB
weight: 5
---

Info:

- https://gtdb.ecogenomic.org/

Tools:

- https://github.com/pirovc/genome_updater, for downloading genomes
- https://github.com/shenwei356/seqkit, for checking sequence files
- https://github.com/shenwei356/rush, for running jobs

Data:

    time genome_updater.sh -d "refseq,genbank" -g "archaea,bacteria" \
        -f "genomic.fna.gz" -o "GTDB_complete" -M "gtdb" -t 12 -m -L curl

    cd GTDB_complete/2024-01-30_19-34-40/


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

Indexing. On a 48-CPU machine, time: 8h:19m:28s, ram: 73 GB, index size: 906 GB.
If you don't have enough memory, please decrease the value of `-b`.

    lexicmap index \
        -I files/ \
        --ref-name-regexp '^(\w{3}_\d{9}\.\d+)' \
        -O gtdb_complete.lmi --log gtdb_complete.lmi.log \
        -b 5000

Files:

    $ du -sh files gtdb_complete.lmi --apparent-size
    413G    files
    907G    gtdb_complete.lmi

    $ dirsize gtdb_complete.lmi
    gtdb_complete.lmi: 905.32 GiB (972,076,048,287)
    542.32 GiB      seeds
    362.99 GiB      genomes
      9.60 MiB      genomes.map.bin
    156.28 KiB      masks.bin
         616 B      info.toml
         168 B      genomes.chunks.bin
