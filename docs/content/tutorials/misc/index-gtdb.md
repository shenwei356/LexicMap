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
- https://github.com/sharkdp/fd, faster `find`, available as `fd-find` in [conda-forge](https://anaconda.org/conda-forge/fd-find)

Data:

    time genome_updater.sh -d "refseq,genbank" -g "archaea,bacteria" \
        -f "genomic.fna.gz" -o "GTDB_complete" -M "gtdb" -t 12 -a -m -L curl

    cd GTDB_complete/2024-01-30_19-34-40/


    # ----------------- check the file integrity -----------------

    genomes=files

    # corrupted files
    # find $genomes -name "*.gz" \
    fd ".gz$" $genomes \
        | rush --eta 'seqkit seq -w 0 {} > /dev/null; if [ $? -ne 0 ]; then echo {}; exit 1; fi' \
            -c -C check-files.rush \
        > failed.txt

    # empty files
    # find $genomes -name "*.gz" -size 0 >> failed.txt
    fd --size 0b $genomes >> failed.txt

    # delete these files
    cat failed.txt | rush '/bin/rm {}'

    # redownload them:
    # run the genome_updater command again, with the flag -i
    
Taxonomy data to limit TaxId in `lexicmap search` since LexicMap v0.7.1.

    # Download gtdb-taxdump files here: https://github.com/shenwei356/gtdb-taxdump
    wget https://github.com/shenwei356/gtdb-taxdump/releases/download/v0.6.0/gtdb-taxdump-R226.tar.gz
    tar -zxvf gtdb-taxdump-R226.tar.gz
    mv gtdb-taxdump-R226 taxdump/
    
    # It has a file mapping assembly accession to TaxId
    ln -s taxdump/taxid.map

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
    gtdb_complete.lmi: 905.34 GiB (972,098,200,328)
    542.34 GiB      seeds
    362.99 GiB      genomes
      9.60 MiB      genomes.map.bin
    156.28 KiB      masks.bin
         616 B      info.toml
         168 B      genomes.chunks.bin
