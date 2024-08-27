---
title: Indexing GenBank+RefSeq
weight: 10
---

Tools:

- https://github.com/pirovc/genome_updater, for downloading genomes
- https://github.com/shenwei356/seqkit, for checking sequence files
- https://github.com/shenwei356/rush, for running jobs
- https://github.com/shenwei356/brename, for batch file renaming

Data:

    time genome_updater.sh -d "refseq,genbank" -g "archaea,bacteria" \
        -f "genomic.fna.gz" -o "genbank" -M "ncbi" -t 12 -m -L curl

    cd genbank/2024-02-15_11-00-51/

    # ----------------- just in case, check the file integrity -----------------
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

    # ----------------- just in case, check the file integrity -----------------


    # creat a new directory containing symbol links to the orginal files
    mkdir gtdb; cd gtdb;
    find ../files -name "*.fna.gz" | rush --eta 'ln -s {}'
    brename -p '^(\w{3}_\d{9}\.\d+).+' -r '$1.fna.gz'
    cd ..

    find gtdb/ -name "*.fna.gz" > files.txt

Indexing. On a 48-CPU machine, time: 54 h, ram: 178 GB, index size: 4.94 TB.
If you don't have enough memory, please decrease the value of `-b`.

    lexicmap index -b 25000 -S -X files.txt -O genbank_refseq.lmi --log genbank_refseq.lmi.log


