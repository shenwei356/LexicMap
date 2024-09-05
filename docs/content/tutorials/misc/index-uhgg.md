---
title: Indexing UHGG
weight: 25
---

Info:

- [Unified Human Gastrointestinal Genome (UHGG) v2.0.2](https://www.ebi.ac.uk/metagenomics/genome-catalogues/human-gut-v2-0-2)
- [A unified catalog of 204,938 reference genomes from the human gut microbiome](https://www.nature.com/articles/s41587-020-0603-3)
- Number of Genomes: 289,232

Tools:

- https://github.com/shenwei356/seqkit, for checking sequence files
- https://github.com/shenwei356/rush, for running jobs

Data:

    # meta data
    wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/genomes-all_metadata.tsv

    # gff url
    sed 1d genomes-all_metadata.tsv | cut -f 20  | sed 's/v2.0/v2.0.2/' | sed -E 's/^ftp/https/' > url.txt

    # download gff files
    mkdir -p files; cd files

    time cat ../url.txt \
        | rush --eta -v 'dir={///%}/{//%}' \
            'mkdir -p {dir}; curl -s -o {dir}/{%} {}' \
            -c -C download.rush -j 12
    cd ..

    # extract sequences from gff files
    find files/ -name "*.gff.gz" \
        | rush --eta \
            'zcat {} | perl -ne "print if \$s; \$s=true if /^##FASTA/" | seqkit seq -w 0 -o {/}/{%:}.fna.gz' \
            -c -C extract.rush


Indexing. On a 48-CPU machine, time: 3 h, ram: 41 GB, index size: 426 GB.
If you don't have enough memory, please decrease the value of `-b`.

    lexicmap index \
        -I files/ \
        -O uhgg.lmi --log uhgg.lmi.log \
        -b 5000

File sizes:

    $ du -sh files/ uhgg.lmi
    658G    files/
    509G    uhgg.lmi

    $ du -sh files/ uhgg.lmi --apparent-size
    425G    files/
    426G    uhgg.lmi

    $ dirsize uhgg.lmi
    uhgg.lmi: 425.15 GiB (456,497,171,291)
    243.47 GiB      seeds
    181.67 GiB      genomes
      6.34 MiB      genomes.map.bin
    312.53 KiB      masks.bin
         330 B      info.toml
