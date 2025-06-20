---
title: Indexing GlobDB
weight: 20
---


Info:

- [GlobDB](https://globdb.org/) , a dereplicated dataset of the species reps of the GTDB, GEM, SPIRE and SMAG datasets a lot.
- https://x.com/daanspeth/status/1822964436950192218


Data:

    # check the latest version here: https://fileshare.lisc.univie.ac.at/globdb/
    
    # download data
    wget https://fileshare.lisc.univie.ac.at/globdb/globdb_r220/globdb_r220_genome_fasta.tar.gz

    tar -zxf globdb_r220_genome_fasta.tar.gz

    # file list
    find globdb_r220_genome_fasta/ -name "*.fa.gz" > files.txt
    
Taxonomy data to limit TaxId in `lexicmap search` since LexicMap v0.7.1.

    wget https://fileshare.lisc.univie.ac.at/globdb/globdb_r220/globdb_r220_tax.tsv
    
    # Create taxdump files with TaxonKit: https://bioinf.shenwei.me/taxonkit/usage/#create-taxdump    
    taxonkit create-taxdump -A 1 -R superkingdom,phylum,class,order,family,genus,species globdb_r220_tax.tsv -O taxdump
    
    # It has a file mapping assembly accession to TaxId
    ln -s taxdump/taxid.map
    
Indexing with LexicMap

    # elapsed time: 3h:40m:38s
    # peak rss: 87.15 GB
    lexicmap index -S -X files.txt -O globdb_r220.lmi --log globdb_r220.lmi -g 50000000

