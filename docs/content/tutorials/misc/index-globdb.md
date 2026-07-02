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
    wget -c https://fileshare.lisc.univie.ac.at/globdb/globdb_r232/globdb_r232_genome_fasta.tar.gz

    tar -zxf globdb_r232_genome_fasta.tar.gz

    # file list
    find globdb_r232_genome_fasta/ -name "*.fa.gz" > files.txt
    
    wc -l files.txt 
    # 346233 files.txt

Indexing with LexicMap

    # elapsed time: 9h:02m:28s
    # peak rss: 89.99 GB
    lexicmap index -S -X files.txt -O globdb_r232.lmi --log globdb_r232.lmi.log -g 50000000 -b 8000 -j 48
    
Taxonomy data for limiting TaxId in `lexicmap search` since LexicMap v0.8.0.

    wget https://fileshare.lisc.univie.ac.at/globdb/globdb_r232/globdb_r232_taxonomy.tsv.gz
    
    # Create taxdump files with TaxonKit: https://bioinf.shenwei.me/taxonkit/usage/#create-taxdump
    # 
    # for older versions with tab-delimited data:
    #    taxonkit create-taxdump -A 1 -R superkingdom,phylum,class,order,family,genus,species globdb_r220_tax.tsv -O taxdump
    #
    taxonkit create-taxdump --gtdb globdb_r232_taxonomy.tsv.gz -O taxdump
    
    # It has a file mapping assembly accession to TaxId (below the species rank)
    ln -s taxdump/taxid.map
    
How does the taxdump data look like?

    $ echo Escherichia coli | taxonkit name2taxid --data-dir taxdump/
    Escherichia coli        599451526
    
    $ echo 599451526 | taxonkit lineage --data-dir taxdump/ -r
    599451526       Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli    species
    
    $ echo 599451526 | taxonkit list --data-dir taxdump/ -nr 
    599451526 [species] Escherichia coli
    1395899945 [no rank] GCF_003697165
    
    $ grep 1395899945 taxdump/taxid.map 
    GCF_003697165   1395899945

    # the number of species
    $ echo 1 | taxonkit list --data-dir taxdump/ -I "" | taxonkit filter -E species --data-dir taxdump/ | wc -l
    346233
    
Search in specific taxa.

    # a species (Escherichia coli)
    lexicmap search -d globdb_r232.lmi -T taxdump/ -G taxdump/taxid.map -t 599451526 b.gene_E_coli_16S.fasta --debug
    
    # a genus (Escherichia)
    # echo Escherichia | taxonkit name2taxid --data-dir taxdump/
    # Escherichia     1028471294
    lexicmap search -d globdb_r232.lmi -T taxdump/ -G taxdump/taxid.map -t 1028471294 b.gene_E_coli_16S.fasta --debug

(last update: 2026-07-02)
