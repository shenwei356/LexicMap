## Demo datasets

### Reference genomes

We choose 15 bacterial genomes for demonstration.

Taxonomy information (NCBI Taxonomy):

    cat taxid.map \
        | taxonkit reformat -I 2 -f '{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}' \
        | csvtk cut -t -f -2 \
        | csvtk add-header -t -n id,superkingdom,phylum,class,order,family,genus,species \
        > taxonomy.tsv

    csvtk pretty -t taxonomy.tsv

    id                superkingdom   phylum           class                 order              family               genus            species
    ---------------   ------------   --------------   -------------------   ----------------   ------------------   --------------   --------------------------
    GCF_000742135.1   Bacteria       Pseudomonadota   Gammaproteobacteria   Enterobacterales   Enterobacteriaceae   Klebsiella       Klebsiella pneumoniae
    GCF_003697165.2   Bacteria       Pseudomonadota   Gammaproteobacteria   Enterobacterales   Enterobacteriaceae   Escherichia      Escherichia coli
    GCF_002949675.1   Bacteria       Pseudomonadota   Gammaproteobacteria   Enterobacterales   Enterobacteriaceae   Shigella         Shigella dysenteriae
    GCF_002950215.1   Bacteria       Pseudomonadota   Gammaproteobacteria   Enterobacterales   Enterobacteriaceae   Shigella         Shigella flexneri
    GCF_000006945.2   Bacteria       Pseudomonadota   Gammaproteobacteria   Enterobacterales   Enterobacteriaceae   Salmonella       Salmonella enterica
    GCF_001544255.1   Bacteria       Bacillota        Bacilli               Lactobacillales    Enterococcaceae      Enterococcus     Enterococcus faecium
    GCF_000392875.1   Bacteria       Bacillota        Bacilli               Lactobacillales    Enterococcaceae      Enterococcus     Enterococcus faecalis
    GCF_001457655.1   Bacteria       Pseudomonadota   Gammaproteobacteria   Pasteurellales     Pasteurellaceae      Haemophilus      Haemophilus influenzae
    GCF_900638025.1   Bacteria       Pseudomonadota   Gammaproteobacteria   Pasteurellales     Pasteurellaceae      Haemophilus      Haemophilus parainfluenzae
    GCF_001027105.1   Bacteria       Bacillota        Bacilli               Bacillales         Staphylococcaceae    Staphylococcus   Staphylococcus aureus
    GCF_006742205.1   Bacteria       Bacillota        Bacilli               Bacillales         Staphylococcaceae    Staphylococcus   Staphylococcus epidermidis
    GCF_001096185.1   Bacteria       Bacillota        Bacilli               Lactobacillales    Streptococcaceae     Streptococcus    Streptococcus pneumoniae
    GCF_000148585.2   Bacteria       Bacillota        Bacilli               Lactobacillales    Streptococcaceae     Streptococcus    Streptococcus mitis
    GCF_009759685.1   Bacteria       Pseudomonadota   Gammaproteobacteria   Moraxellales       Moraxellaceae        Acinetobacter    Acinetobacter baumannii
    GCF_000017205.1   Bacteria       Pseudomonadota   Gammaproteobacteria   Pseudomonadales    Pseudomonadaceae     Pseudomonas      Pseudomonas aeruginosa

Create a file for mapping assembly accessions to species names.

    cat taxid.map \
        | taxonkit reformat -I 2 -f '{s}' \
        | csvtk cut -t -f 1,3 \
        > ass2species.map

### Queries

- A gene sequence: [16S rRNA gene from *Escherichia coli* str. K-12 substr. MG1655](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3).
- A prophage sequence: [Enterobacteria phage P2](https://www.ncbi.nlm.nih.gov/nuccore/NC_001895.1).
- Simulated Oxford Nanopore R10.4.1 long-reads: simulated with [Badread](https://github.com/rrwick/Badread) from the 15 genomes.

        # simulate
        ls refs/*.gz | rush --eta 'badread simulate --reference {} --quantity 1x | seqkit replace -p ".+" -r "{%..}_r{nr}" > {}.fastq'

        # concatenate and remove quality scores to save space
        seqkit fq2fa refs/*.fastq | seqkit shuffle -o q.long-reads.fasta.gz

        # clean
        rm refs/*.fastq

Overview

    $ seqkit stats q.*.fasta q.*.fasta.gz --quiet
    file                   format  type  num_seqs     sum_len  min_len  avg_len  max_len
    q.gene.fasta           FASTA   DNA          1       1,542    1,542    1,542    1,542
    q.prophage.fasta       FASTA   DNA          1      33,593   33,593   33,593   33,593
    q.long-reads.fasta.gz  FASTA   DNA      3,692  54,375,807       67   14,728   90,376


## Building an index

Here we just use the top 5 biggest genomes for computating masks, for saving memory.

    $ lexicmap index -I refs/ -O demo.lmi -n 5
    08:36:24.378 [INFO] removing old output directory: demo.lmi
    08:36:24.381 [INFO] LexicMap v0.3.0
    08:36:24.381 [INFO]   https://github.com/shenwei356/LexicMap
    08:36:24.381 [INFO]
    08:36:24.381 [INFO] checking input files ...
    08:36:24.382 [INFO]   15 input file(s) given
    08:36:24.382 [INFO]
    08:36:24.382 [INFO] --------------------- [ main parameters ] ---------------------
    08:36:24.382 [INFO]
    08:36:24.382 [INFO] input and output:
    08:36:24.382 [INFO]   input directory: refs/
    08:36:24.382 [INFO]     regular expression of input files: (?i)\.(f[aq](st[aq])?|fna)(.gz)?$
    08:36:24.382 [INFO]     *regular expression for extracting reference name from file name: (?i)(.+)\.(f[aq](st[aq])?|fna)(.gz)?$
    08:36:24.382 [INFO]     *regular expressions for filtering out sequences: []
    08:36:24.382 [INFO]   max genome size: 15000000
    08:36:24.382 [INFO]   output directory: demo.lmi
    08:36:24.382 [INFO]
    08:36:24.382 [INFO] k-mer size: 31
    08:36:24.382 [INFO] number of masks: 40000
    08:36:24.382 [INFO] rand seed: 1
    08:36:24.382 [INFO] top N genomes for generating mask: 5
    08:36:24.382 [INFO] prefix extension length: 8
    08:36:24.382 [INFO]
    08:36:24.382 [INFO] seeds data chunks: 16
    08:36:24.382 [INFO] seeds data indexing partitions: 512
    08:36:24.382 [INFO]
    08:36:24.382 [INFO] genome batch size: 10000
    08:36:24.382 [INFO]
    08:36:24.382 [INFO]
    08:36:24.382 [INFO] --------------------- [ generating masks ] ---------------------
    08:36:24.382 [INFO]
    08:36:24.382 [INFO] generating masks from the top 5 out of 15 genomes...
    08:36:24.382 [INFO]
    08:36:24.382 [INFO]   checking genomes sizes of 15 files...
    processed files:  15 / 15 [======================================] ETA: 0s. done
    08:36:24.431 [INFO]     0 genomes longer than 15000000 are filtered out
    08:36:24.431 [INFO]     genome size range in the top 5 files: [4938295, 6588339]
    08:36:24.431 [INFO]
    08:36:24.431 [INFO]   collecting k-mers from 5 files...
    processed files: 5/5
    08:36:31.933 [INFO]
    08:36:31.933 [INFO]   generating masks...
    08:36:31.937 [INFO]     generating 16384 masks covering all 7-bp prefixes...
    processed prefixes: 16384/16384
    08:36:57.657 [INFO]     generating left 23616 masks...
    processed prefixes: 23616/23616
    08:37:29.705 [INFO]
    08:37:29.705 [INFO]   maximum distance between seeds:
    08:37:29.706 [INFO]     GCF_002950215.1: 1249
    08:37:29.707 [INFO]     GCF_000006945.2: 1415
    08:37:29.708 [INFO]     GCF_003697165.2: 1173
    08:37:29.710 [INFO]     GCF_000742135.1: 1243
    08:37:29.711 [INFO]     GCF_000017205.1: 1003
    08:37:29.722 [INFO]
    08:37:29.722 [INFO]   finished generating masks in: 1m5.339634848s
    08:37:29.723 [INFO]
    08:37:29.723 [INFO] --------------------- [ building index ] ---------------------
    08:37:29.970 [INFO]
    08:37:29.970 [INFO]   ------------------------[ batch 0 ]------------------------
    08:37:29.970 [INFO]   building index for batch 0 with 15 files...
    processed files:  15 / 15 [======================================] ETA: 0s. done
    08:37:32.146 [INFO]   writing seeds...
    08:37:32.349 [INFO]   finished writing seeds in 202.743652ms
    08:37:32.349 [INFO]   finished building index for batch 0 in: 2.378872619s
    08:37:32.349 [INFO]
    08:37:32.349 [INFO] finished building LexicMap index from 15 files with 40000 masks in 1m7.971157983s
    08:37:32.349 [INFO] LexicMap index saved: demo.lmi
    08:37:32.349 [INFO]
    08:37:32.349 [INFO] elapsed time: 1m7.971187597s
    08:37:32.349 [INFO]

Overview of index files:

    $ tree demo.lmi/
    demo.lmi/
    ├── genomes
    │   └── batch_0000
    │       ├── genomes.bin
    │       └── genomes.bin.idx
    ├── genomes.map.bin
    ├── info.toml
    ├── masks.bin
    └── seeds
        ├── chunk_000.bin
        ├── chunk_000.bin.idx
        ├── chunk_001.bin
        ├── chunk_001.bin.idx
        ...


    $ dirsize  demo.lmi/
    demo.lmi/: 26.63 MB
      13.39 MB      seeds
      12.93 MB      genomes
     312.53 KB      masks.bin
      375.00 B      genomes.map.bin
      261.00 B      info.toml

## Searching

### A 16S rRNA gene sequence

    $ lexicmap search -d demo.lmi/  q.gene.fasta -o q.gene.fasta.lexicmap.tsv
    09:32:55.551 [INFO] LexicProf v0.3.0
    09:32:55.551 [INFO]   https://github.com/shenwei356/LexicMap
    09:32:55.551 [INFO]
    09:32:55.551 [INFO] checking input files ...
    09:32:55.551 [INFO]   1 input file(s) given
    09:32:55.551 [INFO]
    09:32:55.551 [INFO] loading index: demo.lmi/
    09:32:55.551 [INFO]   reading masks...
    09:32:55.552 [INFO]   reading indexes of seeds (k-mer-value) data...
    09:32:55.555 [INFO]   creating genome reader pools, each batch with 16 readers...
    09:32:55.555 [INFO] index loaded in 4.192051ms
    09:32:55.555 [INFO]
    09:32:55.555 [INFO] searching ...

    09:32:55.596 [INFO]
    09:32:55.596 [INFO] processed queries: 1, speed: 1467.452 queries per minute
    09:32:55.596 [INFO] 100.0000% (1/1) queries matched
    09:32:55.596 [INFO] done searching
    09:32:55.596 [INFO] search results saved to: q.gene.fasta.lexicmap.tsv
    09:32:55.596 [INFO]
    09:32:55.596 [INFO] elapsed time: 45.230604ms
    09:32:55.596 [INFO]

Result preview.
Here we create a `species` column from the genome ID column (`sgnm`) and replace the assemby accessions with species names.

    $ csvtk head -n 28 q.gene.fasta.lexicmap.tsv \
        | csvtk mutate -t -n species -f sgnm \
        | csvtk replace -t -f species -k ass2species.map -p '(.+)' -r '{kv}' \
        | csvtk pretty -t

    query                         qlen   qstart   qend   sgnms   sgnm              seqid           qcovGnm   hsp   qcovHSP   alenHSP   alenFrag   pident   slen      sstart    send      sstr   seeds   species
    ---------------------------   ----   ------   ----   -----   ---------------   -------------   -------   ---   -------   -------   --------   ------   -------   -------   -------   ----   -----   --------------------
    NC_000913.3:4166659-4168200   1542   1        1542   12      GCF_003697165.2   NZ_CP033092.2   100.000   1     99.287    1542      1542       99.287   4903501   3780640   3782181   -      25      Escherichia coli
    NC_000913.3:4166659-4168200   1542   1        1542   12      GCF_003697165.2   NZ_CP033092.2   100.000   2     99.287    1542      1542       99.287   4903501   4551515   4553056   -      25      Escherichia coli
    NC_000913.3:4166659-4168200   1542   1        1542   12      GCF_003697165.2   NZ_CP033092.2   100.000   3     99.287    1542      1542       99.287   4903501   4591684   4593225   -      25      Escherichia coli
    NC_000913.3:4166659-4168200   1542   1        1542   12      GCF_003697165.2   NZ_CP033092.2   100.000   4     99.287    1542      1542       99.287   4903501   458559    460100    +      25      Escherichia coli
    NC_000913.3:4166659-4168200   1542   1        1542   12      GCF_003697165.2   NZ_CP033092.2   100.000   5     99.287    1542      1542       99.287   4903501   4844587   4846128   -      25      Escherichia coli
    NC_000913.3:4166659-4168200   1542   1        1542   12      GCF_003697165.2   NZ_CP033092.2   100.000   6     99.287    1542      1542       99.287   4903501   1285123   1286664   +      25      Escherichia coli
    NC_000913.3:4166659-4168200   1542   1        1542   12      GCF_003697165.2   NZ_CP033092.2   100.000   7     99.092    1542      1542       99.092   4903501   4726193   4727734   -      25      Escherichia coli
    NC_000913.3:4166659-4168200   1542   1        1542   12      GCF_002950215.1   NZ_CP026788.1   100.000   1     99.027    1542      1542       99.027   4659463   3216505   3218046   +      19      Shigella flexneri
    NC_000913.3:4166659-4168200   1542   1        1542   12      GCF_002950215.1   NZ_CP026788.1   100.000   2     98.962    1542      1542       98.962   4659463   3396068   3397609   +      23      Shigella flexneri
    NC_000913.3:4166659-4168200   1542   1        1542   12      GCF_002950215.1   NZ_CP026788.1   100.000   3     98.962    1542      1542       98.962   4659463   3119331   3120872   +      21      Shigella flexneri
    NC_000913.3:4166659-4168200   1542   1        1542   12      GCF_002950215.1   NZ_CP026788.1   100.000   4     98.898    1542      1542       98.898   4659463   4223146   4224687   +      23      Shigella flexneri
    NC_000913.3:4166659-4168200   1542   1        1542   12      GCF_002950215.1   NZ_CP026788.1   100.000   5     98.898    1542      1542       98.898   4659463   3355632   3357173   +      22      Shigella flexneri
    NC_000913.3:4166659-4168200   1542   1        1542   12      GCF_002950215.1   NZ_CP026788.1   100.000   6     98.833    1542      1542       98.833   4659463   3540450   3541991   +      21      Shigella flexneri
    NC_000913.3:4166659-4168200   1542   1        1542   12      GCF_002950215.1   NZ_CP026788.1   100.000   7     98.768    1541      1541       98.832   4659463   2125377   2126917   -      18      Shigella flexneri
    NC_000913.3:4166659-4168200   1542   1        1001   12      GCF_002949675.1   NZ_CP026774.1   100.000   1     96.563    1505      1001       98.937   4395762   2810845   2811845   +      25      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   1039     1542   12      GCF_002949675.1   NZ_CP026774.1   100.000   1     96.563    1505      504        98.937   4395762   2811883   2812386   +      25      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   1        1001   12      GCF_002949675.1   NZ_CP026774.1   100.000   2     96.563    1505      1001       98.937   4395762   2768883   2769883   +      25      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   1039     1542   12      GCF_002949675.1   NZ_CP026774.1   100.000   2     96.563    1505      504        98.937   4395762   2769921   2770424   +      25      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   1        1001   12      GCF_002949675.1   NZ_CP026774.1   100.000   3     96.563    1505      1001       98.937   4395762   2636477   2637477   +      25      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   1039     1542   12      GCF_002949675.1   NZ_CP026774.1   100.000   3     96.563    1505      504        98.937   4395762   2637515   2638018   +      25      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   1        1001   12      GCF_002949675.1   NZ_CP026774.1   100.000   4     96.563    1505      1001       98.937   4395762   3646778   3647778   +      25      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   1039     1542   12      GCF_002949675.1   NZ_CP026774.1   100.000   4     96.563    1505      504        98.937   4395762   3647816   3648319   +      25      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   1        1001   12      GCF_002949675.1   NZ_CP026774.1   100.000   5     96.563    1505      1001       98.937   4395762   3061592   3062592   +      25      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   1039     1542   12      GCF_002949675.1   NZ_CP026774.1   100.000   5     96.563    1505      504        98.937   4395762   3062630   3063133   +      25      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   1        1001   12      GCF_002949675.1   NZ_CP026774.1   100.000   6     96.563    1505      1001       98.937   4395762   2536624   2537624   +      25      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   1039     1542   12      GCF_002949675.1   NZ_CP026774.1   100.000   6     96.563    1505      504        98.937   4395762   2537662   2538165   +      25      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   1        1001   12      GCF_002949675.1   NZ_CP026774.1   100.000   7     96.563    1505      1001       98.937   4395762   1662551   1663551   -      25      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   1039     1542   12      GCF_002949675.1   NZ_CP026774.1   100.000   7     96.563    1505      504        98.937   4395762   1662010   1662513   -      25      Shigella dysenteriae

### A prophage sequence

    $ lexicmap search -d demo.lmi/ q.prophage.fasta -o q.prophage.fasta.lexicmap.tsv

    $ csvtk head -n 20 q.prophage.fasta.lexicmap.tsv \
        | csvtk mutate -t -n species -f sgnm \
        | csvtk replace -t -f species -k ass2species.map -p '(.+)' -r '{kv}' \
        | csvtk pretty -t
    query         qlen    qstart   qend    sgnms   sgnm              seqid           qcovGnm   hsp   qcovHSP   alenHSP   alenFrag   pident   slen      sstart    send      sstr   seeds   species
    -----------   -----   ------   -----   -----   ---------------   -------------   -------   ---   -------   -------   --------   ------   -------   -------   -------   ----   -----   ----------------
    NC_001895.1   33593   3        9357    1       GCF_003697165.2   NZ_CP033092.2   73.822    1     52.422    18922     9357       93.066   4903501   1864413   1873769   +      98      Escherichia coli
    NC_001895.1   33593   17627    22445   1       GCF_003697165.2   NZ_CP033092.2   73.822    1     52.422    18922     4819       93.066   4903501   1882193   1887011   +      98      Escherichia coli
    NC_001895.1   33593   22468    24382   1       GCF_003697165.2   NZ_CP033092.2   73.822    1     52.422    18922     1915       93.066   4903501   1887034   1888948   +      98      Escherichia coli
    NC_001895.1   33593   10592    12485   1       GCF_003697165.2   NZ_CP033092.2   73.822    1     52.422    18922     1894       93.066   4903501   1874130   1876023   +      98      Escherichia coli
    NC_001895.1   33593   15218    15283   1       GCF_003697165.2   NZ_CP033092.2   73.822    1     52.422    18922     66         93.066   4903501   1879476   1879541   +      98      Escherichia coli
    NC_001895.1   33593   17473    17552   1       GCF_003697165.2   NZ_CP033092.2   73.822    1     52.422    18922     80         93.066   4903501   1882043   1882122   +      98      Escherichia coli
    NC_001895.1   33593   14696    15072   1       GCF_003697165.2   NZ_CP033092.2   73.822    1     52.422    18922     377        93.066   4903501   1878954   1879330   +      98      Escherichia coli
    NC_001895.1   33593   12514    12664   1       GCF_003697165.2   NZ_CP033092.2   73.822    1     52.422    18922     151        93.066   4903501   1876052   1876202   +      98      Escherichia coli
    NC_001895.1   33593   10308    10570   1       GCF_003697165.2   NZ_CP033092.2   73.822    1     52.422    18922     263        93.066   4903501   1873846   1874108   +      98      Escherichia coli
    NC_001895.1   33593   24355    27239   1       GCF_003697165.2   NZ_CP033092.2   73.822    2     16.813    5877      2885       96.103   4903501   1853098   1855982   +      40      Escherichia coli
    NC_001895.1   33593   28430    30281   1       GCF_003697165.2   NZ_CP033092.2   73.822    2     16.813    5877      1852       96.103   4903501   1857173   1859024   +      40      Escherichia coli
    NC_001895.1   33593   27262    28401   1       GCF_003697165.2   NZ_CP033092.2   73.822    2     16.813    5877      1140       96.103   4903501   1856005   1857144   +      40      Escherichia coli


### Simulated Oxford Nanopore R10.4.1 long-reads

Here we use the flag `-w/--load-whole-seeds` to accelerate searching.

    $ lexicmap search -d demo.lmi/ q.long-reads.fasta.gz -o q.long-reads.fasta.gz.lexicmap.tsv.gz -w
    10:02:07.729 [INFO] LexicProf v0.3.0
    10:02:07.729 [INFO]   https://github.com/shenwei356/LexicMap
    10:02:07.729 [INFO]
    10:02:07.729 [INFO] checking input files ...
    10:02:07.729 [INFO]   1 input file(s) given
    10:02:07.729 [INFO]
    10:02:07.729 [INFO] loading index: demo.lmi/
    10:02:07.729 [INFO]   reading masks...
    10:02:07.730 [INFO]   reading seeds (k-mer-value) data into memory...
    10:02:07.740 [INFO]   creating genome reader pools, each batch with 16 readers...
    10:02:07.741 [INFO] index loaded in 11.375891ms
    10:02:07.741 [INFO]
    10:02:07.741 [INFO] searching ...
    processed queries: 3584, speed: 6627.022 queries per minute
    10:02:41.167 [INFO]
    10:02:41.167 [INFO] processed queries: 3692, speed: 6627.114 queries per minute
    10:02:41.167 [INFO] 95.8559% (3539/3692) queries matched
    10:02:41.167 [INFO] done searching
    10:02:41.167 [INFO] search results saved to: q.long-reads.fasta.gz.lexicmap.tsv.gz
    10:02:41.171 [INFO]
    10:02:41.171 [INFO] elapsed time: 33.441466803s
    10:02:41.171 [INFO]

Result overview:

    csvtk head -n 26 q.long-reads.fasta.gz.lexicmap.tsv.gz \
        | csvtk mutate -t -n species -f sgnm \
        | csvtk replace -t -f species -k ass2species.map -p '(.+)' -r '{kv}' \
        | csvtk pretty -t

    query                  qlen   qstart   qend   sgnms   sgnm              seqid           qcovGnm   hsp   qcovHSP   alenHSP   alenFrag   pident   slen      sstart    send      sstr   seeds   species
    --------------------   ----   ------   ----   -----   ---------------   -------------   -------   ---   -------   -------   --------   ------   -------   -------   -------   ----   -----   --------------------------
    GCF_003697165.2_r46    2169   446      1666   1       GCF_003697165.2   NZ_CP033092.2   93.730    1     73.721    2033      1257       78.652   4903501   4490226   4491482   +      5       Escherichia coli
    GCF_003697165.2_r46    2169   1794     2023   1       GCF_003697165.2   NZ_CP033092.2   93.730    1     73.721    2033      224        78.652   4903501   4491612   4491835   +      5       Escherichia coli
    GCF_003697165.2_r46    2169   2076     2146   1       GCF_003697165.2   NZ_CP033092.2   93.730    1     73.721    2033      72         78.652   4903501   4491847   4491918   +      5       Escherichia coli
    GCF_003697165.2_r46    2169   1691     1766   1       GCF_003697165.2   NZ_CP033092.2   93.730    1     73.721    2033      76         78.652   4903501   4491509   4491584   +      5       Escherichia coli
    GCF_003697165.2_r46    2169   31       420    1       GCF_003697165.2   NZ_CP033092.2   93.730    1     73.721    2033      404        78.652   4903501   4489794   4490197   +      5       Escherichia coli
    GCF_009759685.1_r164   3132   20       1929   1       GCF_009759685.1   NZ_CP046654.1   98.595    1     91.507    3088      1927       92.811   3980848   1768740   1770666   +      14      Acinetobacter baumannii
    GCF_009759685.1_r164   3132   1969     3121   1       GCF_009759685.1   NZ_CP046654.1   98.595    1     91.507    3088      1161       92.811   3980848   1770695   1771855   +      14      Acinetobacter baumannii
    GCF_000006945.2_r109   3788   37       3278   1       GCF_000006945.2   NC_003197.2     97.888    1     94.483    3708      3243       96.521   4857450   4633813   4637055   -      29      Salmonella enterica
    GCF_000006945.2_r109   3788   3301     3768   1       GCF_000006945.2   NC_003197.2     97.888    1     94.483    3708      465        96.521   4857450   4633323   4633787   -      29      Salmonella enterica
    GCF_900638025.1_r28    6375   3174     5997   1       GCF_900638025.1   NZ_LR134481.1   97.820    1     84.502    6236      2858       86.386   2062405   137873    140730    -      47      Haemophilus parainfluenzae
    GCF_900638025.1_r28    6375   6032     6339   1       GCF_900638025.1   NZ_LR134481.1   97.820    1     84.502    6236      311        86.386   2062405   137524    137834    -      47      Haemophilus parainfluenzae
    GCF_900638025.1_r28    6375   1668     2836   1       GCF_900638025.1   NZ_LR134481.1   97.820    1     84.502    6236      1177       86.386   2062405   141065    142241    -      47      Haemophilus parainfluenzae
    GCF_900638025.1_r28    6375   2863     3152   1       GCF_900638025.1   NZ_LR134481.1   97.820    1     84.502    6236      292        86.386   2062405   140751    141042    -      47      Haemophilus parainfluenzae
    GCF_900638025.1_r28    6375   633      1638   1       GCF_900638025.1   NZ_LR134481.1   97.820    1     84.502    6236      1015       86.386   2062405   142270    143284    -      47      Haemophilus parainfluenzae
    GCF_900638025.1_r28    6375   227      612    1       GCF_900638025.1   NZ_LR134481.1   97.820    1     84.502    6236      387        86.386   2062405   143307    143693    -      47      Haemophilus parainfluenzae
    GCF_900638025.1_r28    6375   6        154    1       GCF_900638025.1   NZ_LR134481.1   97.820    1     84.502    6236      150        86.386   2062405   143787    143936    -      47      Haemophilus parainfluenzae
    GCF_900638025.1_r28    6375   161      206    1       GCF_900638025.1   NZ_LR134481.1   97.820    1     84.502    6236      46         86.386   2062405   143716    143761    -      47      Haemophilus parainfluenzae
    GCF_000006945.2_r8     7258   20       7229   1       GCF_000006945.2   NC_003197.2     99.835    1     95.936    7246      7246       96.094   4857450   4618964   4626209   +      46      Salmonella enterica
    GCF_002949675.1_r249   1937   29       1905   3       GCF_002949675.1   NZ_CP026774.1   97.264    1     89.055    1884      1884       91.561   4395762   3336635   3338518   -      11      Shigella dysenteriae
    GCF_002949675.1_r249   1937   456      1905   3       GCF_002950215.1   NZ_CP026788.1   96.076    1     82.654    1861      1454       86.029   4659463   3959212   3960665   +      5       Shigella flexneri
    GCF_002949675.1_r249   1937   29       425    3       GCF_002950215.1   NZ_CP026788.1   96.076    1     82.654    1861      407        86.029   4659463   3958382   3958788   +      5       Shigella flexneri
    GCF_002949675.1_r249   1937   923      1905   3       GCF_003697165.2   NZ_CP033092.2   82.292    1     72.483    1594      981        88.080   4903501   926776    927756    +      8       Escherichia coli
    GCF_002949675.1_r249   1937   122      350    3       GCF_003697165.2   NZ_CP033092.2   82.292    1     72.483    1594      231        88.080   4903501   925957    926187    +      8       Escherichia coli
    GCF_002949675.1_r249   1937   369      611    3       GCF_003697165.2   NZ_CP033092.2   82.292    1     72.483    1594      244        88.080   4903501   926215    926458    +      8       Escherichia coli
    GCF_002949675.1_r249   1937   740      812    3       GCF_003697165.2   NZ_CP033092.2   82.292    1     72.483    1594      73         88.080   4903501   926590    926662    +      8       Escherichia coli
    GCF_002949675.1_r249   1937   29       93     3       GCF_003697165.2   NZ_CP033092.2   82.292    1     72.483    1594      65         88.080   4903501   925864    925928    +      8       Escherichia coli
