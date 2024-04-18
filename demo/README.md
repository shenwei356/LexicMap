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
Here we create a `species` column from the genome ID column (`sgenome`) and replace the assemby accessions with species names.

    $ csvtk head -n 28 q.gene.fasta.lexicmap.tsv \
        | csvtk mutate -t -n species -f sgenome \
        | csvtk replace -t -f species -k ass2species.map -p '(.+)' -r '{kv}' \
        | csvtk pretty -t

    query                         qlen   qstart   qend   hits   sgenome           sseqid          qcovGnm   hsp   qcovHSP   alenHSP   alenSeg   pident   slen      sstart    send      sstr   seeds   species
    ---------------------------   ----   ------   ----   ----   ---------------   -------------   -------   ---   -------   -------   -------   ------   -------   -------   -------   ----   -----   --------------------
    NC_000913.3:4166659-4168200   1542   1        1542   8      GCF_003697165.2   NZ_CP033092.2   100.000   1     100.000   1542      1542      99.287   4903501   4844587   4846128   -      26      Escherichia coli
    NC_000913.3:4166659-4168200   1542   1        1542   8      GCF_003697165.2   NZ_CP033092.2   100.000   2     100.000   1542      1542      99.287   4903501   4591684   4593225   -      26      Escherichia coli
    NC_000913.3:4166659-4168200   1542   1        1542   8      GCF_003697165.2   NZ_CP033092.2   100.000   3     100.000   1542      1542      99.287   4903501   4551515   4553056   -      26      Escherichia coli
    NC_000913.3:4166659-4168200   1542   1        1542   8      GCF_003697165.2   NZ_CP033092.2   100.000   4     100.000   1542      1542      99.287   4903501   3780640   3782181   -      26      Escherichia coli
    NC_000913.3:4166659-4168200   1542   1        1542   8      GCF_003697165.2   NZ_CP033092.2   100.000   5     100.000   1542      1542      99.287   4903501   458559    460100    +      26      Escherichia coli
    NC_000913.3:4166659-4168200   1542   1        1542   8      GCF_003697165.2   NZ_CP033092.2   100.000   6     100.000   1542      1542      99.287   4903501   1285123   1286664   +      26      Escherichia coli
    NC_000913.3:4166659-4168200   1542   1        1542   8      GCF_003697165.2   NZ_CP033092.2   100.000   7     100.000   1542      1542      99.092   4903501   4726193   4727734   -      26      Escherichia coli
    NC_000913.3:4166659-4168200   1542   1        1542   8      GCF_002950215.1   NZ_CP026788.1   100.000   1     100.000   1542      1542      99.027   4659463   3216505   3218046   +      25      Shigella flexneri
    NC_000913.3:4166659-4168200   1542   1        1542   8      GCF_002950215.1   NZ_CP026788.1   100.000   2     100.000   1542      1542      98.962   4659463   3396068   3397609   +      26      Shigella flexneri
    NC_000913.3:4166659-4168200   1542   1        1542   8      GCF_002950215.1   NZ_CP026788.1   100.000   3     100.000   1542      1542      98.962   4659463   3119331   3120872   +      26      Shigella flexneri
    NC_000913.3:4166659-4168200   1542   1        1542   8      GCF_002950215.1   NZ_CP026788.1   100.000   4     100.000   1542      1542      98.898   4659463   3355632   3357173   +      27      Shigella flexneri
    NC_000913.3:4166659-4168200   1542   1        1542   8      GCF_002950215.1   NZ_CP026788.1   100.000   5     100.000   1542      1542      98.898   4659463   4223146   4224687   +      26      Shigella flexneri
    NC_000913.3:4166659-4168200   1542   1        1542   8      GCF_002950215.1   NZ_CP026788.1   100.000   6     100.000   1542      1542      98.833   4659463   3540450   3541991   +      26      Shigella flexneri
    NC_000913.3:4166659-4168200   1542   1        1542   8      GCF_002950215.1   NZ_CP026788.1   100.000   7     100.000   1542      1542      98.768   4659463   2125377   2126917   -      19      Shigella flexneri
    NC_000913.3:4166659-4168200   1542   1        1001   8      GCF_002949675.1   NZ_CP026774.1   97.601    1     97.601    1505      1001      98.501   4395762   3646778   3647778   +      29      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   1039     1542   8      GCF_002949675.1   NZ_CP026774.1   97.601    1     97.601    1505      504       99.802   4395762   3647816   3648319   +      29      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   1        1001   8      GCF_002949675.1   NZ_CP026774.1   97.601    2     97.601    1505      1001      98.501   4395762   3061592   3062592   +      29      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   1039     1542   8      GCF_002949675.1   NZ_CP026774.1   97.601    2     97.601    1505      504       99.802   4395762   3062630   3063133   +      29      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   1        1001   8      GCF_002949675.1   NZ_CP026774.1   97.601    3     97.601    1505      1001      98.501   4395762   2810845   2811845   +      29      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   1039     1542   8      GCF_002949675.1   NZ_CP026774.1   97.601    3     97.601    1505      504       99.802   4395762   2811883   2812386   +      29      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   1        1001   8      GCF_002949675.1   NZ_CP026774.1   97.601    4     97.601    1505      1001      98.501   4395762   2768883   2769883   +      29      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   1039     1542   8      GCF_002949675.1   NZ_CP026774.1   97.601    4     97.601    1505      504       99.802   4395762   2769921   2770424   +      29      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   1        1001   8      GCF_002949675.1   NZ_CP026774.1   97.601    5     97.601    1505      1001      98.501   4395762   2636477   2637477   +      29      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   1039     1542   8      GCF_002949675.1   NZ_CP026774.1   97.601    5     97.601    1505      504       99.802   4395762   2637515   2638018   +      29      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   1        1001   8      GCF_002949675.1   NZ_CP026774.1   97.601    6     97.601    1505      1001      98.501   4395762   2536624   2537624   +      29      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   1039     1542   8      GCF_002949675.1   NZ_CP026774.1   97.601    6     97.601    1505      504       99.802   4395762   2537662   2538165   +      29      Shigella dysenteriae

### A prophage sequence

    $ lexicmap search -d demo.lmi/ q.prophage.fasta -o q.prophage.fasta.lexicmap.tsv

    $ csvtk head -n 20 q.prophage.fasta.lexicmap.tsv \
        | csvtk mutate -t -n species -f sgenome \
        | csvtk replace -t -f species -k ass2species.map -p '(.+)' -r '{kv}' \
        | csvtk pretty -t
    query         qlen    qstart   qend    hits   sgenome           sseqid          qcovGnm   hsp   qcovHSP   alenHSP   alenSeg   pident    slen      sstart    send      sstr   seeds   species
    -----------   -----   ------   -----   ----   ---------------   -------------   -------   ---   -------   -------   -------   -------   -------   -------   -------   ----   -----   ----------------
    NC_001895.1   33593   1        9357    1      GCF_003697165.2   NZ_CP033092.2   73.289    1     55.878    18771     9357      95.298    4903501   1864411   1873769   +      94      Escherichia coli
    NC_001895.1   33593   17627    22445   1      GCF_003697165.2   NZ_CP033092.2   73.289    1     55.878    18771     4819      93.795    4903501   1882193   1887011   +      94      Escherichia coli
    NC_001895.1   33593   22468    24382   1      GCF_003697165.2   NZ_CP033092.2   73.289    1     55.878    18771     1915      91.593    4903501   1887034   1888948   +      94      Escherichia coli
    NC_001895.1   33593   10592    12485   1      GCF_003697165.2   NZ_CP033092.2   73.289    1     55.878    18771     1894      90.391    4903501   1874130   1876023   +      94      Escherichia coli
    NC_001895.1   33593   15218    15283   1      GCF_003697165.2   NZ_CP033092.2   73.289    1     55.878    18771     66        100.000   4903501   1879476   1879541   +      94      Escherichia coli
    NC_001895.1   33593   17473    17552   1      GCF_003697165.2   NZ_CP033092.2   73.289    1     55.878    18771     80        80.000    4903501   1882043   1882122   +      94      Escherichia coli
    NC_001895.1   33593   14696    15072   1      GCF_003697165.2   NZ_CP033092.2   73.289    1     55.878    18771     377       70.557    4903501   1878954   1879330   +      94      Escherichia coli
    NC_001895.1   33593   10308    10570   1      GCF_003697165.2   NZ_CP033092.2   73.289    1     55.878    18771     263       95.437    4903501   1873846   1874108   +      94      Escherichia coli
    NC_001895.1   33593   24355    27239   1      GCF_003697165.2   NZ_CP033092.2   73.289    2     17.495    5877      2885      96.534    4903501   1853098   1855982   +      41      Escherichia coli
    NC_001895.1   33593   28430    30281   1      GCF_003697165.2   NZ_CP033092.2   73.289    2     17.495    5877      1852      96.058    4903501   1857173   1859024   +      41      Escherichia coli
    NC_001895.1   33593   27262    28401   1      GCF_003697165.2   NZ_CP033092.2   73.289    2     17.495    5877      1140      95.088    4903501   1856005   1857144   +      41      Escherichia coli


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
        | csvtk mutate -t -n species -f sgenome \
        | csvtk replace -t -f species -k ass2species.map -p '(.+)' -r '{kv}' \
        | csvtk pretty -t

    query                  qlen   qstart   qend   hits   sgenome           sseqid              qcovGnm   hsp   qcovHSP   alenHSP   alenSeg   pident   slen      sstart    send      sstr   seeds   species
    --------------------   ----   ------   ----   ----   ---------------   -----------------   -------   ---   -------   -------   -------   ------   -------   -------   -------   ----   -----   --------------------------
    GCF_002950215.1_r182   1718   8        1470   1      GCF_002950215.1   NZ_CP026790.1       93.714    1     93.714    1610      1463      75.598   165702    134670    136167    +      3       Shigella flexneri
    GCF_002950215.1_r182   1718   1561     1707   1      GCF_002950215.1   NZ_CP026790.1       93.714    1     93.714    1610      147       80.272   165702    136263    136412    +      3       Shigella flexneri
    GCF_002950215.1_r182   1718   8        613    1      GCF_002950215.1   NZ_CP026790.1       93.714    2     35.274    606       606       78.713   165702    132258    132873    +      1       Shigella flexneri
    GCF_002950215.1_r182   1718   8        163    1      GCF_002950215.1   NZ_CP026790.1       93.714    2     35.274    606       156       80.769   165702    134670    134825    +      1       Shigella flexneri
    GCF_009759685.1_r164   3132   20       1929   1      GCF_009759685.1   NZ_CP046654.1       97.797    1     97.797    3063      1910      93.298   3980848   1768740   1770666   +      20      Acinetobacter baumannii
    GCF_009759685.1_r164   3132   1969     3121   1      GCF_009759685.1   NZ_CP046654.1       97.797    1     97.797    3063      1153      94.016   3980848   1770695   1771855   +      20      Acinetobacter baumannii
    GCF_003697165.2_r46    2169   446      1666   1      GCF_003697165.2   NZ_CP033092.2       91.655    1     91.655    1988      1221      80.180   4903501   4490226   4491482   +      5       Escherichia coli
    GCF_003697165.2_r46    2169   1794     2023   1      GCF_003697165.2   NZ_CP033092.2       91.655    1     91.655    1988      230       84.348   4903501   4491612   4491835   +      5       Escherichia coli
    GCF_003697165.2_r46    2169   2076     2146   1      GCF_003697165.2   NZ_CP033092.2       91.655    1     91.655    1988      71        91.549   4903501   4491847   4491918   +      5       Escherichia coli
    GCF_003697165.2_r46    2169   1691     1766   1      GCF_003697165.2   NZ_CP033092.2       91.655    1     91.655    1988      76        78.947   4903501   4491509   4491584   +      5       Escherichia coli
    GCF_003697165.2_r46    2169   31       420    1      GCF_003697165.2   NZ_CP033092.2       91.655    1     91.655    1988      390       77.179   4903501   4489794   4490197   +      5       Escherichia coli
    GCF_000006945.2_r109   3788   37       3278   1      GCF_000006945.2   NC_003197.2         97.941    1     97.941    3710      3242      96.730   4857450   4633813   4637055   -      24      Salmonella enterica
    GCF_000006945.2_r109   3788   3301     3768   1      GCF_000006945.2   NC_003197.2         97.941    1     97.941    3710      468       94.658   4857450   4633323   4633787   -      24      Salmonella enterica
    GCF_900638025.1_r28    6375   3174     5997   1      GCF_900638025.1   NZ_LR134481.1       83.247    1     83.247    5307      2824      84.809   2062405   137873    140730    -      45      Haemophilus parainfluenzae
    GCF_900638025.1_r28    6375   6032     6339   1      GCF_900638025.1   NZ_LR134481.1       83.247    1     83.247    5307      308       96.753   2062405   137524    137834    -      45      Haemophilus parainfluenzae
    GCF_900638025.1_r28    6375   1668     2836   1      GCF_900638025.1   NZ_LR134481.1       83.247    1     83.247    5307      1169      88.879   2062405   141065    142241    -      45      Haemophilus parainfluenzae
    GCF_900638025.1_r28    6375   633      1638   1      GCF_900638025.1   NZ_LR134481.1       83.247    1     83.247    5307      1006      87.276   2062405   142270    143284    -      45      Haemophilus parainfluenzae
    GCF_900638025.1_r28    6375   1668     2250   1      GCF_900638025.1   NZ_LR134481.1       83.247    2     9.145     583       583       92.453   2062405   141663    142241    -      1       Haemophilus parainfluenzae
    GCF_000006945.2_r8     7258   20       7229   1      GCF_000006945.2   NC_003197.2         99.339    1     99.339    7210      7210      96.574   4857450   4618964   4626209   +      54      Salmonella enterica
    GCF_000006945.2_r8     7258   20       2148   1      GCF_000006945.2   NC_003197.2         99.339    2     29.333    2129      2129      96.148   4857450   4618964   4621107   +      1       Salmonella enterica
    GCF_001544255.1_r110   9910   17       9910   1      GCF_001544255.1   NZ_BCQD01000005.1   99.839    1     99.839    9894      9894      97.645   191690    155488    165428    +      117     Enterococcus faecium
    GCF_002949675.1_r249   1937   29       1937   3      GCF_002949675.1   NZ_CP026774.1       98.554    1     98.554    1909      1909      92.038   4395762   3336601   3338518   -      19      Shigella dysenteriae
    GCF_002949675.1_r249   1937   456      1937   3      GCF_002950215.1   NZ_CP026788.1       98.554    1     98.554    1909      1482      86.640   4659463   3959212   3960699   +      13      Shigella flexneri
    GCF_002949675.1_r249   1937   29       518    3      GCF_002950215.1   NZ_CP026788.1       98.554    1     98.554    1909      490       83.878   4659463   3958382   3958882   +      13      Shigella flexneri
    GCF_002949675.1_r249   1937   547      611    3      GCF_002950215.1   NZ_CP026788.1       98.554    1     98.554    1909      65        83.077   4659463   3958912   3958976   +      13      Shigella flexneri
    GCF_002949675.1_r249   1937   923      1937   3      GCF_003697165.2   NZ_CP033092.2       83.893    1     83.893    1625      1015      88.768   4903501   926776    927790    +      13      Escherichia coli
