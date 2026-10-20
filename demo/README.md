## Demo datasets

### Reference genomes

We use 15 bacterial genomes for demonstration.

    $ seqkit stats -a refs/* -T | csvtk cut -t -F -f "-Q*\(*,-Avg*" | csvtk pretty -t
    file                         format   type   num_seqs   sum_len   min_len   avg_len     max_len   Q1        Q2        Q3        sum_gap   N50       N50_num   GC(%)   sum_n
    --------------------------   ------   ----   --------   -------   -------   ---------   -------   -------   -------   -------   -------   -------   -------   -----   -----
    refs/GCF_000006945.2.fa.gz   FASTA    DNA    2          4951383   93933     2475691.5   4857450   93933     2475692   4857450   0         4857450   1         52.24   0    
    refs/GCF_000017205.1.fa.gz   FASTA    DNA    1          6588339   6588339   6588339.0   6588339   6588339   6588339   6588339   0         6588339   1         66.45   0    
    refs/GCF_000148585.2.fa.gz   FASTA    DNA    1          1868883   1868883   1868883.0   1868883   1868883   1868883   1868883   0         1868883   1         40.45   0    
    refs/GCF_000392875.1.fa.gz   FASTA    DNA    3          2881400   274762    960466.7    1924212   478594    682426    1303319   0         1924212   1         37.58   4965 
    refs/GCF_000742135.1.fa.gz   FASTA    DNA    5          5545784   16331     1109156.8   5284261   42420     95930     106842    0         5284261   1         56.99   1100 
    refs/GCF_001027105.1.fa.gz   FASTA    DNA    2          2782562   27490     1391281.0   2755072   27490     1391281   2755072   0         2755072   1         32.84   0    
    refs/GCF_001096185.1.fa.gz   FASTA    DNA    24         2117177   366       88215.7     543880    9252      48418     133692    0         163740    4         40.12   4    
    refs/GCF_001457655.1.fa.gz   FASTA    DNA    1          1890645   1890645   1890645.0   1890645   1890645   1890645   1890645   0         1890645   1         38.23   0    
    refs/GCF_001544255.1.fa.gz   FASTA    DNA    38         2484851   563       65390.8     386577    1364      20536     115935    0         191690    5         38.08   0    
    refs/GCF_002949675.1.fa.gz   FASTA    DNA    2          4578459   182697    2289229.5   4395762   182697    2289230   4395762   0         4395762   1         50.98   0    
    refs/GCF_002950215.1.fa.gz   FASTA    DNA    3          4938295   113130    1646098.3   4659463   139416    165702    2412582   0         4659463   1         50.67   0    
    refs/GCF_003697165.2.fa.gz   FASTA    DNA    2          5034834   131333    2517417.0   4903501   131333    2517417   4903501   0         4903501   1         50.64   0    
    refs/GCF_006742205.1.fa.gz   FASTA    DNA    2          2427041   4439      1213520.5   2422602   4439      1213520   2422602   0         2422602   1         32.29   0    
    refs/GCF_009759685.1.fa.gz   FASTA    DNA    2          3990388   9540      1995194.0   3980848   9540      1995194   3980848   0         3980848   1         39.16   0    
    refs/GCF_900638025.1.fa.gz   FASTA    DNA    1          2062405   2062405   2062405.0   2062405   2062405   2062405   2062405   0         2062405   1         39.37   0    

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

Building

    $ lexicmap index -I refs -O demo.lmi
    15:03:45.351 [INFO] LexicMap v0.10.0
    15:03:45.351 [INFO]   https://github.com/shenwei356/LexicMap
    15:03:45.351 [INFO] 
    15:03:45.351 [INFO]  CWD: /home/shenwei/go/src/github.com/shenwei356/LexicMap/demo
    15:03:45.351 [INFO]  CMD: lexicmap index -I refs -O demo.lmi
    15:03:45.351 [INFO] DATE: 2026-10-09
    15:03:45.351 [INFO] 
    15:03:45.360 [INFO] checking input files ...
    15:03:45.360 [INFO]   scanning files from directory: refs
    15:03:45.361 [INFO]   15 input file(s) given
    15:03:45.361 [INFO] 
    15:03:45.361 [INFO] --------------------- [ main parameters ] ---------------------
    15:03:45.361 [INFO] 
    15:03:45.361 [INFO] input and output:
    15:03:45.361 [INFO]   input directory: refs
    15:03:45.361 [INFO]     regular expression of input files: (?i)\.(f[aq](st[aq])?|fna)(\.gz|\.xz|\.zst|\.bz2)?$
    15:03:45.361 [INFO]     *regular expression for extracting reference name from file name: (?i)(.+)\.(f[aq](st[aq])?|fna)(\.gz|\.xz|\.zst|\.bz2)?$
    15:03:45.361 [INFO]     *regular expressions for filtering out sequences: []
    15:03:45.361 [INFO]   min sequence length: 31
    15:03:45.361 [INFO]   max genome size: 20000000
    15:03:45.361 [INFO]   output directory: demo.lmi
    15:03:45.361 [INFO] 
    15:03:45.361 [INFO] mask generation:
    15:03:45.361 [INFO]   k-mer size: 31
    15:03:45.361 [INFO]   number of masks: 20000
    15:03:45.361 [INFO]   rand seed: 1
    15:03:45.361 [INFO] 
    15:03:45.361 [INFO] seed data:
    15:03:45.361 [INFO]   maximum sketching desert length: 100
    15:03:45.361 [INFO]   distance of k-mers to fill deserts: 50
    15:03:45.361 [INFO]   seeds data chunks: 16
    15:03:45.361 [INFO]   seeds data indexing partitions: 4096
    15:03:45.361 [INFO] 
    15:03:45.361 [INFO] general:
    15:03:45.361 [INFO]   genome batch size: 5000
    15:03:45.361 [INFO]   threads: 16
    15:03:45.361 [INFO]   batch merge threads: 8
    15:03:45.361 [INFO] 
    15:03:45.361 [INFO] 
    15:03:45.361 [INFO] --------------------- [ generating masks ] ---------------------
    15:03:45.365 [INFO] 
    15:03:45.365 [INFO] --------------------- [ building index ] ---------------------
    15:03:45.494 [INFO] 
    15:03:45.494 [INFO]   ------------------------[ batch 1/1 ]------------------------
    15:03:45.494 [INFO]   building index for batch 1 with 15 files...
    processed files:  15 / 15 [======================================] ETA: 0s. done
    15:03:46.915 [INFO]   writing seeds...
    15:03:47.042 [INFO]   finished writing seeds in 127.246877ms
    15:03:47.042 [INFO]   finished building index for batch 1 in: 1.548079342s
    15:03:47.043 [INFO] 
    15:03:47.043 [INFO] finished building LexicMap index from 15 files with 20000 masks in 1.691668836s
    15:03:47.043 [INFO] LexicMap index saved: demo.lmi
    15:03:47.043 [INFO] 
    15:03:47.043 [INFO] elapsed time: 1.691699357s
    15:03:47.043 [INFO]
   
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


    $ dirsize demo.lmi/
    demo.lmi/: 78.18 MiB (81,974,001)
     65.08 MiB      seeds
     12.94 MiB      genomes
    156.28 KiB      masks.bin
         921 B      info.toml
         375 B      genomes.map.bin
           0 B      genomes.chunks.bin
           
From v0.10.0, it is optional to create adaptive two-level indexes of seed data.
It can improve seed-matching performance for batch queries with `lexicmap search`.
See the [usage](https://bioinf.shenwei.me/LexicMap/usage/utils/reindex-seeds2) for more details.

    $ lexicmap utils reindex-seeds2 -d demo.lmi/
    15:12:22.910 [INFO] LexicMap v0.10.0
    15:12:22.910 [INFO]   https://github.com/shenwei356/LexicMap
    15:12:22.910 [INFO] 
    15:12:22.910 [INFO]  CWD: /home/shenwei/go/src/github.com/shenwei356/LexicMap/demo
    15:12:22.910 [INFO]  CMD: lexicmap utils reindex-seeds2 -d demo.lmi/
    15:12:22.910 [INFO] DATE: 2026-10-09
    15:12:22.910 [INFO] 
    15:12:22.911 [INFO] recreating adaptive seed indexes with 4096 partitions and a 8.0 KiB block threshold for: demo.lmi/
    processed files:  16 / 16 [======================================] ETA: 0s. done
    15:12:23.098 [INFO] primary-prefix blocks: 81920000 total, 2107085 non-empty, 0 indexed (0.00% of non-empty)
    15:12:23.098 [INFO] idx15 size: 768 B
    15:12:23.099 [INFO] 
    15:12:23.099 [INFO] elapsed time: 188.077603ms
    15:12:23.099 [INFO] 

## Searching

### A 16S rRNA gene sequence

    $ lexicmap search -d demo.lmi/  q.gene.fasta -o q.gene.fasta.lexicmap.tsv
    15:13:46.090 [INFO] LexicMap v0.10.0
    15:13:46.090 [INFO]   https://github.com/shenwei356/LexicMap
    15:13:46.090 [INFO] 
    15:13:46.090 [INFO]  CWD: /home/shenwei/go/src/github.com/shenwei356/LexicMap/demo
    15:13:46.090 [INFO]  CMD: lexicmap search -d demo.lmi/ q.gene.fasta -o q.gene.fasta.lexicmap.tsv
    15:13:46.090 [INFO] DATE: 2026-10-09
    15:13:46.090 [INFO] 
    15:13:46.090 [INFO] checking input files ...
    15:13:46.090 [INFO]   1 input file given: q.gene.fasta
    15:13:46.090 [INFO] 
    15:13:46.090 [INFO] loading index: demo.lmi/
    15:13:46.091 [INFO]   reading masks...
    15:13:46.092 [INFO]   reading indexes of seeds (k-mer-value) data...
    15:13:46.092 [INFO]   creating searcher pools for 16 seed data files, each with 4 searchers...
    15:13:46.389 [INFO]   creating reader pools for 1 genome batches, each with 16 readers...
    15:13:46.389 [INFO] index loaded in 298.870684ms
    15:13:46.389 [INFO] 
    15:13:46.389 [INFO] searching with 16 threads...
    15:13:46.389 [INFO]   maximum number of concurrent queries: 8, force garbage collection for every 64 queries

    15:13:46.424 [INFO] 
    15:13:46.424 [INFO] processed queries: 1, speed: 1709.413 queries per minute
    15:13:46.424 [INFO] 100.0000% (1/1) queries matched
    15:13:46.424 [INFO] done searching
    15:13:46.424 [INFO] search results saved to: q.gene.fasta.lexicmap.tsv
    15:13:46.425 [INFO] 
    15:13:46.425 [INFO] elapsed time: 334.465362ms
    15:13:46.425 [INFO] 

Result preview.
Here we create a `species` column from the genome ID column (`sgenome`) and replace the assembly accessions with species names.

    $ csvtk head -n 21 q.gene.fasta.lexicmap.tsv \
        | csvtk mutate -t -n species -f sgenome \
        | csvtk replace -t -f species -k ass2species.map -p '(.+)' -r '{kv}' \
        | csvtk pretty -t

    query                         qlen   hits   sgenome           sseqid          qcovGnm   cls   hsp   qcovHSP   alenHSP   pident   gaps   qstart   qend   sstart    send      sstr   slen      evalue     bitscore   species             
    ---------------------------   ----   ----   ---------------   -------------   -------   ---   ---   -------   -------   ------   ----   ------   ----   -------   -------   ----   -------   --------   --------   --------------------
    NC_000913.3:4166659-4168200   1542   15     GCF_003697165.2   NZ_CP033092.2   100.000   1     1     100.000   1542      99.805   0      1        1542   458559    460100    +      4903501   0.00e+00   2767       Escherichia coli    
    NC_000913.3:4166659-4168200   1542   15     GCF_003697165.2   NZ_CP033092.2   100.000   2     2     100.000   1542      99.805   0      1        1542   1285123   1286664   +      4903501   0.00e+00   2767       Escherichia coli    
    NC_000913.3:4166659-4168200   1542   15     GCF_003697165.2   NZ_CP033092.2   100.000   3     3     100.000   1542      99.805   0      1        1542   3780640   3782181   -      4903501   0.00e+00   2767       Escherichia coli    
    NC_000913.3:4166659-4168200   1542   15     GCF_003697165.2   NZ_CP033092.2   100.000   4     4     100.000   1542      99.805   0      1        1542   4551515   4553056   -      4903501   0.00e+00   2767       Escherichia coli    
    NC_000913.3:4166659-4168200   1542   15     GCF_003697165.2   NZ_CP033092.2   100.000   5     5     100.000   1542      99.805   0      1        1542   4591684   4593225   -      4903501   0.00e+00   2767       Escherichia coli    
    NC_000913.3:4166659-4168200   1542   15     GCF_003697165.2   NZ_CP033092.2   100.000   6     6     100.000   1542      99.805   0      1        1542   4726193   4727734   -      4903501   0.00e+00   2767       Escherichia coli    
    NC_000913.3:4166659-4168200   1542   15     GCF_003697165.2   NZ_CP033092.2   100.000   7     7     100.000   1542      99.805   0      1        1542   4844587   4846128   -      4903501   0.00e+00   2767       Escherichia coli    
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   1     1     100.000   1542      99.676   0      1        1542   3216505   3218046   +      4659463   0.00e+00   2758       Shigella flexneri   
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   2     2     100.000   1542      99.611   0      1        1542   3119331   3120872   +      4659463   0.00e+00   2755       Shigella flexneri   
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   3     3     100.000   1542      99.611   0      1        1542   3396068   3397609   +      4659463   0.00e+00   2755       Shigella flexneri   
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   4     4     100.000   1542      99.546   0      1        1542   3355632   3357173   +      4659463   0.00e+00   2749       Shigella flexneri   
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   5     5     100.000   1542      99.546   0      1        1542   4223146   4224687   +      4659463   0.00e+00   2749       Shigella flexneri   
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   6     6     100.000   1542      99.546   1      1        1542   2125377   2126917   -      4659463   0.00e+00   2746       Shigella flexneri   
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   7     7     100.000   1542      99.481   0      1        1542   3540450   3541991   +      4659463   0.00e+00   2746       Shigella flexneri   
    NC_000913.3:4166659-4168200   1542   15     GCF_002949675.1   NZ_CP026774.1   100.000   1     1     100.000   1542      99.027   0      1        1542   1662010   1663551   -      4395762   0.00e+00   2713       Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   15     GCF_002949675.1   NZ_CP026774.1   100.000   2     2     100.000   1542      99.027   0      1        1542   2536624   2538165   +      4395762   0.00e+00   2713       Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   15     GCF_002949675.1   NZ_CP026774.1   100.000   3     3     100.000   1542      99.027   0      1        1542   2636477   2638018   +      4395762   0.00e+00   2713       Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   15     GCF_002949675.1   NZ_CP026774.1   100.000   4     4     100.000   1542      99.027   0      1        1542   2768883   2770424   +      4395762   0.00e+00   2713       Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   15     GCF_002949675.1   NZ_CP026774.1   100.000   5     5     100.000   1542      99.027   0      1        1542   2810845   2812386   +      4395762   0.00e+00   2713       Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   15     GCF_002949675.1   NZ_CP026774.1   100.000   6     6     100.000   1542      99.027   0      1        1542   3061592   3063133   +      4395762   0.00e+00   2713       Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   15     GCF_002949675.1   NZ_CP026774.1   100.000   7     7     100.000   1542      99.027   0      1        1542   3646778   3648319   +      4395762   0.00e+00   2713       Shigella dysenteriae

Blast-style format:

```
$ lexicmap search -d demo.lmi/ q.gene.fasta --all \
    | lexicmap utils 2blast --kv-file-genome ass2species.map

Query = NC_000913.3:4166659-4168200
Length = 1542

[Subject genome #1/15] = GCF_003697165.2 Escherichia coli
Query coverage per genome = 100.000%

>NZ_CP033092.2
Length = 4903501

 HSP cluster #1, HSP #1
 Score = 2767 bits, Expect = 0.00e+00
 Query coverage per seq = 100.000%, Aligned length = 1542, Identities = 99.805%, Gaps = 0
 Query range = 1-1542, Subject range = 458559-460100, Strand = Plus/Plus

Query  1       AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAA  60
               ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  458559  AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAA  458618

Query  61      GTCGAACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAA  120
               ||||||||||||||||||| |||||||| |||||||||||||||||||||||||||||||
Sbjct  458619  GTCGAACGGTAACAGGAAGCAGCTTGCTGCTTTGCTGACGAGTGGCGGACGGGTGAGTAA  458678

Query  121     TGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCAT  180
               ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  458679  TGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCAT  458738

Query  181     AACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATG  240
               ||||||||||||||||||||||||||||| ||||||||||||||||||||||||||||||
Sbjct  458739  AACGTCGCAAGACCAAAGAGGGGGACCTTAGGGCCTCTTGCCATCGGATGTGCCCAGATG  458798

Query  241     GGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGA  300
               ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  458799  GGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGA  458858

Query  301     GGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGG  360
               ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  458859  GGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGG  458918

Query  361     GGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCT  420
               ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  458919  GGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCT  458978

Query  421     TCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATT  480
               ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  458979  TCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATT  459038

Query  481     GACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAG  540
               ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  459039  GACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAG  459098

Query  541     GGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCA  600
               ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  459099  GGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCA  459158

Query  601     GATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTC  660
               ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  459159  GATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTC  459218

Query  661     GTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACC  720
               ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  459219  GTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACC  459278

Query  721     GGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCA  780
               ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  459279  GGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCA  459338

Query  781     AACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCC  840
               ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  459339  AACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCC  459398

Query  841     CTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCA  900
               ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  459399  CTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCA  459458

Query  901     AGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAAT  960
               ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  459459  AGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAAT  459518

Query  961     TCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAG  1020
               ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  459519  TCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAG  459578

Query  1021    AATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGA  1080
               ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  459579  AATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGA  459638

Query  1081    AATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGC  1140
               ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  459639  AATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGC  459698

Query  1141    CGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTC  1200
               ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  459699  CGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTC  459758

Query  1201    ATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCG  1260
               ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  459759  ATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCG  459818

Query  1261    ACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAAC  1320
               ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  459819  ACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAAC  459878

Query  1321    TCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGT  1380
               ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  459879  TCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGT  459938

Query  1381    TCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGT  1440
               ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  459939  TCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGT  459998

Query  1441    AGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAA  1500
               ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  459999  AGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAA  460058

Query  1501    CAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTTA  1542
               ||||||||||||||||||||||||||||||||||||||||||
Sbjct  460059  CAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTTA  460100

```


### A prophage sequence

    $ lexicmap search -d demo.lmi/ q.prophage.fasta -o q.prophage.fasta.lexicmap.tsv

    $ csvtk head -n 20 q.prophage.fasta.lexicmap.tsv \
        | csvtk mutate -t -n species -f sgenome \
        | csvtk replace -t -f species -k ass2species.map -p '(.+)' -r '{kv}' \
        | csvtk pretty -t

    query         qlen    hits   sgenome           sseqid          qcovGnm   cls   hsp   qcovHSP   alenHSP   pident   gaps   qstart   qend    sstart    send      sstr   slen      evalue      bitscore   species             
    -----------   -----   ----   ---------------   -------------   -------   ---   ---   -------   -------   ------   ----   ------   -----   -------   -------   ----   -------   ---------   --------   --------------------
    NC_001895.1   33593   3      GCF_003697165.2   NZ_CP033092.2   77.474    1     1     27.890    9371      97.716   2      1        9369    1864411   1873781   +      4903501   0.00e+00    15928      Escherichia coli    
    NC_001895.1   33593   3      GCF_003697165.2   NZ_CP033092.2   77.474    1     2     0.301     101       98.020   0      10308    10408   1873846   1873946   +      4903501   1.72e-43    174        Escherichia coli    
    NC_001895.1   33593   3      GCF_003697165.2   NZ_CP033092.2   77.474    2     3     20.665    6942      96.528   4      17441    24382   1882011   1888948   +      4903501   0.00e+00    11427      Escherichia coli    
    NC_001895.1   33593   3      GCF_003697165.2   NZ_CP033092.2   77.474    3     4     17.685    5941      97.980   0      24355    30295   1853098   1859038   +      4903501   0.00e+00    10174      Escherichia coli    
    NC_001895.1   33593   3      GCF_003697165.2   NZ_CP033092.2   77.474    4     5     8.880     2983      91.720   0      10308    13290   1873846   1876828   +      4903501   0.00e+00    4266       Escherichia coli    
    NC_001895.1   33593   3      GCF_003697165.2   NZ_CP033092.2   77.474    5     6     2.438     820       84.390   1      14540    15358   1878798   1879617   +      4903501   8.15e-261   899        Escherichia coli    
    NC_001895.1   33593   3      GCF_002949675.1   NZ_CP026774.1   0.976     1     1     0.976     331       85.801   3      13919    14246   3704319   3704649   -      4395762   5.95e-106   383        Shigella dysenteriae
    NC_001895.1   33593   3      GCF_002950215.1   NZ_CP026788.1   0.185     1     1     0.185     64        85.938   4      14837    14898   71092     71153     +      4659463   3.94e-11    66         Shigella flexneri   
    NC_001895.1   33593   3      GCF_002950215.1   NZ_CP026788.1   0.185     2     2     0.185     64        85.938   4      14837    14898   4261071   4261132   +      4659463   3.94e-11    66         Shigella flexneri   

### Simulated Oxford Nanopore R10.4.1 long-reads

Here we accelerate searching by
- using the flag `--top-n-genomes 5`, which returns only the top genome matches.
- using the flag `--top-n-chains 1`, which returns only the best match in a genome.

        $ lexicmap search -d demo.lmi/ q.long-reads.fasta.gz -o q.long-reads.fasta.gz.lexicmap.tsv.gz \
            --min-qcov-per-hsp 70 \
            --top-n-genomes 5 \
            --top-n-chains 1
            
        15:16:43.151 [INFO] LexicMap v0.10.0
        15:16:43.151 [INFO]   https://github.com/shenwei356/LexicMap
        15:16:43.151 [INFO] 
        15:16:43.151 [INFO]  CWD: /home/shenwei/go/src/github.com/shenwei356/LexicMap/demo
        15:16:43.151 [INFO]  CMD: lexicmap search -d demo.lmi/ q.long-reads.fasta.gz -o q.long-reads.fasta.gz.lexicmap.tsv.gz --min-qcov-per-hsp 70 --top-n-genomes 5 --top-n-chains 1
        15:16:43.151 [INFO] DATE: 2026-10-09
        15:16:43.151 [INFO] 
        15:16:43.151 [INFO] checking input files ...
        15:16:43.151 [INFO]   1 input file given: q.long-reads.fasta.gz
        15:16:43.151 [INFO] 
        15:16:43.151 [INFO] loading index: demo.lmi/
        15:16:43.152 [INFO]   reading masks...
        15:16:43.152 [INFO]   reading indexes of seeds (k-mer-value) data...
        15:16:43.152 [INFO]   creating searcher pools for 16 seed data files, each with 4 searchers...
        15:16:43.410 [INFO]   creating reader pools for 1 genome batches, each with 16 readers...
        15:16:43.410 [INFO] index loaded in 258.906589ms
        15:16:43.410 [INFO] 
        15:16:43.410 [INFO] searching with 16 threads...
        15:16:43.410 [INFO]   keep the top 5 genomes
        15:16:43.410 [INFO]   keep the top 1 chains
        15:16:43.410 [INFO]   maximum number of concurrent queries: 8, force garbage collection for every 64 queries
        processed queries: 3456, speed: 10623.220 queries per minute
        15:17:04.245 [INFO] 
        15:17:04.245 [INFO] processed queries: 3692, speed: 10632.490 queries per minute
        15:17:04.245 [INFO] 76.2459% (2815/3692) queries matched
        15:17:04.245 [INFO] done searching
        15:17:04.245 [INFO] search results saved to: q.long-reads.fasta.gz.lexicmap.tsv.gz
        15:17:04.251 [INFO] 
        15:17:04.251 [INFO] elapsed time: 21.099691338s
        15:17:04.251 [INFO] 


Result overview:

    $ csvtk head -n 10 q.long-reads.fasta.gz.lexicmap.tsv.gz \
        | csvtk mutate -t -n species -f sgenome \
        | csvtk replace -t -f species -k ass2species.map -p '(.+)' -r '{kv}' \
        | csvtk pretty -t

    query                  qlen    hits   sgenome           sseqid          qcovGnm   cls   hsp   qcovHSP   alenHSP   pident   gaps   qstart   qend    sstart    send      sstr   slen      evalue     bitscore   species                   
    --------------------   -----   ----   ---------------   -------------   -------   ---   ---   -------   -------   ------   ----   ------   -----   -------   -------   ----   -------   --------   --------   --------------------------
    GCF_900638025.1_r28    6375    1      GCF_900638025.1   NZ_LR134481.1   99.953    1     1     99.953    6542      92.831   261    3        6374    137490    143940    -      2062405   0.00e+00   8982       Haemophilus parainfluenzae
    GCF_002949675.1_r183   8176    1      GCF_002949675.1   NZ_CP026774.1   99.682    1     1     99.682    8188      98.107   82     27       8176    4194298   4202441   +      4395762   0.00e+00   13876      Shigella dysenteriae      
    GCF_003697165.2_r46    2169    1      GCF_003697165.2   NZ_CP033092.2   93.177    1     1     93.177    2101      90.243   111    12       2032    4489774   4491843   +      4903501   0.00e+00   2632       Escherichia coli          
    GCF_000006945.2_r8     7258    1      GCF_000006945.2   NC_003197.2     99.724    1     1     99.724    7301      97.630   91     20       7257    4618964   4626236   +      4857450   0.00e+00   12130      Salmonella enterica       
    GCF_000017205.1_r183   14521   1      GCF_000017205.1   NC_009656.1     99.910    1     1     99.910    14685     96.765   269    14       14521   3874715   3889307   +      6588339   0.00e+00   23596      Pseudomonas aeruginosa    
    GCF_009759685.1_r164   3132    1      GCF_009759685.1   NZ_CP046654.1   99.042    1     1     99.042    3152      94.670   86     20       3121    1768740   1771855   +      3980848   0.00e+00   4735       Acinetobacter baumannii   
    GCF_000006945.2_r109   3788    2      GCF_000006945.2   NC_003197.2     99.393    1     1     99.393    3799      97.026   67     10       3774    4633318   4637083   -      4857450   0.00e+00   6145       Salmonella enterica       
    GCF_000006945.2_r109   3788    2      GCF_002950215.1   NZ_CP026788.1   70.301    1     1     70.301    2700      82.296   60     840      3502    2488829   2491505   -      4659463   0.00e+00   2565       Shigella flexneri         
    GCF_001027105.1_r148   20294   1      GCF_001027105.1   NZ_CP011526.1   99.921    1     1     99.921    20481     97.100   307    16       20293   2352020   2372396   +      2755072   0.00e+00   33358      Staphylococcus aureus     
    GCF_002950215.1_r182   1718    1      GCF_002950215.1   NZ_CP026790.1   99.593    1     1     99.593    1781      89.163   95     5        1715    134666    136421    +      165702    0.00e+00   2089       Shigella flexneri         

Blast-style format:

```
$ seqkit seq -g -M 200 q.long-reads.fasta.gz \
    | lexicmap search -d demo.lmi/ -a \
    | csvtk filter2 -t -f '$pident >80 && $pident < 90' \
    | csvtk head -t -n 1 \
    | lexicmap utils 2blast --kv-file-genome ass2species.map

Query = GCF_003697165.2_r40
Length = 186

[Subject genome #1/2] = GCF_002950215.1 Shigella flexneri
Query coverage per genome = 93.548%

>NZ_CP026788.1 
Length = 4659463

 HSP cluster #1, HSP #1
 Score = 212 bits, Expect = 1.18e-54
 Query coverage per seq = 93.548%, Aligned length = 177, Identities = 88.701%, Gaps = 6
 Query range = 13-186, Subject range = 1124816-1124989, Strand = Plus/Plus

Query  13       CGGAAACTGAAACA-CCAGATTCTACGATGATTATGATGATTTA-TGCTTTCTTTACTAA  70
                |||||||||||||| |||||||||| | |||||||||||||||| |||||||||| ||||
Sbjct  1124816  CGGAAACTGAAACAACCAGATTCTATGTTGATTATGATGATTTAATGCTTTCTTTGCTAA  1124875

Query  71       AAAGTAAGCGGCCAAAAAAATGAT-AACACCTGTAATGAGTATCAGAAAAGACACGGTAA  129
                ||    |||||||||||||||||| |||||||||||||||||||||||||||||||||||
Sbjct  1124876  AA--GCAGCGGCCAAAAAAATGATTAACACCTGTAATGAGTATCAGAAAAGACACGGTAA  1124933

Query  130      GAAAACACTCTTTTGGATACCTAGAGTCTGATAAGCGATTATTCTCTCTATGTTACT  186
                 || |||||||||    |||||  |||||||||||||||||||||||| |||| |||
Sbjct  1124934  AAAGACACTCTTTGAAGTACCTGAAGTCTGATAAGCGATTATTCTCTCCATGT-ACT  1124989
```

## Genome similarity searching

A single query

    $ lexicmap genome search -d demo.lmi/ refs/GCF_003697165.2.fa.gz | csvtk pretty -t
    10:15:08.632 [INFO] LexicMap v0.10.0
    10:15:08.632 [INFO]   https://github.com/shenwei356/LexicMap
    10:15:08.632 [INFO] 
    10:15:08.632 [INFO]  CWD: /home/shenwei/go/src/github.com/shenwei356/LexicMap/demo
    10:15:08.632 [INFO]  CMD: lexicmap genome search -d demo.lmi/ refs/GCF_003697165.2.fa.gz
    10:15:08.632 [INFO] DATE: 2026-10-20
    10:15:08.632 [INFO] 
    10:15:08.632 [INFO] checking input files ...
    10:15:08.633 [INFO]   1 input file given: refs/GCF_003697165.2.fa.gz
    10:15:08.633 [INFO] 
    10:15:08.633 [INFO] loading index: demo.lmi/
    10:15:08.633 [INFO]   reading masks...
    10:15:08.634 [INFO]   reading indexes of seeds (k-mer-value) data...
    10:15:08.634 [INFO]   creating searcher pools for 16 seed data files, each with 4 searchers...
    10:15:08.684 [INFO]   creating reader pools for 1 genome batches, each with 16 readers...
    10:15:08.685 [INFO] index loaded in 52.056919ms
    10:15:08.685 [INFO] 
    10:15:08.685 [INFO] searching with 16 threads...
    10:15:08.685 [INFO]   screen genomes with seed data of 1024 masks
    10:15:08.685 [INFO]     keep the top 10 screened genomes, including all ties at the cutoff score
    10:15:08.685 [INFO]   minimum query coverage per HSP: 30.00%
    10:15:08.685 [INFO]   minimum base identity in a HSP segment: 70.00%
    10:15:08.685 [INFO]   maximum evalue: 1.00e-15
    10:15:08.685 [INFO]   maximum number of concurrent queries: 8, force garbage collection for every 4 queries
    10:15:08.685 [INFO]   threads per query: 16

    10:15:12.097 [INFO] 
    10:15:12.097 [INFO] processed queries: 1, speed: 17.583 queries per minute
    10:15:12.097 [INFO] 100.0000% (1/1) queries matched
    10:15:12.097 [INFO] done searching
    10:15:12.098 [INFO] 
    10:15:12.098 [INFO] elapsed time: 3.46511184s
    10:15:12.098 [INFO] 
    query             subject           ANI        qAF       sAF       qcontigs   qsize     scontigs   ssize     sname
    ---------------   ---------------   --------   -------   -------   --------   -------   --------   -------   -----
    GCF_003697165.2   GCF_003697165.2   100.0000   100.000   100.000   2          5034834   2          5034834        
    GCF_003697165.2   GCF_002949675.1   96.8012    67.613    74.353    2          5034834   2          4578459        
    GCF_003697165.2   GCF_002950215.1   96.7220    73.252    74.685    2          5034834   3          4938295        
    GCF_003697165.2   GCF_000006945.2   81.6285    55.366    56.299    2          5034834   2          4951383        
    GCF_003697165.2   GCF_000742135.1   80.5980    46.226    41.967    2          5034834   5          5545784        
    
    # with species name
    $ lexicmap genome search -d demo.lmi/ refs/GCF_003697165.2.fa.gz \
            --taxdump ~/.taxonkit/ --genome2taxid taxid.map --show-species-name --quiet \
        | csvtk pretty -t
        
    query             subject           ANI        qAF       sAF       qcontigs   qsize     scontigs   ssize     sname                
    ---------------   ---------------   --------   -------   -------   --------   -------   --------   -------   ---------------------
    GCF_003697165.2   GCF_003697165.2   100.0000   100.000   100.000   2          5034834   2          5034834   Escherichia coli     
    GCF_003697165.2   GCF_002949675.1   96.8012    67.613    74.353    2          5034834   2          4578459   Shigella dysenteriae 
    GCF_003697165.2   GCF_002950215.1   96.7220    73.252    74.685    2          5034834   3          4938295   Shigella flexneri    
    GCF_003697165.2   GCF_000006945.2   81.6285    55.366    56.299    2          5034834   2          4951383   Salmonella enterica  
    GCF_003697165.2   GCF_000742135.1   80.5980    46.226    41.967    2          5034834   5          5545784   Klebsiella pneumoniae

Self-alignment

    $ lexicmap genome search -d demo.lmi/ refs/* -n 1 --quiet | csvtk pretty -t
    query             subject           ANI        qAF       sAF       qcontigs   qsize     scontigs   ssize     sname
    ---------------   ---------------   --------   -------   -------   --------   -------   --------   -------   -----
    GCF_000148585.2   GCF_000148585.2   100.0000   100.000   100.000   1          1868883   1          1868883        
    GCF_001457655.1   GCF_001457655.1   100.0000   100.000   100.000   1          1890645   1          1890645        
    GCF_001096185.1   GCF_001096185.1   100.0000   100.000   100.000   24         2117177   24         2117177        
    GCF_000392875.1   GCF_000392875.1   100.0000   100.000   100.000   3          2881400   3          2881400        
    GCF_001027105.1   GCF_001027105.1   100.0000   100.000   100.000   2          2782562   2          2782562        
    GCF_001544255.1   GCF_001544255.1   100.0000   100.000   100.000   38         2484851   38         2484851        
    GCF_000006945.2   GCF_000006945.2   100.0000   100.000   100.000   2          4951383   2          4951383        
    GCF_006742205.1   GCF_006742205.1   100.0000   100.000   100.000   2          2427041   2          2427041        
    GCF_000742135.1   GCF_000742135.1   100.0000   100.000   100.000   5          5545784   5          5545784        
    GCF_900638025.1   GCF_900638025.1   100.0000   100.000   100.000   1          2062405   1          2062405        
    GCF_003697165.2   GCF_003697165.2   100.0000   100.000   100.000   2          5034834   2          5034834        
    GCF_009759685.1   GCF_009759685.1   100.0000   100.000   100.000   2          3990388   2          3990388        
    GCF_000017205.1   GCF_000017205.1   100.0000   100.000   100.000   1          6588339   1          6588339        
    GCF_002950215.1   GCF_002950215.1   100.0000   100.000   100.000   3          4938295   3          4938295        
    GCF_002949675.1   GCF_002949675.1   100.0000   100.000   100.000   2          4578459   2          4578459        
    
Pairwise comparison

    $ lexicmap genome compare refs/GCF_003697165.2.fa.gz refs/GCF_002949675.1.fa.gz refs/GCF_002950215.1.fa.gz --quiet \
        | csvtk pretty -t

    genome1           genome2           tANI      ANI1      ANI2      AF1      AF2      ctgs1   size1     ctgs2   size2  
    ---------------   ---------------   -------   -------   -------   ------   ------   -----   -------   -----   -------
    GCF_003697165.2   GCF_002950215.1   73.4098   96.7220   96.8348   73.252   77.978   2       5034834   3       4938295
    GCF_003697165.2   GCF_002949675.1   71.8008   96.8012   96.5448   67.613   81.014   2       5034834   2       4578459
    GCF_002949675.1   GCF_002950215.1   77.0516   97.0923   97.6772   82.972   74.934   2       4578459   3       4938295

