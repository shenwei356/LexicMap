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
    13:01:50.205 [INFO] LexicMap v0.10.1 (3f499bb)
    13:01:50.205 [INFO]   https://github.com/shenwei356/LexicMap
    13:01:50.205 [INFO] 
    13:01:50.205 [INFO]  CWD: /home/shenwei/go/src/github.com/shenwei356/LexicMap/demo
    13:01:50.205 [INFO]  CMD: lexicmap index -I refs -O demo.lmi
    13:01:50.205 [INFO] DATE: 2026-11-06
    13:01:50.205 [INFO] 
    13:01:50.206 [INFO] checking input files ...
    13:01:50.206 [INFO]   scanning files from directory: refs
    13:01:50.206 [INFO]   15 input file(s) given
    13:01:50.206 [INFO] 
    13:01:50.206 [INFO] --------------------- [ main parameters ] ---------------------
    13:01:50.206 [INFO] 
    13:01:50.206 [INFO] input and output:
    13:01:50.206 [INFO]   input directory: refs
    13:01:50.206 [INFO]     regular expression of input files: (?i)\.(f[aq](st[aq])?|fna)(\.gz|\.xz|\.zst|\.bz2)?$
    13:01:50.206 [INFO]     *regular expression for extracting reference name from file name: (?i)(.+)\.(f[aq](st[aq])?|fna)(\.gz|\.xz|\.zst|\.bz2)?$
    13:01:50.206 [INFO]     *regular expressions for filtering out sequences: []
    13:01:50.206 [INFO]   min sequence length: -1
    13:01:50.206 [INFO]   max genome size: 20000000
    13:01:50.206 [INFO]   output directory: demo.lmi
    13:01:50.206 [INFO] 
    13:01:50.206 [INFO] mask generation:
    13:01:50.206 [INFO]   k-mer size: 31
    13:01:50.206 [INFO]   number of masks: 20000
    13:01:50.206 [INFO]   rand seed: 1
    13:01:50.206 [INFO] 
    13:01:50.206 [INFO] seed data:
    13:01:50.206 [INFO]   maximum sketching desert length: 100
    13:01:50.206 [INFO]   distance of k-mers to fill deserts: 50
    13:01:50.206 [INFO]   seeds data chunks: 16
    13:01:50.206 [INFO]   seeds data indexing partitions: 4096
    13:01:50.206 [INFO] 
    13:01:50.206 [INFO] general:
    13:01:50.206 [INFO]   genome batch size: 5000
    13:01:50.206 [INFO]   threads: 16
    13:01:50.206 [INFO]   batch merge threads: 8
    13:01:50.206 [INFO] 
    13:01:50.206 [INFO] 
    13:01:50.206 [INFO] --------------------- [ generating masks ] ---------------------
    13:01:50.211 [INFO] 
    13:01:50.211 [INFO] --------------------- [ building index ] ---------------------
    13:01:50.524 [INFO] 
    13:01:50.524 [INFO]   ------------------------[ batch 1/1 ]------------------------
    13:01:50.524 [INFO]   building index for batch 1 with 15 files...
    processed files:  15 / 15 [======================================] done in 1s
    13:01:51.907 [INFO]   writing seeds...
    13:01:52.054 [INFO]   finished writing seeds in 146.189363ms
    13:01:52.054 [INFO]   finished building index for batch 1 in: 1.52931447s
    13:01:52.054 [INFO] 
    13:01:52.054 [INFO] finished building LexicMap index from 15 files with 20000 masks in 1.848668476s
    13:01:52.054 [INFO] LexicMap index saved: demo.lmi
    13:01:52.054 [INFO] 
    13:01:52.054 [INFO] elapsed time: 1.848704724s
    13:01:52.054 [INFO] 
   
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
    demo.lmi/: 78.17 MiB (81,972,177)
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
    13:04:16.349 [INFO] LexicMap v0.10.1 (3f499bb)
    13:04:16.349 [INFO]   https://github.com/shenwei356/LexicMap
    13:04:16.349 [INFO] 
    13:04:16.349 [INFO]  CWD: /home/shenwei/go/src/github.com/shenwei356/LexicMap/demo
    13:04:16.349 [INFO]  CMD: lexicmap utils reindex-seeds2 -d demo.lmi/
    13:04:16.349 [INFO] DATE: 2026-11-06
    13:04:16.349 [INFO] 
    13:04:16.349 [INFO] recreating adaptive seed indexes with 4096 partitions and a 8.0 KiB block threshold for: demo.lmi
    processed files:  16 / 16 [======================================] done in 0s
    13:04:16.430 [INFO]   idx15 stats:primary-prefix blocks: 81,920,000 total, 2,106,999 non-empty, 0 indexed (0.00% of non-empty)
    13:04:16.430 [INFO]   idx15 size: 768 B
    13:04:16.430 [INFO] 
    13:04:16.430 [INFO] elapsed time: 81.420107ms
    13:04:16.430 [INFO] 

## Searching

### A 16S rRNA gene sequence

    $ lexicmap search -d demo.lmi/  q.gene.fasta -o q.gene.fasta.lexicmap.tsv
    13:04:37.207 [INFO] LexicMap v0.10.1 (3f499bb)
    13:04:37.207 [INFO]   https://github.com/shenwei356/LexicMap
    13:04:37.207 [INFO] 
    13:04:37.207 [INFO]  CWD: /home/shenwei/go/src/github.com/shenwei356/LexicMap/demo
    13:04:37.207 [INFO]  CMD: lexicmap search -d demo.lmi/ q.gene.fasta -o q.gene.fasta.lexicmap.tsv
    13:04:37.207 [INFO] DATE: 2026-11-06
    13:04:37.207 [INFO] 
    13:04:37.207 [INFO] checking input files ...
    13:04:37.207 [INFO]   1 input file given: q.gene.fasta
    13:04:37.207 [INFO] 
    13:04:37.207 [INFO] loading index: demo.lmi
    13:04:37.207 [INFO]   reading masks...
    13:04:37.208 [INFO]   reading indexes of seeds (k-mer-value) data...
    13:04:37.208 [INFO]   creating searcher pools for 16 seed data files, each with 4 searchers...
    13:04:37.549 [INFO]   creating reader pools for 1 genome batches, each with 16 readers...
    13:04:37.549 [INFO] index loaded in 342.537593ms
    13:04:37.550 [INFO] 
    13:04:37.550 [INFO] searching with 16 threads...
    13:04:37.550 [INFO]   maximum number of concurrent queries: 8, force garbage collection for every 64 queries
    13:04:37.550 [INFO]   lazy genome metadata enabled, contig metadata will be read only after pseudoalignment finds a match

    13:04:37.572 [INFO] 
    13:04:37.572 [INFO] processed queries: 1, speed: 2641.147 queries per minute
    13:04:37.572 [INFO] 100.0000% (1/1) queries matched
    13:04:37.572 [INFO] done searching
    13:04:37.572 [INFO] search results saved to: q.gene.fasta.lexicmap.tsv
    13:04:37.573 [INFO] 
    13:04:37.573 [INFO] elapsed time: 365.925336ms
    13:04:37.573 [INFO]

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
    NC_001895.1   33593   3      GCF_003697165.2   NZ_CP033092.2   77.686    1     1     27.890    9371      97.716   2      1        9369    1864411   1873781   +      4903501   0.00e+00    15928      Escherichia coli    
    NC_001895.1   33593   3      GCF_003697165.2   NZ_CP033092.2   77.686    1     2     0.298     100       98.000   0      10308    10407   1873846   1873945   +      4903501   5.94e-43    172        Escherichia coli    
    NC_001895.1   33593   3      GCF_003697165.2   NZ_CP033092.2   77.686    2     3     20.665    6942      96.528   4      17441    24382   1882011   1888948   +      4903501   0.00e+00    11427      Escherichia coli    
    NC_001895.1   33593   3      GCF_003697165.2   NZ_CP033092.2   77.686    3     4     17.685    5941      97.980   0      24355    30295   1853098   1859038   +      4903501   0.00e+00    10174      Escherichia coli    
    NC_001895.1   33593   3      GCF_003697165.2   NZ_CP033092.2   77.686    4     5     8.880     2983      91.720   0      10308    13290   1873846   1876828   +      4903501   0.00e+00    4266       Escherichia coli    
    NC_001895.1   33593   3      GCF_003697165.2   NZ_CP033092.2   77.686    5     6     0.211     71        84.507   0      14213    14283   1878639   1878709   +      4903501   7.15e-15    78         Escherichia coli    
    NC_001895.1   33593   3      GCF_003697165.2   NZ_CP033092.2   77.686    5     7     2.438     820       84.390   1      14540    15358   1878798   1879617   +      4903501   8.15e-261   899        Escherichia coli    
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
            
        13:06:17.741 [INFO] LexicMap v0.10.1 (3f499bb)
        13:06:17.741 [INFO]   https://github.com/shenwei356/LexicMap
        13:06:17.741 [INFO] 
        13:06:17.742 [INFO]  CWD: /home/shenwei/go/src/github.com/shenwei356/LexicMap/demo
        13:06:17.742 [INFO]  CMD: lexicmap search -d demo.lmi/ q.long-reads.fasta.gz -o q.long-reads.fasta.gz.lexicmap.tsv.gz --min-qcov-per-hsp 70 --top-n-genomes 5 --top-n-chains 1
        13:06:17.742 [INFO] DATE: 2026-11-06
        13:06:17.742 [INFO] 
        13:06:17.742 [INFO] checking input files ...
        13:06:17.742 [INFO]   1 input file given: q.long-reads.fasta.gz
        13:06:17.742 [INFO] 
        13:06:17.742 [INFO] loading index: demo.lmi
        13:06:17.742 [INFO]   reading masks...
        13:06:17.743 [INFO]   reading indexes of seeds (k-mer-value) data...
        13:06:17.743 [INFO]   creating searcher pools for 16 seed data files, each with 4 searchers...
        13:06:18.054 [INFO]   creating reader pools for 1 genome batches, each with 16 readers...
        13:06:18.054 [INFO] index loaded in 312.810702ms
        13:06:18.055 [INFO] 
        13:06:18.055 [INFO] searching with 16 threads...
        13:06:18.055 [INFO]   maximum number of concurrent queries: 8, force garbage collection for every 64 queries
        13:06:18.055 [INFO]   keep the top 5 genomes
        13:06:18.055 [INFO]   keep the top 1 chains
        13:06:18.055 [INFO]   lazy genome metadata enabled, contig metadata will be read only after pseudoalignment finds a match
        processed queries: 3456, speed: 11351.217 queries per minute
        13:06:37.555 [INFO] 
        13:06:37.555 [INFO] processed queries: 3692, speed: 11359.810 queries per minute
        13:06:37.555 [INFO] 76.2459% (2815/3692) queries matched
        13:06:37.555 [INFO] done searching
        13:06:37.555 [INFO] search results saved to: q.long-reads.fasta.gz.lexicmap.tsv.gz
        13:06:37.569 [INFO] 
        13:06:37.569 [INFO] elapsed time: 19.827892154s
        13:06:37.569 [INFO] 


Result overview:

    $ csvtk head -n 10 q.long-reads.fasta.gz.lexicmap.tsv.gz \
        | csvtk mutate -t -n species -f sgenome \
        | csvtk replace -t -f species -k ass2species.map -p '(.+)' -r '{kv}' \
        | csvtk pretty -t

    query                  qlen    hits   sgenome           sseqid              qcovGnm   cls   hsp   qcovHSP   alenHSP   pident   gaps   qstart   qend    sstart    send      sstr   slen      evalue     bitscore   species                   
    --------------------   -----   ----   ---------------   -----------------   -------   ---   ---   -------   -------   ------   ----   ------   -----   -------   -------   ----   -------   --------   --------   --------------------------
    GCF_900638025.1_r28    6375    1      GCF_900638025.1   NZ_LR134481.1       99.953    1     1     99.953    6542      92.831   261    3        6374    137490    143940    -      2062405   0.00e+00   8982       Haemophilus parainfluenzae
    GCF_002949675.1_r183   8176    1      GCF_002949675.1   NZ_CP026774.1       99.682    1     1     99.682    8188      98.107   82     27       8176    4194298   4202441   +      4395762   0.00e+00   13876      Shigella dysenteriae      
    GCF_003697165.2_r46    2169    1      GCF_003697165.2   NZ_CP033092.2       93.177    1     1     93.177    2101      90.243   111    12       2032    4489774   4491843   +      4903501   0.00e+00   2632       Escherichia coli          
    GCF_000006945.2_r8     7258    1      GCF_000006945.2   NC_003197.2         99.724    1     1     99.724    7301      97.630   91     20       7257    4618964   4626236   +      4857450   0.00e+00   12130      Salmonella enterica       
    GCF_000017205.1_r183   14521   1      GCF_000017205.1   NC_009656.1         99.910    1     1     99.910    14685     96.765   269    14       14521   3874715   3889307   +      6588339   0.00e+00   23596      Pseudomonas aeruginosa    
    GCF_009759685.1_r164   3132    1      GCF_009759685.1   NZ_CP046654.1       99.042    1     1     99.042    3152      94.670   86     20       3121    1768740   1771855   +      3980848   0.00e+00   4735       Acinetobacter baumannii   
    GCF_001027105.1_r148   20294   1      GCF_001027105.1   NZ_CP011526.1       99.921    1     1     99.921    20481     97.100   307    16       20293   2352020   2372396   +      2755072   0.00e+00   33358      Staphylococcus aureus     
    GCF_001544255.1_r110   9910    1      GCF_001544255.1   NZ_BCQD01000005.1   99.839    1     1     99.839    9983      97.666   131    17       9910    155488    165428    +      191690    0.00e+00   16567      Enterococcus faecium      
    GCF_000006945.2_r109   3788    2      GCF_000006945.2   NC_003197.2         99.393    1     1     99.393    3799      97.026   67     10       3774    4633318   4637083   -      4857450   0.00e+00   6145       Salmonella enterica       
    GCF_000006945.2_r109   3788    2      GCF_002950215.1   NZ_CP026788.1       70.301    1     1     70.301    2700      82.296   60     840      3502    2488829   2491505   -      4659463   0.00e+00   2565       Shigella flexneri       

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
    15:00:46.027 [INFO] LexicMap v0.10.1 (aba7a7c)
    15:00:46.027 [INFO]   https://github.com/shenwei356/LexicMap
    15:00:46.027 [INFO]
    15:00:46.027 [INFO]  CWD: /home/shenwei/go/src/github.com/shenwei356/LexicMap/demo
    15:00:46.027 [INFO]  CMD: lexicmap genome search -d demo.lmi/ refs/GCF_003697165.2.fa.gz
    15:00:46.027 [INFO] DATE: 2026-11-06
    15:00:46.027 [INFO]
    15:00:46.027 [INFO] checking input files ...
    15:00:46.027 [INFO]   1 input file given: refs/GCF_003697165.2.fa.gz
    15:00:46.027 [INFO]
    15:00:46.027 [INFO] loading index: demo.lmi
    15:00:46.027 [INFO]   reading masks...
    15:00:46.029 [INFO]   reading indexes of seeds (k-mer-value) data...
    15:00:46.029 [INFO]   creating searcher pools for 16 seed data files, each with 4 searchers...
    15:00:46.125 [INFO]   creating reader pools for 1 genome batches, each with 16 readers...
    15:00:46.125 [INFO] index loaded in 98.050544ms
    15:00:46.126 [INFO]
    15:00:46.126 [INFO] searching with 16 threads...
    15:00:46.126 [INFO]   screen genomes with seed data of 1024 masks
    15:00:46.126 [INFO]     keep the top 10 screened genomes, including all ties at the cutoff score
    15:00:46.126 [INFO]   minimum query coverage per HSP: 30.00%
    15:00:46.126 [INFO]   minimum base identity in a HSP segment: 70.00%
    15:00:46.126 [INFO]   maximum evalue: 1.00e-15
    15:00:46.126 [INFO]   maximum number of concurrent queries: 8, force garbage collection for every 4 queries
    15:00:46.126 [INFO]   threads per query: 16

    15:00:50.679 [INFO]
    15:00:50.679 [INFO] processed queries: 1, speed: 13.178 queries per minute
    15:00:50.679 [INFO] 100.0000% (1/1) queries matched
    15:00:50.679 [INFO] done searching
    15:00:50.679 [INFO]
    15:00:50.679 [INFO] elapsed time: 4.651703937s
    15:00:50.679 [INFO]

    query             subject           ANI        qAF       sAF       qcontigs   qsize     scontigs   ssize     sname
    ---------------   ---------------   --------   -------   -------   --------   -------   --------   -------   -----
    GCF_003697165.2   GCF_003697165.2   100.0000   100.000   100.000   2          5034834   2          5034834        
    GCF_003697165.2   GCF_002949675.1   96.7896    67.632    74.374    2          5034834   2          4578459        
    GCF_003697165.2   GCF_002950215.1   96.7074    73.272    74.704    2          5034834   3          4938295        
    GCF_003697165.2   GCF_000006945.2   81.5408    55.812    56.753    2          5034834   2          4951383        
    GCF_003697165.2   GCF_000742135.1   80.4634    46.981    42.652    2          5034834   5          5545784        
    
    # with species name
    $ lexicmap genome search -d demo.lmi/ refs/GCF_003697165.2.fa.gz \
            --taxdump ~/.taxonkit/ --genome2taxid taxid.map --show-species-name --quiet \
        | csvtk pretty -t
        
    query             subject           ANI        qAF       sAF       qcontigs   qsize     scontigs   ssize     sname                
    ---------------   ---------------   --------   -------   -------   --------   -------   --------   -------   ---------------------
    GCF_003697165.2   GCF_003697165.2   100.0000   100.000   100.000   2          5034834   2          5034834   Escherichia coli     
    GCF_003697165.2   GCF_002949675.1   96.7896    67.632    74.374    2          5034834   2          4578459   Shigella dysenteriae 
    GCF_003697165.2   GCF_002950215.1   96.7074    73.272    74.704    2          5034834   3          4938295   Shigella flexneri    
    GCF_003697165.2   GCF_000006945.2   81.5408    55.812    56.753    2          5034834   2          4951383   Salmonella enterica  
    GCF_003697165.2   GCF_000742135.1   80.4634    46.981    42.652    2          5034834   5          5545784   Klebsiella pneumoniae


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
    GCF_003697165.2   GCF_002950215.1   73.4218   96.7074   96.8315   73.272   77.989   2       5034834   3       4938295
    GCF_003697165.2   GCF_002949675.1   71.8041   96.7896   96.5454   67.632   81.004   2       5034834   2       4578459
    GCF_002949675.1   GCF_002950215.1   77.0694   97.0916   97.6576   82.979   74.968   2       4578459   3       4938295

