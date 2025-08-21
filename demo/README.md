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

    $ lexicmap index -I refs -O demo.lmi
    20:59:26.445 [INFO] LexicMap v0.7.0
    20:59:26.445 [INFO]   https://github.com/shenwei356/LexicMap
    20:59:26.445 [INFO] 
    20:59:26.445 [INFO] checking input files ...
    20:59:26.445 [INFO]   scanning files from directory: refs
    20:59:26.445 [INFO]   15 input file(s) given
    20:59:26.445 [INFO] 
    20:59:26.445 [INFO] --------------------- [ main parameters ] ---------------------
    20:59:26.445 [INFO] 
    20:59:26.445 [INFO] input and output:
    20:59:26.445 [INFO]   input directory: refs
    20:59:26.445 [INFO]     regular expression of input files: (?i)\.(f[aq](st[aq])?|fna)(\.gz|\.xz|\.zst|\.bz2)?$
    20:59:26.446 [INFO]     *regular expression for extracting reference name from file name: (?i)(.+)\.(f[aq](st[aq])?|fna)(\.gz|\.xz|\.zst|\.bz2)?$
    20:59:26.446 [INFO]     *regular expressions for filtering out sequences: []
    20:59:26.446 [INFO]   min sequence length: 31
    20:59:26.446 [INFO]   max genome size: 15000000
    20:59:26.446 [INFO]   output directory: demo.lmi
    20:59:26.446 [INFO] 
    20:59:26.446 [INFO] mask generation:
    20:59:26.446 [INFO]   k-mer size: 31
    20:59:26.446 [INFO]   number of masks: 20000
    20:59:26.446 [INFO]   rand seed: 1
    20:59:26.446 [INFO] 
    20:59:26.446 [INFO] seed data:
    20:59:26.446 [INFO]   maximum sketching desert length: 100
    20:59:26.446 [INFO]   distance of k-mers to fill deserts: 50
    20:59:26.446 [INFO]   seeds data chunks: 16
    20:59:26.446 [INFO]   seeds data indexing partitions: 4096
    20:59:26.446 [INFO] 
    20:59:26.446 [INFO] general:
    20:59:26.446 [INFO]   genome batch size: 5000
    20:59:26.446 [INFO]   threads: 16
    20:59:26.446 [INFO]   batch merge threads: 8
    20:59:26.446 [INFO] 
    20:59:26.446 [INFO] 
    20:59:26.446 [INFO] --------------------- [ generating masks ] ---------------------
    20:59:26.452 [INFO] 
    20:59:26.453 [INFO] --------------------- [ building index ] ---------------------
    20:59:26.627 [INFO] 
    20:59:26.627 [INFO]   ------------------------[ batch 1/1 ]------------------------
    20:59:26.627 [INFO]   building index for batch 1 with 15 files...
    processed files:  15 / 15 [======================================] ETA: 0s. done
    20:59:30.582 [INFO]   writing seeds...
    20:59:30.688 [INFO]   finished writing seeds in 105.899453ms
    20:59:30.688 [INFO]   finished building index for batch 1 in: 4.061025844s
    20:59:30.688 [INFO] 
    20:59:30.688 [INFO] finished building LexicMap index from 15 files with 20000 masks in 4.254411588s
    20:59:30.689 [INFO] LexicMap index saved: demo.lmi
    20:59:30.689 [INFO] 
    20:59:30.689 [INFO] elapsed time: 4.254439661s
    20:59:30.689 [INFO]

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
    demo.lmi/: 78.36 MiB (82,165,269)
     65.26 MiB      seeds
     12.94 MiB      genomes
    156.28 KiB      masks.bin
         600 B      info.toml
         375 B      genomes.map.bin
           0 B      genomes.chunks.bin

## Searching

### A 16S rRNA gene sequence

    $ lexicmap search -d demo.lmi/  q.gene.fasta -o q.gene.fasta.lexicmap.tsv
    21:00:47.648 [INFO] LexicMap v0.7.0
    21:00:47.648 [INFO]   https://github.com/shenwei356/LexicMap
    21:00:47.648 [INFO] 
    21:00:47.648 [INFO] checking input files ...
    21:00:47.649 [INFO]   1 input file given: q.gene.fasta
    21:00:47.649 [INFO] 
    21:00:47.649 [INFO] loading index: demo.lmi/
    21:00:47.649 [INFO]   reading masks...
    21:00:47.651 [INFO]   reading indexes of seeds (k-mer-value) data...
    21:00:48.568 [INFO]   creating reader pools for 1 genome batches, each with 16 readers...
    21:00:48.568 [INFO] index loaded in 919.499738ms
    21:00:48.568 [INFO] 
    21:00:48.568 [INFO] searching with 16 threads...

    21:00:48.607 [INFO] 
    21:00:48.607 [INFO] processed queries: 1, speed: 1537.363 queries per minute
    21:00:48.607 [INFO] 100.0000% (1/1) queries matched
    21:00:48.607 [INFO] done searching
    21:00:48.607 [INFO] search results saved to: q.gene.fasta.lexicmap.tsv
    21:00:48.607 [INFO] 
    21:00:48.607 [INFO] elapsed time: 958.779224ms
    21:00:48.607 [INFO]

Result preview.
Here we create a `species` column from the genome ID column (`sgenome`) and replace the assemby accessions with species names.

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
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   2     2     100.000   1542      99.546   1      1        1542   2125377   2126917   -      4659463   0.00e+00   2758       Shigella flexneri   
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   3     3     100.000   1542      99.611   0      1        1542   3119331   3120872   +      4659463   0.00e+00   2755       Shigella flexneri   
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   4     4     100.000   1542      99.611   0      1        1542   3396068   3397609   +      4659463   0.00e+00   2755       Shigella flexneri   
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   5     5     100.000   1542      99.546   0      1        1542   3355632   3357173   +      4659463   0.00e+00   2749       Shigella flexneri   
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   6     6     100.000   1542      99.546   0      1        1542   4223146   4224687   +      4659463   0.00e+00   2749       Shigella flexneri   
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
    NC_001895.1   33593   2      GCF_003697165.2   NZ_CP033092.2   77.588    1     1     27.890    9371      97.716   2      1        9369    1864411   1873781   +      4903501   0.00e+00    15953      Escherichia coli    
    NC_001895.1   33593   2      GCF_003697165.2   NZ_CP033092.2   77.588    1     2     0.301     101       98.020   0      10308    10408   1873846   1873946   +      4903501   1.72e-43    174        Escherichia coli    
    NC_001895.1   33593   2      GCF_003697165.2   NZ_CP033092.2   77.588    2     3     20.665    6942      96.528   4      17441    24382   1882011   1888948   +      4903501   0.00e+00    11459      Escherichia coli    
    NC_001895.1   33593   2      GCF_003697165.2   NZ_CP033092.2   77.588    3     4     17.685    5941      97.980   0      24355    30295   1853098   1859038   +      4903501   0.00e+00    10174      Escherichia coli    
    NC_001895.1   33593   2      GCF_003697165.2   NZ_CP033092.2   77.588    4     5     8.993     3021      91.526   0      10308    13328   1873846   1876866   +      4903501   0.00e+00    4295       Escherichia coli    
    NC_001895.1   33593   2      GCF_003697165.2   NZ_CP033092.2   77.588    5     6     2.438     820       84.390   1      14540    15358   1878798   1879617   +      4903501   1.29e-264   911        Escherichia coli    
    NC_001895.1   33593   2      GCF_002949675.1   NZ_CP026774.1   0.976     1     1     0.976     331       85.801   3      13919    14246   3704319   3704649   -      4395762   6.35e-112   403        Shigella dysenteriae


### Simulated Oxford Nanopore R10.4.1 long-reads

Here we use the flag `-w/--load-whole-seeds` to accelerate searching.

    $ lexicmap search -d demo.lmi/ q.long-reads.fasta.gz -o q.long-reads.fasta.gz.lexicmap.tsv.gz -w -q 70
    09:36:23.956 [INFO] LexicMap v0.7.0
    09:36:23.957 [INFO]   https://github.com/shenwei356/LexicMap
    09:36:23.957 [INFO] 
    09:36:23.957 [INFO] checking input files ...
    09:36:23.957 [INFO]   1 input file given: q.long-reads.fasta.gz
    09:36:23.957 [INFO] 
    09:36:23.957 [INFO] loading index: demo.lmi/
    09:36:23.966 [INFO]   reading masks...
    09:36:23.969 [INFO]   reading seeds (k-mer-value) data into memory...
    09:36:24.539 [INFO]   creating reader pools for 1 genome batches, each with 16 readers...
    09:36:24.539 [INFO] index loaded in 583.016377ms
    09:36:24.539 [INFO] 
    09:36:24.539 [INFO] searching with 16 threads...
    processed queries: 3584, speed: 2180.560 queries per minute
    09:38:09.530 [INFO] 
    09:38:09.530 [INFO] processed queries: 3692, speed: 2109.902 queries per minute
    09:38:09.530 [INFO] 76.3543% (2819/3692) queries matched
    09:38:09.530 [INFO] done searching
    09:38:09.530 [INFO] search results saved to: q.long-reads.fasta.gz.lexicmap.tsv.gz
    09:38:09.544 [INFO] 
    09:38:09.544 [INFO] elapsed time: 1m45.587678478s
    09:38:09.544 [INFO] 

Result overview:

    $ csvtk head -n 10 q.long-reads.fasta.gz.lexicmap.tsv.gz \
        | csvtk mutate -t -n species -f sgenome \
        | csvtk replace -t -f species -k ass2species.map -p '(.+)' -r '{kv}' \
        | csvtk pretty -t

    query                  qlen    hits   sgenome           sseqid              qcovGnm   cls   hsp   qcovHSP   alenHSP   pident   gaps   qstart   qend    sstart    send      sstr   slen      evalue     bitscore   species                   
    --------------------   -----   ----   ---------------   -----------------   -------   ---   ---   -------   -------   ------   ----   ------   -----   -------   -------   ----   -------   --------   --------   --------------------------
    GCF_003697165.2_r46    2169    1      GCF_003697165.2   NZ_CP033092.2       93.177    1     1     93.177    2101      90.243   111    12       2032    4489774   4491843   +      4903501   0.00e+00   3699       Escherichia coli          
    GCF_900638025.1_r28    6375    1      GCF_900638025.1   NZ_LR134481.1       99.953    1     1     99.953    6542      92.831   261    3        6374    137490    143940    -      2062405   0.00e+00   11798      Haemophilus parainfluenzae
    GCF_000006945.2_r8     7258    1      GCF_000006945.2   NC_003197.2         99.724    1     1     99.724    7301      97.630   91     20       7257    4618964   4626236   +      4857450   0.00e+00   13135      Salmonella enterica       
    GCF_000006945.2_r109   3788    2      GCF_000006945.2   NC_003197.2         99.393    1     1     99.393    3799      97.026   67     10       3774    4633318   4637083   -      4857450   0.00e+00   6900       Salmonella enterica       
    GCF_000006945.2_r109   3788    2      GCF_002949675.1   NZ_CP026774.1       74.815    1     1     74.815    2871      81.435   64     852      3685    2156464   2159307   +      4395762   0.00e+00   3276       Shigella dysenteriae      
    GCF_001544255.1_r110   9910    1      GCF_001544255.1   NZ_BCQD01000005.1   99.839    1     1     99.839    9983      97.666   131    17       9910    155488    165428    +      191690    0.00e+00   18049      Enterococcus faecium      
    GCF_009759685.1_r164   3132    1      GCF_009759685.1   NZ_CP046654.1       99.042    1     1     99.042    3152      94.670   86     20       3121    1768740   1771855   +      3980848   0.00e+00   5586       Acinetobacter baumannii   
    GCF_000017205.1_r183   14521   1      GCF_000017205.1   NC_009656.1         99.910    1     1     99.910    14685     96.765   269    14       14521   3874715   3889307   +      6588339   0.00e+00   26541      Pseudomonas aeruginosa    
    GCF_000742135.1_r146   19632   1      GCF_000742135.1   NZ_KN046818.1       93.730    1     1     93.730    18695     95.844   412    18       18418   3816815   3835391   +      5284261   0.00e+00   33525      Klebsiella pneumoniae     
    GCF_001027105.1_r148   20294   1      GCF_001027105.1   NZ_CP011526.1       99.921    1     1     99.921    20481     97.100   307    16       20293   2352020   2372396   +      2755072   0.00e+00   36818      Staphylococcus aureus

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
 Score = 279 bits, Expect = 9.66e-75
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
