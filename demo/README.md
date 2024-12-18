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

    19:13:55.369 [INFO] LexicMap v0.5.0 (c66e3a0)
    19:13:55.369 [INFO]   https://github.com/shenwei356/LexicMap
    19:13:55.369 [INFO] 
    19:13:55.369 [INFO] checking input files ...
    19:13:55.369 [INFO]   scanning files from directory: refs/
    19:13:55.370 [INFO]   15 input file(s) given
    19:13:55.370 [INFO] 
    19:13:55.370 [INFO] --------------------- [ main parameters ] ---------------------
    19:13:55.370 [INFO] 
    19:13:55.370 [INFO] input and output:
    19:13:55.370 [INFO]   input directory: refs/
    19:13:55.370 [INFO]     regular expression of input files: (?i)\.(f[aq](st[aq])?|fna)(\.gz|\.xz|\.zst|\.bz2)?$
    19:13:55.370 [INFO]     *regular expression for extracting reference name from file name: (?i)(.+)\.(f[aq](st[aq])?|fna)(\.gz|\.xz|\.zst|\.bz2)?$
    19:13:55.370 [INFO]     *regular expressions for filtering out sequences: []
    19:13:55.370 [INFO]   min sequence length: 31
    19:13:55.370 [INFO]   max genome size: 15000000
    19:13:55.370 [INFO]   output directory: demo.lmi
    19:13:55.370 [INFO] 
    19:13:55.370 [INFO] mask generation:
    19:13:55.370 [INFO]   k-mer size: 31
    19:13:55.370 [INFO]   number of masks: 40000
    19:13:55.370 [INFO]   rand seed: 1
    19:13:55.370 [INFO] 
    19:13:55.370 [INFO] seed data:
    19:13:55.370 [INFO]   maximum sketching desert length: 200
    19:13:55.370 [INFO]   distance of k-mers to fill deserts: 50
    19:13:55.370 [INFO]   seeds data chunks: 16
    19:13:55.370 [INFO]   seeds data indexing partitions: 4096
    19:13:55.370 [INFO] 
    19:13:55.370 [INFO] general:
    19:13:55.370 [INFO]   genome batch size: 5000
    19:13:55.370 [INFO]   threads: 16
    19:13:55.370 [INFO]   batch merge threads: 8
    19:13:55.370 [INFO] 
    19:13:55.370 [INFO] 
    19:13:55.370 [INFO] --------------------- [ generating masks ] ---------------------
    19:13:55.382 [INFO] 
    19:13:55.382 [INFO] --------------------- [ building index ] ---------------------
    19:13:56.018 [INFO] 
    19:13:56.018 [INFO]   ------------------------[ batch 1/1 ]------------------------
    19:13:56.018 [INFO]   building index for batch 1 with 15 files...
    processed files:  15 / 15 [======================================] ETA: 0s. done
    19:14:00.601 [INFO]   writing seeds...
    19:14:00.745 [INFO]   finished writing seeds in 143.225662ms
    19:14:00.745 [INFO]   finished building index for batch 1 in: 4.72683742s
    19:14:00.746 [INFO] 
    19:14:00.746 [INFO] finished building LexicMap index from 15 files with 40000 masks in 5.392303552s
    19:14:00.746 [INFO] LexicMap index saved: demo.lmi
    19:14:00.746 [INFO] 
    19:14:00.746 [INFO] elapsed time: 5.392329816s
    19:14:00.746 [INFO]

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
    demo.lmi/: 73.24 MiB (76,801,302)
     60.00 MiB      seeds
     12.94 MiB      genomes
    312.53 KiB      masks.bin
         563 B      info.toml
         375 B      genomes.map.bin
           0 B      genomes.chunks.bin

## Searching

### A 16S rRNA gene sequence

    $ lexicmap search -d demo.lmi/  q.gene.fasta -o q.gene.fasta.lexicmap.tsv
    19:16:55.757 [INFO] LexicMap v0.5.0 (c66e3a0)
    19:16:55.757 [INFO]   https://github.com/shenwei356/LexicMap
    19:16:55.757 [INFO] 
    19:16:55.757 [INFO] checking input files ...
    19:16:55.757 [INFO]   1 input file given: q.gene.fasta
    19:16:55.757 [INFO] 
    19:16:55.757 [INFO] loading index: demo.lmi/
    19:16:55.758 [INFO]   reading masks...
    19:16:55.762 [INFO]   reading indexes of seeds (k-mer-value) data...
    19:16:58.781 [INFO]   creating genome reader pools, each batch with 16 readers...
    19:16:58.781 [INFO] index loaded in 3.023370768s
    19:16:58.781 [INFO] 
    19:16:58.781 [INFO] searching with 16 threads...

    19:16:58.821 [INFO] 
    19:16:58.821 [INFO] processed queries: 1, speed: 1506.171 queries per minute
    19:16:58.821 [INFO] 100.0000% (1/1) queries matched
    19:16:58.821 [INFO] done searching
    19:16:58.821 [INFO] search results saved to: q.gene.fasta.lexicmap.tsv
    19:16:58.821 [INFO] 
    19:16:58.821 [INFO] elapsed time: 3.063458635s
    19:16:58.821 [INFO]

Result preview.
Here we create a `species` column from the genome ID column (`sgenome`) and replace the assemby accessions with species names.

    $ csvtk head -n 21 q.gene.fasta.lexicmap.tsv \
        | csvtk mutate -t -n species -f sgenome \
        | csvtk replace -t -f species -k ass2species.map -p '(.+)' -r '{kv}' \
        | csvtk pretty -t

    query                         qlen   hits   sgenome           sseqid          qcovGnm   hsp   qcovHSP   alenHSP   pident   gaps   qstart   qend   sstart    send      sstr   slen      species
    ---------------------------   ----   ----   ---------------   -------------   -------   ---   -------   -------   ------   ----   ------   ----   -------   -------   ----   -------   --------------------
    NC_000913.3:4166659-4168200   1542   15     GCF_003697165.2   NZ_CP033092.2   100.000   1     100.000   1542      99.805   0      1        1542   458559    460100    +      4903501   Escherichia coli
    NC_000913.3:4166659-4168200   1542   15     GCF_003697165.2   NZ_CP033092.2   100.000   2     100.000   1542      99.805   0      1        1542   1285123   1286664   +      4903501   Escherichia coli
    NC_000913.3:4166659-4168200   1542   15     GCF_003697165.2   NZ_CP033092.2   100.000   3     100.000   1542      99.805   0      1        1542   3780640   3782181   -      4903501   Escherichia coli
    NC_000913.3:4166659-4168200   1542   15     GCF_003697165.2   NZ_CP033092.2   100.000   4     100.000   1542      99.805   0      1        1542   4551515   4553056   -      4903501   Escherichia coli
    NC_000913.3:4166659-4168200   1542   15     GCF_003697165.2   NZ_CP033092.2   100.000   5     100.000   1542      99.805   0      1        1542   4591684   4593225   -      4903501   Escherichia coli
    NC_000913.3:4166659-4168200   1542   15     GCF_003697165.2   NZ_CP033092.2   100.000   6     100.000   1542      99.805   0      1        1542   4726193   4727734   -      4903501   Escherichia coli
    NC_000913.3:4166659-4168200   1542   15     GCF_003697165.2   NZ_CP033092.2   100.000   7     100.000   1542      99.805   0      1        1542   4844587   4846128   -      4903501   Escherichia coli
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   1     100.000   1542      99.676   0      1        1542   3216505   3218046   +      4659463   Shigella flexneri
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   2     100.000   1542      99.611   0      1        1542   3396068   3397609   +      4659463   Shigella flexneri
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   3     100.000   1542      99.611   0      1        1542   3119331   3120872   +      4659463   Shigella flexneri
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   4     100.000   1542      99.546   0      1        1542   3355632   3357173   +      4659463   Shigella flexneri
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   5     100.000   1542      99.546   0      1        1542   4223146   4224687   +      4659463   Shigella flexneri
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   6     100.000   1542      99.546   1      1        1542   2125377   2126917   -      4659463   Shigella flexneri
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   7     100.000   1542      99.481   0      1        1542   3540450   3541991   +      4659463   Shigella flexneri
    NC_000913.3:4166659-4168200   1542   15     GCF_002949675.1   NZ_CP026774.1   100.000   1     100.000   1542      99.027   0      1        1542   1662010   1663551   -      4395762   Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   15     GCF_002949675.1   NZ_CP026774.1   100.000   2     100.000   1542      99.027   0      1        1542   2536624   2538165   +      4395762   Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   15     GCF_002949675.1   NZ_CP026774.1   100.000   3     100.000   1542      99.027   0      1        1542   2636477   2638018   +      4395762   Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   15     GCF_002949675.1   NZ_CP026774.1   100.000   4     100.000   1542      99.027   0      1        1542   2768883   2770424   +      4395762   Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   15     GCF_002949675.1   NZ_CP026774.1   100.000   5     100.000   1542      99.027   0      1        1542   2810845   2812386   +      4395762   Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   15     GCF_002949675.1   NZ_CP026774.1   100.000   6     100.000   1542      99.027   0      1        1542   3061592   3063133   +      4395762   Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   15     GCF_002949675.1   NZ_CP026774.1   100.000   7     100.000   1542      99.027   0      1        1542   3646778   3648319   +      4395762   Shigella dysenteriae

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

 HSP #1
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

    query         qlen    hits   sgenome           sseqid          qcovGnm   hsp   qcovHSP   alenHSP   pident    gaps   qstart   qend    sstart    send      sstr   slen      species
    -----------   -----   ----   ---------------   -------------   -------   ---   -------   -------   -------   ----   ------   -----   -------   -------   ----   -------   --------------------
    NC_001895.1   33593   2      GCF_003697165.2   NZ_CP033092.2   77.183    1     27.854    9359      97.735    2      1        9357    1864411   1873769   +      4903501   Escherichia coli
    NC_001895.1   33593   2      GCF_003697165.2   NZ_CP033092.2   77.183    2     20.570    6910      96.570    4      17473    24382   1882043   1888948   +      4903501   Escherichia coli
    NC_001895.1   33593   2      GCF_003697165.2   NZ_CP033092.2   77.183    3     17.644    5927      98.043    0      24355    30281   1853098   1859024   +      4903501   Escherichia coli
    NC_001895.1   33593   2      GCF_003697165.2   NZ_CP033092.2   77.183    4     8.844     2971      91.754    0      10308    13278   1873846   1876816   +      4903501   Escherichia coli
    NC_001895.1   33593   2      GCF_003697165.2   NZ_CP033092.2   77.183    5     0.176     59        100.000   0      9299     9357    1873711   1873769   +      4903501   Escherichia coli
    NC_001895.1   33593   2      GCF_003697165.2   NZ_CP033092.2   77.183    6     2.355     791       84.703    0      14543    15333   1878801   1879591   +      4903501   Escherichia coli
    NC_001895.1   33593   2      GCF_002949675.1   NZ_CP026774.1   0.828     1     0.828     281       87.189    3      13919    14196   3704369   3704649   -      4395762   Shigella dysenteriae

### Simulated Oxford Nanopore R10.4.1 long-reads

Here we use the flag `-w/--load-whole-seeds` to accelerate searching.

    $ lexicmap search -d demo.lmi/ q.long-reads.fasta.gz -o q.long-reads.fasta.gz.lexicmap.tsv.gz -w -q 70
    19:17:49.069 [INFO] LexicMap v0.5.0 (c66e3a0)
    19:17:49.069 [INFO]   https://github.com/shenwei356/LexicMap
    19:17:49.069 [INFO] 
    19:17:49.069 [INFO] checking input files ...
    19:17:49.069 [INFO]   1 input file given: q.long-reads.fasta.gz
    19:17:49.069 [INFO] 
    19:17:49.069 [INFO] loading index: demo.lmi/
    19:17:49.069 [INFO]   reading masks...
    19:17:49.073 [INFO]   reading seeds (k-mer-value) data into memory...
    19:17:51.324 [INFO]   creating genome reader pools, each batch with 16 readers...
    19:17:51.325 [INFO] index loaded in 2.256185788s
    19:17:51.325 [INFO] 
    19:17:51.325 [INFO] searching with 16 threads...
    processed queries: 3584, speed: 2235.509 queries per minute
    19:19:33.442 [INFO] 
    19:19:33.442 [INFO] processed queries: 3692, speed: 2169.281 queries per minute
    19:19:33.442 [INFO] 76.3543% (2819/3692) queries matched
    19:19:33.442 [INFO] done searching
    19:19:33.442 [INFO] search results saved to: q.long-reads.fasta.gz.lexicmap.tsv.gz
    19:19:33.449 [INFO] 
    19:19:33.449 [INFO] elapsed time: 1m44.380463612s
    19:19:33.449 [INFO]

Result overview:

    $ csvtk head -n 26 q.long-reads.fasta.gz.lexicmap.tsv.gz \
        | csvtk mutate -t -n species -f sgenome \
        | csvtk replace -t -f species -k ass2species.map -p '(.+)' -r '{kv}' \
        | csvtk pretty -t

    query                  qlen    hits   sgenome           sseqid              qcovGnm   hsp   qcovHSP   alenHSP   pident   gaps   qstart   qend    sstart    send      sstr   slen      species
    --------------------   -----   ----   ---------------   -----------------   -------   ---   -------   -------   ------   ----   ------   -----   -------   -------   ----   -------   --------------------------
    GCF_003697165.2_r46    2169    1      GCF_003697165.2   NZ_CP033092.2       91.886    1     91.886    2072      90.251   109    31       2023    4489794   4491835   +      4903501   Escherichia coli
    GCF_900638025.1_r28    6375    1      GCF_900638025.1   NZ_LR134481.1       99.357    1     99.357    6503      92.849   259    6        6339    137524    143936    -      2062405   Haemophilus parainfluenzae
    GCF_001544255.1_r110   9910    1      GCF_001544255.1   NZ_BCQD01000005.1   99.839    1     99.839    9983      97.666   131    17       9910    155488    165428    +      191690    Enterococcus faecium
    GCF_000006945.2_r8     7258    1      GCF_000006945.2   NC_003197.2         99.339    1     99.339    7273      97.635   90     20       7229    4618964   4626209   +      4857450   Salmonella enterica
    GCF_000006945.2_r109   3788    2      GCF_000006945.2   NC_003197.2         98.522    1     98.522    3764      97.131   63     37       3768    4633323   4637055   -      4857450   Salmonella enterica
    GCF_000006945.2_r109   3788    2      GCF_000742135.1   NZ_KN046818.1       75.422    1     75.422    2942      76.717   156    852      3708    164018    166888    +      5284261   Klebsiella pneumoniae
    GCF_009759685.1_r164   3132    1      GCF_009759685.1   NZ_CP046654.1       99.042    1     99.042    3152      94.670   86     20       3121    1768740   1771855   +      3980848   Acinetobacter baumannii
    GCF_000017205.1_r183   14521   1      GCF_000017205.1   NC_009656.1         99.787    1     99.787    14666     96.782   267    28       14517   3874730   3889304   +      6588339   Pseudomonas aeruginosa
    GCF_001027105.1_r148   20294   1      GCF_001027105.1   NZ_CP011526.1       99.921    1     99.921    20481     97.100   307    16       20293   2352020   2372396   +      2755072   Staphylococcus aureus
    GCF_000742135.1_r146   19632   1      GCF_000742135.1   NZ_KN046818.1       93.689    1     93.689    18687     95.847   412    26       18418   3816823   3835391   +      5284261   Klebsiella pneumoniae
    GCF_001027105.1_r40    28470   1      GCF_001027105.1   NZ_CP011526.1       99.919    1     99.919    28681     97.633   374    24       28470   371992    400532    -      2755072   Staphylococcus aureus
    GCF_000742135.1_r314   17829   1      GCF_000742135.1   NZ_KN046818.1       99.736    1     99.736    17875     98.473   151    23       17804   989908    1007724   -      5284261   Klebsiella pneumoniae
    GCF_000017205.1_r180   14156   1      GCF_000017205.1   NC_009656.1         99.435    1     99.435    14613     88.079   792    47       14122   109205    123562    +      6588339   Pseudomonas aeruginosa
    GCF_002949675.1_r249   1937    4      GCF_002949675.1   NZ_CP026774.1       98.554    1     98.554    1935      95.142   43     29       1937    3336601   3338518   -      4395762   Shigella dysenteriae
    GCF_002949675.1_r249   1937    4      GCF_002950215.1   NZ_CP026788.1       76.510    1     76.510    1499      93.863   28     456      1937    3959212   3960699   +      4659463   Shigella flexneri
    GCF_002949675.1_r249   1937    4      GCF_003697165.2   NZ_CP033092.2       98.554    1     98.554    1944      91.821   52     29       1937    925864    927790    +      4903501   Escherichia coli
    GCF_002949675.1_r249   1937    4      GCF_000006945.2   NC_003197.2         95.044    1     95.044    1886      76.776   89     65       1905    3221659   3223500   -      4857450   Salmonella enterica
    GCF_002949675.1_r183   8176    1      GCF_002949675.1   NZ_CP026774.1       99.682    1     99.682    8188      98.107   82     27       8176    4194298   4202441   +      4395762   Shigella dysenteriae
    GCF_009759685.1_r168   3398    1      GCF_009759685.1   NZ_CP046654.1       98.558    1     98.558    3395      95.523   81     22       3370    3276395   3279754   -      3980848   Acinetobacter baumannii
    GCF_001544255.1_r104   1087    2      GCF_001544255.1   NZ_BCQD01000030.1   95.308    1     95.308    1050      96.571   22     26       1061    9         1050      -      1061      Enterococcus faecium
    GCF_001544255.1_r104   1087    2      GCF_000392875.1   NZ_KB944589.1       95.308    1     95.308    1053      75.783   28     26       1061    649294    650335    -      682426    Enterococcus faecalis
    GCF_001544255.1_r104   1087    2      GCF_000392875.1   NZ_KB944590.1       95.308    2     95.308    1053      75.783   28     26       1061    1763370   1764411   +      1924212   Enterococcus faecalis
    GCF_000006945.2_r43    20355   1      GCF_000006945.2   NC_003197.2         74.576    1     74.576    15396     96.012   336    32       15211   3589949   3605224   +      4857450   Salmonella enterica
    GCF_000392875.1_r181   6365    1      GCF_000392875.1   NZ_KB944590.1       99.042    1     99.042    6547      88.560   359    41       6344    1335036   1341466   +      1924212   Enterococcus faecalis
    GCF_009759685.1_r69    9800    1      GCF_009759685.1   NZ_CP046654.1       99.571    1     99.571    9876      96.223   216    15       9772    121950    131727    +      3980848   Acinetobacter baumannii
    GCF_003697165.2_r248   6741    1      GCF_003697165.2   NZ_CP033092.2       99.733    1     99.733    6828      94.786   176    9        6731    2827221   2833977   -      4903501   Escherichia coli

Blast-style format:

```
# align only one long-read <= 500 bp

$ seqkit seq -g -M 500 q.long-reads.fasta.gz \
    | seqkit head -n 1 \
    | lexicmap search -d demo.lmi/ -a \
    | lexicmap utils 2blast --kv-file-genome ass2species.map

Query = GCF_006742205.1_r100
Length = 431

[Subject genome #1/1] = GCF_006742205.1 Staphylococcus epidermidis
Query coverage per genome = 92.575%

>NZ_AP019721.1
Length = 2422602

 HSP #1
 Query coverage per seq = 92.575%, Aligned length = 402, Identities = 98.507%, Gaps = 4
 Query range = 33-431, Subject range = 1321677-1322077, Strand = Plus/Minus

Query  33       TAAAACGATTGCTAATGAGTCACGTATTTCATCTGGTTCGGTAACTATACCGTCTACTAT  92
                ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1322077  TAAAACGATTGCTAATGAGTCACGTATTTCATCTGGTTCGGTAACTATACCGTCTACTAT  1322018

Query  93       GGACTCAGTGTAACCCTGTAATAAAGAGATTGGCGTACGTAATTCATGTG-TACATTTGC  151
                |||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||
Sbjct  1322017  GGACTCAGTGTAACCCTGTAATAAAGAGATTGGCGTACGTAATTCATGTGATACATTTGC  1321958

Query  152      TATAAAATCTTTTTTCATTTGATCAAGATTATGTTCATTTGTCATATCACAGGATGACCA  211
                |||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||
Sbjct  1321957  TATAAAATCTTTTTTCATTTGATCAAGATTATGTTCATTTGTCATATCAC-GGATGACCA  1321899

Query  212      TGACAATACCACTTCTACCATTTGTTTGAATTCTATCTATATAACTGGAGATAAATACAT  271
                ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1321898  TGACAATACCACTTCTACCATTTGTTTGAATTCTATCTATATAACTGGAGATAAATACAT  1321839

Query  272      AGTACCTTGTATTAATTTCTAATTCTAA-TACTCATTCTGTTGTGATTCAAATGGTGCTT  330
                |||||||||||||||||||||||||||| ||||||||||||||||||||||||| |||||
Sbjct  1321838  AGTACCTTGTATTAATTTCTAATTCTAAATACTCATTCTGTTGTGATTCAAATGTTGCTT  1321779

Query  331      CAATTTGCTGTTCAATAGATTCTTTTGAAAAATCATCAATGTGACGCATAATATAATCAG  390
                |||||||||||||||||||||||||||||||||||||||||||||||||||||| |||||
Sbjct  1321778  CAATTTGCTGTTCAATAGATTCTTTTGAAAAATCATCAATGTGACGCATAATATCATCAG  1321719

Query  391      CCATCTTGTT-GACAATATGATTTCACGTTGATTATTAATGC  431
                |||||||||| |||||||||||||||||||||||||||||||
Sbjct  1321718  CCATCTTGTTTGACAATATGATTTCACGTTGATTATTAATGC  1321677


```
