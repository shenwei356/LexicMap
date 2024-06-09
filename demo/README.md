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


    $ lexicmap index -I refs/ -O demo.lmi
    14:16:59.898 [INFO] removing old output directory: demo.lmi
    14:16:59.899 [INFO] LexicMap v0.4.0
    14:16:59.899 [INFO]   https://github.com/shenwei356/LexicMap
    14:16:59.899 [INFO]
    14:16:59.899 [INFO] checking input files ...
    14:16:59.899 [INFO]   15 input file(s) given
    14:16:59.899 [INFO]
    14:16:59.899 [INFO] --------------------- [ main parameters ] ---------------------
    14:16:59.899 [INFO]
    14:16:59.899 [INFO] input and output:
    14:16:59.899 [INFO]   input directory: refs/
    14:16:59.899 [INFO]     regular expression of input files: (?i)\.(f[aq](st[aq])?|fna)(.gz)?$
    14:16:59.899 [INFO]     *regular expression for extracting reference name from file name: (?i)(.+)\.(f[aq](st[aq])?|fna)(.gz)?$
    14:16:59.899 [INFO]     *regular expressions for filtering out sequences: []
    14:16:59.899 [INFO]   max genome size: 15000000
    14:16:59.899 [INFO]   output directory: demo.lmi
    14:16:59.899 [INFO]
    14:16:59.899 [INFO] k-mer size: 31
    14:16:59.899 [INFO] number of masks: 40000
    14:16:59.899 [INFO] rand seed: 1
    14:16:59.899 [INFO] maximum sketching desert length: 900
    14:16:59.899 [INFO] distance of k-mers to fill deserts: 200
    14:16:59.899 [INFO] prefix length for checking low-complexity in mask generation: 15
    14:16:59.899 [INFO]
    14:16:59.899 [INFO]
    14:16:59.899 [INFO] seeds data chunks: 16
    14:16:59.899 [INFO] seeds data indexing partitions: 512
    14:16:59.899 [INFO] genome batch size: 10000
    14:16:59.899 [INFO]
    14:16:59.899 [INFO]
    14:16:59.899 [INFO] --------------------- [ generating masks ] ---------------------
    14:17:00.187 [INFO]
    14:17:00.187 [INFO] --------------------- [ building index ] ---------------------
    14:17:00.326 [INFO]
    14:17:00.326 [INFO]   ------------------------[ batch 0 ]------------------------
    14:17:00.326 [INFO]   building index for batch 0 with 15 files...
    processed files:  15 / 15 [======================================] ETA: 0s. done
    14:17:01.472 [INFO]   writing seeds...
    14:17:01.689 [INFO]   finished writing seeds in 217.333037ms
    14:17:01.689 [INFO]   finished building index for batch 0 in: 1.362829228s
    14:17:01.690 [INFO]
    14:17:01.690 [INFO] finished building LexicMap index from 15 files with 40000 masks in 1.791529393s
    14:17:01.690 [INFO] LexicMap index saved: demo.lmi
    14:17:01.690 [INFO]
    14:17:01.690 [INFO] elapsed time: 1.791557958s
    14:17:01.690 [INFO]

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
    demo.lmi/: 26.87 MB
      13.64 MB      seeds
      12.93 MB      genomes
     312.53 KB      masks.bin
      375.00 B      genomes.map.bin
      261.00 B      info.toml

## Searching

### A 16S rRNA gene sequence

    $ lexicmap search -d demo.lmi/  q.gene.fasta -o q.gene.fasta.lexicmap.tsv
    09:32:55.551 [INFO] LexicMap v0.4.0
    09:32:55.551 [INFO]   https://github.com/shenwei356/LexicMap
    09:32:55.551 [INFO]
    09:32:55.551 [INFO] checking input files ...
    09:32:55.551 [INFO]   1 input file given: q.gene.fasta
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

### A prophage sequence

    $ lexicmap search -d demo.lmi/ q.prophage.fasta -o q.prophage.fasta.lexicmap.tsv

    $ csvtk head -n 20 q.prophage.fasta.lexicmap.tsv \
        | csvtk mutate -t -n species -f sgenome \
        | csvtk replace -t -f species -k ass2species.map -p '(.+)' -r '{kv}' \
        | csvtk pretty -t

    query         qlen    hits   sgenome           sseqid          qcovGnm   hsp   qcovHSP   alenHSP   pident   gaps   qstart   qend    sstart    send      sstr   slen      species
    -----------   -----   ----   ---------------   -------------   -------   ---   -------   -------   ------   ----   ------   -----   -------   -------   ----   -------   ----------------
    NC_001895.1   33593   1      GCF_003697165.2   NZ_CP033092.2   77.183    1     27.854    9359      97.735   2      1        9357    1864411   1873769   +      4903501   Escherichia coli
    NC_001895.1   33593   1      GCF_003697165.2   NZ_CP033092.2   77.183    2     20.570    6910      96.570   4      17473    24382   1882043   1888948   +      4903501   Escherichia coli
    NC_001895.1   33593   1      GCF_003697165.2   NZ_CP033092.2   77.183    3     8.844     2971      91.754   0      10308    13278   1873846   1876816   +      4903501   Escherichia coli
    NC_001895.1   33593   1      GCF_003697165.2   NZ_CP033092.2   77.183    4     2.355     791       84.703   0      14543    15333   1878801   1879591   +      4903501   Escherichia coli
    NC_001895.1   33593   1      GCF_003697165.2   NZ_CP033092.2   77.183    5     17.644    5927      98.043   0      24355    30281   1853098   1859024   +      4903501   Escherichia coli

### Simulated Oxford Nanopore R10.4.1 long-reads

Here we use the flag `-w/--load-whole-seeds` to accelerate searching.
And `--pseudo-align` is swithced on for using pseudo alignment, which is faster and uses less memory.

    $ lexicmap search -d demo.lmi/ q.long-reads.fasta.gz -o q.long-reads.fasta.gz.lexicmap.tsv.gz -w --pseudo-align -q 70
    20:48:42.090 [INFO] LexicMap v0.4.0
    20:48:42.090 [INFO]   https://github.com/shenwei356/LexicMap
    20:48:42.090 [INFO]
    20:48:42.090 [INFO] checking input files ...
    20:48:42.090 [INFO]   1 input file(s) given
    20:48:42.090 [INFO]
    20:48:42.090 [INFO] loading index: demo.lmi/
    20:48:42.090 [INFO]   reading masks...
    20:48:42.094 [INFO]   reading seeds (k-mer-value) data into memory...
    20:48:42.106 [INFO]   creating genome reader pools, each batch with 16 readers...
    20:48:42.106 [INFO] index loaded in 16.063322ms
    20:48:42.106 [INFO]
    20:48:42.106 [INFO] searching ...
    processed queries: 3584, speed: 2059.159 queries per minute
    20:50:29.538 [INFO]
    20:50:29.538 [INFO] processed queries: 3692, speed: 2061.952 queries per minute
    20:50:29.538 [INFO] 97.0477% (3583/3692) queries matched
    20:50:29.538 [INFO] done searching
    20:50:29.538 [INFO] search results saved to: q.long-reads.fasta.gz.lexicmap.tsv.gz
    20:50:29.542 [INFO]
    20:50:29.542 [INFO] elapsed time: 1m47.45204859s
    20:50:29.542 [INFO]

Result overview:

    csvtk head -n 26 q.long-reads.fasta.gz.lexicmap.tsv.gz \
        | csvtk mutate -t -n species -f sgenome \
        | csvtk replace -t -f species -k ass2species.map -p '(.+)' -r '{kv}' \
        | csvtk pretty -t

    query                  qlen    hits   sgenome           sseqid              qcovGnm   hsp   qcovHSP   alenHSP   pident   gaps   qstart   qend    sstart    send      sstr   slen      species
    --------------------   -----   ----   ---------------   -----------------   -------   ---   -------   -------   ------   ----   ------   -----   -------   -------   ----   -------   --------------------------
    GCF_009759685.1_r164   3132    1      GCF_009759685.1   NZ_CP046654.1       99.042    1     99.042    3102      91.560   -1     20       3121    1768740   1771855   +      3980848   Acinetobacter baumannii
    GCF_002950215.1_r182   1718    1      GCF_002950215.1   NZ_CP026790.1       98.952    1     98.952    1700      70.970   -1     8        1707    134670    136412    +      165702    Shigella flexneri
    GCF_000006945.2_r109   3788    1      GCF_000006945.2   NC_003197.2         98.522    1     98.522    3732      96.009   -1     37       3768    4633323   4637055   -      4857450   Salmonella enterica
    GCF_900638025.1_r28    6375    1      GCF_900638025.1   NZ_LR134481.1       99.357    1     99.357    6334      83.861   -1     6        6339    137524    143936    -      2062405   Haemophilus parainfluenzae
    GCF_000006945.2_r8     7258    1      GCF_000006945.2   NC_003197.2         99.339    1     99.339    7210      95.915   -1     20       7229    4618964   4626209   +      4857450   Salmonella enterica
    GCF_001544255.1_r110   9910    1      GCF_001544255.1   NZ_BCQD01000005.1   99.839    1     99.839    9894      96.972   -1     17       9910    155488    165428    +      191690    Enterococcus faecium
    GCF_003697165.2_r46    2169    1      GCF_003697165.2   NZ_CP033092.2       91.886    1     91.886    1993      74.927   -1     31       2023    4489794   4491835   +      4903501   Escherichia coli
    GCF_000017205.1_r183   14521   1      GCF_000017205.1   NC_009656.1         99.787    1     99.787    14490     94.840   -1     28       14517   3874730   3889304   +      6588339   Pseudomonas aeruginosa
    GCF_002949675.1_r249   1937    3      GCF_002949675.1   NZ_CP026774.1       98.554    1     98.554    1909      91.658   -1     29       1937    3336601   3338518   -      4395762   Shigella dysenteriae
    GCF_002949675.1_r249   1937    3      GCF_002950215.1   NZ_CP026788.1       76.510    1     76.510    1482      86.425   -1     456      1937    3959212   3960699   +      4659463   Shigella flexneri
    GCF_002949675.1_r249   1937    3      GCF_003697165.2   NZ_CP033092.2       98.554    1     98.554    1909      76.492   -1     29       1937    925864    927790    +      4903501   Escherichia coli
    GCF_009759685.1_r168   3398    1      GCF_009759685.1   NZ_CP046654.1       98.558    1     98.558    3349      90.744   -1     22       3370    3276395   3279754   -      3980848   Acinetobacter baumannii
    GCF_001027105.1_r148   20294   1      GCF_001027105.1   NZ_CP011526.1       99.921    1     99.921    20278     95.426   -1     16       20293   2352020   2372396   +      2755072   Staphylococcus aureus
    GCF_000742135.1_r146   19632   1      GCF_000742135.1   NZ_KN046818.1       93.689    1     93.689    18393     92.229   -1     26       18418   3816823   3835391   +      5284261   Klebsiella pneumoniae
    GCF_000006945.2_r43    20355   1      GCF_000006945.2   NC_003197.2         74.576    1     74.576    15180     92.446   -1     32       15211   3589949   3605224   +      4857450   Salmonella enterica
    GCF_001544255.1_r104   1087    1      GCF_001544255.1   NZ_BCQD01000030.1   95.308    1     95.308    1036      95.393   -1     26       1061    9         1050      -      1061      Enterococcus faecium
    GCF_003697165.2_r248   6741    1      GCF_003697165.2   NZ_CP033092.2       99.733    1     99.733    6723      88.560   -1     9        6731    2827221   2833977   -      4903501   Escherichia coli
    GCF_009759685.1_r69    9800    1      GCF_009759685.1   NZ_CP046654.1       99.571    1     99.571    9758      93.465   -1     15       9772    121950    131727    +      3980848   Acinetobacter baumannii
    GCF_001457655.1_r42    8763    1      GCF_001457655.1   NZ_LN831035.1       99.544    1     99.544    8723      89.684   -1     37       8759    637820    646592    +      1890645   Haemophilus influenzae
    GCF_000392875.1_r77    11597   1      GCF_000392875.1   NZ_KB944590.1       91.024    1     91.024    10556     97.640   -1     1015     11570   71860     82452     -      1924212   Enterococcus faecalis
    GCF_000006945.2_r362   20739   1      GCF_000006945.2   NC_003197.2         80.877    1     80.877    16773     77.725   -1     31       16803   237827    254854    -      4857450   Salmonella enterica
    GCF_001027105.1_r40    28470   1      GCF_001027105.1   NZ_CP011526.1       99.919    1     99.919    28447     96.854   -1     24       28470   371992    400532    -      2755072   Staphylococcus aureus
    GCF_002949675.1_r183   8176    1      GCF_002949675.1   NZ_CP026774.1       99.682    1     99.682    8150      97.791   -1     27       8176    4194298   4202441   +      4395762   Shigella dysenteriae
    GCF_002949675.1_r231   2904    3      GCF_002949675.1   NZ_CP026774.1       97.968    1     97.968    2845      94.695   -1     23       2867    960555    963419    -      4395762   Shigella dysenteriae
    GCF_002949675.1_r231   2904    3      GCF_002950215.1   NZ_CP026788.1       94.697    1     94.697    2750      85.498   -1     118      2867    1293581   1296359   -      4659463   Shigella flexneri
    GCF_002949675.1_r231   2904    3      GCF_003697165.2   NZ_CP033092.2       94.697    1     94.697    2750      80.346   -1     118      2867    3015170   3017947   -      4903501   Escherichia coli
