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
    14:16:59.899 [INFO] LexicMap v0.3.0
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
    14:16:59.899 [INFO] prefix for checking low-complexity and choosing k-mers to fill sketching deserts: 15
    14:16:59.899 [INFO] distance of k-mers to fill deserts: 200
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
    demo.lmi/: 26.90 MB
      13.67 MB      seeds
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

    $ csvtk head -n 21 q.gene.fasta.lexicmap.tsv \
        | csvtk mutate -t -n species -f sgenome \
        | csvtk replace -t -f species -k ass2species.map -p '(.+)' -r '{kv}' \
        | csvtk pretty -t

    query                         qlen   hits   sgenome           sseqid          qcovGnm   hsp   qcovHSP   alenHSP   pident   qstart   qend   sstart    send      sstr   slen      seeds   species
    ---------------------------   ----   ----   ---------------   -------------   -------   ---   -------   -------   ------   ------   ----   -------   -------   ----   -------   -----   --------------------
    NC_000913.3:4166659-4168200   1542   15     GCF_003697165.2   NZ_CP033092.2   100.000   1     100.000   1542      99.287   1        1542   458559    460100    +      4903501   23      Escherichia coli
    NC_000913.3:4166659-4168200   1542   15     GCF_003697165.2   NZ_CP033092.2   100.000   2     100.000   1542      99.287   1        1542   1285123   1286664   +      4903501   23      Escherichia coli
    NC_000913.3:4166659-4168200   1542   15     GCF_003697165.2   NZ_CP033092.2   100.000   3     100.000   1542      99.287   1        1542   3780640   3782181   -      4903501   23      Escherichia coli
    NC_000913.3:4166659-4168200   1542   15     GCF_003697165.2   NZ_CP033092.2   100.000   4     100.000   1542      99.287   1        1542   4551515   4553056   -      4903501   23      Escherichia coli
    NC_000913.3:4166659-4168200   1542   15     GCF_003697165.2   NZ_CP033092.2   100.000   5     100.000   1542      99.287   1        1542   4591684   4593225   -      4903501   23      Escherichia coli
    NC_000913.3:4166659-4168200   1542   15     GCF_003697165.2   NZ_CP033092.2   100.000   6     100.000   1542      99.287   1        1542   4844587   4846128   -      4903501   23      Escherichia coli
    NC_000913.3:4166659-4168200   1542   15     GCF_003697165.2   NZ_CP033092.2   100.000   7     100.000   1542      99.092   1        1542   4726193   4727734   -      4903501   23      Escherichia coli
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   1     100.000   1542      99.027   1        1542   3216505   3218046   +      4659463   22      Shigella flexneri
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   2     100.000   1542      98.962   1        1542   3396068   3397609   +      4659463   24      Shigella flexneri
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   3     100.000   1542      98.962   1        1542   3119331   3120872   +      4659463   23      Shigella flexneri
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   4     100.000   1542      98.898   1        1542   3355632   3357173   +      4659463   24      Shigella flexneri
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   5     100.000   1542      98.898   1        1542   4223146   4224687   +      4659463   24      Shigella flexneri
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   6     100.000   1542      98.833   1        1542   3540450   3541991   +      4659463   23      Shigella flexneri
    NC_000913.3:4166659-4168200   1542   15     GCF_002950215.1   NZ_CP026788.1   100.000   7     100.000   1542      98.768   1        1542   2125377   2126917   -      4659463   18      Shigella flexneri
    NC_000913.3:4166659-4168200   1542   15     GCF_002949675.1   NZ_CP026774.1   97.601    1     97.601    1505      98.937   1        1542   1662010   1663551   -      4395762   30      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   15     GCF_002949675.1   NZ_CP026774.1   97.601    2     97.601    1505      98.937   1        1542   2536624   2538165   +      4395762   30      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   15     GCF_002949675.1   NZ_CP026774.1   97.601    3     97.601    1505      98.937   1        1542   2636477   2638018   +      4395762   30      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   15     GCF_002949675.1   NZ_CP026774.1   97.601    4     97.601    1505      98.937   1        1542   2768883   2770424   +      4395762   30      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   15     GCF_002949675.1   NZ_CP026774.1   97.601    5     97.601    1505      98.937   1        1542   2810845   2812386   +      4395762   30      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   15     GCF_002949675.1   NZ_CP026774.1   97.601    6     97.601    1505      98.937   1        1542   3061592   3063133   +      4395762   30      Shigella dysenteriae
    NC_000913.3:4166659-4168200   1542   15     GCF_002949675.1   NZ_CP026774.1   97.601    7     97.601    1505      98.937   1        1542   3646778   3648319   +      4395762   30      Shigella dysenteriae

### A prophage sequence

    $ lexicmap search -d demo.lmi/ q.prophage.fasta -o q.prophage.fasta.lexicmap.tsv

    $ csvtk head -n 20 q.prophage.fasta.lexicmap.tsv \
        | csvtk mutate -t -n species -f sgenome \
        | csvtk replace -t -f species -k ass2species.map -p '(.+)' -r '{kv}' \
        | csvtk pretty -t
    query         qlen    hits   sgenome           sseqid          qcovGnm   hsp   qcovHSP   alenHSP   pident   qstart   qend    sstart    send      sstr   slen      seeds   species
    -----------   -----   ----   ---------------   -------------   -------   ---   -------   -------   ------   ------   -----   -------   -------   ----   -------   -----   ----------------
    NC_001895.1   33593   1      GCF_003697165.2   NZ_CP033092.2   73.289    1     55.878    18771     93.495   1        24382   1864411   1888948   +      4903501   81      Escherichia coli
    NC_001895.1   33593   1      GCF_003697165.2   NZ_CP033092.2   73.289    2     17.495    5877      96.103   24355    30281   1853098   1859024   +      4903501   52      Escherichia coli


### Simulated Oxford Nanopore R10.4.1 long-reads

Here we use the flag `-w/--load-whole-seeds` to accelerate searching.

    $ lexicmap search -d demo.lmi/ q.long-reads.fasta.gz -o q.long-reads.fasta.gz.lexicmap.tsv.gz -w
    13:32:54.624 [INFO] LexicProf v0.3.0
    13:32:54.624 [INFO]   https://github.com/shenwei356/LexicMap
    13:32:54.624 [INFO]
    13:32:54.624 [INFO] checking input files ...
    13:32:54.624 [INFO]   1 input file(s) given
    13:32:54.624 [INFO]
    13:32:54.624 [INFO] loading index: demo.lmi/
    13:32:54.624 [INFO]   reading masks...
    13:32:54.628 [INFO]   reading seeds (k-mer-value) data into memory...
    13:32:54.639 [INFO]   creating genome reader pools, each batch with 16 readers...
    13:32:54.639 [INFO] index loaded in 15.420708ms
    13:32:54.639 [INFO]
    13:32:54.639 [INFO] searching ...
    processed queries: 3584, speed: 12460.410 queries per minute
    13:33:12.528 [INFO]
    13:33:12.528 [INFO] processed queries: 3692, speed: 12383.331 queries per minute
    13:33:12.528 [INFO] 97.0748% (3584/3692) queries matched
    13:33:12.528 [INFO] done searching
    13:33:12.528 [INFO] search results saved to: q.long-reads.fasta.gz.lexicmap.tsv.gz
    13:33:12.532 [INFO]
    13:33:12.532 [INFO] elapsed time: 17.908286401s
    13:33:12.532 [INFO]

Result overview:

    csvtk head -n 26 q.long-reads.fasta.gz.lexicmap.tsv.gz \
        | csvtk mutate -t -n species -f sgenome \
        | csvtk replace -t -f species -k ass2species.map -p '(.+)' -r '{kv}' \
        | csvtk pretty -t

    query                  qlen    hits   sgenome           sseqid              qcovGnm   hsp   qcovHSP   alenHSP   pident   qstart   qend    sstart    send      sstr   slen      seeds   species
    --------------------   -----   ----   ---------------   -----------------   -------   ---   -------   -------   ------   ------   -----   -------   -------   ----   -------   -----   --------------------------
    GCF_003697165.2_r46    2169    2      GCF_003697165.2   NZ_CP033092.2       91.655    1     91.655    1988      80.181   31       2146    4489794   4491918   +      4903501   7       Escherichia coli
    GCF_003697165.2_r46    2169    2      GCF_002950215.1   NZ_CP026788.1       23.559    1     23.559    511       72.994   51       664     3495342   3495975   +      4659463   2       Shigella flexneri
    GCF_002950215.1_r182   1718    2      GCF_002950215.1   NZ_CP026790.1       93.714    1     93.714    1610      75.528   8        1707    134670    136412    +      165702    2       Shigella flexneri
    GCF_002950215.1_r182   1718    2      GCF_002949675.1   NZ_CP026774.1       28.231    1     28.231    485       69.691   132      616     3617642   3618136   -      4395762   1       Shigella dysenteriae
    GCF_002950215.1_r182   1718    2      GCF_002949675.1   NZ_CP026774.1       28.231    2     28.231    485       69.691   132      616     4179234   4179728   +      4395762   1       Shigella dysenteriae
    GCF_002950215.1_r182   1718    2      GCF_002949675.1   NZ_CP026774.1       28.231    3     28.231    485       69.691   132      616     4008855   4009349   +      4395762   2       Shigella dysenteriae
    GCF_002950215.1_r182   1718    2      GCF_002949675.1   NZ_CP026774.1       28.231    4     28.231    485       69.691   132      616     1515091   1515585   -      4395762   1       Shigella dysenteriae
    GCF_002950215.1_r182   1718    2      GCF_002949675.1   NZ_CP026774.1       28.231    5     28.231    485       69.485   132      616     3425441   3425935   -      4395762   2       Shigella dysenteriae
    GCF_002950215.1_r182   1718    2      GCF_002949675.1   NZ_CP026774.1       28.231    6     28.231    485       69.485   132      616     1488485   1488979   +      4395762   1       Shigella dysenteriae
    GCF_002950215.1_r182   1718    2      GCF_002949675.1   NZ_CP026774.1       28.231    7     28.231    485       69.485   132      616     1390293   1390787   +      4395762   2       Shigella dysenteriae
    GCF_002950215.1_r182   1718    2      GCF_002949675.1   NZ_CP026774.1       28.231    8     28.231    485       68.660   132      616     2885716   2886210   +      4395762   1       Shigella dysenteriae
    GCF_002950215.1_r182   1718    2      GCF_002949675.1   NZ_CP026775.1       28.231    9     28.231    485       68.247   132      616     159236    159730    -      182697    1       Shigella dysenteriae
    GCF_002950215.1_r182   1718    2      GCF_002949675.1   NZ_CP026774.1       28.231    10    28.231    485       67.835   132      616     544073    544567    +      4395762   2       Shigella dysenteriae
    GCF_002950215.1_r182   1718    2      GCF_002949675.1   NZ_CP026774.1       28.231    11    28.231    485       65.361   132      616     3713569   3714063   +      4395762   1       Shigella dysenteriae
    GCF_009759685.1_r164   3132    1      GCF_009759685.1   NZ_CP046654.1       97.797    1     97.797    3063      93.144   20       3121    1768740   1771855   +      3980848   17      Acinetobacter baumannii
    GCF_900638025.1_r28    6375    1      GCF_900638025.1   NZ_LR134481.1       96.188    1     96.188    6132      86.954   6        6339    137524    143936    -      2062405   44      Haemophilus parainfluenzae
    GCF_000006945.2_r109   3788    2      GCF_000006945.2   NC_003197.2         97.941    1     97.941    3710      96.604   37       3768    4633323   4637055   -      4857450   15      Salmonella enterica
    GCF_000006945.2_r109   3788    2      GCF_000742135.1   NZ_KN046818.1       21.753    1     21.753    824       66.141   1057     3410    164225    166586    +      5284261   2       Klebsiella pneumoniae
    GCF_000006945.2_r8     7258    2      GCF_000006945.2   NC_003197.2         99.339    1     99.339    7210      96.394   20       7229    4618964   4626209   +      4857450   59      Salmonella enterica
    GCF_000006945.2_r8     7258    2      GCF_000006945.2   NC_003197.2         99.339    2     29.044    2108      95.778   20       2127    4618964   4621086   +      4857450   1       Salmonella enterica
    GCF_000006945.2_r8     7258    2      GCF_002949675.1   NZ_CP026774.1       14.329    1     14.329    1040      71.731   4722     7048    2167923   2170245   -      4395762   1       Shigella dysenteriae
    GCF_001544255.1_r110   9910    2      GCF_001544255.1   NZ_BCQD01000005.1   99.839    1     99.839    9894      97.433   17       9910    155488    165428    +      191690    117     Enterococcus faecium
    GCF_001544255.1_r110   9910    2      GCF_000392875.1   NZ_KB944590.1       0.585     1     0.585     58        65.517   6219     6276    837470    837526    +      1924212   1       Enterococcus faecalis
    GCF_000017205.1_r183   14521   1      GCF_000017205.1   NC_009656.1         99.325    1     99.325    14423     95.840   28       14517   3874730   3889304   +      6588339   125     Pseudomonas aeruginosa
    GCF_002949675.1_r249   1937    5      GCF_002949675.1   NZ_CP026774.1       98.554    1     98.554    1909      92.090   29       1937    3336601   3338518   -      4395762   16      Shigella dysenteriae
    GCF_002949675.1_r249   1937    5      GCF_002950215.1   NZ_CP026788.1       97.006    1     97.006    1879      86.748   29       1937    3958382   3960699   +      4659463   10      Shigella flexneri
