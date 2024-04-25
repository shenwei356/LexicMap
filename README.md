## <a href="https://bioinf.shenwei.me/LexicMap"><img src="logo.svg" width="30"/></a> LexicMap: efficient sequence alignment against millions of prokaryotic genomes​

LexicMap is a sequence alignment tool aiming to efficiently query gene/plasmid/virus/long-read sequences against up to millions of prokaryotic genomes.

For example, **querying a 51.5-kb plasmid in all 2,340,672 Genbank+Refseq prokaryotic genomes takes only 5 minutes and 2 seconds with 13.7 GB RAM and 48 CPUs, with 17,822 genome hits returned**.
By contrast, BLASTN is unable to run with the same dataset on common servers because it requires >2000 GB RAM. See [performance](#performance).

LexicMap uses a modified [LexicHash](https://doi.org/10.1093/bioinformatics/btad652) algorithm, which supports variable-length substring matching rather than classical fixed-length k-mers matching, to compute seeds for sequence alignment and uses multiple-level storage for fast and low-memory quering of seeds data. See [algorithm overview](#algorithm-overview).

LexicMap is easy to [install](http://bioinf.shenwei.me/LexicMap/installation/) (a binary file with no dependencies) and use ([tutorials](http://bioinf.shenwei.me/LexicMap/tutorials/index/) and [usages](http://bioinf.shenwei.me/LexicMap/usage/lexicmap/)).

More documents: http://bioinf.shenwei.me/LexicMap.

## Quick start

Building an index (see the tutorial of [building an index](http://bioinf.shenwei.me/LexicMap/tutorials/index/)).

    # From a directory with multiple genome files
    lexicmap index -I genomes/ -O db.lmi

    # From a file list with one file per line
    lexicmap index -X files.txt -O db.lmi

Querying (see the tutorial of [searching](http://bioinf.shenwei.me/LexicMap/tutorials/search/)).

    # For short queries like genes or long reads, returning top N hits.
    lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
        --min-qcov-per-genome 70 --min-match-pident 70 --min-qcov-per-hsp 70 --top-n-genomes 1000

    # For longer queries like plasmids, returning all hits.
    lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
        --min-qcov-per-genome 50 --min-match-pident 70 --min-qcov-per-hsp 0  --top-n-genomes 0


Sample output (queries are a few Nanopore Q20 reads). See [output format details](http://bioinf.shenwei.me/LexicMap/tutorials/search/#output).

    query                qlen   qstart   qend   hits   sgenome           sseqid              qcovGnm   hsp   qcovHSP   alenHSP   alenSeg   pident    slen      sstart    send      sstr   seeds   species
    ------------------   ----   ------   ----   ----   ---------------   -----------------   -------   ---   -------   -------   -------   -------   -------   -------   -------   ----   -----   -----------------------------
    ERR5396170.1000017   516    27       514    1      GCF_013394085.1   NZ_CP040910.1       94.574    1     94.574    488       488       100.000   1887974   293509    293996    +      3       Limosilactobacillus fermentum
    ERR5396170.1000047   960    24       812    1      GCF_001027105.1   NZ_CP011526.1       90.521    1     90.521    869       789       89.480    2755072   2204718   2205520   -      6       Staphylococcus aureus
    ERR5396170.1000047   960    881      960    1      GCF_001027105.1   NZ_CP011526.1       90.521    1     90.521    869       80        100.000   2755072   2204568   2204647   -      6       Staphylococcus aureus
    ERR5396170.1000016   740    71       733    1      GCF_013394085.1   NZ_CP040910.1       89.595    1     89.595    663       663       98.492    1887974   13515     14177     +      12      Limosilactobacillus fermentum
    ERR5396170.1000000   698    53       650    1      GCF_001457615.1   NZ_LN831024.1       85.673    1     85.673    598       598       97.324    6316979   4452083   4452685   +      4       Pseudomonas aeruginosa
    ERR5396170.1000005   2516   38       2510   5      GCF_000006945.2   NC_003197.2         98.291    1     98.291    2473      2473      99.151    4857450   3198806   3201283   +      15      Salmonella enterica
    ERR5396170.1000005   2516   38       2497   5      GCF_008692785.1   NZ_VXJV01000001.1   97.774    1     97.774    2460      2460      96.098    797633    423400    425865    +      14      Salmonella enterica
    ERR5396170.1000005   2516   40       2510   5      GCA_900478215.1   LS483478.1          98.211    1     98.211    2471      2471      95.548    4624613   785866    788342    -      13      Salmonella enterica
    ERR5396170.1000005   2516   1350     2497   5      GCF_008692845.1   NZ_VXJW01000004.1   86.765    1     86.765    2183      1148      95.557    366711    6705      7855      +      12      Salmonella enterica
    ERR5396170.1000005   2516   634      1309   5      GCF_008692845.1   NZ_VXJW01000004.1   86.765    1     86.765    2183      676       89.053    366711    5991      6664      +      12      Salmonella enterica
    ERR5396170.1000005   2516   387      608    5      GCF_008692845.1   NZ_VXJW01000004.1   86.765    1     86.765    2183      222       85.135    366711    5745      5965      +      12      Salmonella enterica
    ERR5396170.1000005   2516   69       205    5      GCF_008692845.1   NZ_VXJW01000004.1   86.765    1     86.765    2183      137       83.212    366711    5426      5563      +      12      Salmonella enterica
    ERR5396170.1000005   2516   1830     2263   5      GCF_000252995.1   NC_015761.1         78.378    1     78.378    1972      434       97.696    4460105   2898281   2898717   +      7       Salmonella bongori
    ERR5396170.1000005   2516   2307     2497   5      GCF_000252995.1   NC_015761.1         78.378    1     78.378    1972      191       76.963    4460105   2898761   2898951   +      7       Salmonella bongori
    ERR5396170.1000005   2516   415      938    5      GCF_000252995.1   NC_015761.1         78.378    1     78.378    1972      524       86.641    4460105   2896865   2897391   +      7       Salmonella bongori
    ERR5396170.1000005   2516   1113     1807   5      GCF_000252995.1   NC_015761.1         78.378    1     78.378    1972      695       71.511    4460105   2897564   2898258   +      7       Salmonella bongori
    ERR5396170.1000005   2516   961      1088   5      GCF_000252995.1   NC_015761.1         78.378    1     78.378    1972      128       85.156    4460105   2897414   2897541   +      7       Salmonella bongori

    Note: the column `species` is added by mapping genome ID (column `sgenome`) to taxonomic information.

Learn more [tutorials](http://bioinf.shenwei.me/LexicMap/tutorials/index/) and [usages](http://bioinf.shenwei.me/LexicMap/usage/lexicmap/).

## Performance

|dataset          |genomes  |query          |query_len|genome_hits|time    |RAM    |
|:----------------|--------:|:--------------|--------:|----------:|-------:|------:|
|GTDB repr        |85,205   |a MutL gene    |1,956 bp |2          |0.9 s   |460 MB |
|                 |         |a 16S rRNA gene|1,542 bp |13,466     |4.0 s   |765 MB |
|                 |         |a plasmid      |51,466 bp|2          |1.1 s   |752 MB |
|GTDB complete    |402,538  |a MutL gene    |1,956 bp |268        |3.8 s   |544 MB |
|                 |         |a 16S rRNA gene|1,542 bp |169,480    |2 m 14 s|2.9 GB |
|                 |         |a plasmid      |51,466 bp|3,649      |56 s    |2.9 GB |
|Genbank+RefSeq   |2,340,672|a MutL gene    |1,956 bp |817        |10.0 s  |2.3 GB |
|                 |         |a 16S rRNA gene|1,542 bp |1,148,049  |5 m 34 s|11.8 GB|
|                 |         |a plasmid      |51,466 bp|19,265     |3 m 32 s|15.7 GB|
|AllTheBacteria HQ|1,858,610|a MutL gene    |1,956 bp |404        |18.7 s  |2.4 GB |
|                 |         |a 16S rRNA gene|1,542 bp |1,193,874  |13 m 8 s|9.4 GB |
|                 |         |a plasmid      |51,466 bp|10,954     |5 m 25 s|9.7 GB |

Notes:
- All files are stored on a server with HDD disks.
- Tests are performed in a single cluster node with 48 CPU cores (Intel Xeon Gold 6336Y CPU @ 2.40 GHz).
- Index building parameters: `-k 31 -m 40000`. Genome batch size: `-b 10000` for GTDB datasets, `-b 131072` for others.
- Searching parameters: `--top-n-genomes 0 --min-qcov-per-genome 50 --min-match-pident 70 --min-qcov-per-hsp 0`.

## Installation

LexicMap is implemented in [Go](https://go.dev/) programming language,
executable binary files **for most popular operating systems** are freely available
in [release](https://github.com/shenwei356/lexicmap/releases) page.

Or install with `conda`:

    conda install -c bioconda lexicmap

## Algorithm overview

<img src="overview.svg" alt="LexicMap overview" width="800"/>

## Related projects

- High-performance [LexicHash](https://github.com/shenwei356/LexicHash) computation in Go.

## Support

Please [open an issue](https://github.com/shenwei356/LexicMap/issues) to report bugs,
propose new functions or ask for help.

## License

[MIT License](https://github.com/shenwei356/LexicMap/blob/master/LICENSE)

