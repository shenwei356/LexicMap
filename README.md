## LexicMap

LexicMap: efficient sequence alignment against millions of microbial genomes​.

## Quick start

Building an index

    # from a directory
    lexicmap index -I genomes/ -O db.lmi --force

    # from a file list
    lexicmap index -X files.txt -S -O db.lmi --force

Querying

    lexicmap search -d db.lmi query.fasta -o query.fasta.txt \
        --min-aligned-fraction 60 --min-identity 60

Sample output (queries are a few Nanopure Q20 reads).
The column `species` is added by mapping genome ID to taxonomic information

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
     query                qlen   refs   ref               seqid                afrac     ident      tlen    tstart      tend   strand   seeds   species
    --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
     ERR5396170.1000019   1560    2     GCF_000227465.1   NC_016047.1         92.885    92.547   4207222   2215650   2217212     -          5   Bacillus spizizenii
     ERR5396170.1000019   1560    2     GCF_000332645.1   NZ_AMXN01000009.1   92.628    79.862    199646     46019     47578     +          1   Bacillus inaquosorum
     ERR5396170.1000017    516    1     GCF_013394085.1   NZ_CP040910.1       94.574   100.000   1887974    293482    293998     +          2   Limosilactobacillus fermentum
     ERR5396170.1000016    740    1     GCF_013394085.1   NZ_CP040910.1       89.595    98.492   1887974     13445     14184     +         10   Limosilactobacillus fermentum
     ERR5396170.1000011   2747    4     GCF_000006945.2   NC_003197.2         98.216    99.481   4857450   4627934   4630682     +         12   Salmonella enterica
     ERR5396170.1000011   2747    4     GCA_900478215.1   LS483478.1          95.013    89.157   4624613   4044803   4047572     -          6   Salmonella enterica
     ERR5396170.1000011   2747    4     GCF_008692845.1   NZ_VXJW01000007.1   95.704    87.067    285234    181513    184261     +          6   Salmonella enterica
     ERR5396170.1000011   2747    4     GCF_008692785.1   NZ_VXJV01000010.1   77.321    80.650    167301     89095     91851     -          2   Salmonella enterica
     ERR5396170.1000000    698    1     GCF_001457615.1   NZ_LN831024.1       86.390    96.517   6316979   4452036   4452733     +          2   Pseudomonas aeruginosa
     ERR5396170.1000005   2516    5     GCF_000006945.2   NC_003197.2         98.490    98.951   4857450   3198773   3201289     +         10   Salmonella enterica
     ERR5396170.1000005   2516    5     GCF_008692785.1   NZ_VXJV01000001.1   98.013    95.864    797633    423368    425884     +          5   Salmonella enterica
     ERR5396170.1000005   2516    5     GCA_900478215.1   LS483478.1          98.450    95.317   4624613    785863    788377     -          7   Salmonella enterica
     ERR5396170.1000005   2516    5     GCF_008692845.1   NZ_VXJW01000004.1   87.599    91.742    366711      5358      7871     +          6   Salmonella enterica
     ERR5396170.1000005   2516    5     GCF_000252995.1   NC_015761.1         88.593    80.485   4460105   2896454   2898967     +          8   Salmonella bongori
    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

## Installation

LexicMap is implemented in [Go](https://go.dev/) programming language,
executable binary files **for most popular operating systems** are freely available
in [release](https://github.com/shenwei356/lexicmap/releases) page.

Or install with `conda`:

    conda install -c bioconda lexicmap

## Related projects

- High-performance [LexicHash](https://github.com/shenwei356/LexicHash) computation in Go.

## Support

Please [open an issue](https://github.com/shenwei356/LexicMap/issues) to report bugs,
propose new functions or ask for help.

## License

[MIT License](https://github.com/shenwei356/LexicMap/blob/master/LICENSE)

