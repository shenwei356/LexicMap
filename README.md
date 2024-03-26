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

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
     query                qlen   qstart   qend   refs   ref               seqid                   afrac    ident      tlen    tstart      tend   strand   seeds   species
    --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
     ERR5396170.1000020   5146      137   5155    6     GCF_000307025.1   NC_018584.1            96.619   95.012   2951805    866862    871860     -         16   Listeria monocytogenes
     ERR5396170.1000020   5146      142   5160    6     GCF_900187225.1   NZ_LT906436.1          91.061   86.044   2864663    833240    838237     -         11   Listeria monocytogenes
     ERR5396170.1000020   5146      132   5111    6     GCF_013282665.1   NZ_CP054041.1          89.992   82.552   2788056   1974929   1979886     -          8   Listeria monocytogenes
     ERR5396170.1000020   5146       99   5268    6     GCA_000183865.1   CM001047.1             82.744   74.049   2884551    813289    818277     -          7   Listeria marthii
     ERR5396170.1000020   5146      133   5113    6     GCF_014229095.1   NZ_JAARZC010000001.1   83.230   72.052   1397946    543086    548044     -          2   Listeria cossartiae
     ERR5396170.1000020   5146      113   5132    6     GCF_014229645.1   NZ_JAATOD010000001.1   82.880   72.192   1428095    563087    568080     -          5   Listeria swaminathanii
     ERR5396170.1000029    480       38    474    1     GCF_001027105.1   NZ_CP011526.1          91.042   95.423   2755072    821077    821514     +          2   Staphylococcus aureus
     ERR5396170.1000006    796       28    796    3     GCF_013394085.1   NZ_CP040910.1          96.608   97.269   1887974   1138941   1139706     +          1   Limosilactobacillus fermentum
     ERR5396170.1000006    796       36    796    3     GCF_013394085.1   NZ_CP040910.1          95.603   97.372   1887974     32649     33406     +          4   Limosilactobacillus fermentum
     ERR5396170.1000006    796       36    796    3     GCF_013394085.1   NZ_CP040910.1          95.603   97.372   1887974    134468    135225     -          4   Limosilactobacillus fermentum
     ERR5396170.1000006    796       36    796    3     GCF_013394085.1   NZ_CP040910.1          95.603   97.240   1887974   1768938   1769695     +          2   Limosilactobacillus fermentum
     ERR5396170.1000006    796       36    796    3     GCF_013394085.1   NZ_CP040910.1          95.603   97.109   1887974    242012    242769     -          3   Limosilactobacillus fermentum
     ERR5396170.1000006    796       36    796    3     GCF_009663775.1   NZ_RDBR01000008.1      91.709   86.438     52610     21391     22148     -          3   Lactobacillus sp. 0.1XD8-4
     ERR5396170.1000006    796       65    796    3     GCF_001591685.1   NZ_BCVJ01000102.1      88.065   88.445      1933       435      1165     +          4   Ligilactobacillus murinus
     ERR5396170.1000021   2020       31   2011    1     GCF_013394085.1   NZ_CP040910.1          98.069   99.849   1887974   1050019   1052001     +         22   Limosilactobacillus fermentum
     ERR5396170.1000032   5782      229   5786    3     GCF_000307025.1   NC_018584.1            95.278   98.094   2951805   1260709   1266255     -         17   Listeria monocytogenes
     ERR5396170.1000032   5782      240   4784    3     GCF_900187225.1   NZ_LT906436.1          73.193   86.011   2864663   1225249   1229767     -          8   Listeria monocytogenes
     ERR5396170.1000032   5782      241   4801    3     GCF_013282665.1   NZ_CP054041.1          70.425   82.195   2788056   2352018   2356571     -          4   Listeria monocytogenes
     ERR5396170.1000000    698       48    650    1     GCF_001457615.1   NZ_LN831024.1          86.390   96.517   6316979   4452088   4452685     +          2   Pseudomonas aeruginosa
    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


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

