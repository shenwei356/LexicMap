## LexicMap

LexicMap: efficient sequence alignment against millions of microbial genomes​.

## Quick start

Building an index

    # from a directory
    lexicmap index -I genomes/ -O db.lmi --force \
        -k 31 --masks 10000

    # from a file list
    lexicmap index -X files.txt -S -O db.lmi --force
        -k 31 --masks 10000

Querying

    lexicmap search -d db.lmi query.fasta -o query.fasta.txt \
        --min-aligned-fraction 60 --min-identity 60

Sample output

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    query     qlen   refs   ref               seqid           afrac    ident       tlen    tstart      tend   strand   seeds
    --------------------------------------------------------------------------------------------------------------------------
    read_1     333    1     GCF_013394085.1   NZ_CP040910.1   86.486   87.847   1887974    196495    196827     -        2
    read_2     578    1     GCF_013394085.1   NZ_CP040910.1   87.024   98.211   1887974   1709172   1709749     +        3
    read_3     297    1     GCF_013394085.1   NZ_CP040910.1   73.401   99.083   1887974   1879531   1879827     +        1
    read_4     586    1     GCF_013394085.1   NZ_CP040910.1   84.642   89.718   1887974   1624549   1625134     +        1
    read_5     749    1     GCF_000392875.1   NZ_KB944589.1   89.987   91.691    682426    342117    342866     +        2
    read_6     699    1     GCF_013394085.1   NZ_CP040910.1   79.542   94.784   1887974     38079     38776     +        3
    read_7    1779    1     GCF_013394085.1   NZ_CP040910.1   96.684   93.372   1887974   1622969   1624753     -        5
    read_8    3139    3     GCF_000307025.1   NC_018584.1     82.606   93.868   2951805    126520    129674     -        6
    read_9    3139    3     GCF_900187225.1   NZ_LT906436.1   80.917   83.504   2864663    112898    116051     -        5
    read_10   3139    3     GCF_013282665.1   NZ_CP054041.1   75.852   77.908   2788056   1295185   1298323     -        2
    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

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

