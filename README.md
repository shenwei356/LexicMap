## LexicMap

LexicMap: efficient sequence alignment against millions of microbial genomes​.

## Quick start

Building index

    # from a directory
    lexicmap index -I genomes/ -O db.lmi --force \
        -k 31 --masks 10000

    # from a file list
    lexicmap index -X files.txt -S -O db.lmi --force
        -k 31 --masks 10000

Querying

    lexicmap search -d db.lmi query.fasta -o query.fasta.txt \
        --min-aligned-fraction 60 --min-identity 60

## Installation

LexicMap is implemented in [Go](https://go.dev/) programming language,
executable binary files **for most popular operating systems** are freely available
in [release](https://github.com/shenwei356/lexicmap/releases) page.

Or install with `conda`:

    conda install -c bioconda lexicmap

## Support

Please [open an issue](https://github.com/shenwei356/LexicMap/issues) to report bugs,
propose new functions or ask for help.

## License

[MIT License](https://github.com/shenwei356/LexicMap/blob/master/LICENSE)

