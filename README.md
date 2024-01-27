## LexicMap

LexicMap: efficient sequence alignment against millions of microbial genomes​.

## Quick start

    # building index
    lexicmap index -I genomes/ -O db.lmi --force

    # query
    lexicmap map -d db.lmi query.fasta -o query.fasta.txt

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

