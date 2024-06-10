## <a href="https://bioinf.shenwei.me/LexicMap"><img src="logo.svg" width="30"/></a> LexicMap: efficient sequence alignment against millions of prokaryotic genomes​

[![Latest Version](https://img.shields.io/github/release/shenwei356/LexicMap.svg?style=flat?maxAge=86400)](https://github.com/shenwei356/LexicMap/releases)
[![Anaconda Cloud](https://anaconda.org/bioconda/lexicmap/badges/version.svg)](https://anaconda.org/bioconda/lexicmap)
[![Cross-platform](https://img.shields.io/badge/platform-any-ec2eb4.svg?style=flat)](http://bioinf.shenwei.me/LexicMap/installation/)

LexicMap is a nucleotide sequence alignment tool for efficiently querying gene, plasmid, viral, or long-read sequences against up to millions of prokaryotic genomes.

**Motivation**: Alignment against a database of genomes is a fundamental operation in bioinformatics, popularised by BLAST. However, given the increasing rate at which genomes are sequenced, existing tools struggle to scale. Current tools either attempt full alignment but face challenges of high memory consumption and slow speeds, or they fall back on k-mer indexing, without information of where matches occur in the genome.

**Results**: In LexicMap, a [modified version](https://github.com/shenwei356/lexichash) of the sequence sketching method [LexicHash](https://doi.org/10.1093/bioinformatics/btad652) is adopted to compute alignment seeds.
A multi-level index enables fast and low-memory variable-length seed matching and alignment on a single server
at the scale of millions of genomes (see [algorithm overview](#algorithm-overview)),
successfully indexing and searching both RefSeq+GenBank, and the [AllTheBacteria](https://www.biorxiv.org/content/10.1101/2024.03.08.584059v1) datasets (2.3 and 1.9 million genomes respectively).
Running at this scale has previously only been achieved by [Phylign](https://github.com/karel-brinda/Phylign) (previously called mof-search).

For example, for searching in all 2,340,672 Genbank+Refseq prokaryotic genomes, *BLASTN is unable to run with this dataset on common servers because it requires >2000 GB RAM*.  (see [performance](#performance)). **With LexicMap**,
- Querying with a **1.3-kb marker gene** took **4 seconds** with **2.12 GB RAM** and 48 CPUs, with **16,788 genome hits** returned.
- Querying with a **52.8-kb plasmid** took **4 minutes** with **22.46 GB RAM** and 48 CPUs, with **495,915 genome hits** returned.
- Querying with a **1.5-kb 16S rRNA gene** took **5.5 minutes** with **17.25 GB RAM** and 48 CPUs, with **1,894,943 genome hits** returned.

LexicMap is easy to [install](http://bioinf.shenwei.me/LexicMap/installation/) (binary files with no dependencies are provided for most common platforms, ) and use ([tutorials](http://bioinf.shenwei.me/LexicMap/tutorials/index/) and [usages](http://bioinf.shenwei.me/LexicMap/usage/lexicmap/)).

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
        --min-qcov-per-hsp 70 --min-qcov-per-genome 70  --top-n-genomes 1000

    # For longer queries like plasmids, returning all hits.
    lexicmap search -d db.lmi query.fasta -o query.fasta.lexicmap.tsv \
        --min-qcov-per-hsp 0  --min-qcov-per-genome 0   --top-n-genomes 0


    # Extracting similar sequences for a query gene.
      # search matches with query coverage >= 90%
      lexicmap search -d gtdb_complete.lmi/ b.gene_E_faecalis_SecY.fasta --all -o results.tsv \
          --min-qcov-per-hsp 90

      # extract matched sequences as FASTA format
      sed 1d results.tsv | awk -F'\t' '{print ">"$5":"$14"-"$15":"$16"\n"$20;}' | seqkit seq -g > results.fasta


Sample output (queries are a few Nanopore Q20 reads). See [output format details](http://bioinf.shenwei.me/LexicMap/tutorials/search/#output-format).

    query                qlen   hits   sgenome           sseqid          qcovGnm   hsp   qcovHSP   alenHSP   pident    gaps   qstart   qend   sstart    send      sstr   slen
    ------------------   ----   ----   ---------------   -------------   -------   ---   -------   -------   -------   ----   ------   ----   -------   -------   ----   -------
    ERR5396170.1000016   740    1      GCF_013394085.1   NZ_CP040910.1   89.595    1     89.595    663       99.246    0      71       733    13515     14177     +      1887974
    ERR5396170.1000000   698    1      GCF_001457615.1   NZ_LN831024.1   85.673    1     85.673    603       98.010    5      53       650    4452083   4452685   +      6316979
    ERR5396170.1000017   516    1      GCF_013394085.1   NZ_CP040910.1   94.574    1     94.574    489       99.591    2      27       514    293509    293996    +      1887974
    ERR5396170.1000012   848    1      GCF_013394085.1   NZ_CP040910.1   95.165    1     95.165    811       97.411    7      22       828    190329    191136    -      1887974
    ERR5396170.1000038   1615   1      GCA_000183865.1   CM001047.1      64.706    1     60.000    973       95.889    13     365      1333   88793     89756     -      2884551
    ERR5396170.1000038   1615   1      GCA_000183865.1   CM001047.1      64.706    2     4.706     76        98.684    0      266      341    89817     89892     -      2884551
    ERR5396170.1000036   1159   1      GCF_013394085.1   NZ_CP040910.1   95.427    1     95.427    1107      99.729    1      32       1137   1400097   1401203   +      1887974
    ERR5396170.1000031   814    4      GCF_013394085.1   NZ_CP040910.1   86.486    1     86.486    707       99.151    3      104      807    242235    242941    -      1887974
    ERR5396170.1000031   814    4      GCF_013394085.1   NZ_CP040910.1   86.486    2     86.486    707       98.444    3      104      807    1138777   1139483   +      1887974
    ERR5396170.1000031   814    4      GCF_013394085.1   NZ_CP040910.1   86.486    3     84.152    688       98.983    4      104      788    154620    155306    -      1887974
    ERR5396170.1000031   814    4      GCF_013394085.1   NZ_CP040910.1   86.486    4     84.029    687       99.127    3      104      787    32477     33163     +      1887974
    ERR5396170.1000031   814    4      GCF_013394085.1   NZ_CP040910.1   86.486    5     72.727    595       98.992    3      104      695    1280183   1280777   +      1887974
    ERR5396170.1000031   814    4      GCF_013394085.1   NZ_CP040910.1   86.486    6     11.671    95        100.000   0      693      787    1282480   1282574   +      1887974
    ERR5396170.1000031   814    4      GCF_013394085.1   NZ_CP040910.1   86.486    7     82.064    671       99.106    3      120      787    1768782   1769452   +      1887974


CIGAR string, aligned query and subject sequences can be outputted as extra columns via the flag `-a/--all`.

    # export blast-like alignment text
    lexicmap search -d db.lmi query.fasta --all --quiet \
        | sed 1d | awk -F'\t' '{print ">"$1":"$12"-"$13" vs "$5":"$14"-"$15":"$16" pident:"$10" gaps:"$11"\n"$18"\n"$19"\n"$21"\n"$20"\n";}'

    >NC_000913.3:4166659-4168200:5-120 vs CAMDMN010000161.1:25-140:- pident:87.069 gaps:0
    14M1X25M1X24M1X1M1X2M1X1M4X8M4X2M1X2M1X22M
    TGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAA
    |||||||||||||| ||||||||||||||||||||||||| |||||||||||||||||||||||| | || |    ||||||||    || || ||||||||||||||||||||||
    TGAAGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGGCAGCATGGTCTAGCTTGCTAGACTGATGGCGAGTGGCGGACGGGTGAGTAA


Learn more [tutorials](http://bioinf.shenwei.me/LexicMap/tutorials/index/) and [usages](http://bioinf.shenwei.me/LexicMap/usage/lexicmap/).

## Performance

See [performance](https://bioinf.shenwei.me/LexicMap/introduction/#performance).

## Installation

LexicMap is implemented in [Go](https://go.dev/) programming language,
executable binary files **for most popular operating systems** are freely available
in [release page](https://github.com/shenwei356/lexicmap/releases).

Or install with `conda`:

    conda install -c bioconda lexicmap

## Algorithm overview

<img src="overview.svg" alt="LexicMap overview" width="800"/>

## Related projects

- High-performance [LexicHash](https://github.com/shenwei356/LexicHash) computation in Go.
- [Wavefront alignment algorithm (WFA) in Golang](https://github.com/shenwei356/wfa).

## Support

Please [open an issue](https://github.com/shenwei356/LexicMap/issues) to report bugs,
propose new functions or ask for help.

## License

[MIT License](https://github.com/shenwei356/LexicMap/blob/master/LICENSE)

