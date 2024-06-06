## <a href="https://bioinf.shenwei.me/LexicMap"><img src="logo.svg" width="30"/></a> LexicMap: efficient sequence alignment against millions of prokaryotic genomes​

[![Latest Version](https://img.shields.io/github/release/shenwei356/LexicMap.svg?style=flat?maxAge=86400)](https://github.com/shenwei356/LexicMap/releases)
[![Anaconda Cloud](https://anaconda.org/bioconda/lexicmap/badges/version.svg)](https://anaconda.org/bioconda/lexicmap)
[![Cross-platform](https://img.shields.io/badge/platform-any-ec2eb4.svg?style=flat)](http://bioinf.shenwei.me/LexicMap/installation/)

LexicMap is a nucleotide sequence alignment tool for efficiently querying gene, plasmid, viral, or long-read sequences against up to millions of prokaryotic genomes.

**Motivation**:
Alignment against a database of genomes is a fundamental operation in bioinformatics, popularised by BLAST.
However, given the increasing rate at which genomes are sequenced, existing tools struggle to scale.
Current tools either attempt full alignment but face challenges of high memory consumption and slow speeds,
or they fall back on k-mer indexing, without information of where matches occur in the genome.

**Results**:
In LexicMap, a [modified version](https://github.com/shenwei356/lexichash) of the sequence sketching method
[LexicHash](https://doi.org/10.1093/bioinformatics/btad652) is adopted to compute alignment seeds.
A multi-level index enables fast and low-memory variable-length seed matching and alignment on a single server
at the scale of millions of genomes (see [algorithm overview](#algorithm-overview)),
successfully indexing and searching both RefSeq+GenBank,
and the [AllTheBacteria](https://www.biorxiv.org/content/10.1101/2024.03.08.584059v1) datasets (2.3 and 1.9 million genomes respectively).
Running at this scale has previously only been achieved by [Phylign](https://github.com/karel-brinda/Phylign) (previously called mof-search).

For example, **querying a 52.8-kb plasmid in all 2,340,672 Genbank+Refseq prokaryotic genomes takes only 4 minutes and 8 seconds with 15 GB RAM and 48 CPUs, with 494,860 genome hits returned**.
In contrast, BLASTN is unable to run with the same dataset on common servers because it requires >2000 GB RAM. See [performance](#performance).

**LexicMap is easy to [install](http://bioinf.shenwei.me/LexicMap/installation/)** (a binary file with no dependencies) **and use** ([tutorials](http://bioinf.shenwei.me/LexicMap/tutorials/index/) and [usages](http://bioinf.shenwei.me/LexicMap/usage/lexicmap/)).

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
      sed 1d results.tsv | awk '{print ">"$5":"$13"-"$14":"$15"\n"$19;}' > results.fasta


Sample output (queries are a few Nanopore Q20 reads). See [output format details](http://bioinf.shenwei.me/LexicMap/tutorials/search/#output).

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


Matched query and subject sequences can be outputted as extra two columns via the flag `-a/--all`.

Learn more [tutorials](http://bioinf.shenwei.me/LexicMap/tutorials/index/) and [usages](http://bioinf.shenwei.me/LexicMap/usage/lexicmap/).

## Performance

### Indexing

|dataset          |genomes  |gzip_size|tool    |db_size|time     |RAM   |
|:----------------|--------:|--------:|:-------|------:|--------:|-----:|
|GTDB complete    |402,538  |578 GB   |LexicMap|510 GB |3 h 26 m |35 GB |
|                 |         |         |Blastn  |360 GB |3 h 11 m |718 MB|
|Genbank+RefSeq   |2,340,672|3.5 TB   |LexicMap|2.91 TB|16 h 40 m|79 GB |
|                 |         |         |Blastn  |2.15 TB|14 h 4 m |4.3 GB|
|AllTheBacteria HQ|1,858,610|3.1 TB   |LexicMap|2.32 TB|10 h 48 m|43 GB |
|                 |         |         |Blastn  |1.76 TB|14 h 3 m |2.9 GB|
|                 |         |         |Phylign |248 GB |/        |/     |

### Searching

Blastn failed to run as it requires >2000GB RAM for Genbank+RefSeq and AllTheBacteria datasets.
Phylign only has the index for AllTheBacteria HQ dataset.


GTDB complete (402,538 genomes):

|query          |query_len|tool    |genome_hits|time     |RAM   |
|:--------------|--------:|:-------|----------:|--------:|-----:|
|a marker gene  |1,299 bp |LexicMap|3,588      |1.4 s    |598 MB|
|               |         |Blastn  |7,121      |36 m 11 s|351 GB|
|a 16S rRNA gene|1,542 bp |LexicMap|294,285    |1 m 19 s |3.0 GB|
|               |         |Blastn  |301,197    |39 m 13 s|378 GB|
|a plasmid      |52,830 bp|LexicMap|58,930     |40 s     |3.2 GB|
|               |         |Blastn  |69,311     |37 m 42 s|365 GB|

AllTheBacteria HQ (1,858,610 genomes):

|query          |query_len|tool           |genome_hits|time     |RAM    |
|:--------------|--------:|:--------------|----------:|--------:|------:|
|a marker gene  |1,299 bp |LexicMap       |10,837     |11.3 s   |1.1 GB |
|               |         |Phylign_local  |7,937      |29 m 13 s|78.5 GB|
|               |         |Phylign_cluster|7,937      |32 m 26 s|/      |
|a 16S rRNA gene|1,542 bp |LexicMap       |1,853,846  |12 m 31 s|9.7 GB |
|               |         |Phylign_local  |1,032,948  |1 h 58 m |73.3 GB|
|               |         |Phylign_cluster|1,032,948  |83 m 36 s|/      |
|a plasmid      |52,830 bp|LexicMap       |427,112    |4 m 13 s |12.9 GB|
|               |         |Phylign_local  |46,822     |38 m 55 s|76.6 GB|
|               |         |Phylign_cluster|46,822     |36 m 30s |/      |

Genbank+RefSeq (2,340,672 genomes):

|query          |query_len|tool    |genome_hits|time    |RAM    |
|:--------------|--------:|:-------|----------:|-------:|------:|
|a marker gene  |1,299 bp |LexicMap|16,556     |6.2 s   |1.3 GB |
|a 16S rRNA gene|1,542 bp |LexicMap|1,875,260  |8 m 29 s|10.8 GB|
|a plasmid      |52,830 bp|LexicMap|494,860    |4 m 08 s|15.0 GB|

Notes:
- All files are stored on a server with HDD disks. No files are cached in memory.
- Tests are performed in a single cluster node with 48 CPU cores (Intel Xeon Gold 6336Y CPU @ 2.40 GHz).
- LexicMap index building parameters: `-k 31 -m 40000`. Genome batch size: `-b 10000` for GTDB datasets, `-b 50000` for others.
- Main searching parameters:
    - LexicMap v0.3.0: `--threads 48 --top-n-genomes 0 --min-qcov-per-genome 0 --min-qcov-per-hsp 0 --min-match-pident 50`.
    - Blastn v2.15.0+: `-num_threads 48 -max_target_seqs 10000000`.
    - Phylign (AllTheBacteria fork 9fc65e6): `threads: 48, cobs_kmer_thres: 0.33, nb_best_hits: 5000000, max_ram_gb: 100`; For cluster: maximum number of slurm jobs is `100`.

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
- [Wavefront alignment algorithm (WFA) in Golang](https://github.com/shenwei356/wfa).

## Support

Please [open an issue](https://github.com/shenwei356/LexicMap/issues) to report bugs,
propose new functions or ask for help.

## License

[MIT License](https://github.com/shenwei356/LexicMap/blob/master/LICENSE)

