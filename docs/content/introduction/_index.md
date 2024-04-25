---
title: Introduction
weight: 10
---

LexicMap is a sequence alignment tool aiming to efficiently query gene or plasmid sequences against up to millions of prokaryotic genomes.

For example, **querying a 51.5-kb plasmid in all 2,340,672 Genbank+Refseq prokaryotic genomes takes only 5 minutes and 2 seconds with 13.7 GB RAM and 48 CPUs, with 17,822 genome hits returned**.
By contrast, BLASTN is unable to run with the same dataset on common servers because it requires >2000 GB RAM. See [performance](#performance).

LexicMap uses a modified [LexicHash](https://doi.org/10.1093/bioinformatics/btad652) algorithm, which supports variable-length substring matching rather than classical fixed-length k-mers matching, to compute seeds for sequence alignment and uses multiple-level storage for fast and low-memory quering of seeds data. See [algorithm overview](#algorithm-overview).

LexicMap is also very easy to [install](http://bioinf.shenwei.me/LexicMap/installation/) (a binary file with no dependencies) and use ([tutorials](http://bioinf.shenwei.me/LexicMap/tutorials/index/) and [usages](http://bioinf.shenwei.me/LexicMap/usage/lexicmap/)).


## Performance

{{< tabs "uniqueid" >}}
{{< tab "GTDB repr" >}}

### Index information

|dataset          |genomes  |gzip_size|db_size|indexing_time|indexing_RAM|
|:----------------|--------:|--------:|------:|------------:|-----------:|
|GTDB repr        |85,205   |75 GB    |110 GB |1 h 30 m     |38 GB       |


### Query performance

|query          |query_len|genome_hits|time     |RAM    |
|:--------------|--------:|----------:|--------:|------:|
|a MutL gene    |1,956 bp |2          |3.2 s    |436 MB |
|a 16S rRNA gene|1,542 bp |4,374      |38.5 s   |747 MB |
|a plasmid      |51,466 bp|0          |13.0 s   |768 MB |

{{< /tab >}}

{{< tab "GTDB complete" >}}


### Index information

|dataset          |genomes  |gzip_size|db_size|indexing_time|indexing_RAM|
|:----------------|--------:|--------:|------:|------------:|-----------:|
|GTDB complete    |402,538  |578 GB   |510 GB |3 h 26 m     |35 GB       |

### Query performance

|query          |query_len|genome_hits|time     |RAM    |
|:--------------|--------:|----------:|--------:|------:|
|a MutL gene    |1,956 bp |268        |2.8 s    |571 MB |
|a 16S rRNA gene|1,542 bp |107,557    |3 m 38 s |2.6 GB |
|a plasmid      |51,466 bp|3,220      |56.2 s   |3.0 GB |

{{< /tab>}}


{{< tab "Genbank+RefSeq" >}}

### Index information

|dataset          |genomes  |gzip_size|db_size |indexing_time|indexing_RAM|
|:----------------|--------:|--------:|-------:|------------:|-----------:|
|Genbank+RefSeq   |2,340,672|3.5 TB   |2.91 TB |16 h 40 m    |79 GB       |

### Query performance

|query          |query_len|genome_hits|time     |RAM    |
|:--------------|--------:|----------:|--------:|------:|
|a MutL gene    |1,956 bp |817        |6.0 s    |1.4 GB |
|a 16S rRNA gene|1,542 bp |832,161    |18 m 58 s|8.3 GB |
|a plasmid      |51,466 bp|17,822     |5 m 02 s |13.7 GB|

{{< /tab>}}


{{< tab "AllTheBacteria HQ" >}}

### Index information

|dataset          |genomes  |gzip_size|db_size |indexing_time|indexing_RAM|
|:----------------|--------:|--------:|-------:|------------:|-----------:|
|AllTheBacteria HQ|1,858,610|3.1 TB   |2.32 TB |10 h 48 m    |43 GB       |

### Query performance

|query          |query_len|genome_hits|time     |RAM    |
|:--------------|--------:|----------:|--------:|------:|
|a MutL gene    |1,956 bp |404        |4.7 s    |1.1 GB |
|a 16S rRNA gene|1,542 bp |1,031,705  |17 m 54 s|8.4 GB |
|a plasmid      |51,466 bp|10,897     |4 m 07 s |10.8 GB|



{{< /tab>}}


{{< /tabs >}}


Notes:
- All files are stored on a server with HDD disks. No files are cached in memory.
- Tests are performed in a single cluster node with 48 CPU cores (Intel Xeon Gold 6336Y CPU @ 2.40 GHz).
- Index building parameters: `-k 31 -m 40000`. Genome batch size: `-b 10000` for GTDB datasets, `-b 50000` for others.
- Searching parameters: `--top-n-genomes 0 --min-qcov-per-genome 50 --min-match-pident 70 --min-qcov-per-hsp 0`.

## Algorithm overview

<img src="/LexicMap/overview.svg" alt="LexicMap overview" width="900"/>


## Related projects

- High-performance [LexicHash](https://github.com/shenwei356/LexicHash) computation in Go.

## Support

Please [open an issue](https://github.com/shenwei356/LexicMap/issues) to report bugs,
propose new functions or ask for help.

## License

[MIT License](https://github.com/shenwei356/LexicMap/blob/master/LICENSE)
