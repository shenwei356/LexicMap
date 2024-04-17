---
title: Introduction
weight: 10
---

LexicMap is a sequence alignment tool aiming to efficiently query gene or plasmid sequences against up to millions of prokaryotic genomes.

For example, **querying a 51.5-kb plasmid in all 2,340,672 Genbank+Refseq prokaryotic genomes takes only 3 minutes and 32 seconds with 15.7 GB RAM and 48 CPUs, with 19,265 genome hits returned**.
By contrast, BLASTN is unable to run with the same dataset on common servers because it requires >2000 GB RAM. See [performance](#performance).

LexicMap uses a modified [LexicHash](https://doi.org/10.1093/bioinformatics/btad652) algorithm, which supports variable-length substring matching rather than classical fixed-length k-mers matching, to compute seeds for sequence alignment and uses multiple-level storage for fast and low-memory quering of seeds data. See [algorithm overview](#algorithm-overview).

LexicMap is also very easy to [install](http://bioinf.shenwei.me/LexicMap/installation/) (a binary file with no dependencies) and use ([tutorials](http://bioinf.shenwei.me/LexicMap/tutorials/index/) and [usages](http://bioinf.shenwei.me/LexicMap/usage/lexicmap/)).


## Performance

{{< tabs "uniqueid" >}}
{{< tab "GTDB repr" >}}

### Index information

|dataset          |genomes  |gzip_size|db_size|indexing_time|indexing_RAM|
|:----------------|--------:|--------:|------:|------------:|-----------:|
|GTDB repr        |85,205   |75 GB    |110 GB |1 h 25 m     |49 GB       |


### Query performance

|query          |query_len|genome_hits|time    |RAM    |
|:--------------|--------:|----------:|-------:|------:|
|a MutL gene    |1,956 bp |2          |0.9 s   |460 MB |
|a 16S rRNA gene|1,542 bp |13,466     |4.0 s   |765 MB |
|a plasmid      |51,466 bp|2          |1.1 s   |752 MB |
{{< /tab >}}

{{< tab "GTDB complete" >}}


### Index information

|dataset          |genomes  |gzip_size|db_size|indexing_time|indexing_RAM|
|:----------------|--------:|--------:|------:|------------:|-----------:|
|GTDB complete    |402,538  |578 GB   |507 GB |5 h 18 m     |46 GB       |

### Query performance

|query          |query_len|genome_hits|time    |RAM    |
|:--------------|--------:|----------:|-------:|------:|
|a MutL gene    |1,956 bp |268        |3.8 s   |544 MB |
|a 16S rRNA gene|1,542 bp |169,480    |2 m 14 s|2.9 GB |
|a plasmid      |51,466 bp|3,649      |56 s    |2.9 GB |

{{< /tab>}}


{{< tab "Genbank+RefSeq" >}}

### Index information

|dataset          |genomes  |gzip_size|db_size|indexing_time|indexing_RAM|
|:----------------|--------:|--------:|------:|------------:|-----------:|
|Genbank+RefSeq   |2,340,672|3.5 TB   |2.9 TB |31 h 19 m    |148 GB      |

### Query performance

|query          |query_len|genome_hits|time    |RAM    |
|:--------------|--------:|----------:|-------:|------:|
|a MutL gene    |1,956 bp |817        |10.0 s  |2.3 GB |
|a 16S rRNA gene|1,542 bp |1,148,049  |5 m 34 s|11.8 GB|
|a plasmid      |51,466 bp|19,265     |3 m 32 s|15.7 GB|

{{< /tab>}}


{{< tab "AllTheBacteria HQ" >}}

### Index information

|dataset          |genomes  |gzip_size|db_size|indexing_time|indexing_RAM|
|:----------------|--------:|--------:|------:|------------:|-----------:|
|AllTheBacteria HQ|1,858,610|3.1 TB   |2.4 TB |16 h 24 m    |70 GB       |

### Query performance

|query          |query_len|genome_hits|time    |RAM    |
|:--------------|--------:|----------:|-------:|------:|
|a MutL gene    |1,956 bp |404        |18.7 s  |2.4 GB |
|a 16S rRNA gene|1,542 bp |1,193,874  |13 m 8 s|9.4 GB |
|a plasmid      |51,466 bp|10,954     |5 m 25 s|9.7 GB |

{{< /tab>}}


{{< /tabs >}}


Notes:
- All files are stored on a server with HDD disks.
- Tests are performed in a single cluster node with 48 CPU cores (Intel Xeon Gold 6336Y CPU @ 2.40 GHz).
- Index building parameters: `-k 31 -m 40000`. Genome batch size: `-b 10000` for GTDB datasets, `-b 131072` for others.

## Algorithm overview

<img src="/LexicMap/overview.svg" alt="LexicMap overview" width="900"/>


## Related projects

- High-performance [LexicHash](https://github.com/shenwei356/LexicHash) computation in Go.

## Support

Please [open an issue](https://github.com/shenwei356/LexicMap/issues) to report bugs,
propose new functions or ask for help.

## License

[MIT License](https://github.com/shenwei356/LexicMap/blob/master/LICENSE)
