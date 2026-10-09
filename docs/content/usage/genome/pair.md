---
title: lexicmap genome pair
linkTitle: pair
weight: 20
---

## Usage

```plain
$ lexicmap utils 2blast -h
Find similar genome pairs in the index

Output format:
  Tab-delimited format with 7 columns.

    1.  genome1,    Genome 1.
    2.  genome2,    Genome 2.
    3.  minPrefix,  Minimum common prefix length between two seeds (-p/--min-prefix).
    4.  fracMasks,  Fraction of masks with seeds sharing a common prefix >= minPrefix.
    5.  nMasks,     Number   of masks with seeds sharing a common prefix >= minPrefix.
    6.  sumPrefix,  Total length of common prefixes, with only the best seed pair
                    (with the longest common prefix) considered for each mask.
    7.  avgPrefix,  Average prefix length (sumPrefix / nMasks).

Limitations:
  1. Not suitable for very large genome sets (e.g., > 1 million genomes) due to memory constraints.
  2. Genome sets with many highly similar genomes may result in a large number of genome pairs, 
     which can require substantial time and memory.
  3. Genomes stored in multiple chunks are not evaluated as a whole.

Notes on pair-table memory:
  1. For indexes with no more than --max-dense-genomes genomes, compact dense pair
     tables are used. They are faster and use less memory when many genome pairs
     match, but their size grows quadratically with the number of genomes. Set the
     value to 0 to force the sparse-map implementation.
  2. Raising --max-dense-genomes enables dense tables for a larger index and can
     greatly increase baseline memory. Lowering it can save memory for indexes in
     which relatively few genome pairs match. At the default limit of 20,000
     genomes, the global tables need about 1.14 GiB and each mask worker needs
     about 215 MiB; the default 768 MiB budget therefore permits up to three mask
     workers. These values exclude k-mer windows and final sorting results.
  3. --dense-mask-memory limits only the temporary tables used by concurrent mask
     workers; it is not a total RSS limit. One table is allocated up front for each
     permitted worker. Raising the value may increase parallelism, but speed usually
     has diminishing returns and can even regress because workers compete for memory
     bandwidth and CPU cache, while completed masks are merged into the global dense
     table by one collector goroutine. The global pair statistics, k-mer windows, and
     final sorting results are additional memory.

Notes on probabilistic pruning:
  1. -f/--min-mask-fraction is the final reporting threshold. In contrast,
     -s/--prob-threshold controls an early heuristic that discards pairs unlikely
     to reach that threshold. Higher values prune more aggressively; 0 disables
     the heuristic and computes exact pair statistics.
  2. In sparse-map mode, pruning can substantially reduce the number of map entries
     and therefore memory use. In dense-table mode, all pair slots are allocated
     up front, so pruning does not reduce the main table memory or skip pair
     comparisons within a mask. It can still save some time by avoiding updates
     to discarded pairs and by keeping the active result set smaller.
  3. The default -s 0.001 is retained for compatibility and sparse-map workloads.
     With -s 0, results are independent of mask completion order, so changing
     --dense-mask-memory affects only memory and parallelism. With -s > 0, changing
     concurrency can alter the completion order and may affect marginal pairs
     considered by the pruning heuristic.

Usage:
  lexicmap genome pair [flags] -d <index path> [-o out.tsv.gz]

Flags:
      --dense-mask-memory int     ► Memory budget in MiB for concurrent temporary dense-mask tables
                                  (not total RSS). Bigger values permit more mask workers but may have
                                  diminishing returns because of memory-bandwidth contention and serial
                                  global merging. See the notes above. (default 768)
  -h, --help                      help for pair
  -d, --index string              ► Index directory created by "lexicmap index".
  -m, --masks int                 ► Only use seed data of N masks. It should be 0 (for all masks in
                                  the index) or power of 4 (needs to be >= 64, e.g., 64, 256, 1024,
                                  4096, 16384). (default 1024)
      --max-dense-genomes int     ► Maximum number of genomes for using compact dense pair tables
                                  (maximum: 92682; 0 disables them). Bigger values can greatly increase
                                  memory usage. (default 20000)
  -f, --min-mask-fraction float   ► Minimum fraction of masks that must match for a genome pair to be
                                  reported. (default 0.25)
  -p, --min-prefix int            ► Minimum prefix length between k-mers captured by a mask. (default 21)
  -o, --out-file string           ► Out file, supports and recommends a ".gz" suffix ("-" for stdout).
                                  (default "-")
  -s, --prob-threshold float      ► Probability threshold for early pruning (higher = more aggressive;
                                  0 disables pruning and computes exact pair statistics). (default 0.001)

Global Flags:
  -X, --infile-list string   ► File of input file list (one file per line). If given, they are
                             appended to files from CLI arguments.
      --log string           ► Log file.
      --quiet                ► Do not print any verbose information. But you can write them to a file
                             with --log.
  -j, --threads int          ► Number of CPU cores to use. By default, it uses all available cores.
                             (default 16)
```

Examples:

    $ lexicmap genome pair -d demo.lmi/ \
        | csvtk pretty -t
    
    genome1           genome2           minPrefix   fracMasks   nMasks   sumPrefix   avgPrefix
    ---------------   ---------------   ---------   ---------   ------   ---------   ---------
    GCF_002949675.1   GCF_002950215.1   21          0.6436      659      20178       30.62    
    GCF_003697165.2   GCF_002950215.1   21          0.6045      619      18834       30.43    
    GCF_003697165.2   GCF_002949675.1   21          0.5479      561      17081       30.45    
