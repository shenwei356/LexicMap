---
title: FAQs
weight: 60
---
## Table of contents

{{< toc format=html >}}

## Does LexicMap support short reads?

LexicMap is mainly designed for sequence alignment with a small number of queries (gene/plasmid/virus/phage sequences) longer than 100 bp by default.

If you just want to search long (>1kb) queries for highly similar (>95%) targets, you can build an index with a bigger `-D/--seed-max-desert` (default 100) and `-d/--seed-in-desert-dist` (default 50), e.g.,

    --seed-max-desert 300 --seed-in-desert-dist 150

Bigger values decrease the search sensitivity for distant targets, speed up the indexing
speed, decrease the indexing memory occupation and decrease the index size. While the
alignment speed is almost not affected.


## Does LexicMap support fungi genomes?

Yes. LexicMap mainly supports small genomes including prokaryotic, viral, and plasmid genomes.
**Fungi can also be supported, just remember to increase the value of `-g/--max-genome` when running `lexicmap index`,
which is used to skip genomes larger than 15Mb by default**.

```
  -g, --max-genome int            ► Maximum genome size. Extremely large genomes (e.g., non-isolate
                                  assemblies from Genbank) will be skipped. (default 15000000)
```

Maximum genome size is about 268 Mb (268,435,456). More precisely:

    $total_bases + ($num_contigs - 1) * 1000 <= 268,435,456

as we concatenate contigs with 1000-bp intervals of N’s to reduce the sequence scale to index.

For big and complex genomes, like the human genome (chr1 is ~248 Mb) which has many repetitive sequences, LexicMap would be slow to align.


## How's the hardware requirement?

- For index building. See details [hardware requirement](https://bioinf.shenwei.me/LexicMap/tutorials/index/#hardware-requirements).
- For seaching. See details [hardware requirement](https://bioinf.shenwei.me/LexicMap/tutorials/search/#hardware-requirements).



## How to resume the indexing as Slurm job time limit is almost reached while lexicmap index is still in the merging step?

Use [lexicmap utils remerge](https://bioinf.shenwei.me/LexicMap/usage/utils/remerge/) (available since v0.5.0), which reruns the merging step for an unfinished index.

> When to use this command?
> - Only one thread is used for merging indexes, which happens when there are
>  a lot (>200 batches) of batches (`$inpu_files / --batch-size`) and the value
>  of `--max-open-files` is not big enough.
> - The Slurm/PBS job time limit is almost reached and the merging step won't be finished before that.
> - Disk quota is reached in the merging step.

So you can stop the indexing command by press `Ctrl` + `C` (**make sure it is in the merging step**, see example below), and run `lexicmap utils remerge -d index.lmi`,
where `index.lmi` is the output index directory in `lexicmap index`.

Optionally, you might set bigger values of
flag `--max-open-files` and `-J/--seed-data-threads` if you have hundreds of thousands of input genomes or have set
a small batch size with `-b/--batch-size`. E.g.,

    22:54:24.420 [INFO] merging 297 indexes...
    22:54:24.455 [INFO]   [round 1]
    22:54:24.455 [INFO]     batch 1/1, merging 297 indexes to xxx.lmi.tmp/r1_b1 with 1 threads...

There's only one thread was used for seed data merging, it would take a long time.
So we can set a larger `--max-open-files`, e.g., `4096`,
and it would allow `4096 / (297+2) = 13.7` threads for merging, let's set `--seed-data-threads 12`.

    # specify the maximum open files per process
    ulimit -n 4096

    lexicmap utils remerge -d index.lmi --max-open-files 4096 --seed-data-threads 12


## Can I extract the matched sequences?

Yes, `lexicmap search` has a flag

```
  -a, --all                            ► Output more columns, e.g., matched sequences. Use this if you
                                       want to output blast-style format with "lexicmap utils 2blast".
```

to output CIGAR string, aligned query and subject sequences.

```
20. cigar,    CIGAR string of the alignment.                      (optional with -a/--all)
21. qseq,     Aligned part of query sequence.                     (optional with -a/--all)
22. sseq,     Aligned part of subject sequence.                   (optional with -a/--all)
23. align,    Alignment text ("|" and " ") between qseq and sseq. (optional with -a/--all)
```

An example:

    # Extracting similar sequences for a query gene.

    # search matches with query coverage >= 90%
    lexicmap search -d gtdb_complete.lmi/ b.gene_E_faecalis_SecY.fasta -o results.tsv \
        --min-qcov-per-hsp 90 --all

    # extract matched sequences as FASTA format
    sed 1d results.tsv | awk -F'\t' '{print ">"$5":"$14"-"$15":"$16"\n"$22;}' \
        | seqkit seq -g > results.fasta

    seqkit head -n 1 results.fasta | head -n 3
    >NZ_JALSCK010000007.1:39224-40522:-
    TTGTTCAAGCTATTAAAGAACGCCTTTAAAGTCAAAGACATTAGATCAAAAATCTTATTT
    ACAGTTTTAATCTTGTTTGTATTTCGCCTAGGTGCGCACATTACTGTGCCCGGGGTGAAT


And `lexicmap util 2blast` can help to convert the tabular format to Blast-style format,
see [examples](https://bioinf.shenwei.me/LexicMap/usage/utils/2blast/#examples).

## How can I extract the upstream and downstream flanking sequences of matched regions?

[lexicmap utils subseq](https://bioinf.shenwei.me/LexicMap/usage/utils/subseq/)
can extract subsequencess via genome ID, sequence ID and positions.
So you can use these information from the search result and expand the region positions to extract flanking sequences.



## Why isn't the pident 100% when aligning with a sequence from the reference genomes?

It happens if there are some degenerate bases (e.g., `N`) in the query sequence.
In the indexing step, all degenerate bases are converted to their lexicographic first bases. E.g., `N` is converted to `A`.
While for the query sequences, we don't convert them.


## Why is LexicMap slow for batch searching?

LexicMap is mainly designed for sequence alignment with a small number of queries against a database with a huge number (up to 17 million) of genomes.

There are some ways to improve the search speed of `lexicmap search`: 
http://bioinf.shenwei.me/LexicMap/tutorials/search/#improving-searching-speed

{{< button relref="/usage/search"  >}}Click{{< /button >}}  to read more detail of the usage.

