---
title: kmers
weight: 10
---

```plain
$ lexicmap utils kmers -h
View k-mers captured by the masks

Attention:
  1. Mask index (column mask) is 1-based.
  2. Prefix means the length of shared prefix between a k-mer and the mask.
  3. K-mer positions (column pos) are 1-based.
     For reference genomes with multiple sequences, the sequences were
     concatenated to a single sequence with intervals of N's.
  4. Reversed means if the k-mer is reversed for suffix matching.

Usage:
  lexicmap utils kmers [flags] -d <index path> [-m <mask index>] [-o out.tsv.gz]

Flags:
  -h, --help              help for kmers
  -d, --index string      ► Index directory created by "lexicmap index".
  -m, --mask int          ► View k-mers captured by Xth mask. (0 for all) (default 1)
  -f, --only-forward      ► Only output forward k-mers.
  -o, --out-file string   ► Out file, supports and recommends a ".gz" suffix ("-" for stdout).
                          (default "-")

Global Flags:
  -X, --infile-list string   ► File of input file list (one file per line). If given, they are
                             appended to files from CLI arguments.
      --log string           ► Log file.
      --quiet                ► Do not print any verbose information. But you can write them to a file
                             with --log.
  -j, --threads int          ► Number of CPU cores to use. By default, it uses all available cores.
                             (default 16)
```

## Examples

1. The default output is captured k-mers of the first mask.

        $ lexicmap utils kmers --quiet -d demo.lmi/ | head -n 20 | csvtk pretty -t
        mask   kmer                              prefix   number   ref               pos       strand   reversed
        ----   -------------------------------   ------   ------   ---------------   -------   ------   --------
        1      AAAAAAAAACAAACATTTGCGGCGGGGCCAT   8        1        GCF_000742135.1   2043044   +        no
        1      AAAAAAAAACCAGAAATCACACGCCAACTCC   8        1        GCF_002949675.1   1345415   +        yes
        1      AAAAAAAAACGATTATCCTCAATTAATTTCT   8        1        GCF_000392875.1   814251    +        no
        1      AAAAAAAAACGCTTCTACATCGAGCAGCGAG   8        1        GCF_001457655.1   941619    +        yes
        1      AAAAAAAAACGCTTTGTAACTCGATTGATAG   8        1        GCF_009759685.1   997945    +        yes
        1      AAAAAAAAACTGCTGTCCCTGGTCCGTCAGG   8        1        GCF_002950215.1   4262890   -        yes
        1      AAAAAAAAAGATTTGATTTTTTTCATTAATA   8        1        GCF_000392875.1   766998    -        yes
        1      AAAAAAAAAGCATTTTTTCGATCTCTTTACG   8        1        GCF_000392875.1   1623731   +        yes
        1      AAAAAAAAAGTTTCCGGGACACTACCTAACC   8        1        GCF_000017205.1   5804200   -        yes
        1      AAAAAAAAATTATTTTGCTAATCAATAGGTC   8        1        GCF_000006945.2   4886411   -        yes
        1      AAAAAAAACAAAGAATTATTACACAACATTC   8        1        GCF_003697165.2   4055655   +        yes
        1      AAAAAAAACACGGACTTATTGAAATCGTATT   8        1        GCF_000392875.1   746746    +        yes
        1      AAAAAAAACCAACTTTGAAAAAAGTAATGTA   8        1        GCF_000148585.2   917529    -        yes
        1      AAAAAAAACCATATTATGTCCGATCCTCACA   8        1        GCF_000392875.1   1060650   +        yes
        1      AAAAAAAACCCGCCGAAGCGGGTTTTTTTAT   8        1        GCF_000742135.1   1612499   +        no
        1      AAAAAAAACCTAATGGTAAATAACGTTTTGG   8        1        GCF_006742205.1   2346818   +        yes
        1      AAAAAAAACGAAAAACGGTAACACGGGAATT   8        1        GCF_001544255.1   1605298   +        yes
        1      AAAAAAAACGACTCCAGAGAGATCATCGTAT   8        1        GCF_000392875.1   1279686   +        yes
        1      AAAAAAAACGAGTCATTTCCCCTACTGAACC   8        2        GCF_002949675.1   2284659   -        yes

    Only forward k-mers.

        $ lexicmap utils kmers --quiet -d demo.lmi/ -f | head -n 20 | csvtk pretty -t
        mask   kmer                              prefix   number   ref               pos       strand   reversed
        ----   -------------------------------   ------   ------   ---------------   -------   ------   --------
        1      AAAAAAAAACAAACATTTGCGGCGGGGCCAT   8        1        GCF_000742135.1   2043044   +        no
        1      AAAAAAAAACGATTATCCTCAATTAATTTCT   8        1        GCF_000392875.1   814251    +        no
        1      AAAAAAAACCCGCCGAAGCGGGTTTTTTTAT   8        1        GCF_000742135.1   1612499   +        no
        1      AAAAAAAACGGTTCAGCTGACCAGCCAGCTG   8        1        GCF_002950215.1   401140    +        no
        1      AAAAAAAAGATATTGAAGTTAAAGTAATTTG   9        1        GCF_000742135.1   3038258   +        no
        1      AAAAAAAAGCCCACGAACCGGGGGCAATATC   9        1        GCF_002950215.1   3578394   +        no
        1      AAAAAAAAGCCCCGCCGAAGCGGGGCTTTTT   9        1        GCF_000017205.1   5110420   +        no
        1      AAAAAAAAGGATTATAACAAAATTTTGTCAT   9        1        GCF_001544255.1   426716    +        no
        1      AAAAAAAAGTAATTGCAGCTATTATTGGGAC   10       1        GCF_001027105.1   437272    +        no
        1      AAAAAAAAGTATTAAGCAACTGACTAAAAGT   10       1        GCF_006742205.1   1841209   +        no
        1      AAAAAAAAGTCACAATTATTGGTGCCGGTTT   13       1        GCF_000392875.1   1508457   -        no
        1      AAAAAAAAGTCATCAAGGATTATTTGAGTTA   12       1        GCF_001457655.1   1847867   +        no
        1      AAAAAAAAGTCATCGCTTTATCTGTCAGTAT   12       1        GCF_001544255.1   156689    -        no
        1      AAAAAAAAGTCATCTTCGGATGGCTTTTTTA   12       1        GCF_000148585.2   1363150   -        no
        1      AAAAAAAAGTCCATCCTGCAGCATAAAATAA   11       1        GCF_000742135.1   4671015   +        no
        1      AAAAAAAAGTCCCTGCTGTTTGCCCAGTCCT   11       1        GCF_000006945.2   3796      -        no
        1      AAAAAAAAGTCCGCTGATAAGGCTTGAAAAG   11       3        GCF_002949675.1   2356807   +        no
        1      AAAAAAAAGTCCGCTGATAAGGCTTGAAAAG   11       3        GCF_002950215.1   3051946   +        no
        1      AAAAAAAAGTCCGCTGATAAGGCTTGAAAAG   11       3        GCF_003697165.2   16156     +        no


1. Specify the mask.

        $ lexicmap utils kmers --quiet -d demo.lmi/ --mask 12345 | head -n 20 | csvtk pretty -t
        mask    kmer                              prefix   number   ref               pos       strand   reversed
        -----   -------------------------------   ------   ------   ---------------   -------   ------   --------
        12345   CATTAGTAAAAACCAACTTAGTTACGACACG   8        1        GCF_001027105.1   1823411   +        no
        12345   CATTAGTAAAACATTTTGAACCTGTGATTGA   8        1        GCF_006742205.1   1192019   +        no
        12345   CATTAGTAAAAGTCGTTTGGTAAAGCGATTA   8        1        GCF_001027105.1   1334989   +        yes
        12345   CATTAGTAAACGTACAAAACTATTGGTTAGA   8        1        GCF_001027105.1   2037559   +        yes
        12345   CATTAGTAAATCCAGGAATCCTAACCGACGA   8        1        GCF_001027105.1   963152    +        yes
        12345   CATTAGTAACGCGTACGAAACCGTAGTAAGT   8        1        GCF_001027105.1   1958187   +        yes
        12345   CATTAGTAAGTTGTCGGTCTAACGCGGATTA   8        1        GCF_002950215.1   2882180   +        yes
        12345   CATTAGTACATTCAAGTATTATTCATTAAAC   8        1        GCF_009759685.1   665376    +        yes
        12345   CATTAGTACCGATAGGACATCATGAACACAA   8        1        GCF_002950215.1   4677222   +        yes
        12345   CATTAGTACCTTCATCGCTATCCCATTAGGC   8        1        GCF_000006945.2   92542     +        yes
        12345   CATTAGTACGTGTCCCGCAAAGAGAAAGAAC   8        1        GCF_000006945.2   3412102   +        yes
        12345   CATTAGTAGAAAAATACAAAGGCATTTATGA   11       1        GCF_900638025.1   665985    -        no
        12345   CATTAGTAGAAAATTGATAATCTAAGAGTTC   11       1        GCF_002950215.1   2940281   +        no
        12345   CATTAGTAGAAATGGGCAAAGAATAGGAAAA   11       1        GCF_000148585.2   81286     +        no
        12345   CATTAGTAGAAGAAATTGCAGCAAGTATTAA   14       1        GCF_001027105.1   621160    +        no
        12345   CATTAGTAGAAGAACTGAAGTTAGTGCCTAT   14       1        GCF_001096185.1   2113047   +        no
        12345   CATTAGTAGAAGAAGACCAAGCACGACGCAT   15       1        GCF_000392875.1   891723    +        no
        12345   CATTAGTAGAAGAGTTGTTCGTCAGTTACGG   13       1        GCF_001544255.1   831068    -        no
        12345   CATTAGTAGAAGATTTAGTGGCAAGCTCAAT   13       1        GCF_001457655.1   1280653   +        no

    "reversed" means means if the k-mer is reversed for suffix matching.
    E.g., `CATTAGTAAAAGTCGTTTGGTAAAGCGATTA` is reversed, so you need to reverse it before searching in the genome.


        $ seqkit locate -p $(echo CATTAGTAAAAGTCGTTTGGTAAAGCGATTA | rev) refs/GCF_001027105.1.fa.gz -M | csvtk pretty -t
        seqID           patternName                       pattern                           strand   start     end
        -------------   -------------------------------   -------------------------------   ------   -------   -------
        NZ_CP011526.1   ATTAGCGAAATGGTTTGCTGAAAATGATTAC   ATTAGCGAAATGGTTTGCTGAAAATGATTAC   +        1334989   1335019


1. For all masks. The result might be very big, therefore, writing to gzip format is recommended.


        $ lexicmap utils kmers -d demo.lmi/ --mask 0 -o kmers.tsv.gz

        $ zcat kmers.tsv.gz | csvtk freq -t -f mask -nr | head -n 10
        mask    frequency
        24088   322
        15814   295
        13923   293
        27102   291
        13922   282
        15967   281
        10001   280
        15986   272
        16440   269

    a faster way

        seq 1 $(lexicmap utils masks -d demo.lmi/ --quiet | wc -l) \
            | rush --eta 'echo -e {}"\t"$(lexicmap utils kmers -d demo.lmi/ -m {} -f --quiet | csvtk nrow)' \
            | csvtk add-header -t -n mask,seeds \
            | csvtk sort -t -k seeds:nr \
            | head -n 10


1. Lengths of shared prefixes between probes and captured k-mers.

        zcat kmers.tsv.gz \
          | csvtk grep -t -f reversed -p  no \
          | csvtk plot hist -t -f prefix -o prefix.hist.png \
              --xlab "length of common prefixes between captured k-mers and masks"


    <img src="/LexicMap/prefix.hist.png" alt="" width="400"/>

The output (TSV format) is formatted with [csvtk pretty](https://github.com/shenwei356/csvtk).
