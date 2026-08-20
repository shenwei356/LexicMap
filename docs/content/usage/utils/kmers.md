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
        1      AAAAAAAAACATTTTTTTCAAACACCTCCTT   7        1        GCF_001544255.1   109698    +        yes     
        1      AAAAAAAAACGAAAAAGATTTTCCCTCATAC   7        1        GCF_000392875.1   2088530   +        yes     
        1      AAAAAAAAACGCTTCTACATCGAGCAGCGAG   7        1        GCF_001457655.1   941619    +        yes     
        1      AAAAAAAAACGTATCCCTCTTTATTACTTAT   7        1        GCF_000006945.2   3392260   -        yes     
        1      AAAAAAAAACTAGGGTTAAATGCCTTATGTT   7        1        GCF_009759685.1   442423    +        yes     
        1      AAAAAAAAACTCAGCTAGACGTGGCCGTCGT   7        1        GCF_000017205.1   1108195   +        yes     
        1      AAAAAAAAAGATTTGATTTTTTTCATTAATA   7        1        GCF_000392875.1   766998    -        yes     
        1      AAAAAAAAAGGGGTTAGCGCGCAATTGCGGA   7        1        GCF_002950215.1   337801    -        yes     
        1      AAAAAAAAATCTATTTTAAAACCTAATCACG   7        1        GCF_000392875.1   2201506   +        yes     
        1      AAAAAAAAATGTCACAACAGCCCAACCTCCA   7        1        GCF_000392875.1   860216    +        yes     
        1      AAAAAAAACAAAAACTAGTTCGAGTGCCGAA   7        1        GCF_000006945.2   1587885   -        yes     
        1      AAAAAAAACACCGACGGGGAGTCTCCTCTTT   7        1        GCF_000742135.1   4064602   -        yes     
        1      AAAAAAAACCATATTATGTCCGATCCTCACA   7        1        GCF_000392875.1   1060650   +        yes     
        1      AAAAAAAACGAAAAACGGTAACACGGGAATT   7        1        GCF_001544255.1   1605298   +        yes     
        1      AAAAAAAACGACGCAGAAAACGACATTGCGA   7        1        GCF_003697165.2   564733    +        yes     
        1      AAAAAAAACGACTCCAGAGAGATCATCGTAT   7        1        GCF_000392875.1   1279686   +        yes     
        1      AAAAAAAACGAGCGATTGGTTGCATTAAGGA   7        1        GCF_002949675.1   3914985   -        yes     
        1      AAAAAAAACGAGCGCTCGGTTGCATTAAGGA   7        2        GCF_003697165.2   1514669   -        yes     
        1      AAAAAAAACGAGCGCTCGGTTGCATTAAGGA   7        2        GCF_002949675.1   2061956   -        yes     

    Only forward k-mers.

        $ lexicmap utils kmers --quiet -d demo.lmi/ -f | head -n 20 | csvtk pretty -t
        mask   kmer                              prefix   number   ref               pos       strand   reversed
        ----   -------------------------------   ------   ------   ---------------   -------   ------   --------
        1      AAAAAAATAAAAACTTAGTTGTCCCATAACA   8        1        GCF_000392875.1   1044207   -        no      
        1      AAAAAAATAAAATTCCTTCAACCTTTCACTA   8        1        GCF_009759685.1   1863608   -        no      
        1      AAAAAAATAAATCTGCGATGGCTGTTGATGG   8        1        GCF_002950215.1   462416    +        no      
        1      AAAAAAATAAATTTTAATTTCTTCATTAACC   8        1        GCF_001096185.1   843683    +        no      
        1      AAAAAAATAACTCAATGAGGTTATGGGCATG   8        1        GCF_000742135.1   4160317   -        no      
        1      AAAAAAATAACTCCATTTCCTCCTGTTTCTG   8        1        GCF_001544255.1   1605589   -        no      
        1      AAAAAAATAACTGCTTTACTCTTTGCTCTTT   8        1        GCF_009759685.1   2134145   +        no      
        1      AAAAAAATAAGAACACAAAAAAGGTATCTAG   8        1        GCF_001544255.1   1050935   +        no      
        1      AAAAAAATAAGAAGGTAGCACCAATAACTTT   8        1        GCF_900638025.1   137037    -        no      
        1      AAAAAAATAATATTGTCCATTCTCCTAGCAA   8        1        GCF_001544255.1   173045    -        no      
        1      AAAAAAATAATTATATATGTGAACCAAGAAC   8        1        GCF_000148585.2   355059    +        no      
        1      AAAAAAATACATGGATTTATATTTCCAACGC   8        1        GCF_001457655.1   146298    +        no      
        1      AAAAAAATACTTGTTTGATTCTGTATTACGT   8        1        GCF_000392875.1   493472    +        no      
        1      AAAAAAATAGAAAATGAGTCAACACCACTAT   8        1        GCF_006742205.1   1365300   +        no      
        1      AAAAAAATAGAATTATATCGTGAACGTTTTG   8        1        GCF_009759685.1   2234982   +        no      
        1      AAAAAAATAGAGGATTAAATGCTAATTCATA   8        1        GCF_001457655.1   671915    +        no      
        1      AAAAAAATAGCGAATTTTCCAACGACAAAAG   8        1        GCF_002950215.1   650464    -        no      
        1      AAAAAAATAGCGATTTTTCCAACGACAAAAG   8        1        GCF_002949675.1   188498    -        no      
        1      AAAAAAATAGTATAAATCCGCCATATAAAAT   8        1        GCF_001457655.1   1222761   -        no      

1. Specify the mask.

        $ lexicmap utils kmers --quiet -d demo.lmi/ --mask 12345 | head -n 20 | csvtk pretty -t
        mask    kmer                              prefix   number   ref               pos       strand   reversed
        -----   -------------------------------   ------   ------   ---------------   -------   ------   --------
        12345   GCTGCACAAAGTACGATTACGATGCAAGCCC   8        1        GCF_002949675.1   716651    +        no      
        12345   GCTGCACAAAGTGGGCATGCTGGAGCACGCC   8        1        GCF_000742135.1   2084867   -        no      
        12345   GCTGCACAACAGGCTGCGGCTGGTGTTGCGG   8        1        GCF_000742135.1   4128289   -        no      
        12345   GCTGCACAACCAGGCAGAAAAAATAATGGGA   8        1        GCF_002950215.1   3009005   -        no      
        12345   GCTGCACAACGATTAGAAAAAATGGGGTACG   8        1        GCF_001544255.1   2041481   -        no      
        12345   GCTGCACAACTATCCCAATGCCGAGGTGGAA   8        1        GCF_000017205.1   5101754   +        no      
        12345   GCTGCACAAGCGCTCGGTTTAGAGCAAACAC   8        1        GCF_009759685.1   1232954   -        no      
        12345   GCTGCACAAGGGGCCACTTTCGTACATCGTC   8        1        GCF_000742135.1   3888020   +        yes     
        12345   GCTGCACAAGTACCTGCTGGCCTACGCCTCG   8        1        GCF_000017205.1   1166094   +        no      
        12345   GCTGCACAAGTCATCGATGTTAAAGACTTAG   8        1        GCF_009759685.1   514250    +        no      
        12345   GCTGCACAAGTTGCAAAACAGCTGATTAAGG   8        1        GCF_000392875.1   908172    +        no      
        12345   GCTGCACAATATCGATTTGAACATTGCTCAG   8        1        GCF_003697165.2   3212441   +        no      
        12345   GCTGCACAATATTTCATAATGACTTACGGCA   8        1        GCF_002950215.1   3443237   +        no      
        12345   GCTGCACAATCAGGAGTAGGCCACTTAATCA   8        1        GCF_009759685.1   3500297   -        no      
        12345   GCTGCACAATCCGCTGGGCTGGGTGCTCAAC   8        1        GCF_000742135.1   1083211   -        no      
        12345   GCTGCACAATCGCCAGCCCCAGCCCTGTGCC   8        1        GCF_000006945.2   3658390   +        no      
        12345   GCTGCACAATTACCACGTGAATTATTTGAAG   8        1        GCF_900638025.1   304434    -        no      
        12345   GCTGCACAATTGCCAGCCCTAATCCCGTGCC   8        1        GCF_002950215.1   2671971   +        no      
        12345   GCTGCACAATTTGTGAAAGAAGCGAAAGCAT   8        1        GCF_001457655.1   1175308   +        no      

    "reversed" means means if the k-mer is reversed for suffix matching.
    E.g., `GCTGCACAAGGGGCCACTTTCGTACATCGTC` is reversed, so you need to reverse it before searching in the genome.


        $ seqkit locate -p $(echo GCTGCACAAGGGGCCACTTTCGTACATCGTC | rev) refs/GCF_000742135.1.fa.gz -M | csvtk pretty -t
        seqID           patternName                       pattern                           strand   start     end    
        -------------   -------------------------------   -------------------------------   ------   -------   -------
        NZ_KN046818.1   CTGCTACATGCTTTCACCGGGGAACACGTCG   CTGCTACATGCTTTCACCGGGGAACACGTCG   +        3888020   3888050


1. For all masks. The result might be very big, therefore, writing to gzip format is recommended.


        $ lexicmap utils kmers -d demo.lmi/ --mask 0 -o kmers.tsv.gz

        $ zcat kmers.tsv.gz | csvtk freq -t -f mask -nr | head -n 10
        mask    frequency
        8206    716
        11723   646
        4974    643
        12043   619
        14979   617
        16230   612
        18      603
        3724    586
        19979   583

1. Lengths of shared prefixes between probes and captured k-mers.

        zcat kmers.tsv.gz \
            | csvtk grep -t -f reversed -p  no \
            | csvtk plot hist -t -f prefix -o prefix.hist.png \
                --xlab "length of common prefixes between captured k-mers and masks"


    <img src="/LexicMap/prefix.hist.png" alt="" width="400"/>

The output (TSV format) is formatted with [csvtk pretty](https://github.com/shenwei356/csvtk).
