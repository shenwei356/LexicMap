---
title: 2sam
weight: 0.5
---

## Usage

```plain
$ lexicmap utils 2sam -h
Convert the default search output to SAM format

Input:
   - Output file of 'lexicmap search' with the flag -a/--all.
   - Do not support STDIN.

Output:
   - Clipped regions in SEQ are represented as N's,
     as the input contains only the aligned portion of the sequences.
   - Different from SAM files produced by Minimap2,
     'X' (mismatch) in CIGAR is not converted to 'M' (match).
   - NM (edit distance) and AS (alignment score) fields are produced.

Usage:
  lexicmap utils 2sam [flags] 

Flags:
  -b, --buffer-size string      ► Size of buffer, supported unit: K, M, G. You need increase the value
                                when "bufio.Scanner: token too long" error reported (default "20M")
  -c, --concat-sgenome-sseqid   ► Concatenate sgenome and sseqid to make sure the reference sequence
                                names are distinct.
  -h, --help                    help for 2sam
  -o, --out-file string         ► Out file, supports and recommends a ".gz" suffix ("-" for stdout).
                                (default "-")
  -s, --separater string        ► Separater between sgenome and sseqid (default "~")

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

Index with demo data.

```
lexicmap index -I refs/ -O demo.lmi/
```

Search with a 16S rRNA gene.

```
# use -a to output more fields.
# use -n to limit matched genomes. But in practice, the value can't be too small.

lexicmap search -d demo.lmi/ q.gene.fasta -o q.gene.fasta.lexicmap_all.tsv -a -n 2
```

Preview

```text
$ cat q.gene.fasta.lexicmap_all.tsv \
    | csvtk cut -t -f query,sgenome,sseqid,qcovHSP,alenHSP,pident,bitscore,cigar \
    | csvtk pretty -t -S 3line
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 query                         sgenome           sseqid          qcovHSP   alenHSP   pident   bitscore   cigar                                                    
──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 NC_000913.3:4166659-4168200   GCF_003697165.2   NZ_CP033092.2   100.000   1542      99.805   2767       79M1X8M1X120M1X1332M                                     
 NC_000913.3:4166659-4168200   GCF_003697165.2   NZ_CP033092.2   100.000   1542      99.805   2767       79M1X8M1X120M1X1332M                                     
 NC_000913.3:4166659-4168200   GCF_003697165.2   NZ_CP033092.2   100.000   1542      99.805   2767       79M1X8M1X120M1X1332M                                     
 NC_000913.3:4166659-4168200   GCF_003697165.2   NZ_CP033092.2   100.000   1542      99.805   2767       79M1X8M1X120M1X1332M                                     
 NC_000913.3:4166659-4168200   GCF_003697165.2   NZ_CP033092.2   100.000   1542      99.805   2767       79M1X8M1X120M1X1332M                                     
 NC_000913.3:4166659-4168200   GCF_003697165.2   NZ_CP033092.2   100.000   1542      99.805   2767       75M1X3M1X8M1X1453M                                       
 NC_000913.3:4166659-4168200   GCF_003697165.2   NZ_CP033092.2   100.000   1542      99.805   2767       79M1X8M1X120M1X1332M                                     
 NC_000913.3:4166659-4168200   GCF_002949675.1   NZ_CP026774.1   100.000   1542      99.027   2713       75M1X3M1X8M1X658M1X253M1X3M2X2M2X6M2X2M2X14M1X327M1X176M 
 NC_000913.3:4166659-4168200   GCF_002949675.1   NZ_CP026774.1   100.000   1542      99.027   2713       75M1X3M1X8M1X658M1X253M1X3M2X2M2X6M2X2M2X14M1X327M1X176M 
 NC_000913.3:4166659-4168200   GCF_002949675.1   NZ_CP026774.1   100.000   1542      99.027   2713       75M1X3M1X8M1X658M1X253M1X3M2X2M2X6M2X2M2X14M1X327M1X176M 
 NC_000913.3:4166659-4168200   GCF_002949675.1   NZ_CP026774.1   100.000   1542      99.027   2713       75M1X3M1X8M1X658M1X253M1X3M2X2M2X6M2X2M2X14M1X327M1X176M 
 NC_000913.3:4166659-4168200   GCF_002949675.1   NZ_CP026774.1   100.000   1542      99.027   2713       75M1X3M1X8M1X658M1X253M1X3M2X2M2X6M2X2M2X14M1X327M1X176M 
 NC_000913.3:4166659-4168200   GCF_002949675.1   NZ_CP026774.1   100.000   1542      99.027   2713       75M1X3M1X8M1X658M1X253M1X3M2X2M2X6M2X2M2X14M1X327M1X176M 
 NC_000913.3:4166659-4168200   GCF_002949675.1   NZ_CP026774.1   100.000   1542      99.027   2713       75M1X3M1X8M1X658M1X253M1X3M2X2M2X6M2X2M2X14M1X327M1X176M 
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

2sam

```text
$ lexicmap utils 2sam q.gene.fasta.lexicmap_all.tsv 
10:25:03.274 [INFO] round 1/2: extracting subject sequence IDs and lengths ...
10:25:03.274 [INFO]   elapsed time: 188.502µs
10:25:03.274 [INFO] round 2/2: converting to SAM format ...
10:25:03.274 [INFO]   elapsed time: 75.289µs
@HD     VN:1.6  SO:unsorted     GO:query
@SQ     SN:NZ_CP033092.2        LN:4903501
@SQ     SN:NZ_CP026774.1        LN:4395762
@PG     ID:lexicmap     PN:lexicmap     VN:0.9.0
NC_000913.3:4166659-4168200     0       NZ_CP033092.2   458559  0       79M1X8M1X120M1X1332M    *       0       1542    AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTTA     *       NM:i:3  AS:i:3067
NC_000913.3:4166659-4168200     256     NZ_CP033092.2   1285123 0       79M1X8M1X120M1X1332M    *       0       1542    *       *       NM:i:3  AS:i:3067
NC_000913.3:4166659-4168200     272     NZ_CP033092.2   3780640 0       79M1X8M1X120M1X1332M    *       0       1542    *       *       NM:i:3  AS:i:3067
NC_000913.3:4166659-4168200     272     NZ_CP033092.2   4551515 0       79M1X8M1X120M1X1332M    *       0       1542    *       *       NM:i:3  AS:i:3067
NC_000913.3:4166659-4168200     272     NZ_CP033092.2   4591684 0       79M1X8M1X120M1X1332M    *       0       1542    *       *       NM:i:3  AS:i:3067
NC_000913.3:4166659-4168200     272     NZ_CP033092.2   4726193 0       75M1X3M1X8M1X1453M      *       0       1542    *       *       NM:i:3  AS:i:3067
NC_000913.3:4166659-4168200     272     NZ_CP033092.2   4844587 0       79M1X8M1X120M1X1332M    *       0       1542    *       *       NM:i:3  AS:i:3067
NC_000913.3:4166659-4168200     272     NZ_CP026774.1   1662010 0       75M1X3M1X8M1X658M1X253M1X3M2X2M2X6M2X2M2X14M1X327M1X176M        *       0       1542    *       * NM:i:15  AS:i:3007
NC_000913.3:4166659-4168200     256     NZ_CP026774.1   2536624 0       75M1X3M1X8M1X658M1X253M1X3M2X2M2X6M2X2M2X14M1X327M1X176M        *       0       1542    *       * NM:i:15  AS:i:3007
NC_000913.3:4166659-4168200     256     NZ_CP026774.1   2636477 0       75M1X3M1X8M1X658M1X253M1X3M2X2M2X6M2X2M2X14M1X327M1X176M        *       0       1542    *       * NM:i:15  AS:i:3007
NC_000913.3:4166659-4168200     256     NZ_CP026774.1   2768883 0       75M1X3M1X8M1X658M1X253M1X3M2X2M2X6M2X2M2X14M1X327M1X176M        *       0       1542    *       * NM:i:15  AS:i:3007
NC_000913.3:4166659-4168200     256     NZ_CP026774.1   2810845 0       75M1X3M1X8M1X658M1X253M1X3M2X2M2X6M2X2M2X14M1X327M1X176M        *       0       1542    *       * NM:i:15  AS:i:3007
NC_000913.3:4166659-4168200     256     NZ_CP026774.1   3061592 0       75M1X3M1X8M1X658M1X253M1X3M2X2M2X6M2X2M2X14M1X327M1X176M        *       0       1542    *       * NM:i:15  AS:i:3007
NC_000913.3:4166659-4168200     256     NZ_CP026774.1   3646778 0       75M1X3M1X8M1X658M1X253M1X3M2X2M2X6M2X2M2X14M1X327M1X176M        *       0       1542    *       * NM:i:15  AS:i:3007
```

Validate with samtools

```
$ lexicmap utils 2sam q.gene.fasta.lexicmap_all.tsv | samtool view -bS > q.gene.fasta.lexicmap_all.tsv.bam
```
