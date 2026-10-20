---
title: lexicmap utils 2sam
linkTitle: 2sam
weight: 0.5
---

## Usage

```plain
$ lexicmap utils 2sam -h
Convert the default search output to SAM format

Input:
   - Output file of 'lexicmap search' with the flag -a/--all.
     Results with --show-avg-qual are also supported.
   - Do not support STDIN.

Output:
   - Clipped regions in SEQ are represented as N's,
     as the input contains only the aligned portion of the sequences.
   - Different from SAM files produced by Minimap2,
     'X' (mismatch) in CIGAR is not converted to 'M' (match).
   - NM (edit distance) and AS (alignment score) fields are produced.
   - Custom integer length tags are emitted in this order:
     ZA:i: alignment length (alenHSP, without the optional quality suffix).
     ZQ:i: query alignment span (qend - qstart + 1).
     ZR:i: reference alignment span (send - sstart + 1).
     Both spans exclude clipping; ZQ includes insertions, ZR includes deletions.
   - Reverse-strand alignments use reverse-complemented SEQ and reversed
     CIGAR operation order, including soft clipping.
   - Records are unpaired, with RNEXT=*, PNEXT=0, and TLEN=0.

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

lexicmap search -d demo.lmi/ q.gene.fasta -o q.gene.fasta.lexicmap_top-2-genomes_all.tsv -a -n 2
```

Preview

```text
$ cat q.gene.fasta.lexicmap_top-2-genomes_all.tsv \
    | csvtk cut -t -f query,sgenome,sseqid,qcovHSP,alenHSP,pident,bitscore,cigar \
    | csvtk pretty -t -S 3line

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 query                         sgenome           sseqid          qcovHSP   alenHSP   pident   bitscore   cigar                               
─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 NC_000913.3:4166659-4168200   GCF_003697165.2   NZ_CP033092.2   100.000   1542      99.805   2767       79M1X8M1X120M1X1332M                
 NC_000913.3:4166659-4168200   GCF_003697165.2   NZ_CP033092.2   100.000   1542      99.805   2767       79M1X8M1X120M1X1332M                
 NC_000913.3:4166659-4168200   GCF_003697165.2   NZ_CP033092.2   100.000   1542      99.805   2767       79M1X8M1X120M1X1332M                
 NC_000913.3:4166659-4168200   GCF_003697165.2   NZ_CP033092.2   100.000   1542      99.805   2767       79M1X8M1X120M1X1332M                
 NC_000913.3:4166659-4168200   GCF_003697165.2   NZ_CP033092.2   100.000   1542      99.805   2767       79M1X8M1X120M1X1332M                
 NC_000913.3:4166659-4168200   GCF_003697165.2   NZ_CP033092.2   100.000   1542      99.805   2767       75M1X3M1X8M1X1453M                  
 NC_000913.3:4166659-4168200   GCF_003697165.2   NZ_CP033092.2   100.000   1542      99.805   2767       79M1X8M1X120M1X1332M                
 NC_000913.3:4166659-4168200   GCF_002950215.1   NZ_CP026788.1   100.000   1542      99.676   2758       79M1X8M2X2M1X801M1X647M             
 NC_000913.3:4166659-4168200   GCF_002950215.1   NZ_CP026788.1   100.000   1542      99.611   2755       79M1X8M2X2M1X380M1X373M1X694M       
 NC_000913.3:4166659-4168200   GCF_002950215.1   NZ_CP026788.1   100.000   1542      99.611   2755       79M1X8M2X2M1X717M1X36M1X694M        
 NC_000913.3:4166659-4168200   GCF_002950215.1   NZ_CP026788.1   100.000   1542      99.546   2749       79M1X8M2X2M1X380M1X336M1X36M1X694M  
 NC_000913.3:4166659-4168200   GCF_002950215.1   NZ_CP026788.1   100.000   1542      99.546   2749       79M1X8M2X2M1X322M1X394M1X36M1X694M  
 NC_000913.3:4166659-4168200   GCF_002950215.1   NZ_CP026788.1   100.000   1542      99.546   2746       79M1I1X7M1X3M1X156M1X2M1X19M1X1269M 
 NC_000913.3:4166659-4168200   GCF_002950215.1   NZ_CP026788.1   100.000   1542      99.481   2746       78M2X8M2X2M1X717M1X36M1X46M1X647M   
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

2sam

```text
$ lexicmap utils 2sam q.gene.fasta.lexicmap_top-2-genomes_all.tsv
10:25:03.274 [INFO] round 1/2: extracting subject sequence IDs and lengths ...
10:25:03.274 [INFO]   elapsed time: 188.502µs
10:25:03.274 [INFO] round 2/2: converting to SAM format ...
10:25:03.274 [INFO]   elapsed time: 75.289µs
@HD     VN:1.6  SO:unsorted     GO:query
@SQ     SN:NZ_CP033092.2        LN:4903501
@SQ     SN:NZ_CP026788.1        LN:4659463
@PG     ID:lexicmap     PN:lexicmap     VN:0.10.0
NC_000913.3:4166659-4168200	0	NZ_CP033092.2	458559	0	79M1X8M1X120M1X1332M	*	0	0	AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTTA	*	NM:i:3	AS:i:3067	ZA:i:1542	ZQ:i:1542	ZR:i:1542
NC_000913.3:4166659-4168200	256	NZ_CP033092.2	1285123	0	79M1X8M1X120M1X1332M	*	0	0	*	*	NM:i:3	AS:i:3067	ZA:i:1542	ZQ:i:1542	ZR:i:1542
NC_000913.3:4166659-4168200	272	NZ_CP033092.2	3780640	0	1332M1X120M1X8M1X79M	*	0	0	*	*	NM:i:3	AS:i:3067	ZA:i:1542	ZQ:i:1542	ZR:i:1542
NC_000913.3:4166659-4168200	272	NZ_CP033092.2	4551515	0	1332M1X120M1X8M1X79M	*	0	0	*	*	NM:i:3	AS:i:3067	ZA:i:1542	ZQ:i:1542	ZR:i:1542
NC_000913.3:4166659-4168200	272	NZ_CP033092.2	4591684	0	1332M1X120M1X8M1X79M	*	0	0	*	*	NM:i:3	AS:i:3067	ZA:i:1542	ZQ:i:1542	ZR:i:1542
NC_000913.3:4166659-4168200	272	NZ_CP033092.2	4726193	0	1453M1X8M1X3M1X75M	*	0	0	*	*	NM:i:3	AS:i:3067	ZA:i:1542	ZQ:i:1542	ZR:i:1542
NC_000913.3:4166659-4168200	272	NZ_CP033092.2	4844587	0	1332M1X120M1X8M1X79M	*	0	0	*	*	NM:i:3	AS:i:3067	ZA:i:1542	ZQ:i:1542	ZR:i:1542
NC_000913.3:4166659-4168200	256	NZ_CP026788.1	3216505	0	79M1X8M2X2M1X801M1X647M	*	0	0	*	*	NM:i:5	AS:i:3057	ZA:i:1542	ZQ:i:1542	ZR:i:1542
NC_000913.3:4166659-4168200	256	NZ_CP026788.1	3119331	0	79M1X8M2X2M1X380M1X373M1X694M	*	0	0	*	*	NM:i:6	AS:i:3053	ZA:i:1542	ZQ:i:1542	ZR:i:1542
NC_000913.3:4166659-4168200	256	NZ_CP026788.1	3396068	0	79M1X8M2X2M1X717M1X36M1X694M	*	0	0	*	*	NM:i:6	AS:i:3053	ZA:i:1542	ZQ:i:1542	ZR:i:1542
NC_000913.3:4166659-4168200	256	NZ_CP026788.1	3355632	0	79M1X8M2X2M1X380M1X336M1X36M1X694M	*	0	0	*	*	NM:i:7	AS:i:3047	ZA:i:1542	ZQ:i:1542	ZR:i:1542
NC_000913.3:4166659-4168200	256	NZ_CP026788.1	4223146	0	79M1X8M2X2M1X322M1X394M1X36M1X694M	*	0	0	*	*	NM:i:7	AS:i:3047	ZA:i:1542	ZQ:i:1542	ZR:i:1542
NC_000913.3:4166659-4168200	272	NZ_CP026788.1	2125377	0	1269M1X19M1X2M1X156M1X3M1X7M1X1I79M	*	0	0	*	*	NM:i:7	AS:i:3043	ZA:i:1542	ZQ:i:1542	ZR:i:1541
NC_000913.3:4166659-4168200	256	NZ_CP026788.1	3540450	0	78M2X8M2X2M1X717M1X36M1X46M1X647M	*	0	0	*	*	NM:i:8	AS:i:3043	ZA:i:1542	ZQ:i:1542	ZR:i:1542
```

Validate with samtools

```
$ lexicmap utils 2sam q.gene.fasta.lexicmap_top-2-genomes_all.tsv | samtools view -bS > q.gene.fasta.lexicmap_all.tsv.bam
```
