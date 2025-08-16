---
title: edit-genome-ids
weight: 70
---

## Usage

```plain
Edit genome IDs in the index via a regular expression

Use cases:
  In the 'lexicmap index' command, users might forget to use the flag
  -N/--ref-name-regexp to extract the genome ID from the sequence file.
  A genome file from NCBI looks like:

    GCF_009818595.1_ASM981859v1_genomic.fna.gz

  In this case, the genome ID would be GCF_009818595.1_ASM981859v1_genomic,
  which is too long. So we can use this command to extract the assembly
  accession via:

    lexicmap utils edit-genome-ids -d t.lmi/ -p '^(\w{3}_\d{9}\.\d+).*' -r '$1'

Tips:
  - A backup file (genomes.map.bin.bak) will be created on the first run.

Usage:
  lexicmap utils edit-genome-ids [flags] 

Flags:
  -h, --help                 help for edit-genome-ids
  -d, --index string         ► Index directory created by "lexicmap index".
  -p, --pattern string       ► Search regular expression".
  -r, --replacement string   ► Replacement. Supporting capture variables.  e.g. $1 represents the text
                             of the first submatch. ATTENTION: for *nix OS, use SINGLE quote NOT double
                             quotes or use the \ escape character.

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

Suppose that we have one genome.

    $ ls GCF_009818595.1_ASM981859v1_genomic.fna.gz 
    GCF_009818595.1_ASM981859v1_genomic.fna.gz

During the indexing, one might use the default parameters and **forgot to set `-N/--ref-name-regexp` to extract the genome ID from the sequence file**.

    $ lexicmap index GCF_009818595.1_ASM981859v1_genomic.fna.gz -O t.lmi
    
So the genome ids in the index are something like this.

    $ lexicmap utils genomes -d t.lmi/ | csvtk pretty -Ht
    ref                                   chunked
    GCF_009818595.1_ASM981859v1_genomic 
    
Well, those IDs will be in the search result, which are too long and not convenient for downstream analyis.

    $ lexicmap search -d t.lmi/ b.gene_E_coli_16S.fasta --quiet \
        | csvtk cut -t -f 1-11 \
        | csvtk pretty -t
      
    query                         qlen   hits   sgenome                               sseqid          qcovGnm   cls   hsp   qcovHSP   alenHSP   pident
    ---------------------------   ----   ----   -----------------------------------   -------------   -------   ---   ---   -------   -------   ------
    NC_000913.3:4166659-4168200   1542   1      GCF_009818595.1_ASM981859v1_genomic   NZ_OX216966.1   99.676    1     1     99.676    1562      77.337
    NC_000913.3:4166659-4168200   1542   1      GCF_009818595.1_ASM981859v1_genomic   NZ_OX216966.1   99.676    2     2     99.676    1562      77.209
    NC_000913.3:4166659-4168200   1542   1      GCF_009818595.1_ASM981859v1_genomic   NZ_OX216966.1   99.676    3     3     99.676    1562      77.145
    NC_000913.3:4166659-4168200   1542   1      GCF_009818595.1_ASM981859v1_genomic   NZ_OX216966.1   99.676    4     4     99.676    1562      77.145
    NC_000913.3:4166659-4168200   1542   1      GCF_009818595.1_ASM981859v1_genomic   NZ_OX216966.1   99.676    5     5     99.676    1562      77.145
    NC_000913.3:4166659-4168200   1542   1      GCF_009818595.1_ASM981859v1_genomic   NZ_OX216966.1   99.676    6     6     99.676    1562      77.145
    NC_000913.3:4166659-4168200   1542   1      GCF_009818595.1_ASM981859v1_genomic   NZ_OX216966.1   99.676    7     7     99.676    1562      77.145
    NC_000913.3:4166659-4168200   1542   1      GCF_009818595.1_ASM981859v1_genomic   NZ_OX216966.1   99.676    8     8     99.676    1562      77.081

Luckily, we can use this command to fix this.

    $ lexicmap utils edit-genome-ids -d t.lmi/ -p '^(\w{3}_\d{9}\.\d+).*' -r '$1'
    15:51:29.584 [INFO] 1 of 1 genome IDs are changed
    
Check again:

    $ lexicmap utils genomes -d t.lmi/ | csvtk pretty -Ht
    ref               chunked
    GCF_009818595.1        
    
    $ lexicmap search -d t.lmi/ b.gene_E_coli_16S.fasta --quiet \
        | csvtk cut -t -f 1-11 \
        | csvtk pretty -t
    query                         qlen   hits   sgenome           sseqid          qcovGnm   cls   hsp   qcovHSP   alenHSP   pident
    ---------------------------   ----   ----   ---------------   -------------   -------   ---   ---   -------   -------   ------
    NC_000913.3:4166659-4168200   1542   1      GCF_009818595.1   NZ_OX216966.1   99.676    1     1     99.676    1562      77.337
    NC_000913.3:4166659-4168200   1542   1      GCF_009818595.1   NZ_OX216966.1   99.676    2     2     99.676    1562      77.209
    NC_000913.3:4166659-4168200   1542   1      GCF_009818595.1   NZ_OX216966.1   99.676    3     3     99.676    1562      77.145
    NC_000913.3:4166659-4168200   1542   1      GCF_009818595.1   NZ_OX216966.1   99.676    4     4     99.676    1562      77.145
    NC_000913.3:4166659-4168200   1542   1      GCF_009818595.1   NZ_OX216966.1   99.676    5     5     99.676    1562      77.145
    NC_000913.3:4166659-4168200   1542   1      GCF_009818595.1   NZ_OX216966.1   99.676    6     6     99.676    1562      77.145
    NC_000913.3:4166659-4168200   1542   1      GCF_009818595.1   NZ_OX216966.1   99.676    7     7     99.676    1562      77.145
    NC_000913.3:4166659-4168200   1542   1      GCF_009818595.1   NZ_OX216966.1   99.676    8     8     99.676    1562      77.081
