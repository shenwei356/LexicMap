---
title: genome-details
weight: 21
---

## Usage

```plain
Extract or view genome details in the index

On the first run, this command will extract genome details and save them to 'genomes.details.bin'.
If the file exists, it will be read directly and details will be printed.

Output format:
  Tab-delimited format.

    ref             genome id
    genome_size     genome size (sum of all genome chunks)
    chunks          the number of genome chunks
    chunk           nth genome chunk
    chunk_size      genome (chunk) size
    seqs            the number of sequences in the genome (chunk)
    seqsizes        comma-separated sequence sizes in the genome (chunk)    (optional with -e/--extra)
    seqids          comma-separated sequence ids in the genome (chunk)      (optional with -e/--extra)
                    only available when the genome details file is created with the -i/--save-seqids flag.

  Note that genome chunks are created when a genome is too large, and the chunk size is determined by the
  "-g/--max-genome" parameter in "lexicmap index". If a genome is not chunked, it will be treated as one chunk.

Usage:
  lexicmap utils genome-details [flags] 

Flags:
  -e, --extra             ► Show extra columns, including seqsizes and seqids.
  -h, --help              help for genome-details
  -d, --index string      ► Index directory created by "lexicmap index".
  -o, --out-file string   ► Out file, supports the ".gz" suffix ("-" for stdout). (default "-")
  -i, --save-seqids       ► Extract and save sequence ids. This will increase the file size

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

Suppose that we have a few fungal genomes. 

```plain
$ ls t.fungi/*.fna.gz | seqkit stats -X - 
processed files:  3 / 3 [======================================] ETA: 0s. done
file                            format  type  num_seqs     sum_len  min_len      avg_len     max_len
t.fungi/GCF_000003525.1.fna.gz  FASTA   DNA         27  66,608,683    2,055  2,466,988.3  10,302,040
t.fungi/GCF_000146465.1.fna.gz  FASTA   DNA         11   2,216,898  160,332    201,536.2     236,244
t.fungi/GCF_001481775.2.fna.gz  FASTA   DNA        172  37,908,251    1,023    220,396.8   1,830,400
```

We can see that two of them are larger than 15 Mb.
During indexing, they will be chunked.

```plain
$ lexicmap index -I t.fungi -O t.fungi.lmi --force
...
22:10:36.329 [INFO] --------------------- [ building index ] ---------------------
22:10:36.588 [INFO] 
22:10:36.588 [INFO]   ------------------------[ batch 1/1 ]------------------------
22:10:36.588 [INFO]   building index for batch 1 with 3 files...
processed files:  0 / 3 [--------------------------------------] ETA: 0s
22:10:36.798 [WARN]   splitting a big genome into 3 chunks: t.fungi/GCF_001481775.2.fna.gz
processed files:  0 / 3 [--------------------------------------] ETA: 0s
22:10:36.920 [WARN]   splitting a big genome into 6 chunks: t.fungi/GCF_000003525.1.fna.gz
processed files:  3 / 3 [======================================] ETA: 0s. done
...
```
    
Check which genomes are chunked.

```plain
$ lexicmap utils genomes -d t.fungi.lmi/  | csvtk pretty -t
ref               chunked
---------------   -------
GCF_000146465.1          
GCF_000003525.1   yes    
GCF_000003525.1   yes    
GCF_001481775.2   yes    
GCF_000003525.1   yes    
GCF_000003525.1   yes    
GCF_000003525.1   yes    
GCF_000003525.1   yes    
GCF_001481775.2   yes    
GCF_001481775.2   yes
```
    
This command can show more details

```plain
$ lexicmap utils genome-details -d t.fungi.lmi/ -i -o t.txt
22:11:48.239 [INFO] extracting genome details and saving to t.fungi.lmi/genomes.details.bin
22:11:48.239 [INFO]   reading index info file
22:11:48.239 [INFO]   creating reader pools for 1 genome batches, each with 1 reader(s)...
22:11:48.239 [INFO]   reading genome chunk data files
22:11:48.240 [INFO]   elapsed time: 369.671µs
22:11:48.240 [INFO] 
22:11:48.240 [INFO] reading genome details from t.fungi.lmi/genomes.details.bin
22:11:48.240 [INFO]   elapsed time: 32.407µs

$ dirsize t.fungi.lmi/
t.fungi.lmi/: 134.47 MiB (141,005,063)
108.80 MiB      seeds
 25.51 MiB      genomes
156.28 KiB      masks.bin
  4.28 KiB      genomes.details.bin
     921 B      info.toml
     250 B      genomes.map.bin
      88 B      genomes.chunks.bin

$ csvtk pretty -t t.txt
ref               genome_size   chunks   chunk   chunk_size   seqs
---------------   -----------   ------   -----   ----------   ----
GCF_000146465.1   2216898       1        1       2216898      11  
GCF_000003525.1   66608683      6        1       4140423      15  
GCF_000003525.1   66608683      6        2       8458534      1   
GCF_000003525.1   66608683      6        3       10302040     1   
GCF_000003525.1   66608683      6        4       14887486     3   
GCF_000003525.1   66608683      6        5       14159393     2   
GCF_000003525.1   66608683      6        6       14660807     5   
GCF_001481775.2   37908251      3        1       8428770      50  
GCF_001481775.2   37908251      3        2       14930313     33  
GCF_001481775.2   37908251      3        3       14549168     89
```

show more

```plain
$ lexicmap utils genome-details -d t.fungi.lmi/ -e -i | csvtk pretty -t -W 30 -x ,
22:22:14.035 [INFO] reading genome details from t.fungi.lmi/genomes.details.bin
22:22:14.035 [INFO]   elapsed time: 102.876µs
ref               genome_size   chunks   chunk   chunk_size   seqs   seqsizes                         seqids                        
---------------   -----------   ------   -----   ----------   ----   ------------------------------   ------------------------------
GCF_000146465.1   2216898       1        1       2216898      11     160332,175776,176815,193740,     NC_014415.1,NC_014416.1,      
                                                                     196642,198217,205935,204910,     NC_014417.1,NC_014418.1,      
                                                                     233397,234890,236244             NC_014419.1,NC_014420.1,      
                                                                                                      NC_014421.1,NC_014422.1,      
                                                                                                      NC_014423.1,NC_014424.1,      
                                                                                                      NC_014425.1                   
GCF_000003525.1   66608683      6        1       4140423      15     1266001,1197428,1165216,         NW_025544786.1,NW_025544787.1,
                                                                     291280,80333,27336,21986,        NW_025544788.1,NW_025544789.1,
                                                                     14585,14672,11961,4442,2714,     NW_025544790.1,NW_025544791.1,
                                                                     2055,7105,33309                  NW_025544792.1,NW_025544793.1,
                                                                                                      NW_025544794.1,NW_025544795.1,
                                                                                                      NW_025544796.1,NW_025544797.1,
                                                                                                      NW_025544798.1,NW_025544800.1,
                                                                                                      NW_025544799.1                
GCF_000003525.1   66608683      6        2       8458534      1      8458534                          NW_025544775.1                
GCF_000003525.1   66608683      6        3       10302040     1      10302040                         NW_025544774.1                
GCF_000003525.1   66608683      6        4       14887486     3      5552162,5526002,3809322          NW_025544778.1,NW_025544779.1,
                                                                                                      NW_025544780.1                
GCF_000003525.1   66608683      6        5       14159393     2      7673272,6486121                  NW_025544776.1,NW_025544777.1 
GCF_000003525.1   66608683      6        6       14660807     5      3470708,3420641,3384717,         NW_025544781.1,NW_025544782.1,
                                                                     2312387,2072354                  NW_025544783.1,NW_025544784.1,
                                                                                                      NW_025544785.1                
GCF_001481775.2   37908251      3        1       8428770      50     206451,184594,182966,181086,     NW_020171391.1,NW_020171392.1,
                                                                     179346,177996,1098913,170513,    NW_020171393.1,NW_020171394.1,
                                                                     165395,160443,144741,135584,     NW_020171395.1,NW_020171396.1,
                                                                     133725,128829,125597,114733,     NW_020171343.1,NW_020171397.1,
                                                                     111993,1077425,107249,99142,     NW_020171398.1,NW_020171399.1,
                                                                     89022,87285,86882,85792,80156,   NW_020171400.1,NW_020171401.1,
                                                                     77599,76596,75850,1061107,       NW_020171402.1,NW_020171403.1,
                                                                     75732,62243,60583,59792,59162,   NW_020171404.1,NW_020171405.1,
                                                                     56578,47222,43337,41586,35694,   NW_020171406.1,NW_020171344.1,
                                                                     1046289,32390,29939,26384,       NW_020171407.1,NW_020171408.1,
                                                                     26131,21558,20447,19681,19381,   NW_020171409.1,NW_020171410.1,
                                                                     19283,18348                      NW_020171411.1,NW_020171412.1,
                                                                                                      NW_020171413.1,NW_020171414.1,
                                                                                                      NW_020171415.1,NW_020171416.1,
                                                                                                      NW_020171345.1,NW_020171417.1,
                                                                                                      NW_020171418.1,NW_020171419.1,
                                                                                                      NW_020171420.1,NW_020171421.1,
                                                                                                      NW_020171422.1,NW_020171423.1,
                                                                                                      NW_020171424.1,NW_020171425.1,
                                                                                                      NW_020171426.1,NW_020171346.1,
                                                                                                      NW_020171427.1,NW_020171428.1,
                                                                                                      NW_020171429.1,NW_020171430.1,
                                                                                                      NW_020171431.1,NW_020171432.1,
                                                                                                      NW_020171433.1,NW_020171434.1,
                                                                                                      NW_020171435.1,NW_020171436.1 
GCF_001481775.2   37908251      3        2       14930313     33     571579,543409,541099,534031,     NW_020171361.1,NW_020171362.1,
                                                                     496733,486107,1534703,440006,    NW_020171363.1,NW_020171364.1,
                                                                     416161,387184,375046,373527,     NW_020171365.1,NW_020171366.1,
                                                                     366839,353419,347232,346848,     NW_020171340.1,NW_020171367.1,
                                                                     342359,1499607,319338,306888,    NW_020171368.1,NW_020171369.1,
                                                                     302879,278528,272581,253825,     NW_020171370.1,NW_020171371.1,
                                                                     250542,241653,240401,229821,     NW_020171372.1,NW_020171373.1,
                                                                     1402731,226051,222589,215524,    NW_020171374.1,NW_020171375.1,
                                                                     211073                           NW_020171376.1,NW_020171341.1,
                                                                                                      NW_020171377.1,NW_020171378.1,
                                                                                                      NW_020171379.1,NW_020171380.1,
                                                                                                      NW_020171381.1,NW_020171382.1,
                                                                                                      NW_020171383.1,NW_020171384.1,
                                                                                                      NW_020171385.1,NW_020171386.1,
                                                                                                      NW_020171342.1,NW_020171387.1,
                                                                                                      NW_020171388.1,NW_020171389.1,
                                                                                                      NW_020171390.1                
GCF_001481775.2   37908251      3        3       14549168     89     1830400,934585,17109,17138,      NW_020171338.1,NW_020171347.1,
                                                                     16741,16533,14459,13906,13683,   NW_020171438.1,NW_020171437.1,
                                                                     12824,12397,12016,880928,        NW_020171439.1,NW_020171440.1,
                                                                     11499,10620,10563,9521,9285,     NW_020171441.1,NW_020171442.1,
                                                                     8649,7997,7831,7373,7395,        NW_020171443.1,NW_020171444.1,
                                                                     851809,7209,6934,6731,6518,      NW_020171445.1,NW_020171446.1,
                                                                     6515,6508,5567,4885,4820,        NW_020171348.1,NW_020171447.1,
                                                                     830279,4776,4334,4266,4099,      NW_020171448.1,NW_020171449.1,
                                                                     3951,3860,3777,3439,3390,3379,   NW_020171450.1,NW_020171451.1,
                                                                     829916,3334,3221,3074,3076,      NW_020171452.1,NW_020171453.1,
                                                                     2914,2792,2778,2643,2512,2429,   NW_020171454.1,NW_020171456.1,
                                                                     826040,2407,2263,2252,2104,      NW_020171455.1,NW_020171349.1,
                                                                     2061,1954,1894,1889,1880,1872,   NW_020171457.1,NW_020171458.1,
                                                                     780056,1791,1711,1711,1662,      NW_020171459.1,NW_020171460.1,
                                                                     1653,1646,1640,1509,1480,1328,   NW_020171461.1,NW_020171462.1,
                                                                     702779,1310,1230,1068,1023,      NW_020171463.1,NW_020171464.1,
                                                                     697391,687281,1710302,676739,    NW_020171465.1,NW_020171350.1,
                                                                     673579,632956,601520             NW_020171466.1,NW_020171467.1,
                                                                                                      NW_020171468.1,NW_020171469.1,
                                                                                                      NW_020171470.1,NW_020171471.1,
                                                                                                      NW_020171472.1,NW_020171473.1,
                                                                                                      NW_020171474.1,NW_020171475.1,
                                                                                                      NW_020171351.1,NW_020171476.1,
                                                                                                      NW_020171477.1,NW_020171479.1,
                                                                                                      NW_020171478.1,NW_020171480.1,
                                                                                                      NW_020171481.1,NW_020171482.1,
                                                                                                      NW_020171483.1,NW_020171484.1,
                                                                                                      NW_020171485.1,NW_020171352.1,
                                                                                                      NW_020171486.1,NW_020171487.1,
                                                                                                      NW_020171488.1,NW_020171489.1,
                                                                                                      NW_020171490.1,NW_020171491.1,
                                                                                                      NW_020171492.1,NW_020171493.1,
                                                                                                      NW_020171494.1,NW_020171495.1,
                                                                                                      NW_020171353.1,NW_020171496.1,
                                                                                                      NW_020171497.1,NW_020171498.1,
                                                                                                      NW_020171499.1,NW_020171500.1,
                                                                                                      NW_020171501.1,NW_020171502.1,
                                                                                                      NW_020171503.1,NW_020171504.1,
                                                                                                      NW_020171505.1,NW_020171354.1,
                                                                                                      NW_020171506.1,NW_020171507.1,
                                                                                                      NW_020171508.1,NW_020171509.1,
                                                                                                      NW_020171355.1,NW_020171356.1,
                                                                                                      NW_020171339.1,NW_020171357.1,
                                                                                                      NW_020171358.1,NW_020171359.1,
                                                                                                      NW_020171360.1                
```
