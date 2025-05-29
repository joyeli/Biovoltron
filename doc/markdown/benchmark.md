# Benchmark {#benchmark}


[TOC]

## Burrow Wheeler Aligner

The elapsed time of mapping sequencing data with 50× coverage 
(2×180GB) is less than 2 hours under 80 cores and ~18 GB memory 
usage and has comparable performance to [bwa-mem2]( https://github.com/bwa-mem2/bwa-mem2).

### Performance
PrecisionFDA [Truth Challenge]( https://precision.fda.gov/challenges/truth) 
benchmark versus [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)
(Dataset [hs37d5](https://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz)):

![](https://i.imgur.com/1T35IWl.png)

### CPU profile
![](https://i.imgur.com/ypMLT1J.png)
