# Tutorials {#tutorials}

[TOC]

Welcome to the Biovoltron tutorial. Explore these tutorials to gain a deeper understanding of Biovoltron's capabilities in various bioinformatics tasks.
Each tutorial provides practical guidance, code examples, and best practices tailored to specific bioinformatics tasks. Click on the links to begin your journey with the tutorials of your choice.

## Applications {#tutorials-applications}
High-level bioinformatics applications built using Biovoltron's modules.

1. [Burrow-Wheeler Aligner](applications/burrow_wheeler_aligner.md) <br/>
   A paired-end short-read aligner optimized for hs37d5 using an FM-index and SSE-accelerated Smith-Waterman.
2. [HaplotypeCaller](applications/haplotypecaller.md) <br/>
   A tool for calling SNPs and indels via local re-assembly of haplotypes.
3. [Adapter Trimmer (EARRINGS)](applications/adapter_trimmer.md) <br/>
   A fast, multithreaded, and SIMD-optimized adapter trimmer for single-end and paired-end reads.
4. [Tailor Pipeline](applications/tailor_pipeline.md) <br/>
   A specialized pipeline for miRNA analysis, including alignment, non-template addition detection, and trimming.

## File Input/Output {#tutorials-fileio}
Modules for reading and writing common bioinformatics file formats.

1. [Core I/O Concepts](file_io/core.md) <br/>
   Learn about the fundamental `Record` and `Header` structures for file parsing.
2. [FASTA Format](file_io/fasta.md) <br/>
   Handle nucleotide or protein sequences stored in FASTA format.
3. [FASTQ Format](file_io/fastq.md) <br/>
   Parse and manipulate sequence data and quality scores from FASTQ files.
4. [SAM/BAM Format](file_io/sam.md) <br/>
   Process Sequence Alignment/Map (SAM) and its binary version (BAM) files.
5. [VCF Format](file_io/vcf.md) <br/>
   Manage and analyze genomic variations stored in Variant Call Format (VCF).
6. [BED Format](file_io/bed.md) <br/>
   Handle genomic annotations stored in Browser Extensible Data (BED) files.
7. [BED6+ Format](file_io/bedsix.md) <br/>
   Work with 6-column BED files containing gene type and name annotations.
8. [GFF/GTF Format](file_io/gff.md) <br/>
   Parse and manipulate data from General Feature Format (GFF) files.
9. [WIG Format](file_io/wig.md) <br/>
   Represent genome-wide signal data using the Wiggle (WIG) format.
10. [CIGAR Strings](file_io/cigar.md) <br/>
    Understand and manipulate CIGAR strings representing sequence alignments.
11. [Reference Genome Format](file_io/reference.md) <br/>
    Work with `ReferenceRecord` to store, save, and load reference genomes efficiently.

## Containers {#tutorials-container}
Specialized, efficient data structures for bioinformatics data.

1. [X-bit Vector (`DibitVector`, `QuadbitVector`)](container/xbit_vector.md) <br/>
   Use `DibitVector` and `QuadbitVector` for space-efficient storage of DNA sequences or other data requiring 2 or 4 bits per value.

## Algorithms {#tutorials-algo}
Core algorithms for sequence analysis, assembly, and annotation.

### Alignment {#tutorials-algo-align}
Algorithms for aligning sequences to a reference or to each other.

1. [FM-Index Exact Matching](algo/align/exact_match/fm_index.md) <br/>
   Use a fast, compressed FM-Index for rapid and precise exact sequence matching.
2. [Pair HMM Alignment](algo/align/inexact_match/pairhmm.md) <br/>
   Perform probabilistic alignment of reads to haplotypes using a Pair Hidden Markov Model. An [AVX-accelerated version](algo/align/inexact_match/pairhmm_avx.md) is also available.
3. [Smith-Waterman Local Alignment](algo/align/inexact_match/smithwaterman.md) <br/>
   Learn the fundamentals of Smith-Waterman for optimal local sequence alignment. Accelerated versions for [AVX2](algo/align/inexact_match/smithwaterman_avx.md) and [SSE](algo/align/inexact_match/smithwaterman_sse.md) are available.
4. [SIMD-accelerated Alignment](algo/align/inexact_match/simd_alignment.md) <br/>
   Explore a generic SIMD-optimized engine for local and global sequence alignment.
5. [Partial Order Alignment (SPOA)](algo/align/spoa.md) <br/>
   Generate consensus sequences from a set of related sequences using Partial Order Alignment.
6. [Wavefront Partial Order Alignment (WFPOA)](algo/align/wfpoa.md) <br/>
   Utilize the Wavefront algorithm for efficient Partial Order Alignment and consensus generation.
7. [Mapping Quality (MAPQ) Calculation](algo/align/mapq.md) <br/>
   Understand and calculate mapping quality scores for single-end and paired-end alignments.
8. [Tailor for miRNA Analysis](algo/align/tailor/tailor.md) <br/>
   A specialized tool for analyzing miRNA non-template additions and trimming events.

### Assembly {#tutorials-algo-assemble}
Algorithms for constructing larger sequences from smaller fragments.

1. [De Bruijn Graph Assembler](algo/assemble/assembler.md) <br/>
   Assemble haplotypes or discover adapters from reads using a de Bruijn graph-based approach.

### Annotation {#tutorials-algo-annotate}
Algorithms for annotating genomic features.

1. [Genomic Interval Annotator](algo/annotate/annotator.md) <br/>
   Use an Interval Tree to efficiently find and annotate overlapping genomic features.

### Suffix Sorting {#tutorials-algo-sort}
A collection of algorithms for constructing suffix arrays.

1. [Suffix Sorter Algorithms](algo/suffix_sorter/core/suffix_sorter.md) <br/>
   An overview of suffix sorting concepts in Biovoltron.
2. [SAIS-based Sorters](algo/sort/sais_sorter.md) <br/>
   Explore various Suffix-Array-Induced-Sorting (SAIS) based algorithms, including [pSAIS](algo/sort/psais_sorter.md), [k-prefix pSAIS](algo/sort/kpsais_sorter.md), [k-ISS1](algo/sort/kiss1_sorter.md), and [k-ISS2](algo/sort/kiss2_sorter.md).
3. [Stable Sorter](algo/sort/stable_sorter.md) <br/>
   A straightforward, parallelized suffix sorter using `std::stable_sort`.

## Pipes {#tutorials-pipe}
Tools for chaining multiple bioinformatics operations into a single pipeline.

1. [Algorithm Pipe](pipe/algo_pipe.md) <br/>
   Learn how to use the `|` operator to chain together building an index, aligning reads, and calling variants in a single expression.

## Utilities {#tutorials-utility}
Fundamental data structures and helper functions used throughout the library.

1. [DNA Representations (`dna4`, `dna5`)](utility/dna4and5.md) <br/>
   Work with DNA sequences using 4-base (A,C,G,T) or 5-base (A,C,G,T,N) alphabets.
2. [Integer String (`istring`)](utility/istring.md) <br/>
   Use integer-based strings for efficient storage and manipulation of DNA sequences.
3. [Genomic Interval](utility/interval.md) <br/>
   Define and operate on genomic regions with chromosome, start, end, and strand.
4. [Haplotype Representation](utility/haplotype.md) <br/>
   Understand the data structure for representing a single haplotype sequence and its associated variants.
5. [Variant Representation](utility/variant.md) <br/>
   Learn about the data structure for storing genomic variants, including location, alleles, and genotype information.
6. [Genotype Representation](utility/genotype.md) <br/>
   Work with genotype data using a simple `std::pair` structure.
7. [Read Utilities](utility/read.md) <br/>
   Utilities for filtering and clipping sequence reads based on various criteria.
8. [Range Utilities](utility/range.md) <br/>
   Helper functions for working with C++20 ranges.
9. [Archive Utilities (`gzstream`, `Serializer`)](utility/archive/serializer.md) <br/>
   Tools for data serialization and handling gzipped streams, including a tutorial for [gzstream](utility/archive/gzstream.md).
10. [Expression Analysis](utility/expression/mirna/mirna_exp.md) <br/>
    Data structures and functions for miRNA expression quantification and [normalization](utility/expression/normalization.md).
11. [hs37d5 Reference Utilities](utility/ref/hs37d5.md) <br/>
    Constants and helper functions specific to the hs37d5 human reference genome.
12. [Thread Pool](utility/threadpool.md) <br/>
    A utility for managing a pool of threads for parallel execution of tasks.

## Math {#tutorials-math}
Mathematical functions tailored for bioinformatics calculations.

1. [Math Utilities](math/math_utils.md) <br/>
   A collection of mathematical functions for tasks like quality score conversion, factorial calculations, and logarithmic operations.