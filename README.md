# XY_RNAseq
Aligning RNA-sequencing reads to a sex chromosome informed reference genome. 
The effects of expression estimates when aligning RNAseq reads to a reference genome that is sensitive to the sex chromosome complement of the reference genome compared to using a default reference genome that does not account for the sex of the sample.

## Background
The human X and Y chromosomes share an evolutionary origin and sequence homology, including regions of 100% identity, that can confound estimates of transcript abundance if sex chromosome complement is not taken into account when aligning reads. Here we present XY_RNAseq, a pipeline for aligning RNAseq reads to an apportiate sex chromosome complement reference to infer differential expression between male XY and female XX RNAseq samples. 

## Preprint

Please see our preprint for more information:

*Aligning RNA-Seq reads to a sex chromosome informed reference genome increases ability to detect sex differences in gene expression* 2019. Olney KC;Brotman SM; Valverde-Vesling V; Andrews J; Wilson MA. 

If you use XY_RNAseq or discuss/correct for bias in mapping on the sex chromosomes, please cite this preprint.

## Contents:
1. SexChromosomeComplementRefernceGenomes - Download and create sex chromosome complement reference genomes and gene annotation files
2. DifferentialExpression - Suggested pipeline for processing raw RNAseq fastq files for differential expression with HISAT and STAR read aligners followed by limma/voom for differential expression analysis
