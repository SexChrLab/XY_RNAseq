# XY_RNAseq
Aligning RNA-sequencing reads to a sex chromosome informed reference genome. 
The effects of expression estimates when aligning RNAseq reads to a reference genome that is sensitive to the sex chromosome complement of the reference genome compared to using a default reference genome that does not account for the sex of the sample.

## Background
The human X and Y chromosomes share an evolutionary origin and sequence homology, including regions of 100% identity, that can confound estimates of transcript abundance if sex chromosome complement is not taken into account when aligning reads. Here we present XY_RNAseq, a pipeline for aligning RNAseq reads to an apportiate sex chromosome complement reference to infer differential expression between male XY and female XX RNAseq samples. 

## Preprint

Please see our preprint for more information:

*Aligning RNA-Seq reads to a sex chromosome informed reference genome increases ability to detect sex differences in gene expression* 2019. Olney KC; Brotman SM; Valverde-Vesling V; Andrews J; Wilson MA. 

If you use XY_RNAseq or discuss/correct for bias in mapping on the sex chromosomes, please cite this preprint.

## Contents:
1. Sex chromosome complement refernce genomes - Download and create sex chromosome complement reference genomes and gene annotation files

2. Differential expression pipeline - Suggested pipeline for processing raw RNAseq fastq files for differential expression with HISAT and STAR read aligners followed by limma/voom for differential expression analysis

### 1. Sex chromosome complement refernce genomes:
Creating sex chromosome complement reference genomes for genetic male XY and genetic female XX samples. 
Run the following command to download the ensembl GRCh38 human reference genome and create two sex chromosome complement reference genomes. 

### Sex chromosome complement reference genomes 
--- | --- |
Y-masked | for aligning genetic female XX individuals |
YPARs-masked | for aligning genetic male XX individuals |
default | the default ensemble GRCh38 human refernce genome will include sequences for all chromosomes including the sex chromosomes, X and Y |

To download and create the sex chromosome complement reference genomes run the following script:

  'sh downloadGRCh38HumanReferenceGenomeCreateSexChrComplement.sh'

Note: The GRCh38 human reference genome is large and will take several hours to download depending on processing and memory available. 


### 2. Differential expression pipeline:
We have put together a pipeline for inferring differential expression between males XY and females XX using two read aligners, HISAT and STAR, and limma/voom for computing differential expression. These tools are publicly available and we ask that if you use this pipeline to cite the tools used:

### publicly available tools used in this analysis
Tool | usage | citation
--- | --- |  ---
Trimmomatic | |
HISAT | |
STAR | |
Limma/voom | |



## Group Members
Name | email | github ID
--- | --- |  ---
Kimberly Olney | olneykimberly@gmail.com | @olneykimberly
Sarah Brotman |
Valeria Valverde-Vesling |
Jocelyn Andrews |
Melissa A. Wilson | melissa.wilsonsayres@asu.edu | @mwilsonsayres

