## Differential expression work flow

## Contents:
1. Download RNAseq transcriptome data 
2. convert RNAseq SRA files to be FASTQ format
3. Create FASTQC reports
4. Trim for quality
5. Create FASTQC and MULTIQC reports on pre and post trimmed fastq files
6. Align to the reference genomes using STAR
7. Align to the reference genomes using HISAT2
8. Generate stats on initial BAM files
9. Sort BAM files 
10. Mark duplicates
11. Add or replace read groups 
12. Index BAM files
13. Generate stats on read group bam files
14. get raw transcript counts 
15. Create chromosome CSV file
16. Create phenotype TSV file
17. Create counts TSV file
18. Differential expression using LimmaVoom

## 1. Download Data
The Genotyping-Tissue Expression (GTEx) Project was initiated to give researchers a resource to analyze RNAseq data among human individuals across multiple tissues (GTEx Consortium 2015, 2013). GTEx Project includes 544 recently deceased donors over 53 tissues types for 8,555 total samples, with multiple tissues collected per individual. RNA was performed using a non-strand-specific with a poly-A selection using Illumina TrueSeq and resulted in an average of 50 million 76 base pairs (bp) paired-end reads per sample. 

In sra (sequence read archive, known as short-read archive) format. Include GEO accession number 
`wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRR1850937.sra`

## 2. Convert sra to fastq
Files will need to be converted from sra to fastq for downstream analysis. FASTQ format is a text-based format for storing both a biological sequence (usually nucleotide sequence) and its corresponding quality scores. Both the sequence letter and quality score are each encoded with a single ASCII character for brevity.

`fastq-dump sampleID.sra`

- fastq-dump - program part of the SRA toolkit that converts SRA files to fastq format
- sampleID.sra - path and name of sample.sra that you would like to convert to fastq format 

## 3. Create fastqc reports
Fastqc reads raw sequence data from high throughput sequencers and runs a set of quality checks to produce a report. Best reports are those whose "per base sequence quality" are included in the green area of the graph & kmer content is good or average.

`fastqc sampleID.fastq`

- fastqc - Babraham bioinformatics program that that checks for quality of reads 
- sampleID.fastq - path and name of sampleID in fastq format, may also be in fastq.gz format

## 4. Trim for quality 
Since it has been found that raw untrimmed data leads to errors in read-mapping (Del Fabbro et al. 2013), we tested the effects of trimming versus no-trimming on read abundance, and 
approach that scans reads in the 5’-3’ direction and calculates the average quality of a group of 4 bases, read groups on the 3’-end whose quality scores were lower than the phred score 30 were removed.

The current trimming steps are:
- ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
- SLIDINGWINDOW: Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.
- LEADING: Cut bases off the start of a read, if below a threshold quality
- TRAILING: Cut bases off the end of a read, if below a threshold quality
- CROP: Cut the read to a specified length
- HEADCROP: Cut the specified number of bases from the start of the read
- MINLEN: Drop the read if it is below a specified length
- TOPHRED33: Convert quality scores to Phred-33
- TOPHRED64: Convert quality scores to Phred-64

It works with FASTQ (using phred + 33 or phred + 64 quality scores, depending on the Illumina pipeline used), either uncompressed or gzipp'ed FASTQ. Use of gzip format is determined based on the .gz extension.

The parameters selected were 
- slidingwindow:4:30 leading10 trailing25 minlen40 phred33

`java -jar trimmomatic-0.36.jar PE -phred33 sampleID_input_1.fastq sampleID_input_2.fastq sampleID_output_1_paired.fastq sampleID_output_1_unpaired.fastq sampleID_output_2_paired.fastq sampleID_output_2_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:10 TRAILING:25 SLIDINGWINDOW:4:30 MINLEN:40`

- java - indicates that this is a java program and will require java in order to run
- -jar - jar file to follow
- trimmomatic-0.36.jar - tool that will trim the raw fastq files
- PE - PE is for pair end reads. If single end then SE. 
- -phred33 - Using phred + 33 or phred + 64 quality scores, depending on the Illumina pipeline used, either uncompressed or gzipp'ed FASTQ 
- sampleID_input.fastq - sampeID in fastq format
- sampleID_output.fastq - sampleID output file. Use a descriptive name such as sampleID_minlen50_sliding430_leading30_trailing40.fq
- ILLUMINACLIP:TruSeq3-PE:2:30:10 - Remove Illumina adapters provided in the TruSeq3-PE.fa file (provided). Initially Trimmomatic will look for seed matches (16 bases) allowing maximally 2 mismatches. These seeds will be extended and clipped if in the case of paired end reads a score of 30 is reached (about 50 bases), or in the case of single ended reads a score of 10, (about 17 bases).
- LEADING:10 - Cut bases off the start of a read, if below a threshold quality of 10
- TRAILING:25 - Cut bases off the end of a read, if below a threshold quality of 25
- SLIDINGWINDOW:4:30 - Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 30
- MINLEN:40 - Drop the read if it is below a specified length of 40
- adapters - add pathway to adapters directory

## 5. Create fastqc reports on trimmed files
Fastqc reads sequence data from high throughput sequencers and runs a set of quality checks to produce a report. Best reports are those whose "per base sequence quality" are included in the green area of the graph & kmer content is good or average. This command will create two outputs: an .html file & an .zip file. Will output sampleID_fastqc.html and sampleID_fastqc.zip files

`fastqc sampleID_output_1_paired.fastq sampleID_output_2_paired.fastq`

- fastqc - Babraham bioinformatics program that that checks for quality of reads 
- sampleID.fastq - path and name of sampleID in fastq format, may also be in fastq.gz format

# MULTIQC
`cd multiqc sampleID1_output_1_paired_fastqc/ sampleID1_output_2_paired_fastqc/ sampleID2_output_1_paired_fastqc/ etc...`

## 6. Aligning to the reference genomes using STAR
STAR read aligner is a 2 pass process. The user supplies the genome files generated in the pervious step (generate genome indexes), as well as the RNA-seq reads (sequences) in the form of FASTA or FASTQ files. STAR maps the reads to the genome, and writes several output files, such as alignments (SAM/BAM), mapping summary statistics, splice junctions, unmapped reads, signal (wiggle) tracks etc. Mapping is controlled by a variety of input parameters (options). STAR highly recommends using --sjdbGTFfile which specifies the path to the file with annotated transcripts in the standard GTF format. Where STAR will extract splice junctions from this file and use them to greatly improve accuracy of the mapping. While this is optional, and STAR can be run without annotations, using annotations is highly recommended whenever they are available.However this option should not be included for projects that include hybrids, as this might cause a bias towards the reference. 

Align male samples to the default genome and the YPARs_masked genome 

`STAR --genomeDir gencode.GRCh38.p7_YPARsMasked/ --sjdbGTFfile gencode.GRCh38.p7_YPARsMasked/gencode.v25.chr_patch_hapl_scaff.annotation.gtf --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --readFilesIn sampleID_1_trim_paired.fastq sampleID_2_trim_paired.fastq --outSAMtype BAM Unsorted --outFileNamePrefix sampleID_STAR_M_std_YPARsMasked. --runThreadN 4`

`STAR --genomeDir gencode.GRCh38.p7_default/ --sjdbGTFfile gencode.GRCh38.p7_default/gencode.v25.chr_patch_hapl_scaff.annotation.gtf --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --readFilesIn sampleID_1_trim_paired.fastq sampleID_2_trim_paired.fastq --outSAMtype BAM Unsorted --outFileNamePrefix sampleID_STAR_M_std_default. --runThreadN 4`

Align female samples to the default genome and the Y_masked genome 

`STAR --genomeDir gencode.GRCh38.p7_Ymasked/ --sjdbGTFfile gencode.GRCh38.p7_Ymasked/gencode.v25.chr_patch_hapl_scaff.annotation.gtf --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --readFilesIn sampleID_1_trim_paired.fastq sampleID_2_trim_paired.fastq --outSAMtype BAM Unsorted --outFileNamePrefix sampleID_STAR_F_std_Ymasked. --runThreadN 4`

`STAR --genomeDir gencode.GRCh38.p7_default/ --sjdbGTFfile gencode.GRCh38.p7_default/gencode.v25.chr_patch_hapl_scaff.annotation.gtf --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --readFilesIn sampleID_1_trim_paired.fastq sampleID_2_trim_paired.fastq --outSAMtype BAM Unsorted --outFileNamePrefix sampleID_STAR_F_std_default. --runThreadN 4`

- STAR - STAR read aligner package
- --genomeDir - define where the genome is location, star_genome 
- Project/refrence_genome - path and directory to reference genome
- --genomeLoad - mode of shared memory usage for the genome files
- LoadAndKeep - load genome into shared and keep it in memory after run
- --sjdbGTFfile - path to the GTF file with annotations
- --readFilesIn - if pair end reads, include path to both reads with a " " space inbetween _1 _2, /geuvadis_fastq/sampleID_1.fastq /home/kcolney/map_geuvadis/geuvadis_fastq/sampleID_2.fastq 
- sampleID_1.fastq - name and path to sample in fastq format. If paired end samples include both pairs and separate with a space i.e (sample_1.fastq sample_2.fastq)
- --outSAMtype - indicate which output format, BAM unsorted
- BAM Unsorted - output unsorted Aligned.out.bam file. The paired ends of an alignment are always adjacent, and multiple alignments of a read are adjacent as well. This ”unsorted” file can be directly used with downstream software such as HTseq, without the need of name sorting. The order of the reads will match that of the input FASTQ(A) files only if one thread is used
- --sjdbFileChrStartEnd	- path to the pass_1.SJ.out.tab files made in the first pass
- sampleID1_pass1.SJ.out.tab - 4 columns separated by tabs: Chr \tab Start \tab End \tab Strand=+/-/. Here Start and End are first and last bases of the introns (1-based chromosome coordinates). This file can be used in addition to the --sjdbGTFfile, in which case STAR will extract junctions from both files.
- sampleID2_pass1.SJ.out.tab - List all the samples from the first pass or all the samples in a group (i.e population, cases and controls, hybrids, males and females)
- --outFileNamePrefix - define the sample id prefix, sampleID_pass1. (bam, will be added by the STAR program)
- --runThreadN - for computing purpose allocate the number of threads, 14 


## 7. Align to the reference genomes using HISAT2
using -q to specify reads are fastq, --phred33 to indicate that input qualities are ASCII chars equal to the Phred+33 encoding which is used by the GTEx Illumina processing pipeline. HISAT2 parameters -p 8 launched 8 number of parallel search threads which increased alignment throughput by approximately a multiple of the number of threads, and finally -x followed by the basename of the index for the reference genome being either Def, Y-masked or YPARs-masked.

Align male samples to the default genome and the YPARs_masked genome 

`hisat2 -q --phred33 -p 8 -x GRCh38_YPARsMasked_reference_HISAT2 -s no -1 sampleID_1_trimpaired.fastq -2 sampleID_2_trimpaired.fastq -b sampleID_stdtrim_HISAT_YPARsMasked.sam`

`hisat2 -q --phred33 -p 8 -x GRCh38_default_reference_HISAT2 -s no -1 sampleID_1_trimpaired.fastq -2 sampleID_2_trimpaired.fastq -b sampleID_stdtrim_HISAT_default.sam`

Align female samples to the default genome and the Y_masked genome 

`hisat2 -q --phred33 -p 8 -x GRCh38_Ymasked_reference_HISAT2 -s no -1 sampleID_1_trim_paired.fastq -2 sampleID_2_trim_paired.fastq -b sampleID_stdtrim_HISAT_Ymasked.sam`

`hisat2 -q --phred33 -p 8 -x GRCh38_default_reference_HISAT2 -s no -1 sampleID_1_trim_paired.fastq -2 sampleID_2_trim_paired.fastq -b sampleID_stdtrim_HISAT_default.sam`

- hisat2 - HISAT2 read aligner package
- -q -
- --phred33
- -p
- 8
- -x
- -s
- no
- -1
- sampleID_1.fastq - name and path to forward paired sample in fastq format
- -2
- sampleID_2.fastq - name and path to reverse paired sample in fastq format
- -b - output in bam format

All post alignment processing described above was completed for the brain cortex, and whole blood tissues that were aligned to both the default genome and to the reference genome informed on the sex chromosome complement of the subject.

## 8. Generate stats on initial BAM files
`bamtools stats -in sampleID.bam > sampleID.txt`

- bamtools - package 
- stats	- command to get general - alignment statistics
- -in - indicates input file
- sampleID.bam	- path and name to bam file 
- ">" - directs output     
- sampleID_pass2.txt - indicated output file name

Will print basic statistics from input BAM file(s)
- Total reads:       
- Mapped reads:      
- Forward strand:    
- Reverse strand:    
- Failed QC:         
- Duplicates:        
- Paired-end reads:
- 'Proper-pairs':    
- Both pairs mapped: 
- Read 1:            
- Read 2:            
- Singletons:   

## 9. Sort BAM files
For each sample, sort the BAM file because BAM files are compressed. Sorting helps to give a better compression ratio because similar sequences are grouped together. An appropriate @HD-SO sort order header tag will be added or an existing one updated if necessary.

`bamtools sort -in sampleID.bam -out sampleID.sorted.bam`

- bamtools - package 
- sort - command to add header tags
- -in - indicates input file
- sampleID.bam	- path and name to bam file 
- -out - indicates output file
- sampleID.sorted.bam - output file

## 10. Mark duplicates
Mark duplicates: "Flags" where the duplicate reads are

`java -Xmx8g -jar picard.jar MarkDuplicates INPUT=sampleID.sorted.bam OUTPUT=sampleID.sorted.markdup.bam METRICS_FILE=sampleID.markdup.picardMetrics.txt REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT`

- java - program
- -Xmx8g - declares memory 
- picard.jar - path to picard jar file 
- MarkDuplicates - command to create sequence dictionary 
- INPUT= - path to input file, sorted bam file per sample	
- OUTPUT= - path and name or output file .markdup to indicate this file will contain duplicates that have been marked
- METRICS_FILE=	- file to write duplication metrics to save as sampleID.markdup.picardMetrics.txt
- REMOVE_DUPLICATES=false - If true do not write duplicates to the output file instead of writing them with appropriate flags set. Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false}
- ASSUME_SORTED=true - BAM files are sorted because we sorted them in step 6
- VALIDATION_STRINGENCY=LENIENT	- setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.

## 11. Add or replace read groups
For each sample, add a read group to the mark duplicate BAM files (a read group is a "tag" such as a sample ID)

`java -Xmx8g -jar picard.jar AddOrReplaceReadGroups INPUT=sampleID.sorted.markdup.bam OUTPUT=sampleID.sorted.markdup.addReadGr.bam RGLB=sampleID RGPL=machineUsed RGPU=laneUsed RGSM=sampleName RGCN=location RGDS=species VALIDATION_STRINGENCY=LENIENT`

- java - program called 
- Xmx8g - declares memory 
- picard.jar - path to picard jar file 
- AddOrReplaceReadGroups - Replaces all read groups in the INPUT file with a single new read group and assigns all reads to this read group in the OUTPUT BAM
- INPUT= - path to input file, sorted bam file per sample	
- OUTPUT= - path and name or output file .markdup to indicate this file will contain duplicates that have been marked
- RGLB=	- Read Group Library Required (sampleID)
- RGID=	- Read Groupsample ID
- RGPL=	- Read Group platform (e.g. illumina, solid) Required 
- RGPU=	- Read Group platform unit (eg. run barcode) Required	(laneUsed)
- RGSM=	- Read Group sample name Required (sampleName or sampleID)
- RGCN= - Read Group sequencing center name Default value: null (i.e. ASU)
- RGDS=	- Read Group description Default value: null (speciesName)
- VALIDATION_STRINGENCY=LENIENT

## 12. Index BAM files
For each sample index the processed BAM files that are sorted, have marked duplicates, and have read groups added. These will be used to identify callable loci. Indexing is used to "sort" by chromosome and region. Output will be sampleID.sorted.markdup.addReadGr.bam.bai

`bamtools index -in sampleID.sorted.markdup.addReadGr.bam`

- bamtools - package 
- index	- Generates index for BAM file
- -in - indicates input file
- sampleID.sorted.markdup.addReadGrbam - output file

## 13. Generate stats on final processed bam files
For each sample get the read stats for the remove duplicates and add read groups BAM files. Statistics on the BAM files should be the same as before the previous step when read groups were modified. Compare stat results for each sample to the markdup.bam file (sanity check: is there the same number of reads as the original BAM file?)

`bamtools stats -in sampleID.sorted.markdup.addReadGr.bam`

- bamtools - package 
- stats	- command to get general - alignment statistics
- -in - indicates input file
- sampleID.bam	- path and name to bam file 
- ">" - directs output     
- sampleID_pass2.txt - indicated output file name

## 14. Get raw transcript counts
For each sample get the read stats for the remove duplicates and add read groups BAM files. Statistics on the BAM files should be the same as before the previous step when read groups were modified. Compare stat results for each sample to the markdup.bam file (sanity check: is there the same number of reads as the original BAM file?)

`featureCounts -T 8 --primary -p -s 0 -t gene -g gene_name -a Homo_sapiens.GRCh38.89.gtf -o sampleID_FC.txt sampleID.sorted.markdup.addReadGr.bam`

- featureCounts - package 
- -T 8- command to get general - alignment statistics
- --primary
- -p 
- -s 0    
- -t gene
- -g gene_name
- -a Homo_sapiens.GRCh38.89.gtf


## 15. Create gene chromosome CSV file

create a file containing each gene of interest, and which chromosome(S) it is located on. There will be one file total. 

- gene	Chr
- DDX11L1	1
- WASH7P	1
- MIR6859-1	1
- MIR1302-2HG	1


## 16. Create phenotype CSV file

create a file containing each sampleID, sex of the sample, genome the sample is aligned to, and which aligner was used. Thre will be a different file for each tissue used. 

- sampleID	sex	genome	aligner
- SRRID  male  default HI
- SRRID  male  SS  HI
- SRRID  male  default STAR
- SRRID  male  SS  STAR
- SRRID  female  default HI
- SRRID  female  SS  HI
- SRRID  female  default STAR
- SRRID  female  SS   STAR

## 17. Create counts TSV file

Create a file containing each sampleID.sorted.markdup.addReadGr.bam file for a specific tissue. Use all sample for that specific tissue including male, female, STAR, HISAT, whole genome, and sex specific. there will be a different file for each tissue. 

- sampleID1 sampleID2 sampleID3 SampleID4
- 5 10  4 3
- 0 0 0 0
- 3 0 6 8
- 8 0 6 4

## 18. Differential expression using LimmaVoom
Designed to assign mapped reads or fragments from pair-end genomic features from genes, exons, and promoters, featureCounts with the limma/voom (Law et al. 2014) differential expression pipeline is highly rated as one of the best-performing pipelines for the analyses of RNAseq data (SEQC/MAQC-III Consortium 2014) and was therefore chosen for our analysis.  

A gene-level information file associated with the rows of the counts matrix was created using the Homo_sapiens.GRCh38.89.gtf gene annotation file, which was used in the subread featureCounts to generate the gene count data, the gene-level.csv file contains unique gene ids for each row and the corresponding chromosome location of the gene. The gene order is the same in both the annotation Homo_sapiens.GRCh38.89.gtf and the DGEList gene-level.csv data object. 

The limma/voom vebayesfit is an empirical Bayes moderation that takes information from all genes to obtain a more precise estimate of gene-wise variability and is recommended for RNAseq analysis (Law et al. 2014). 

LimmaVoom uses an R script for  differential expression analysis. The R scripts used in this analysis are located in the R folder. 
