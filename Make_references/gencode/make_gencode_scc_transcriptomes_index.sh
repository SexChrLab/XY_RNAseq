#!/bin/bash
#SBATCH --job-name=HISAT_BrainCortex_snakemake # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=END,FAIL           # notifications for job done & fail
#SBATCH --mail-user=kcolney@asu.edu # send-to address
#SBATCH -n 12
#SBATCH -t 96:00:00
##SBATCH --qos=mwilsons

source activate XY_RNA-Seq
newgrp wilsonlab

# Y masked gencode transcriptome
grep "^>" <(cat GRCh38.p12.genome.XXonly.fa) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt

cat gencode.v29.transcripts_Ymasked_XX.fa GRCh38.p12.genome.XXonly.fa > gentrome.fa.gz

salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i gencode_salmon_index_XXonly --gencode

rm decoys.txt
rm gentrome.fa.gz

# YPARs masked gencode transcriptome
grep "^>" <(cat GRCh38.p12.genome.XYonly.fa) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt

cat gencode.v29.transcripts_YPARs_masked_XY.fa GRCh38.p12.genome.XYonly.fa > gentrome.fa.gz


