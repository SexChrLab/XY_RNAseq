#!/bin/bash
#################
# This script will download and rreate the following GRCh38 human reference genomes:
#	YPARs-masked GRCh38 1-22,MT,X,YPARsN,nonchromosomal regions
#	Y-Masked GRCh38 1-22,MT,X,YN,nonchromosomal regions
#	
# Will need access to the internet
# The use of following bash commands: 
#	wget 
#	cat 
#	sed
################

# Download Ensembl GRCh38 reference genome by chromosome
#------------------------------------------------------------------
# NOTE: Ensembl reference genome has the PARs masked in the Y chromosome (https://useast.ensembl.org/info/genome/genebuild/human_PARS.html)

# Repeat for all chromosomes, including autosomes 1-22, X, Y, MT, and nonchromosomal
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.2.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.3.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.4.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.5.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.6.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.7.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.8.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.9.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.10.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.11.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.12.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.14.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.16.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.17.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.18.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.19.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.nonchromosomal.fa.gz


# YPARs-masked
cat Homo_sapiens.GRCh38.dna.chromosome.1.fa Homo_sapiens.GRCh38.dna.chromosome.2.fa Homo_sapiens.GRCh38.dna.chromosome.3.fa Homo_sapiens.GRCh38.dna.chromosome.4.fa Homo_sapiens.GRCh38.dna.chromosome.5.fa Homo_sapiens.GRCh38.dna.chromosome.6.fa Homo_sapiens.GRCh38.dna.chromosome.7.fa Homo_sapiens.GRCh38.dna.chromosome.8.fa Homo_sapiens.GRCh38.dna.chromosome.9.fa Homo_sapiens.GRCh38.dna.chromosome.10.fa Homo_sapiens.GRCh38.dna.chromosome.11.fa Homo_sapiens.GRCh38.dna.chromosome.12.fa Homo_sapiens.GRCh38.dna.chromosome.13.fa Homo_sapiens.GRCh38.dna.chromosome.14.fa Homo_sapiens.GRCh38.dna.chromosome.15.fa Homo_sapiens.GRCh38.dna.chromosome.16.fa Homo_sapiens.GRCh38.dna.chromosome.17.fa Homo_sapiens.GRCh38.dna.chromosome.18.fa Homo_sapiens.GRCh38.dna.chromosome.19.fa Homo_sapiens.GRCh38.dna.chromosome.20.fa Homo_sapiens.GRCh38.dna.chromosome.21.fa Homo_sapiens.GRCh38.dna.chromosome.22.fa Homo_sapiens.GRCh38.dna.chromosome.X.fa Homo_sapiens.GRCh38.dna.chromosome.Y.fa Homo_sapiens.GRCh38.dna.chromosome.MT.fa Homo_sapiens.GRCh38.dna.nonchromosomal.fa > GRCh38_YPARs_masked_reference.fa

# hard-mask all of the nucleotides by changing to N for the Y chromosome
sed -e 's/[ATCG]/N/g' Homo_sapiens.GRCh38.dna.chromosome.Y.fa > Homo_sapiens.GRCh38.dna.chromosome.YMASKED.fa

# Ymasked
cat Homo_sapiens.GRCh38.dna.chromosome.1.fa Homo_sapiens.GRCh38.dna.chromosome.2.fa Homo_sapiens.GRCh38.dna.chromosome.3.fa Homo_sapiens.GRCh38.dna.chromosome.4.fa Homo_sapiens.GRCh38.dna.chromosome.5.fa Homo_sapiens.GRCh38.dna.chromosome.6.fa Homo_sapiens.GRCh38.dna.chromosome.7.fa Homo_sapiens.GRCh38.dna.chromosome.8.fa Homo_sapiens.GRCh38.dna.chromosome.9.fa Homo_sapiens.GRCh38.dna.chromosome.10.fa Homo_sapiens.GRCh38.dna.chromosome.11.fa Homo_sapiens.GRCh38.dna.chromosome.12.fa Homo_sapiens.GRCh38.dna.chromosome.13.fa Homo_sapiens.GRCh38.dna.chromosome.14.fa Homo_sapiens.GRCh38.dna.chromosome.15.fa Homo_sapiens.GRCh38.dna.chromosome.16.fa Homo_sapiens.GRCh38.dna.chromosome.17.fa Homo_sapiens.GRCh38.dna.chromosome.18.fa Homo_sapiens.GRCh38.dna.chromosome.19.fa Homo_sapiens.GRCh38.dna.chromosome.20.fa Homo_sapiens.GRCh38.dna.chromosome.21.fa Homo_sapiens.GRCh38.dna.chromosome.22.fa Homo_sapiens.GRCh38.dna.chromosome.X.fa Homo_sapiens.GRCh38.dna_rmY.chromosome.YMASKED.fa Homo_sapiens.GRCh38.dna.chromosome.MT.fa Homo_sapiens.GRCh38.dna.nonchromosomal.fa > GRCh38_Ymasked_reference.fa

################

