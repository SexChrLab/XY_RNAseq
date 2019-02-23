###################################################################
# Create reference genomes
# 		Default GRCh38 1-22,MT,X,Y,nonchromosomal regions
#		Y-Masked GRCh38 1-22,MT,X,YN,nonchromosomal regions
#		YPARs-masked GRCh38 1-22,MT,X,YPARsN,nonchromosomal regions
#
###################################################################

# Download Ensembl GRCh38 reference genome by chromosome
#------------------------------------------------------------------
# NOTE: Ensembl reference genome has the PARs masked in the Y chromosome (https://useast.ensembl.org/info/genome/genebuild/human_PARS.html)

# Repeat for all chromosomes, including X, Y, MT, and nonchromosomal
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz

# Download individual regions of the Y chromosome
#------------------------------------------------------------------
# http://rest.ensembl.org/documentation/info/sequence_region (NOTE: Can use whichever programming language preferred)

# EXAMPLE: getYchr1-10000.py

	import requests, sys
 
	server = "http://rest.ensembl.org"
	ext = "/sequence/region/human/Y:1..10000?coord_system_version=GRCh38"
 
	r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"})
 
	if not r.ok:
	  r.raise_for_status()
	  sys.exit()
 
 
	print (r.text)

# Get the Y chromosome regions and remove headers
#------------------------------------------------------------------
# NOTE: Only can pull regions of less than 10,000,000 base pairs
# NOTE: PAR regions for Ensembl are defined here: https://useast.ensembl.org/info/genome/genebuild/human_PARS.html

python getYchr1-10000.py | sed -e '1d' > Ychr1-10000_noHeader.fa #PAR1 region
python getYchr10001-2781479.py | sed -e '1d' > Ychr10001-2781479_noHeader.fa #PAR1 region
python getYchr2781480-12781479.py | sed -e '1d' > Ychr2781480-12781479_noHeader.fa
python getYchr12781480-22781479.py | sed -e '1d' > Ychr12781480-22781479_noHeader.fa
python getYchr22781480-32781479.py | sed -e '1d' > Ychr22781480-32781479_noHeader.fa
python getYchr32781480-42781479.py | sed -e '1d' > Ychr32781480-42781479_noHeader.fa
python getYchr42781480-52781479.py | sed -e '1d' > Ychr42781480-52781479_noHeader.fa
python getYchr52781480-56887902.py | sed -e '1d' > Ychr52781480-56887902_noHeader.fa
python getYchr56887903-57217415.py | sed -e '1d' > Ychr56887903-57217415_noHeader.fa #PAR2 region
python getYchr57217416-57227415.py | sed -e '1d' > Ychr57217416-57227415_noHeader.fa #PAR2 region

# Create the Y chromosome (including PARs)
#------------------------------------------------------------------

# first cat all files together (this will leave you with line breaks for each file)
cat Ychr1-10000_noHeader.fa Ychr10001-2781479_noHeader.fa Ychr2781480-12781479_noHeader.fa Ychr12781480-22781479_noHeader.fa Ychr22781480-32781479_noHeader.fa Ychr32781480-42781479_noHeader.fa Ychr42781480-52781479_noHeader.fa Ychr52781480-56887902_noHeader.fa Ychr56887903-57217415_noHeader.fa Ychr57217416-57227415_noHeader.fa > YchrWhole.fa

# Turn file into single line to remove whitespace
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' YchrWhole.fa > YchrWhole_single.fa

# return to multiple lines of 60 for fasta format
fold -w 60 YchrWhole_single.fa > Homo_sapiens.GRCh38.dna.chromosome.Y.fa

# Create the Y chromosome with PARs hard masked (YPARs-masked)
#------------------------------------------------------------------

# hard-mask all of the nucleotides by changing to N for both PAR regions
sed -e 's/[ATCG]/N/g' Ychr10001-2781479_noHeader.fa > PAR1_masked.fa
sed -e 's/[ATCG]/N/g' Ychr56887903-57217415_noHeader.fa > PAR2_masked.fa

# cat all files together with the PARs masked (this will leave you with line breaks for each file)
cat Ychr1-10000_noHeader.fa PAR1_masked.fa Ychr2781480-12781479_noHeader.fa Ychr12781480-22781479_noHeader.fa Ychr22781480-32781479_noHeader.fa Ychr32781480-42781479_noHeader.fa Ychr42781480-52781479_noHeader.fa Ychr52781480-56887902_noHeader.fa PAR2_masked.fa Ychr57217416-57227415_noHeader.fa > YchrPARs-masked.fa

# Turn file into single line to remove whitespace
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' YchrPARs-masked.fa > YchrPARs-masked_single.fa

# return to multiple lines of 60 for fasta format
fold -w 60 YchrPARs-masked_single.fa > Homo_sapiens.GRCh38.dna_rmYPARs.chromosome.Y.fa

# Create the Y chromosome hard masked (Ymasked)
#------------------------------------------------------------------

# hard-mask all of the nucleotides by changing to N for the Y chromosome
sed -e 's/[ATCG]/N/g' Homo_sapiens.GRCh38.dna.chromosome.Y.fa > Homo_sapiens.GRCh38.dna_rmY.chromosome.Y.fa

# Create the reference genomes
#------------------------------------------------------------------

# Default Genome
cat Homo_sapiens.GRCh38.dna.chromosome.1.fa Homo_sapiens.GRCh38.dna.chromosome.2.fa Homo_sapiens.GRCh38.dna.chromosome.3.fa Homo_sapiens.GRCh38.dna.chromosome.4.fa Homo_sapiens.GRCh38.dna.chromosome.5.fa Homo_sapiens.GRCh38.dna.chromosome.6.fa Homo_sapiens.GRCh38.dna.chromosome.7.fa Homo_sapiens.GRCh38.dna.chromosome.8.fa Homo_sapiens.GRCh38.dna.chromosome.9.fa Homo_sapiens.GRCh38.dna.chromosome.10.fa Homo_sapiens.GRCh38.dna.chromosome.11.fa Homo_sapiens.GRCh38.dna.chromosome.12.fa Homo_sapiens.GRCh38.dna.chromosome.13.fa Homo_sapiens.GRCh38.dna.chromosome.14.fa Homo_sapiens.GRCh38.dna.chromosome.15.fa Homo_sapiens.GRCh38.dna.chromosome.16.fa Homo_sapiens.GRCh38.dna.chromosome.17.fa Homo_sapiens.GRCh38.dna.chromosome.18.fa Homo_sapiens.GRCh38.dna.chromosome.19.fa Homo_sapiens.GRCh38.dna.chromosome.20.fa Homo_sapiens.GRCh38.dna.chromosome.21.fa Homo_sapiens.GRCh38.dna.chromosome.22.fa Homo_sapiens.GRCh38.dna.chromosome.X.fa Homo_sapiens.GRCh38.dna.chromosome.Y.fa Homo_sapiens.GRCh38.dna.chromosome.MT.fa Homo_sapiens.GRCh38.dna.nonchromosomal.fa > GRCh38_wholeGenome_reference.fa

# YPARs-masked
cat Homo_sapiens.GRCh38.dna.chromosome.1.fa Homo_sapiens.GRCh38.dna.chromosome.2.fa Homo_sapiens.GRCh38.dna.chromosome.3.fa Homo_sapiens.GRCh38.dna.chromosome.4.fa Homo_sapiens.GRCh38.dna.chromosome.5.fa Homo_sapiens.GRCh38.dna.chromosome.6.fa Homo_sapiens.GRCh38.dna.chromosome.7.fa Homo_sapiens.GRCh38.dna.chromosome.8.fa Homo_sapiens.GRCh38.dna.chromosome.9.fa Homo_sapiens.GRCh38.dna.chromosome.10.fa Homo_sapiens.GRCh38.dna.chromosome.11.fa Homo_sapiens.GRCh38.dna.chromosome.12.fa Homo_sapiens.GRCh38.dna.chromosome.13.fa Homo_sapiens.GRCh38.dna.chromosome.14.fa Homo_sapiens.GRCh38.dna.chromosome.15.fa Homo_sapiens.GRCh38.dna.chromosome.16.fa Homo_sapiens.GRCh38.dna.chromosome.17.fa Homo_sapiens.GRCh38.dna.chromosome.18.fa Homo_sapiens.GRCh38.dna.chromosome.19.fa Homo_sapiens.GRCh38.dna.chromosome.20.fa Homo_sapiens.GRCh38.dna.chromosome.21.fa Homo_sapiens.GRCh38.dna.chromosome.22.fa Homo_sapiens.GRCh38.dna.chromosome.X.fa Homo_sapiens.GRCh38.dna_rmYPARs.chromosome.Y.fa Homo_sapiens.GRCh38.dna.chromosome.MT.fa Homo_sapiens.GRCh38.dna.nonchromosomal.fa > GRCh38_YPARs_masked_reference.fa

# Ymasked
cat Homo_sapiens.GRCh38.dna.chromosome.1.fa Homo_sapiens.GRCh38.dna.chromosome.2.fa Homo_sapiens.GRCh38.dna.chromosome.3.fa Homo_sapiens.GRCh38.dna.chromosome.4.fa Homo_sapiens.GRCh38.dna.chromosome.5.fa Homo_sapiens.GRCh38.dna.chromosome.6.fa Homo_sapiens.GRCh38.dna.chromosome.7.fa Homo_sapiens.GRCh38.dna.chromosome.8.fa Homo_sapiens.GRCh38.dna.chromosome.9.fa Homo_sapiens.GRCh38.dna.chromosome.10.fa Homo_sapiens.GRCh38.dna.chromosome.11.fa Homo_sapiens.GRCh38.dna.chromosome.12.fa Homo_sapiens.GRCh38.dna.chromosome.13.fa Homo_sapiens.GRCh38.dna.chromosome.14.fa Homo_sapiens.GRCh38.dna.chromosome.15.fa Homo_sapiens.GRCh38.dna.chromosome.16.fa Homo_sapiens.GRCh38.dna.chromosome.17.fa Homo_sapiens.GRCh38.dna.chromosome.18.fa Homo_sapiens.GRCh38.dna.chromosome.19.fa Homo_sapiens.GRCh38.dna.chromosome.20.fa Homo_sapiens.GRCh38.dna.chromosome.21.fa Homo_sapiens.GRCh38.dna.chromosome.22.fa Homo_sapiens.GRCh38.dna.chromosome.X.fa Homo_sapiens.GRCh38.dna_rmY.chromosome.Y.fa Homo_sapiens.GRCh38.dna.chromosome.MT.fa Homo_sapiens.GRCh38.dna.nonchromosomal.fa > GRCh38_Ymasked_reference.fa
