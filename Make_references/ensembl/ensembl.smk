rule ensembl_fasta:
    """Download ensembl fasta files"""
    output:
        fastas=["Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz",
"Homo_sapiens.GRCh38.dna.chromosome.2.fa.gz",
"Homo_sapiens.GRCh38.dna.chromosome.3.fa.gz",
"Homo_sapiens.GRCh38.dna.chromosome.4.fa.gz",
"Homo_sapiens.GRCh38.dna.chromosome.5.fa.gz",
"Homo_sapiens.GRCh38.dna.chromosome.6.fa.gz",
"Homo_sapiens.GRCh38.dna.chromosome.7.fa.gz",
"Homo_sapiens.GRCh38.dna.chromosome.8.fa.gz",
"Homo_sapiens.GRCh38.dna.chromosome.9.fa.gz",
"Homo_sapiens.GRCh38.dna.chromosome.10.fa.gz",
"Homo_sapiens.GRCh38.dna.chromosome.11.fa.gz",
"Homo_sapiens.GRCh38.dna.chromosome.12.fa.gz",
"Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz",
"Homo_sapiens.GRCh38.dna.chromosome.14.fa.gz",
"Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz",
"Homo_sapiens.GRCh38.dna.chromosome.16.fa.gz",
"Homo_sapiens.GRCh38.dna.chromosome.17.fa.gz",
"Homo_sapiens.GRCh38.dna.chromosome.18.fa.gz",
"Homo_sapiens.GRCh38.dna.chromosome.19.fa.gz",
"Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz",
"Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz",
"Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz",
"Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz",
"Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz",
"Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz",
"Homo_sapiens.GRCh38.dna.nonchromosomal.fa.gz"
]
    shell:
        """
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.2.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.3.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.4.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.5.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.6.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.7.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.8.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.9.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.10.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.11.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.12.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.14.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.16.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.17.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.18.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.19.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz;
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.nonchromosomal.fa.gz
"""

# YPARs-masked
rule ypars_masked_reference:
    """Make the ypars masked reference"""
    input: rules.ensembl_fasta.output.fastas
    output: "GRCh38_YPARs_masked_reference.fa"
    shell:
        "cat Homo_sapiens.GRCh38.dna.chromosome.1.fa Homo_sapiens.GRCh38.dna.chromosome.2.fa Homo_sapiens.GRCh38.dna.chromosome.3.fa Homo_sapiens.GRCh38.dna.chromosome.4.fa Homo_sapiens.GRCh38.dna.chromosome.5.fa Homo_sapiens.GRCh38.dna.chromosome.6.fa Homo_sapiens.GRCh38.dna.chromosome.7.fa Homo_sapiens.GRCh38.dna.chromosome.8.fa Homo_sapiens.GRCh38.dna.chromosome.9.fa Homo_sapiens.GRCh38.dna.chromosome.10.fa Homo_sapiens.GRCh38.dna.chromosome.11.fa Homo_sapiens.GRCh38.dna.chromosome.12.fa Homo_sapiens.GRCh38.dna.chromosome.13.fa Homo_sapiens.GRCh38.dna.chromosome.14.fa Homo_sapiens.GRCh38.dna.chromosome.15.fa Homo_sapiens.GRCh38.dna.chromosome.16.fa Homo_sapiens.GRCh38.dna.chromosome.17.fa Homo_sapiens.GRCh38.dna.chromosome.18.fa Homo_sapiens.GRCh38.dna.chromosome.19.fa Homo_sapiens.GRCh38.dna.chromosome.20.fa Homo_sapiens.GRCh38.dna.chromosome.21.fa Homo_sapiens.GRCh38.dna.chromosome.22.fa Homo_sapiens.GRCh38.dna.chromosome.X.fa Homo_sapiens.GRCh38.dna.chromosome.Y.fa Homo_sapiens.GRCh38.dna.chromosome.MT.fa Homo_sapiens.GRCh38.dna.nonchromosomal.fa > GRCh38_YPARs_masked_reference.fa"

rule ymasked_mask:
    """Mask the Y chrom # hard-mask all of the nucleotides by changing to N for the Y chromosome"""
    input: "Homo_sapiens.GRCh38.dna.chromosome.Y.fa"
    output: fa="Homo_sapiens.GRCh38.dna.chromosome.YMASKED.fa"
    shell:
        """sed -e 's/[ATCG]/N/g' Homo_sapiens.GRCh38.dna.chromosome.Y.fa > Homo_sapiens.GRCh38.dna.chromosome.YMASKED.fa"""

rule ymasked_cat:
    """Ymasked"""
    input: rules.ymasked_mask.output.fa, rules.ensembl_fasta.output.fastas
    output: "GRCh38_Ymasked_reference.fa"
    shell:
        """cat Homo_sapiens.GRCh38.dna.chromosome.1.fa Homo_sapiens.GRCh38.dna.chromosome.2.fa Homo_sapiens.GRCh38.dna.chromosome.3.fa Homo_sapiens.GRCh38.dna.chromosome.4.fa Homo_sapiens.GRCh38.dna.chromosome.5.fa Homo_sapiens.GRCh38.dna.chromosome.6.fa Homo_sapiens.GRCh38.dna.chromosome.7.fa Homo_sapiens.GRCh38.dna.chromosome.8.fa Homo_sapiens.GRCh38.dna.chromosome.9.fa Homo_sapiens.GRCh38.dna.chromosome.10.fa Homo_sapiens.GRCh38.dna.chromosome.11.fa Homo_sapiens.GRCh38.dna.chromosome.12.fa Homo_sapiens.GRCh38.dna.chromosome.13.fa Homo_sapiens.GRCh38.dna.chromosome.14.fa Homo_sapiens.GRCh38.dna.chromosome.15.fa Homo_sapiens.GRCh38.dna.chromosome.16.fa Homo_sapiens.GRCh38.dna.chromosome.17.fa Homo_sapiens.GRCh38.dna.chromosome.18.fa Homo_sapiens.GRCh38.dna.chromosome.19.fa Homo_sapiens.GRCh38.dna.chromosome.20.fa Homo_sapiens.GRCh38.dna.chromosome.21.fa Homo_sapiens.GRCh38.dna.chromosome.22.fa Homo_sapiens.GRCh38.dna.chromosome.X.fa Homo_sapiens.GRCh38.dna.chromosome.YMASKED.fa Homo_sapiens.GRCh38.dna.chromosome.MT.fa Homo_sapiens.GRCh38.dna.nonchromosomal.fa > GRCh38_Ymasked_reference.fa"""

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


import requests, sys

def getchr(chr_text, start, end):
    """# Download individual regions of the Y chromosome
#------------------------------------------------------------------
# http://rest.ensembl.org/documentation/info/sequence_region (NOTE: Can use whichever programming language preferred)

# EXAMPLE: getYchr1-10000.py
"""
    server = "http://rest.ensembl.org"
    ext = f"/sequence/region/human/{chr}:{start}..{end}?coord_system_version=GRCh38"
    r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"})

    if not r.ok:
        r.raise_for_status()
	    sys.exit()
    return (r.text)

def savetofile(text, filename):
    open(filename, 'w').write(text)

# Get the Y chromosome regions and remove headers
#------------------------------------------------------------------
# NOTE: Only can pull regions of less than 10,000,000 base pairs
# NOTE: PAR regions for Ensembl are defined here: https://useast.ensembl.org/info/genome/genebuild/human_PARS.html

region_config = [("Y", 1, 10000, "Ychr1-10000_noHeader.fa"),
                 ("Y", 10001, 2781479, "Ychr10001-2781479_noHeader.fa"), #PAR1 region
                  ("Y", 2781480, 12781479, "Ychr2781480-12781479_noHeader.fa"),
                  ("Y", 12781480, 22781479, "Ychr12781480-22781479_noHeader.fa"),
                  ("Y", 22781480, 32781479, "Ychr22781480-32781479_noHeader.fa"),
                  ("Y", 32781480, 42781479, "Ychr32781480-42781479_noHeader.fa"),
                  ("Y", 42781480, 52781479, "Ychr42781480-52781479_noHeader.fa"),
                  ("Y", 52781480, 56887902, "Ychr52781480-56887902_noHeader.fa"),
                  ("Y", 56887903, 57217415, "Ychr56887903-57217415_noHeader.fa"), #PAR2 region
                  ("Y", 57217416, 57227415, "Ychr57217416-57227415_noHeader.fa")] #PAR2 region

rule generate_region_files_with_header:
    """python getYchr1-10000.py | sed -e '1d' > Ychr1-10000_noHeader.fa #PAR1 region
python getYchr10001-2781479.py | sed -e '1d' > Ychr10001-2781479_noHeader.fa #PAR1 region
python getYchr2781480-12781479.py | sed -e '1d' > Ychr2781480-12781479_noHeader.fa
python getYchr12781480-22781479.py | sed -e '1d' > Ychr12781480-22781479_noHeader.fa
python getYchr22781480-32781479.py | sed -e '1d' > Ychr22781480-32781479_noHeader.fa
python getYchr32781480-42781479.py | sed -e '1d' > Ychr32781480-42781479_noHeader.fa
python getYchr42781480-52781479.py | sed -e '1d' > Ychr42781480-52781479_noHeader.fa
python getYchr52781480-56887902.py | sed -e '1d' > Ychr52781480-56887902_noHeader.fa
python getYchr56887903-57217415.py | sed -e '1d' > Ychr56887903-57217415_noHeader.fa #PAR2 region
python getYchr57217416-57227415.py | sed -e '1d' > Ychr57217416-57227415_noHeader.fa #PAR2 region"""
    output:
        actually_with_headers=[filename+"_butactually_With_header" for chr_text,start,end,filename in region_config]
    run:
        for chr_text,start,end,filename in region_config
            text = getchr(chr_text, start, end)
            savetofile(text, filename+"_butactually_with_header")

rule generate_region_files_without_header:
    """generate the region files without header"""
    input:
        rules.generate_region_files_with_header.output.actually_with_headers
    output:
        without_headers=[filename for chr_text,start,end,filename in region_config]
    shell:
        ";".join([f"cat {filename}_butactually_with_header | sed -e '1d' > {filename}" for chr_text,start,end,filename in region_config])

rule create_y_chromosome:
    """# Create the Y chromosome (including PARs)
#------------------------------------------------------------------

# first cat all files together (this will leave you with line breaks for each file)"""
    input: rules.generate_region_files_without_header.output.without_headers
    output: "YchrWhole.fa"
    shell:
        "cat Ychr1-10000_noHeader.fa Ychr10001-2781479_noHeader.fa Ychr2781480-12781479_noHeader.fa Ychr12781480-22781479_noHeader.fa Ychr22781480-32781479_noHeader.fa Ychr32781480-42781479_noHeader.fa Ychr42781480-52781479_noHeader.fa Ychr52781480-56887902_noHeader.fa Ychr56887903-57217415_noHeader.fa Ychr57217416-57227415_noHeader.fa > YchrWhole.fa"

rule single_line:
    """# Turn file into single line to remove whitespace"""
    input: "YchrWhole.fa"
    output: "YchrWhole_single.fa"
    shell:
        """awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' YchrWhole.fa > YchrWhole_single.fa"""

rule fold:
    """# return to multiple lines of 60 for fasta format"""
    input: "YchrWhole_single.fa"
    output: "Homo_sapiens.GRCh38.dna.chromosome.Y.fa"
    shell: """fold -w 60 YchrWhole_single.fa > Homo_sapiens.GRCh38.dna.chromosome.Y.fa"""


# Create the Y chromosome with PARs hard masked (YPARs-masked)
#------------------------------------------------------------------

rule mask_pars:
    """# hard-mask all of the nucleotides by changing to N for both PAR regions"""
    input:
        "Ychr10001-2781479_noHeader.fa",
        "Ychr56887903-57217415_noHeader.fa"
    output:
        "PAR1_masked.fa",
        "PAR2_masked.fa"
    shell:
        """sed -e 's/[ATCG]/N/g' Ychr10001-2781479_noHeader.fa > PAR1_masked.fa
sed -e 's/[ATCG]/N/g' Ychr56887903-57217415_noHeader.fa > PAR2_masked.fa"""

rule cat_all_files:
    """
# cat all files together with the PARs masked (this will leave you with line breaks for each file)"""
    input: ["Ychr1-10000_noHeader.fa", "PAR1_masked.fa", "Ychr2781480-12781479_noHeader.fa", "Ychr12781480-22781479_noHeader.fa", "Ychr22781480-32781479_noHeader.fa", "Ychr32781480-42781479_noHeader.fa", "Ychr42781480-52781479_noHeader.fa", "Ychr52781480-56887902_noHeader.fa", "PAR2_masked.fa", "Ychr57217416-57227415_noHeader.fa"]
    output: "YchrPARs-masked.fa"
    shell:
        """cat Ychr1-10000_noHeader.fa PAR1_masked.fa Ychr2781480-12781479_noHeader.fa Ychr12781480-22781479_noHeader.fa Ychr22781480-32781479_noHeader.fa Ychr32781480-42781479_noHeader.fa Ychr42781480-52781479_noHeader.fa Ychr52781480-56887902_noHeader.fa PAR2_masked.fa Ychr57217416-57227415_noHeader.fa > YchrPARs-masked.fa"""

rule single_line:
    """# Turn file into single line to remove whitespace"""
    input: "YchrPARs-masked.fa"
    output: "YchrPARs-masked_single.fa"
    shell: """awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' YchrPARs-masked.fa > YchrPARs-masked_single.fa"""

rule fold_ychrmasked:
    """# return to multiple lines of 60 for fasta format"""
    input: "YchrPARs-masked_single.fa"
    output: "Homo_sapiens.GRCh38.dna_rmYPARs.chromosome.Y.fa"
    shell: """fold -w 60 YchrPARs-masked_single.fa > Homo_sapiens.GRCh38.dna_rmYPARs.chromosome.Y.fa"""

# Create the Y chromosome hard masked (Ymasked)
#------------------------------------------------------------------

rule hardmask_ychr:
    """# hard-mask all of the nucleotides by changing to N for the Y chromosome"""
    input: "Homo_sapiens.GRCh38.dna.chromosome.Y.fa"
    output: "Homo_sapiens.GRCh38.dna_rmY.chromosome.Y.fa"
    shell:
        """sed -e 's/[ATCG]/N/g' Homo_sapiens.GRCh38.dna.chromosome.Y.fa > Homo_sapiens.GRCh38.dna_rmY.chromosome.Y.fa"""

# Create the reference genomes
#------------------------------------------------------------------

rule create_reference_genomes:
    """# Default Genome"""
    input: ["Homo_sapiens.GRCh38.dna.chromosome.1.fa", "Homo_sapiens.GRCh38.dna.chromosome.2.fa", "Homo_sapiens.GRCh38.dna.chromosome.3.fa", "Homo_sapiens.GRCh38.dna.chromosome.4.fa", "Homo_sapiens.GRCh38.dna.chromosome.5.fa", "Homo_sapiens.GRCh38.dna.chromosome.6.fa", "Homo_sapiens.GRCh38.dna.chromosome.7.fa", "Homo_sapiens.GRCh38.dna.chromosome.8.fa", "Homo_sapiens.GRCh38.dna.chromosome.9.fa", "Homo_sapiens.GRCh38.dna.chromosome.10.fa", "Homo_sapiens.GRCh38.dna.chromosome.11.fa", "Homo_sapiens.GRCh38.dna.chromosome.12.fa", "Homo_sapiens.GRCh38.dna.chromosome.13.fa", "Homo_sapiens.GRCh38.dna.chromosome.14.fa", "Homo_sapiens.GRCh38.dna.chromosome.15.fa", "Homo_sapiens.GRCh38.dna.chromosome.16.fa", "Homo_sapiens.GRCh38.dna.chromosome.17.fa", "Homo_sapiens.GRCh38.dna.chromosome.18.fa", "Homo_sapiens.GRCh38.dna.chromosome.19.fa", "Homo_sapiens.GRCh38.dna.chromosome.20.fa", "Homo_sapiens.GRCh38.dna.chromosome.21.fa", "Homo_sapiens.GRCh38.dna.chromosome.22.fa", "Homo_sapiens.GRCh38.dna.chromosome.X.fa", "Homo_sapiens.GRCh38.dna.chromosome.Y.fa", "Homo_sapiens.GRCh38.dna.chromosome.MT.fa", "Homo_sapiens.GRCh38.dna.nonchromosomal.fa"]
    output: "GRCh38_wholeGenome_reference.fa"
    shell:
        """cat Homo_sapiens.GRCh38.dna.chromosome.1.fa Homo_sapiens.GRCh38.dna.chromosome.2.fa Homo_sapiens.GRCh38.dna.chromosome.3.fa Homo_sapiens.GRCh38.dna.chromosome.4.fa Homo_sapiens.GRCh38.dna.chromosome.5.fa Homo_sapiens.GRCh38.dna.chromosome.6.fa Homo_sapiens.GRCh38.dna.chromosome.7.fa Homo_sapiens.GRCh38.dna.chromosome.8.fa Homo_sapiens.GRCh38.dna.chromosome.9.fa Homo_sapiens.GRCh38.dna.chromosome.10.fa Homo_sapiens.GRCh38.dna.chromosome.11.fa Homo_sapiens.GRCh38.dna.chromosome.12.fa Homo_sapiens.GRCh38.dna.chromosome.13.fa Homo_sapiens.GRCh38.dna.chromosome.14.fa Homo_sapiens.GRCh38.dna.chromosome.15.fa Homo_sapiens.GRCh38.dna.chromosome.16.fa Homo_sapiens.GRCh38.dna.chromosome.17.fa Homo_sapiens.GRCh38.dna.chromosome.18.fa Homo_sapiens.GRCh38.dna.chromosome.19.fa Homo_sapiens.GRCh38.dna.chromosome.20.fa Homo_sapiens.GRCh38.dna.chromosome.21.fa Homo_sapiens.GRCh38.dna.chromosome.22.fa Homo_sapiens.GRCh38.dna.chromosome.X.fa Homo_sapiens.GRCh38.dna.chromosome.Y.fa Homo_sapiens.GRCh38.dna.chromosome.MT.fa Homo_sapiens.GRCh38.dna.nonchromosomal.fa > GRCh38_wholeGenome_reference.fa"""

rule create_reference_genome_ypars_masked:
    """# YPARs-masked"""
    input: ["Homo_sapiens.GRCh38.dna.chromosome.1.fa", "Homo_sapiens.GRCh38.dna.chromosome.2.fa", "Homo_sapiens.GRCh38.dna.chromosome.3.fa", "Homo_sapiens.GRCh38.dna.chromosome.4.fa", "Homo_sapiens.GRCh38.dna.chromosome.5.fa", "Homo_sapiens.GRCh38.dna.chromosome.6.fa", "Homo_sapiens.GRCh38.dna.chromosome.7.fa", "Homo_sapiens.GRCh38.dna.chromosome.8.fa", "Homo_sapiens.GRCh38.dna.chromosome.9.fa", "Homo_sapiens.GRCh38.dna.chromosome.10.fa", "Homo_sapiens.GRCh38.dna.chromosome.11.fa", "Homo_sapiens.GRCh38.dna.chromosome.12.fa", "Homo_sapiens.GRCh38.dna.chromosome.13.fa", "Homo_sapiens.GRCh38.dna.chromosome.14.fa", "Homo_sapiens.GRCh38.dna.chromosome.15.fa", "Homo_sapiens.GRCh38.dna.chromosome.16.fa", "Homo_sapiens.GRCh38.dna.chromosome.17.fa", "Homo_sapiens.GRCh38.dna.chromosome.18.fa", "Homo_sapiens.GRCh38.dna.chromosome.19.fa", "Homo_sapiens.GRCh38.dna.chromosome.20.fa", "Homo_sapiens.GRCh38.dna.chromosome.21.fa", "Homo_sapiens.GRCh38.dna.chromosome.22.fa", "Homo_sapiens.GRCh38.dna.chromosome.X.fa", "Homo_sapiens.GRCh38.dna_rmYPARs.chromosome.Y.fa", "Homo_sapiens.GRCh38.dna.chromosome.MT.fa", "Homo_sapiens.GRCh38.dna.nonchromosomal.fa"]
    output: "GRCh38_YPARs_masked_reference.fa"
    shell:
        """cat Homo_sapiens.GRCh38.dna.chromosome.1.fa Homo_sapiens.GRCh38.dna.chromosome.2.fa Homo_sapiens.GRCh38.dna.chromosome.3.fa Homo_sapiens.GRCh38.dna.chromosome.4.fa Homo_sapiens.GRCh38.dna.chromosome.5.fa Homo_sapiens.GRCh38.dna.chromosome.6.fa Homo_sapiens.GRCh38.dna.chromosome.7.fa Homo_sapiens.GRCh38.dna.chromosome.8.fa Homo_sapiens.GRCh38.dna.chromosome.9.fa Homo_sapiens.GRCh38.dna.chromosome.10.fa Homo_sapiens.GRCh38.dna.chromosome.11.fa Homo_sapiens.GRCh38.dna.chromosome.12.fa Homo_sapiens.GRCh38.dna.chromosome.13.fa Homo_sapiens.GRCh38.dna.chromosome.14.fa Homo_sapiens.GRCh38.dna.chromosome.15.fa Homo_sapiens.GRCh38.dna.chromosome.16.fa Homo_sapiens.GRCh38.dna.chromosome.17.fa Homo_sapiens.GRCh38.dna.chromosome.18.fa Homo_sapiens.GRCh38.dna.chromosome.19.fa Homo_sapiens.GRCh38.dna.chromosome.20.fa Homo_sapiens.GRCh38.dna.chromosome.21.fa Homo_sapiens.GRCh38.dna.chromosome.22.fa Homo_sapiens.GRCh38.dna.chromosome.X.fa Homo_sapiens.GRCh38.dna_rmYPARs.chromosome.Y.fa Homo_sapiens.GRCh38.dna.chromosome.MT.fa Homo_sapiens.GRCh38.dna.nonchromosomal.fa > GRCh38_YPARs_masked_reference.fa"""

rule create_reference_genome_y_masked:
    """# Ymasked"""
    input: ["Homo_sapiens.GRCh38.dna.chromosome.1.fa", "Homo_sapiens.GRCh38.dna.chromosome.2.fa", "Homo_sapiens.GRCh38.dna.chromosome.3.fa", "Homo_sapiens.GRCh38.dna.chromosome.4.fa", "Homo_sapiens.GRCh38.dna.chromosome.5.fa", "Homo_sapiens.GRCh38.dna.chromosome.6.fa", "Homo_sapiens.GRCh38.dna.chromosome.7.fa", "Homo_sapiens.GRCh38.dna.chromosome.8.fa", "Homo_sapiens.GRCh38.dna.chromosome.9.fa", "Homo_sapiens.GRCh38.dna.chromosome.10.fa", "Homo_sapiens.GRCh38.dna.chromosome.11.fa", "Homo_sapiens.GRCh38.dna.chromosome.12.fa", "Homo_sapiens.GRCh38.dna.chromosome.13.fa", "Homo_sapiens.GRCh38.dna.chromosome.14.fa", "Homo_sapiens.GRCh38.dna.chromosome.15.fa", "Homo_sapiens.GRCh38.dna.chromosome.16.fa", "Homo_sapiens.GRCh38.dna.chromosome.17.fa", "Homo_sapiens.GRCh38.dna.chromosome.18.fa", "Homo_sapiens.GRCh38.dna.chromosome.19.fa", "Homo_sapiens.GRCh38.dna.chromosome.20.fa", "Homo_sapiens.GRCh38.dna.chromosome.21.fa", "Homo_sapiens.GRCh38.dna.chromosome.22.fa", "Homo_sapiens.GRCh38.dna.chromosome.X.fa", "Homo_sapiens.GRCh38.dna_rmY.chromosome.Y.fa", "Homo_sapiens.GRCh38.dna.chromosome.MT.fa", "Homo_sapiens.GRCh38.dna.nonchromosomal.fa"]
    output: "GRCh38_Ymasked_reference.fa"
    shell:
        """cat Homo_sapiens.GRCh38.dna.chromosome.1.fa Homo_sapiens.GRCh38.dna.chromosome.2.fa Homo_sapiens.GRCh38.dna.chromosome.3.fa Homo_sapiens.GRCh38.dna.chromosome.4.fa Homo_sapiens.GRCh38.dna.chromosome.5.fa Homo_sapiens.GRCh38.dna.chromosome.6.fa Homo_sapiens.GRCh38.dna.chromosome.7.fa Homo_sapiens.GRCh38.dna.chromosome.8.fa Homo_sapiens.GRCh38.dna.chromosome.9.fa Homo_sapiens.GRCh38.dna.chromosome.10.fa Homo_sapiens.GRCh38.dna.chromosome.11.fa Homo_sapiens.GRCh38.dna.chromosome.12.fa Homo_sapiens.GRCh38.dna.chromosome.13.fa Homo_sapiens.GRCh38.dna.chromosome.14.fa Homo_sapiens.GRCh38.dna.chromosome.15.fa Homo_sapiens.GRCh38.dna.chromosome.16.fa Homo_sapiens.GRCh38.dna.chromosome.17.fa Homo_sapiens.GRCh38.dna.chromosome.18.fa Homo_sapiens.GRCh38.dna.chromosome.19.fa Homo_sapiens.GRCh38.dna.chromosome.20.fa Homo_sapiens.GRCh38.dna.chromosome.21.fa Homo_sapiens.GRCh38.dna.chromosome.22.fa Homo_sapiens.GRCh38.dna.chromosome.X.fa Homo_sapiens.GRCh38.dna_rmY.chromosome.Y.fa Homo_sapiens.GRCh38.dna.chromosome.MT.fa Homo_sapiens.GRCh38.dna.nonchromosomal.fa > GRCh38_Ymasked_reference.fa"""
