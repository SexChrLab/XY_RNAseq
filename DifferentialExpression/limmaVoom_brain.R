#------------------------------------------------------------------------------------
# Single factor DE analysis using limma/voom
#
#   brain 
#       
#
#------------------------------------------------------------------------------------
# To install packages use biocLite from bioconductor source
#source("http://bioconductor.org/biocLite.R")
#biocLite
#install.packages("BiocUpgrade")
#install.packages("Glimma")
#biocLite("Glimma")

#library(Glimma)
#library(UncertainInterval)
library(limma)
library(edgeR)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(devtools)
library(gplots) 
library(VennDiagram)

#--------------
# set working directory 
#--------------
setwd("~/brain/")

#-----------
# Read in count, genes, and phenotype data
#-----------
counts <- read.delim("brain_std_counts.tsv", header=FALSE, sep="\t")
genes <- read.csv("genes_chr.csv", header=TRUE)
chrX_genes <- read.csv("chrX_startEnd.txt", header=TRUE, sep="\t")
genes_chrX <- subset(genes, Chr=="X")
X_gene<-(chrX_genes$gene)
pheno <- read.csv("brain_pheno.csv", header=TRUE, sep=",")

# PAR is PAR1, XTR, and PAR2
PAR <- read.delim("~/01_code/geneLists/chrX_geneLists/X_PARsXTR_bybase.txt", header=TRUE, sep="\t")
PAR<-data.frame(PAR)
names(PAR)[1] <- "gene"
PAR_gene <- (PAR$gene)

PAR1 <- read.delim("~/01_code/geneLists/chrX_geneLists/X_PAR1_bybase.txt", header=FALSE, sep="\t")
PAR1<-data.frame(PAR1)
names(PAR1)[1] <- "gene"
PAR1_gene <- (PAR1$gene)

PAR2 <- read.delim("~/01_code/geneLists/chrX_geneLists/X_PAR2_bybase.txt", header=FALSE, sep="\t")
PAR2<-data.frame(PAR2)
names(PAR2)[1] <- "gene"
PAR2_gene <- (PAR2$gene)

XTR <- read.delim("~/01_code/geneLists/chrX_geneLists/X_XTR_bybase.txt", header=FALSE, sep="\t")
XTR<-data.frame(XTR)
names(XTR)[1] <- "gene"
XTR_gene <- (XTR$gene)

XDR <- read.delim("~/01_code/geneLists/chrX_geneLists/X_XDR.txt", header=FALSE, sep="\t")
XDR<-data.frame(XDR)
names(XDR)[1] <- "gene"
XDR_gene <- (XDR$gene)

XAR <- read.delim("~/01_code/geneLists/chrX_geneLists/X_XAR.txt", header=FALSE, sep="\t")
XAR<-data.frame(XAR)
names(XAR)[1] <- "gene"
XAR_gene <- (XAR$gene)

XCR <- read.delim("~/01_code/geneLists/chrX_geneLists/X_XCR.txt", header=FALSE, sep="\t")
XCR<-data.frame(XCR)
names(XCR)[1] <- "gene"
XCR_gene <- (XCR$gene)

XYHOM <- read.delim("~/01_code/geneLists/XY_homologousGenes.txt", header=FALSE, sep="\t")
XYHOM <-data.frame(XYHOM)
names(XYHOM)[1] <- "gene"
XYHOM_gene <- (XYHOM$gene)

X_5Mb <- read.delim("~/01_code/geneLists/X_5Mb.txt", header=FALSE, sep="\t")
X_5Mb <-data.frame(X_5Mb)
names(X_5Mb)[1] <- "gene"
X_5Mb_gene <- (X_5Mb$gene)

#-----------
# create a DGElist 
#-----------
dge <- DGEList(counts=counts, genes=genes)
dim(dge)

#-----------
# organize sample information 
#-----------
samplenames <- (pheno$sampleid) # sample names are not unique because a sample may belong to multiple groups
colnames(dge) <- samplenames

# create groups for the samples  
sex <- factor(pheno$sex, levels=c("female", "male"))
genome <- factor(pheno$genome, levels=c("default", "SS"))
aligner<- factor(pheno$aligner, levels=c("STAR", "HI"))

# add groups to samples in dge list 
dge$samples$sex <- sex
dge$samples$genome <- genome
dge$samples$aligner <- aligner
#dge$samples
dge$genes <- genes

#-----------
# data pre-processing
#-----------
cpm <- cpm(dge) #log2-transformed counts per gene per million mapped reads (cpm).
lcpm <- cpm(dge, log=FALSE) # log transformting the data 
# filter to only include counts of 1 or greater 
keep<- rowSums(cpm>1)>=0
table(keep)
dge <- dge[keep,, keep.lib.sizes=FALSE]
dim(dge)

#lcpm <- cpm(dge, log=TRUE)
png(filename ="figures/CLUSTER_SEX_GENOME_ALIGNER_dim23_brain.png",width=6,height=3,units="in",res=1200)
par(mfrow=c(1,3))
col.sex <- sex
levels(col.sex) <-  brewer.pal(nlevels(col.sex), "Set1")
col.sex <- as.character(col.sex)

col.genome <- genome
levels(col.genome) <-  brewer.pal(nlevels(col.genome), "Set2")
col.genome <- as.character(col.genome)

col.aligner <- aligner
levels(col.aligner) <-  brewer.pal(nlevels(col.aligner), "Set1")
col.aligner <- as.character(col.aligner)

plotMDS(lcpm, labels=sex, col=col.sex, dim.plot = c(2,3))
title(main="A.sexs")
plotMDS(lcpm, labels=genome, col=col.genome, dim.plot = c(2,3))
title(main="B.genome")
plotMDS(lcpm, labels=aligner, col=col.aligner, dim.plot = c(2,3))
title(main="C.aligner")
dev.off()
dev.off()

#lcpm <- cpm(dge, log=TRUE)
png(filename ="figures/CLUSTER_SEX_GENOME_ALIGNER_dim12_brain.png",width=6,height=3,units="in",res=1200)
par(mfrow=c(1,3))
col.sex <- sex
levels(col.sex) <-  brewer.pal(nlevels(col.sex), "Set1")
col.sex <- as.character(col.sex)

col.genome <- genome
levels(col.genome) <-  brewer.pal(nlevels(col.genome), "Set2")
col.genome <- as.character(col.genome)

col.aligner <- aligner
levels(col.aligner) <-  brewer.pal(nlevels(col.aligner), "Set1")
col.aligner <- as.character(col.aligner)

plotMDS(lcpm, labels=sex, col=col.sex)
title(main="A.sexs")
plotMDS(lcpm, labels=genome, col=col.genome)
title(main="B.genome")
plotMDS(lcpm, labels=aligner, col=col.aligner)
title(main="C.aligner")
dev.off()
dev.off()

#-------------------
# model matrix
#-------------------

c<-model.matrix(~sex+sex:genome+genome:aligner+aligner)
design<-model.matrix(~0+sex:genome:aligner)#:aligner+aligner)
colnames(design) <- make.names(colnames(design))

#-------------------
# Running limma voom
#-------------------
v <- voom(dge, design, normalize="quantile")

# Assigning transformed, quantile normalized expression values to a variable
transformed <- v$E

# Fitting the linear model
fit <- lmFit(v,design)

# contrast matrix all 
contr.matrix <- makeContrasts(
  MvsF_def_STAR = sexmale.genomedefault.alignerSTAR - sexfemale.genomedefault.alignerSTAR,
  MvsF_def_HI =sexmale.genomedefault.alignerHI - sexfemale.genomedefault.alignerHI,
  MvsF_SS_STAR = sexmale.genomeSS.alignerSTAR - sexfemale.genomeSS.alignerSTAR,
  MvsF_SS_HI = sexmale.genomeSS.alignerHI- sexfemale.genomeSS.alignerHI,
  FvsF_def_SS_STAR = sexfemale.genomedefault.alignerSTAR - sexfemale.genomeSS.alignerSTAR, 
  FvsF_def_SS_HISAT = sexfemale.genomedefault.alignerHI - sexfemale.genomeSS.alignerHI,
  MvsM_def_SS_STAR = sexmale.genomedefault.alignerSTAR - sexmale.genomeSS.alignerSTAR, 
  MvsM_def_SS_HISAT = sexmale.genomedefault.alignerHI - sexmale.genomeSS.alignerHI,
  levels = colnames(design))
contr.matrix

# Running contrast analysis
vfit <- contrasts.fit(fit, contrasts=contr.matrix)
summary(decideTests(vfit))

# Getting genes with logFC >2 and adjusted P-value <0.01 
vebayesfit <- eBayes(vfit)
MvsF_def_STAR_vebayesfit <- topTable(vebayesfit, coef=1, n=Inf, p.value=0.01, lfc=1) #MvsF_def_STAR
MvsF_def_HI_vebayesfit <- topTable(vebayesfit, coef=2, n=Inf, p.value=0.01, lfc=1) #MvsF_def_HI
MvsF_SS_STAR_vebayesfit <- topTable(vebayesfit, coef=3, n=Inf, p.value=0.01, lfc=1) #MvsF_SS_STAR
MvsF_SS_HI_vebayesfit <- topTable(vebayesfit, coef=4, n=Inf, p.value=0.01, lfc=1) #MvsF_SS_HI
FvsF_def_SS_STAR_vebayesfit <- topTable(vebayesfit, coef=5, n=Inf, p.value=0.01, lfc=1) #FvsF_def_SS_STAR
FvsF_def_SS_HISAT_vebayesfit <- topTable(vebayesfit, coef=6, n=Inf, p.value=0.01, lfc=1) #FvsF_def_SS_HISAT
MvsM_def_SS_STAR_vebayesfit <- topTable(vebayesfit, coef=7, n=Inf, p.value=0.01, lfc=1) #FvsF_def_SS_STAR
MvsM_def_SS_HISAT_vebayesfit <- topTable(vebayesfit, coef=8, n=Inf, p.value=0.01, lfc=1) #FvsF_def_SS_HISAT

# Number of selected significant genes
ngenes_MvsF_def_STAR_vebayesfit <- nrow(MvsF_def_STAR_vebayesfit)
ngenes_MvsF_def_STAR_vebayesfit
ngenes_MvsF_def_HI_vebayesfit <- nrow(MvsF_def_HI_vebayesfit)
ngenes_MvsF_def_HI_vebayesfit
ngenes_MvsF_SS_STAR_vebayesfit <- nrow(MvsF_SS_STAR_vebayesfit)
ngenes_MvsF_SS_STAR_vebayesfit
ngenes_MvsF_SS_HI_vebayesfit <- nrow(MvsF_SS_HI_vebayesfit)
ngenes_MvsF_SS_HI_vebayesfit
ngenes_FvsF_def_SS_STAR_vebayesfit <- nrow(FvsF_def_SS_STAR_vebayesfit)
ngenes_FvsF_def_SS_STAR_vebayesfit
ngenes_FvsF_def_SS_HISAT_vebayesfit <- nrow(FvsF_def_SS_HISAT_vebayesfit)
ngenes_FvsF_def_SS_HISAT_vebayesfit
ngenes_MvsM_def_SS_STAR_vebayesfit <- nrow(MvsM_def_SS_STAR_vebayesfit)
ngenes_MvsM_def_SS_STAR_vebayesfit
ngenes_MvsM_def_SS_HISAT_vebayesfit <- nrow(MvsM_def_SS_HISAT_vebayesfit)
ngenes_MvsM_def_SS_HISAT_vebayesfit

# Writing results into a comma-delimited file
#write.table(ngenes_MvsF_def_STAR_vebayesfit, file="genelists/genes_MvsF_def_STAR_vebayesfit_brain.txt",sep="\t")
#write.table(ngenes_MvsF_def_HI_vebayesfit, file="genelists/genes_MvsF_def_HI_vebayesfit_brain.txt",sep="\t")
#write.table(ngenes_MvsF_SS_STAR_vebayesfit, file="genelists/genes_MvsF_SS_STAR_vebayesfit_brain.txt",sep="\t")
#write.table(ngenes_MvsF_SS_HI_vebayesfit, file="genelists/genes_MvsF_SS_HI_vebayesfit_brain.txt",sep="\t")
#write.table(ngenes_FvsF_def_SS_STAR_vebayesfit, file="genelists/genes_FvsF_def_SS_STAR_vebayesfit_brain.txt",sep="\t")
#write.table(ngenes_FvsF_def_SS_HISAT_vebayesfit, file="genelists/genes_FvsF_def_SS_HISAT_vebayesfit_brain.txt",sep="\t")
#write.table(ngenes_MvsM_def_SS_STAR_vebayesfit, file="genelists/genes_MvsM_def_SS_STAR_vebayesfit_brain.txt",sep="\t")
#write.table(ngenes_MvsM_def_SS_HISAT_vebayesfit, file="genelists/genes_MvsM_def_SS_HISAT_vebayesfit_brain.txt",sep="\t")

# Decide tests
dt <- decideTests(vebayesfit)
summary(decideTests(vebayesfit))


vebayesfit <- eBayes(vfit)
MvsF_def_STAR_vebayesfit_lfc1 <- topTable(vebayesfit, coef=1, n=Inf, lfc=1) #MvsF_def_STAR
MvsF_def_HI_vebayesfit_lfc1 <- topTable(vebayesfit, coef=2, n=Inf, lfc=1) #MvsF_def_HI
MvsF_SS_STAR_vebayesfit_lfc1 <- topTable(vebayesfit, coef=3, n=Inf, lfc=1) #MvsF_SS_STAR
MvsF_SS_HI_vebayesfit_lfc1 <- topTable(vebayesfit, coef=4, n=Inf, lfc=1) #MvsF_SS_HI
FvsF_def_SS_STAR_vebayesfit_lfc1 <- topTable(vebayesfit, coef=5, n=Inf, lfc=1) #FvsF_def_SS_STAR
FvsF_def_SS_HISAT_vebayesfit_lfc1 <- topTable(vebayesfit, coef=6, n=Inf, lfc=1) #FvsF_def_SS_HISAT
MvsM_def_SS_STAR_vebayesfit_lfc1 <- topTable(vebayesfit, coef=7, n=Inf, lfc=1) #FvsF_def_SS_STAR
MvsM_def_SS_HISAT_vebayesfit_lfc1 <- topTable(vebayesfit, coef=8, n=Inf, lfc=1) #FvsF_def_SS_HISAT

# Number of selected significant genes
ngenes_MvsF_def_STAR_vebayesfit_lfc1 <- nrow(MvsF_def_STAR_vebayesfit_lfc1)
ngenes_MvsF_def_STAR_vebayesfit_lfc1
ngenes_MvsF_def_HI_vebayesfit_lfc1 <- nrow(MvsF_def_HI_vebayesfit_lfc1)
ngenes_MvsF_def_HI_vebayesfit_lfc1
ngenes_MvsF_SS_STAR_vebayesfit_lfc1 <- nrow(MvsF_SS_STAR_vebayesfit_lfc1)
ngenes_MvsF_SS_STAR_vebayesfit_lfc1
ngenes_MvsF_SS_HI_vebayesfit_lfc1 <- nrow(MvsF_SS_HI_vebayesfit_lfc1)
ngenes_MvsF_SS_HI_vebayesfit_lfc1
ngenes_FvsF_def_SS_STAR_vebayesfit_lfc1 <- nrow(FvsF_def_SS_STAR_vebayesfit_lfc1)
ngenes_FvsF_def_SS_STAR_vebayesfit_lfc1
ngenes_FvsF_def_SS_HISAT_vebayesfit_lfc1 <- nrow(FvsF_def_SS_HISAT_vebayesfit_lfc1)
ngenes_FvsF_def_SS_HISAT_vebayesfit_lfc1
ngenes_MvsM_def_SS_STAR_vebayesfit_lfc1 <- nrow(MvsM_def_SS_STAR_vebayesfit_lfc1)
ngenes_MvsM_def_SS_STAR_vebayesfit_lfc1
ngenes_MvsM_def_SS_HISAT_vebayesfit_lfc1 <- nrow(MvsM_def_SS_HISAT_vebayesfit_lfc1)
ngenes_MvsM_def_SS_HISAT_vebayesfit_lfc1

# Writing results into a comma-delimited file
#write.table(MvsF_def_STAR_vebayesfit, file="genelists/genes_MvsF_def_STAR_vebayesfit_lfc1_brain.txt",sep="\t")
#write.table(MvsF_def_HI_vebayesfit, file="genelists/genes_MvsF_def_HI_vebayesfit_lfc1_brain.txt",sep="\t")
#write.table(MvsF_SS_STAR_vebayesfit, file="genelists/genes_MvsF_SS_STAR_vebayesfit_lfc1_brain.txt",sep="\t")
#write.table(MvsF_SS_HI_vebayesfit, file="genelists/genes_MvsF_SS_HI_vebayesfit_lfc1_brain.txt",sep="\t")
#write.table(FvsF_def_SS_STAR_vebayesfit, file="genelists/genes_FvsF_def_SS_STAR_vebayesfit_lfc1_brain.txt",sep="\t")
#write.table(FvsF_def_SS_HISAT_vebayesfit, file="genelists/genes_FvsF_def_SS_HISAT_vebayesfit_lfc1_brain.txt",sep="\t")
#write.table(MvsM_def_SS_STAR_vebayesfit, file="genelists/genes_MvsM_def_SS_STAR_vebayesfit_lfc1_brain.txt",sep="\t")
#write.table(MvsM_def_SS_HISAT_vebayesfit, file="genelists/genes_MvsM_def_SS_HISAT_vebayesfit_lfc1_brain.txt",sep="\t")

# with only p.value = 0.01
vebayesfit_pva01 <- eBayes(vfit)
MvsF_def_STAR_vebayesfit_pval01 <- topTable(vebayesfit, coef=1, n=Inf, p.value=0.01) #MvsF_def_STAR
MvsF_def_HI_vebayesfit_pval01 <- topTable(vebayesfit, coef=2, n=Inf, p.value=0.01) #MvsF_def_HI
MvsF_SS_STAR_vebayesfit_pval01 <- topTable(vebayesfit, coef=3, n=Inf, p.value=0.01) #MvsF_SS_STAR
MvsF_SS_HI_vebayesfit_pval01 <- topTable(vebayesfit, coef=4, n=Inf, p.value=0.01) #MvsF_SS_HI
FvsF_def_SS_STAR_vebayesfit_pval01 <- topTable(vebayesfit, coef=5, n=Inf, p.value=0.01) #FvsF_def_SS_STAR
FvsF_def_SS_HISAT_vebayesfit_pval01 <- topTable(vebayesfit, coef=6, n=Inf, p.value=0.01) #FvsF_def_SS_HISAT
MvsM_def_SS_STAR_vebayesfit_pval01 <- topTable(vebayesfit, coef=7, n=Inf, p.value=0.01) #FvsF_def_SS_STAR
MvsM_def_SS_HISAT_vebayesfit_pval01 <- topTable(vebayesfit, coef=8, n=Inf, p.value=0.01) #FvsF_def_SS_HISAT

# Number of selected significant genes
ngenes_MvsF_def_STAR_vebayesfit_pval01 <- nrow(MvsF_def_STAR_vebayesfit_pval01)
ngenes_MvsF_def_STAR_vebayesfit_pval01
ngenes_MvsF_def_HI_vebayesfit_pval01 <- nrow(MvsF_def_HI_vebayesfit_pval01)
ngenes_MvsF_def_HI_vebayesfit_pval01
ngenes_MvsF_SS_STAR_vebayesfit_pval01 <- nrow(MvsF_SS_STAR_vebayesfit_pval01)
ngenes_MvsF_SS_STAR_vebayesfit_pval01
ngenes_MvsF_SS_HI_vebayesfit_pval01 <- nrow(MvsF_SS_HI_vebayesfit_pval01)
ngenes_MvsF_SS_HI_vebayesfit_pval01
ngenes_FvsF_def_SS_STAR_vebayesfit_pval01 <- nrow(FvsF_def_SS_STAR_vebayesfit_pval01)
ngenes_FvsF_def_SS_STAR_vebayesfit_pval01
ngenes_FvsF_def_SS_HISAT_vebayesfit_pval01 <- nrow(FvsF_def_SS_HISAT_vebayesfit_pval01)
ngenes_FvsF_def_SS_HISAT_vebayesfit_pval01
ngenes_MvsM_def_SS_STAR_vebayesfit_pval01 <- nrow(MvsM_def_SS_STAR_vebayesfit_pval01)
ngenes_MvsM_def_SS_STAR_vebayesfit_pval01
ngenes_MvsM_def_SS_HISAT_vebayesfit_pval01 <- nrow(MvsM_def_SS_HISAT_vebayesfit_pval01)
ngenes_MvsM_def_SS_HISAT_vebayesfit_pval01

# Writing results into a comma-delimited file
#write.table(MvsF_def_STAR_vebayesfit_pval01, file="genelists/genes_MvsF_def_STAR_vebayesfitpval01_brain.txt",sep="\t")
#write.table(MvsF_def_HI_vebayesfit_pval01, file="genelists/genes_MvsF_def_HI_vebayesfitpval01_brain.txt",sep="\t")
#write.table(MvsF_SS_STAR_vebayesfit_pval01, file="genelists/genes_MvsF_SS_STAR_vebayesfitpval01_brain.txt",sep="\t")
#write.table(MvsF_SS_HI_vebayesfit_pval01, file="genelists/genes_MvsF_SS_HI_vebayesfitpval01_brain.txt",sep="\t")
#write.table(FvsF_def_SS_STAR_vebayesfit_pval01, file="genelists/genes_FvsF_def_SS_STAR_vebayesfitpval01_brain.txt",sep="\t")
#write.table(FvsF_def_SS_HISAT_vebayesfit_pval01, file="genelists/genes_FvsF_def_SS_HISAT_vebayesfitpval01_brain.txt",sep="\t")
#write.table(MvsM_def_SS_STAR_vebayesfit_pval01, file="genelists/genes_MvsM_def_SS_STAR_vebayesfitpval01_brain.txt",sep="\t")
#write.table(MvsM_def_SS_HISAT_vebayesfit_pval01, file="genelists/genes_MvsM_def_SS_HISAT_vebayesfitpval01_brain.txt",sep="\t")

#------------------------------------------------------------------------
# M vs F def HI
#-------------------------------------------------------------------------

#--------------------------------
# Volcano plot with X and Y linked genes highlighted
#--------------------------------

# Assigning color to highlight genes that have an absolute fold change > 2 
# and a p-value < Bonferroni cut-off
MvsF_def_HI <- topTable(vebayesfit, n=Inf, coef=2)
write.table(MvsF_def_HI, "genelists/MvsF_def_HI.txt", sep="\t", quote = FALSE,row.names=FALSE)

P.Value <- c(MvsF_def_HI$P.Value)
logFC <- c(MvsF_def_HI$logFC)
adj.P.Val <- c(MvsF_def_HI$adj.P.Val)

df <- data.frame(P.Value, logFC, adj.P.Val, MvsF_def_HI$gene, MvsF_def_HI$Chr)

df.MB <- subset(df, adj.P.Val < 0.05 & MvsF_def_HI$Chr != "Y") #define male-biased, black
df.MB <- cbind(df.MB, rep(1, nrow(df.MB)))
colnames(df.MB)[6] <- "Color"

df.N <- subset(df, adj.P.Val >= 0.05) #define non-significant, grey
df.N<- cbind(df.N, rep(2, nrow(df.N)))
colnames(df.N)[6] <- "Color"

df.FB <- subset(df,  adj.P.Val < 0.05 & MvsF_def_HI$Chr != "X") #define female-biased, black
df.FB <- cbind(df.FB, rep(3, nrow(df.FB)))
colnames(df.FB)[6] <- "Color"

df.X <- subset(df, adj.P.Val < 0.05  & MvsF_def_HI$Chr == "X")#| MvsF_def_HI$Chr ="X;X" & adj.P.Val < 0.05) #define X-chromosomal, red
df.X <- cbind(df.X, rep(4, nrow(df.X)))
colnames(df.X)[6] <- "Color"

df.Y <- subset(df, adj.P.Val < 0.05  & MvsF_def_HI$Chr == "Y") #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(5, nrow(df.Y)))
colnames(df.Y)[6] <- "Color"

df.t <- rbind(df.MB, df.N, df.FB, df.X, df.Y)
df.t$Color <- as.factor(df.t$Color)
colnames(df.t)[4] <- "gene"
colnames(df.t)[5] <- "chr"

DEG_MvsF_def_HI_brain <- subset(df.t, adj.P.Val < 0.05)

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=Color )) +
  geom_point(alpha = 0.5, size = 1.75) +
  theme( legend.position = "none") +
  xlim(c(-10, 10)) + ylim(c(0, 17)) +
  scale_color_manual(values = c("black", "azure3", "black", "red", "blue")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
png(filename ="figures/MvsF_def_HI_VOLCANOplotFC_brain.png")
p
dev.off()
dev.off()
# Getting a subset of genes to label
subset_data_Y <- subset(df.t, adj.P.Val<0.05 & (chr=="Y") & (Color=="5"))
subset_data_X <- subset(df.t, adj.P.Val<0.05 & (chr == "X") & (Color=="4"))
subset_data <- rbind(subset_data_Y,subset_data_X)

png(filename ="figures/MvsF_def_HI_VOLCANOplotFC_brain.png",width=6,height=6,units="in",res=1200)
p+geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene, color=Color)
)+geom_hline(yintercept = 2, colour="#990000", linetype="dashed"
) + geom_vline(xintercept = 2, colour="#990000", linetype="dashed"
) + geom_vline(xintercept = -2, colour="#990000", linetype="dashed")
dev.off()
dev.off()

#--------------------------------
# log(CPM) male to female
#--------------------------------
# Getting the sex of each sample from the metada
males <- pheno[pheno$sex == 'male' & pheno$genome == 'default' & pheno$aligner == 'HI',]
males <- as.vector(males["sampleid"])
males <- as.vector(unlist(males))
females <- pheno[pheno$sex == 'female' & pheno$genome == 'default' & pheno$aligner == 'HI',]
females <- as.vector(females["sampleid"])
females <- as.vector(unlist(females))

# Subsetting male and female data
transformed_males <- as.data.frame(transformed)[males]
transformed_females <- as.data.frame(transformed)[females]

# Adding columns with mean expression values
transformed_males$mean = apply(transformed_males,1,mean,na.rm=TRUE)
transformed_females$mean = apply(transformed_females,1,mean,na.rm=TRUE)

# New dataframe with gene names and male & female mean expression values
mean_exp <- data.frame(matrix(nrow=nrow(transformed)))
mean_exp$gene <- v$genes$gene
mean_exp$chr <- v$genes$Chr
mean_exp$female <- transformed_females$mean
mean_exp$male <- transformed_males$mean
mean_exp_MvsF_def_HI <- mean_exp 

# Plotting with ggplot
# Assigning color to highlight X and Y linked genes

# Y chromsomes significant 
Y_sig_genes<-subset(df.Y, select=c(MvsF_def_HI.gene))
names(Y_sig_genes)[1] <- "gene"
Y_sig_genes <- (Y_sig_genes$gene) 
mean_exp.y<- subset(mean_exp, gene %in% Y_sig_genes)
mean_exp.y <- subset(mean_exp.y, female >= 0.00001 | male >= 0.00001)
mean_exp.y <- cbind(mean_exp.y, rep(1, nrow(mean_exp.y)))
colnames(mean_exp.y)[6] <- "Color"

# X chromsomes significant 
X_sig_genes<-subset(df.X, select=c(MvsF_def_HI.gene))
names(X_sig_genes)[1] <- "gene"
X_sig_genes <- (X_sig_genes$gene) 
mean_exp.x<- subset(mean_exp, gene %in% X_sig_genes)
mean_exp.x <- subset(mean_exp.x, female >= 0.00001 | male >= 0.00001)
mean_exp.x <- cbind(mean_exp.x, rep(2, nrow(mean_exp.x)))
colnames(mean_exp.x)[6] <- "Color"

# autosomes 
mean_exp.a <- subset(mean_exp, chr != "X" & chr != "Y" & chr !="X;X")
mean_exp.a<- cbind(mean_exp.a, rep(3, nrow(mean_exp.a)))
colnames(mean_exp.a)[6] <- "Color"

# female bias 
FB_sig_genes<-subset(df.FB, select=c(MvsF_def_HI.gene))
names(FB_sig_genes)[1] <- "gene"
FB_sig_genes <- (FB_sig_genes$gene) 
mean_exp.FB <- subset(mean_exp, gene %in% FB_sig_genes)
mean_exp.FB <- subset(mean_exp.FB, female >= 0.00001 | male >= 0.00001)
mean_exp.FB <- cbind(mean_exp.FB, rep(4, nrow(mean_exp.FB)))
colnames(mean_exp.FB)[6] <- "Color"

# male bias 
MB_sig_genes<-subset(df.MB, select=c(MvsF_def_HI.gene))
names(MB_sig_genes)[1] <- "gene"
MB_sig_genes <- (MB_sig_genes$gene) 
mean_exp.MB <- subset(mean_exp, gene %in% MB_sig_genes)
mean_exp.MB <- subset(mean_exp.MB, female >= 0.00001 | male >= 0.00001)
mean_exp.MB <- cbind(mean_exp.MB, rep(5, nrow(mean_exp.MB)))
colnames(mean_exp.MB)[6] <- "Color"

# Combine together 
mean_exp.c <- rbind(mean_exp.a, mean_exp.FB, mean_exp.MB, mean_exp.y, mean_exp.x)
mean_exp.c$Color <- as.factor(mean_exp.c$Color)
subset(mean_exp.c.df, gene=="AWAT1")
mean_exp.c_minCPM<-subset(mean_exp.c, female >= 0.00001 | male >= 0.00001 )
expressedIn_brain_def_HISAT <- mean_exp.c_minCPM
#write.table(mean_exp.c_minCPM, "genelists/expressedIn_brain_def_HISAT.txt", sep="\t", quote = FALSE,row.names=FALSE)

# if gene in
mean_exp.c.df<-data.frame(mean_exp.c_minCPM)
mean_exp.c.df_gene <- (mean_exp.c.df$gene)
DEG_CPM1_MvsF_def_HI_brain<- subset(DEG_MvsF_def_HI_brain, gene %in% mean_exp.c.df_gene)
DEG_CPM1_MvsF_def_HI_brain <- DEG_MvsF_def_HI_brain[-c(6)]
DEG_CPM1_MvsF_def_HI_brain<-DEG_CPM1_MvsF_def_HI_brain[duplicated(DEG_CPM1_MvsF_def_HI_brain), ]
DEG_CPM1_MvsF_def_HI_brain<- subset(DEG_MvsF_def_HI_brain, gene %in% mean_exp.c.df_gene)
write.table(DEG_CPM1_MvsF_def_HI_brain, "genelists/DEG_CPM1_MvsF_def_HISAT_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)
DEG_CPM1_MvsF_def_HI_maleBias_brain<-subset(DEG_CPM1_MvsF_def_HI_brain, logFC>=0)
DEG_CPM1_MvsF_def_HI_maleBias_brain<-DEG_CPM1_MvsF_def_HI_maleBias_brain[c(4)]
write.table(DEG_CPM1_MvsF_def_HI_maleBias_brain, "genelists/DEG_CPM1_MvsF_def_HI_maleBias_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)
DEG_CPM1_MvsF_def_HI_femaleBias_brain<-subset(DEG_CPM1_MvsF_def_HI_brain, logFC<0)
DEG_CPM1_MvsF_def_HI_femaleBias_brain<-DEG_CPM1_MvsF_def_HI_femaleBias_brain[c(4)]
write.table(DEG_CPM1_MvsF_def_HI_femaleBias_brain, "genelists/DEG_CPM1_MvsF_def_HI_femaleBias_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)

# $E component of the output object contains a "numeric matrix of normalized expression values 
# on the log2 scale": negative values indicate low levels of (normalized) expression, <1 read.
mean_exp_plot <- ggplot(data = mean_exp.c, aes(x = female, y = male, color=Color )) +
  geom_point(alpha = 0.5, size = 1.75) +
  theme( legend.position = "none") +
  xlim(c(-6, 17)) + ylim(c(-6, 17)) +
  scale_color_manual(values = c("blue", "red", "azure3", "black", "black")) +
  labs(x="Female log(CPM)", y="Male log(CPM)") +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
#png(filename ="figures/MvsF_def_HI_MEANEXPFC.png",width=6,height=6,units="in",res=1200)
mean_exp_plot
mean_exp_plot + geom_abline(intercept = 0)
dev.off()
dev.off()


# Getting a subset of genes to label
subset_data <- subset(mean_exp.c, abs(female - male) > 5 | abs(male - female) > 5)
subset_data <- mean_exp.c
subset_data<- subset(subset_data, (chr == "X" | chr=="Y") & (Color=="1" | Color=="2"))

# Plot with gene labels and regression
png(filename ="figures/MvsF_def_HI_MEANEXP_XYFC_brain.png",width=6,height=6,units="in",res=1200)
mean_exp_plot+coord_cartesian(ylim=c(-6,15))+coord_cartesian(xlim=c(-6,15))
mean_exp_plot+geom_text_repel(data=subset_data, aes(x = female, y = male, label=gene, color=Color))+
  geom_abline(intercept = 0, slope = 1, color="gray", linetype="dashed", size=1.5)+
  theme(legend.position="none")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_y_continuous(breaks=seq(-7,15,2))+ 
  scale_x_continuous(breaks=seq(-7,15,2)) 
dev.off()
dev.off()

exp_FemalesMales_def_HISAT<-cbind(transformed_females, transformed_males, mean_exp_MvsF_def_HI)
geneticSEXgenes_def_HISAT<- subset(exp_FemalesMales_def_HISAT, gene=="XIST" | gene=="ZFY" |gene=="ZFX" | gene=="SRY" | gene=="PCDH11X" |gene=="PCDH11Y" | gene=="USP9X" |gene=="USP9Y" | gene=="UTY" |gene=="DDX3X" |gene=="DDX3Y"| gene=="KDM6A")
geneticSEXgenes_def_HISAT<-subset(geneticSEXgenes_def_HISAT, select = -c(mean,mean.1,matrix.nrow...nrow.transformed..,female,male))
write.table(geneticSEXgenes_def_HISAT, "genelists/geneticSEXgenes_def_HI_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)
exp_FemalesMales_def_HISAT$CPM_female <-exp(exp_FemalesMales_def_HISAT$female)
exp_FemalesMales_def_HISAT$CPM_male <-exp(exp_FemalesMales_def_HISAT$male)
write.table(exp_FemalesMales_def_HISAT, "genelists/exp_FemalesMales_def_HISAT_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)

# Assigning color to STARghlight genes that have an absolute fold change > 2 
# and a p-value < Bonferroni cut-off
MvsF_def_STAR <- topTable(vebayesfit, n=Inf, coef=1)
write.table(MvsF_def_STAR, "genelists/MvsF_def_STAR.txt", sep="\t", quote = FALSE,row.names=FALSE)

P.Value <- c(MvsF_def_STAR$P.Value)
logFC <- c(MvsF_def_STAR$logFC)
adj.P.Val <- c(MvsF_def_STAR$adj.P.Val)

df <- data.frame(P.Value, logFC, adj.P.Val, MvsF_def_STAR$gene, MvsF_def_STAR$Chr)

df.MB <- subset(df, adj.P.Val < 0.05 & MvsF_def_STAR$Chr != "Y") #define male-biased, black
df.MB <- cbind(df.MB, rep(1, nrow(df.MB)))
colnames(df.MB)[6] <- "Color"

df.N <- subset(df, adj.P.Val >= 0.05) #define non-significant, grey
df.N<- cbind(df.N, rep(2, nrow(df.N)))
colnames(df.N)[6] <- "Color"

df.FB <- subset(df,  adj.P.Val < 0.05 & MvsF_def_STAR$Chr != "X") #define female-biased, black
df.FB <- cbind(df.FB, rep(3, nrow(df.FB)))
colnames(df.FB)[6] <- "Color"

df.X <- subset(df, adj.P.Val < 0.05  & MvsF_def_STAR$Chr == "X") #define X-chromosomal, red
df.X <- cbind(df.X, rep(4, nrow(df.X)))
colnames(df.X)[6] <- "Color"

df.Y <- subset(df, adj.P.Val < 0.05  & MvsF_def_STAR$Chr == "Y") #define Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(5, nrow(df.Y)))
colnames(df.Y)[6] <- "Color"

df.t <- rbind(df.MB, df.N, df.FB, df.X, df.Y)
df.t$Color <- as.factor(df.t$Color)
colnames(df.t)[4] <- "gene"
colnames(df.t)[5] <- "chr"
DEG_MvsF_def_STAR_brain <- subset(df.t, adj.P.Val < 0.05)

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=Color )) +
  geom_point(alpha = 0.5, size = 1.75) +
  theme( legend.position = "none") +
  xlim(c(-10, 10)) + ylim(c(0, 17)) +
  scale_color_manual(values = c("black", "azure3", "black", "red", "blue")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
png(filename ="figures/MvsF_def_STAR_VOLCANOplotFC_brain.png")
p
dev.off()
dev.off()
# Getting a subset of genes to label
subset_data_Y <- subset(df.t, adj.P.Val<0.05 & (chr=="Y") & (Color=="5"))
subset_data_X <- subset(df.t, adj.P.Val<0.05 & (chr == "X") & (Color=="4"))
subset_data <- rbind(subset_data_Y,subset_data_X)

# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
png(filename ="figures/MvsF_def_STAR_VOLCANOplotFC_brain.png",width=6,height=6,units="in",res=1200)
p+geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene, color=Color)
)+geom_hline(yintercept = 2, colour="#990000", linetype="dashed"
) + geom_vline(xintercept = 2, colour="#990000", linetype="dashed"
) + geom_vline(xintercept = -2, colour="#990000", linetype="dashed")
dev.off()
dev.off()


#--------------------------------
# log(CPM) male to female
#--------------------------------
# Getting the sex of each sample from the metada
males <- pheno[pheno$sex == 'male' & pheno$genome == 'default' & pheno$aligner == 'STAR',]
males <- as.vector(males["sampleid"])
males <- as.vector(unlist(males))
females <- pheno[pheno$sex == 'female' & pheno$genome == 'default' & pheno$aligner == 'STAR',]
females <- as.vector(females["sampleid"])
females <- as.vector(unlist(females))

# Subsetting male and female data
transformed_males <- as.data.frame(transformed)[males]
transformed_females <- as.data.frame(transformed)[females]

# Adding columns with mean expression values
transformed_males$mean = apply(transformed_males,1,mean,na.rm=TRUE)
transformed_females$mean = apply(transformed_females,1,mean,na.rm=TRUE)
#------------------------------------------------------------------------
# M vs F def STAR
#-------------------------------------------------------------------------

#--------------------------------
# Volcano plot with X and Y linked genes highlighted
#--------------------------------

# New dataframe with gene names and male & female mean expression values
mean_exp <- data.frame(matrix(nrow=nrow(transformed)))
mean_exp$gene <- v$genes$gene
mean_exp$chr <- v$genes$Chr
mean_exp$female <- transformed_females$mean
mean_exp$male <- transformed_males$mean
mean_exp_MvsF_def_STAR <- mean_exp
# Plotting with ggplot
# Assigning color to STARghlight X and Y linked genes

# Y chromsomes significant 
Y_sig_genes<-subset(df.Y, select=c(MvsF_def_STAR.gene))
names(Y_sig_genes)[1] <- "gene"
Y_sig_genes <- (Y_sig_genes$gene) 
mean_exp.y<- subset(mean_exp, gene %in% Y_sig_genes)
mean_exp.y <- subset(mean_exp.y, female >= 0.00001 | male >= 0.00001)
mean_exp.y <- cbind(mean_exp.y, rep(1, nrow(mean_exp.y)))
colnames(mean_exp.y)[6] <- "Color"

# X chromsomes significant 
X_sig_genes<-subset(df.X, select=c(MvsF_def_STAR.gene))
names(X_sig_genes)[1] <- "gene"
X_sig_genes <- (X_sig_genes$gene) 
mean_exp.x<- subset(mean_exp, gene %in% X_sig_genes)
mean_exp.x <- subset(mean_exp.x, female >= 0.00001 | male >= 0.00001)
mean_exp.x <- cbind(mean_exp.x, rep(2, nrow(mean_exp.x)))
colnames(mean_exp.x)[6] <- "Color"

# autosomes 
mean_exp.a <- subset(mean_exp, chr != "X" & chr != "Y" & chr !="X;X")
mean_exp.a<- cbind(mean_exp.a, rep(3, nrow(mean_exp.a)))
colnames(mean_exp.a)[6] <- "Color"

# female bias 
FB_sig_genes<-subset(df.FB, select=c(MvsF_def_STAR.gene))
names(FB_sig_genes)[1] <- "gene"
FB_sig_genes <- (FB_sig_genes$gene) 
mean_exp.FB <- subset(mean_exp, gene %in% FB_sig_genes)
mean_exp.FB <- subset(mean_exp.FB, female >= 0.00001 | male >= 0.00001)
mean_exp.FB <- cbind(mean_exp.FB, rep(4, nrow(mean_exp.FB)))
colnames(mean_exp.FB)[6] <- "Color"

# male bias 
MB_sig_genes<-subset(df.MB, select=c(MvsF_def_STAR.gene))
names(MB_sig_genes)[1] <- "gene"
MB_sig_genes <- (MB_sig_genes$gene) 
mean_exp.MB <- subset(mean_exp, gene %in% MB_sig_genes)
mean_exp.MB <- subset(mean_exp.MB, female >= 0.00001 | male >= 0.00001)
mean_exp.MB <- cbind(mean_exp.MB, rep(5, nrow(mean_exp.MB)))
colnames(mean_exp.MB)[6] <- "Color"

# Combine together 
mean_exp.c <- rbind(mean_exp.a, mean_exp.FB, mean_exp.MB, mean_exp.y, mean_exp.x)
mean_exp.c$Color <- as.factor(mean_exp.c$Color)
mean_exp.c_minCPM<-subset(mean_exp.c, female >= 0.00001 | male >=0.00001 )
expressedIn_brain_def_STAR <- mean_exp.c_minCPM
#write.table(mean_exp.c, "genelists/expressedIn_brain_def_STAR.txt", sep="\t", quote = FALSE,row.names=FALSE)

# if gene in
mean_exp.c.df<-data.frame(mean_exp.c_minCPM)
mean_exp.c.df_gene <- (mean_exp.c.df$gene)
DEG_CPM1_MvsF_def_STAR_brain <- subset(DEG_MvsF_def_STAR_brain, gene %in% mean_exp.c.df_gene)
DEG_CPM1_MvsF_def_STAR_brain <- DEG_MvsF_def_STAR_brain[-c(6)]
DEG_CPM1_MvsF_def_STAR_brain<-DEG_CPM1_MvsF_def_STAR_brain[duplicated(DEG_CPM1_MvsF_def_STAR_brain), ]
DEG_CPM1_MvsF_def_STAR_brain <- subset(DEG_MvsF_def_STAR_brain, gene %in% mean_exp.c.df_gene)
write.table(DEG_CPM1_MvsF_def_STAR_brain, "genelists/DEG_CPM1_MvsF_def_STAR_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)
DEG_CPM1_MvsF_def_STAR_maleBias_brain<-subset(DEG_CPM1_MvsF_def_STAR_brain, logFC>=0)
DEG_CPM1_MvsF_def_STAR_maleBias_brain<-DEG_CPM1_MvsF_def_STAR_maleBias_brain[c(4)]
write.table(DEG_CPM1_MvsF_def_STAR_maleBias_brain, "genelists/DEG_CPM1_MvsF_def_STAR_maleBias_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)
DEG_CPM1_MvsF_def_STAR_femaleBias_brain<-subset(DEG_CPM1_MvsF_def_STAR_brain, logFC<0)
DEG_CPM1_MvsF_def_STAR_femaleBias_brain<-DEG_CPM1_MvsF_def_STAR_femaleBias_brain[c(4)]
write.table(DEG_CPM1_MvsF_def_STAR_femaleBias_brain, "genelists/DEG_CPM1_MvsF_def_STAR_femaleBias_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)

# $E component of the output object contains a "numeric matrix of normalized expression values 
# on the log2 scale": negative values indicate low levels of (normalized) expression, <1 read.
mean_exp_plot <- ggplot(data = mean_exp.c, aes(x = female, y = male, color=Color )) +
  geom_point(alpha = 0.5, size = 1.75) +
  theme( legend.position = "none") +
  xlim(c(-6, 17)) + ylim(c(-6, 17)) +
  scale_color_manual(values = c("blue", "red", "azure3", "black", "black")) +
  labs(x="Female log(CPM)", y="Male log(CPM)") +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
#png(filename ="figures/MvsF_def_STAR_MEANEXPFC.png")
mean_exp_plot
dev.off()
dev.off()


# Getting a subset of genes to label
subset_data <- subset(mean_exp.c, abs(female - male) > 5 | abs(male - female) > 5)
subset_data <- mean_exp.c
subset_data<- subset(subset_data, (chr == "X" | chr=="Y") & (Color=="1" | Color=="2"))

# Plot with gene labels and regression
png(filename ="figures/MvsF_def_STAR_MEANEXP_XYFC_brain.png",width=6,height=6,units="in",res=1200)
mean_exp_plot+coord_cartesian(ylim=c(-6,15))+coord_cartesian(xlim=c(-6,15))
mean_exp_plot+geom_text_repel(data=subset_data, aes(x = female, y = male, label=gene, color=Color))+
  geom_abline(intercept = 0, slope = 1, color="gray", linetype="dashed", size=1.5)+
  theme(legend.position="none")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_y_continuous(breaks=seq(-7,15,2))+ 
  scale_x_continuous(breaks=seq(-7,15,2)) 
dev.off()
dev.off()

exp_FemalesMales_def_STAR<-cbind(transformed_females, transformed_males, mean_exp_MvsF_def_STAR)
geneticSEXgenes_def_STAR<- subset(exp_FemalesMales_def_STAR, gene=="XIST" | gene=="ZFY" |gene=="ZFX" | gene=="SRY" | gene=="PCDH11X" |gene=="PCDH11Y" | gene=="USP9X" |gene=="USP9Y" | gene=="UTY" |gene=="DDX3X" |gene=="DDX3Y"| gene=="KDM6A")
geneticSEXgenes_def_STAR<-subset(geneticSEXgenes_def_STAR, select = -c(mean,mean.1,matrix.nrow...nrow.transformed..,female,male))
write.table(geneticSEXgenes_def_STAR, "genelists/geneticSEXgenes_def_STAR_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)
exp_FemalesMales_def_STAR$CPM_female <-exp(exp_FemalesMales_def_STAR$female)
exp_FemalesMales_def_STAR$CPM_male <-exp(exp_FemalesMales_def_STAR$male)
write.table(exp_FemalesMales_def_STAR, "genelists/exp_FemalesMales_def_STAR_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)

#------------------------------------------------------------------------
# M vs F SS STAR
#-------------------------------------------------------------------------

#--------------------------------
# Volcano plot with X and Y linked genes highlighted
#--------------------------------

# Assigning color to STARghlight genes that have an absolute fold change > 2 
# and a p-value < Bonferroni cut-off
MvsF_SS_STAR <- topTable(vebayesfit, n=Inf, coef=3)
write.table(MvsF_SS_STAR, "genelists/MvsF_SS_STAR.txt", sep="\t", quote = FALSE,row.names=FALSE)

P.Value <- c(MvsF_SS_STAR$P.Value)
logFC <- c(MvsF_SS_STAR$logFC)
adj.P.Val <- c(MvsF_SS_STAR$adj.P.Val)

df <- data.frame(P.Value, logFC, adj.P.Val, MvsF_SS_STAR$gene, MvsF_SS_STAR$Chr)

df.MB <- subset(df, adj.P.Val < 0.05 & MvsF_SS_STAR$Chr != "Y") #SSine male-biased, black
df.MB <- cbind(df.MB, rep(1, nrow(df.MB)))
colnames(df.MB)[6] <- "Color"

df.N <- subset(df, adj.P.Val >= 0.05) #SSine non-significant, grey
df.N<- cbind(df.N, rep(2, nrow(df.N)))
colnames(df.N)[6] <- "Color"

df.FB <- subset(df,  adj.P.Val < 0.05 & MvsF_SS_STAR$Chr != "X") #SSine female-biased, black
df.FB <- cbind(df.FB, rep(3, nrow(df.FB)))
colnames(df.FB)[6] <- "Color"

df.X <- subset(df, adj.P.Val < 0.05  & MvsF_SS_STAR$Chr == "X") #SSine X-chromosomal, red
df.X <- cbind(df.X, rep(4, nrow(df.X)))
colnames(df.X)[6] <- "Color"

df.Y <- subset(df, adj.P.Val < 0.05  & MvsF_SS_STAR$Chr == "Y") #SSine Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(5, nrow(df.Y)))
colnames(df.Y)[6] <- "Color"

df.t <- rbind(df.MB, df.N, df.FB, df.X, df.Y)
df.t$Color <- as.factor(df.t$Color)
colnames(df.t)[4] <- "gene"
colnames(df.t)[5] <- "chr"
DEG_MvsF_SS_STAR_brain <- subset(df.t, adj.P.Val < 0.05)

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=Color )) +
  geom_point(alpha = 0.5, size = 1.75) +
  theme( legend.position = "none") +
  xlim(c(-10, 10)) + ylim(c(0, 17)) +
  scale_color_manual(values = c("black", "azure3", "black", "red", "blue")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
png(filename ="figures/MvsF_SS_STAR_VOLCANOplotFC_brain.png")
p
dev.off()
dev.off()
# Getting a subset of genes to label
subset_data_Y <- subset(df.t, adj.P.Val<0.05 & (chr=="Y") & (Color=="5"))
subset_data_X <- subset(df.t, adj.P.Val<0.05 & (chr == "X") & (Color=="4"))
subset_data <- rbind(subset_data_Y,subset_data_X)

# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
png(filename ="figures/MvsF_SS_STAR_VOLCANOplotFC_brain.png",width=6,height=6,units="in",res=1200)
p+geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene, color=Color)
)+geom_hline(yintercept = 2, colour="#990000", linetype="dashed"
) + geom_vline(xintercept = 2, colour="#990000", linetype="dashed"
) + geom_vline(xintercept = -2, colour="#990000", linetype="dashed")
dev.off()
dev.off()


#--------------------------------
# log(CPM) male to female
#--------------------------------
# Getting the sex of each sample from the metada
males <- pheno[pheno$sex == 'male' & pheno$genome == 'SS' & pheno$aligner == 'STAR',]
males <- as.vector(males["sampleid"])
males <- as.vector(unlist(males))
females <- pheno[pheno$sex == 'female' & pheno$genome == 'SS' & pheno$aligner == 'STAR',]
females <- as.vector(females["sampleid"])
females <- as.vector(unlist(females))

# Subsetting male and female data
transformed_males <- as.data.frame(transformed)[males]
transformed_females <- as.data.frame(transformed)[females]

# Adding columns with mean expression values
transformed_males$mean = apply(transformed_males,1,mean,na.rm=TRUE)
transformed_females$mean = apply(transformed_females,1,mean,na.rm=TRUE)

# New dataframe with gene names and male & female mean expression values
mean_exp <- data.frame(matrix(nrow=nrow(transformed)))
mean_exp$gene <- v$genes$gene
mean_exp$chr <- v$genes$Chr
mean_exp$female <- transformed_females$mean
mean_exp$male <- transformed_males$mean
mean_exp_MvsF_SS_STAR <- mean_exp

# Plotting with ggplot
# Assigning color to STARghlight X and Y linked genes

# Y chromsomes significant 
Y_sig_genes<-subset(df.Y, select=c(MvsF_SS_STAR.gene))
names(Y_sig_genes)[1] <- "gene"
Y_sig_genes <- (Y_sig_genes$gene) 
mean_exp.y<- subset(mean_exp, gene %in% Y_sig_genes)
mean_exp.y <- subset(mean_exp.y, female >= 0.00001 | male >= 0.00001)
mean_exp.y <- cbind(mean_exp.y, rep(1, nrow(mean_exp.y)))
colnames(mean_exp.y)[6] <- "Color"

# X chromsomes significant 
X_sig_genes<-subset(df.X, select=c(MvsF_SS_STAR.gene))
names(X_sig_genes)[1] <- "gene"
X_sig_genes <- (X_sig_genes$gene) 
mean_exp.x<- subset(mean_exp, gene %in% X_sig_genes)
mean_exp.x <- subset(mean_exp.x, female >= 0.00001 | male >= 0.00001)
mean_exp.x <- cbind(mean_exp.x, rep(2, nrow(mean_exp.x)))
colnames(mean_exp.x)[6] <- "Color"

# autosomes 
mean_exp.a <- subset(mean_exp, chr != "X" & chr != "Y" & chr !="X;X")
mean_exp.a<- cbind(mean_exp.a, rep(3, nrow(mean_exp.a)))
colnames(mean_exp.a)[6] <- "Color"

# female bias 
FB_sig_genes<-subset(df.FB, select=c(MvsF_SS_STAR.gene))
names(FB_sig_genes)[1] <- "gene"
FB_sig_genes <- (FB_sig_genes$gene) 
mean_exp.FB <- subset(mean_exp, gene %in% FB_sig_genes)
mean_exp.FB <- subset(mean_exp.FB, female >= 0.00001 | male >= 0.00001)
mean_exp.FB <- cbind(mean_exp.FB, rep(4, nrow(mean_exp.FB)))
colnames(mean_exp.FB)[6] <- "Color"

# male bias 
MB_sig_genes<-subset(df.MB, select=c(MvsF_SS_STAR.gene))
names(MB_sig_genes)[1] <- "gene"
MB_sig_genes <- (MB_sig_genes$gene) 
mean_exp.MB <- subset(mean_exp, gene %in% MB_sig_genes)
mean_exp.MB <- subset(mean_exp.MB, female >= 0.00001 | male >= 0.00001)
mean_exp.MB <- cbind(mean_exp.MB, rep(5, nrow(mean_exp.MB)))
colnames(mean_exp.MB)[6] <- "Color"


# Combine together 
mean_exp.c <- rbind(mean_exp.a, mean_exp.FB, mean_exp.MB, mean_exp.y, mean_exp.x)
mean_exp.c$Color <- as.factor(mean_exp.c$Color)
mean_exp.c_minCPM<-subset(mean_exp.c, female >= 0.00001 | male >=0.00001)
expressedIn_brain_SS_STAR <-mean_exp.c_minCPM
#write.table(mean_exp.c, "genelists/expressedIn_brain_SS_STAR.txt", sep="\t", quote = FALSE,row.names=FALSE)

# if gene in
mean_exp.c.df<-data.frame(mean_exp.c_minCPM)
mean_exp.c.df_gene <- (mean_exp.c.df$gene)
DEG_CPM1_MvsF_SS_STAR_brain <- subset(DEG_MvsF_SS_STAR_brain, gene %in% mean_exp.c.df_gene)
DEG_CPM1_MvsF_SS_STAR_brain <- DEG_MvsF_SS_STAR_brain[-c(6)]
DEG_CPM1_MvsF_SS_STAR_brain<-DEG_CPM1_MvsF_SS_STAR_brain[duplicated(DEG_CPM1_MvsF_SS_STAR_brain), ]
DEG_CPM1_MvsF_SS_STAR_brain <- subset(DEG_MvsF_SS_STAR_brain, gene %in% mean_exp.c.df_gene)
write.table(DEG_CPM1_MvsF_SS_STAR_brain, "genelists/DEG_CPM1_MvsF_SS_STAR_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)
DEG_CPM1_MvsF_SS_STAR_maleBias_brain<-subset(DEG_CPM1_MvsF_SS_STAR_brain, logFC>=0)
DEG_CPM1_MvsF_SS_STAR_maleBias_brain<-DEG_CPM1_MvsF_SS_STAR_maleBias_brain[c(4)]
write.table(DEG_CPM1_MvsF_SS_STAR_maleBias_brain, "genelists/DEG_CPM1_MvsF_SS_STAR_maleBias_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)
DEG_CPM1_MvsF_SS_STAR_femaleBias_brain<-subset(DEG_CPM1_MvsF_SS_STAR_brain, logFC<0)
DEG_CPM1_MvsF_SS_STAR_femaleBias_brain<-DEG_CPM1_MvsF_SS_STAR_femaleBias_brain[c(4)]
write.table(DEG_CPM1_MvsF_SS_STAR_femaleBias_brain, "genelists/DEG_CPM1_MvsF_SS_STAR_femaleBias_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)


# $E component of the output object contains a "numeric matrix of normalized expression values 
# on the log2 scale": negative values indicate low levels of (normalized) expression, <1 read.
mean_exp_plot <- ggplot(data = mean_exp.c, aes(x = female, y = male, color=Color )) +
  geom_point(alpha = 0.5, size = 1.75) +
  theme( legend.position = "none") +
  xlim(c(-7, 17)) + ylim(c(-7, 17)) +
  scale_color_manual(values = c("blue", "red", "azure3", "black", "black")) +
  labs(x="Female log(CPM)", y="Male log(CPM)") +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
#png(filename ="figures/MvsF_SS_STAR_MEANEXPFC.png")
mean_exp_plot
dev.off()
dev.off()


# Getting a subset of genes to label
subset_data <- subset(mean_exp.c, abs(female - male) > 5 | abs(male - female) > 5)
subset_data <- mean_exp.c
subset_data<- subset(subset_data, (chr == "X" | chr=="Y") & (Color=="1" | Color=="2"))

# Plot with gene labels and regression
png(filename ="figures/MvsF_SS_STAR_MEANEXP_XYFC_brain.png",width=6,height=6,units="in",res=1200)
mean_exp_plot+coord_cartesian(ylim=c(-6,15))+coord_cartesian(xlim=c(-6,15))
mean_exp_plot+geom_text_repel(data=subset_data, aes(x = female, y = male, label=gene, color=Color))+
  geom_abline(intercept = 0, slope = 1, color="gray", linetype="dashed", size=1.5)+
  theme(legend.position="none")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_y_continuous(breaks=seq(-7,15,2))+ 
  scale_x_continuous(breaks=seq(-7,15,2)) 
dev.off()
dev.off()

exp_FemalesMales_SS_STAR<-cbind(transformed_females, transformed_males, mean_exp_MvsF_SS_STAR)
geneticSEXgenes_SS_STAR <- subset(exp_FemalesMales_SS_STAR, gene=="XIST" | gene=="ZFY" |gene=="ZFX" | gene=="SRY" | gene=="PCDH11X" |gene=="PCDH11Y" | gene=="USP9X" |gene=="USP9Y" | gene=="UTY" |gene=="DDX3X" |gene=="DDX3Y"| gene=="KDM6A")
geneticSEXgenes_SS_STAR<-subset(geneticSEXgenes_SS_STAR, select = -c(mean,mean.1,matrix.nrow...nrow.transformed..,female,male))
write.table(geneticSEXgenes_SS_STAR, "genelists/geneticSEXgenes_SS_STAR_brain.txt", sep="\t")
exp_FemalesMales_SS_STAR$CPM_female <-exp(exp_FemalesMales_SS_STAR$female)
exp_FemalesMales_SS_STAR$CPM_male <-exp(exp_FemalesMales_SS_STAR$male)
write.table(exp_FemalesMales_SS_STAR, "genelists/exp_FemalesMales_SS_STAR_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)

#------------------------------------------------------------------------
# M vs F SS HI
#-------------------------------------------------------------------------

#--------------------------------
# Volcano plot with X and Y linked genes highlighted
#--------------------------------

# Assigning color to HIghlight genes that have an absolute fold change > 2 
# and a p-value < Bonferroni cut-off
MvsF_SS_HI <- topTable(vebayesfit, n=Inf, coef=3)
write.table(MvsF_SS_HI, "genelists/MvsF_SS_HI.txt", sep="\t", quote = FALSE,row.names=FALSE)

P.Value <- c(MvsF_SS_HI$P.Value)
logFC <- c(MvsF_SS_HI$logFC)
adj.P.Val <- c(MvsF_SS_HI$adj.P.Val)

df <- data.frame(P.Value, logFC, adj.P.Val, MvsF_SS_HI$gene, MvsF_SS_HI$Chr)

df.MB <- subset(df, adj.P.Val < 0.05 & MvsF_SS_HI$Chr != "Y") #SSine male-biased, black
df.MB <- cbind(df.MB, rep(1, nrow(df.MB)))
colnames(df.MB)[6] <- "Color"

df.N <- subset(df, adj.P.Val >= 0.05) #SSine non-significant, grey
df.N<- cbind(df.N, rep(2, nrow(df.N)))
colnames(df.N)[6] <- "Color"

df.FB <- subset(df,  adj.P.Val < 0.05 & MvsF_SS_HI$Chr != "X") #SSine female-biased, black
df.FB <- cbind(df.FB, rep(3, nrow(df.FB)))
colnames(df.FB)[6] <- "Color"

df.X <- subset(df, adj.P.Val < 0.05  & MvsF_SS_HI$Chr == "X") #SSine X-chromosomal, red
df.X <- cbind(df.X, rep(4, nrow(df.X)))
colnames(df.X)[6] <- "Color"

df.Y <- subset(df, adj.P.Val < 0.05  & MvsF_SS_HI$Chr == "Y") #SSine Y-chromosomal, blue
df.Y <- cbind(df.Y, rep(5, nrow(df.Y)))
colnames(df.Y)[6] <- "Color"

df.t <- rbind(df.MB, df.N, df.FB, df.X, df.Y)
df.t$Color <- as.factor(df.t$Color)
colnames(df.t)[4] <- "gene"
colnames(df.t)[5] <- "chr"
DEG_MvsF_SS_HI_brain <- subset(df.t, adj.P.Val < 0.05)

##Construct the plot object
p <- ggplot(data = df.t, aes(x = logFC, y = -log10(P.Value), color=Color )) +
  geom_point(alpha = 0.5, size = 1.75) +
  theme( legend.position = "none") +
  xlim(c(-10, 10)) + ylim(c(0, 17)) +
  scale_color_manual(values = c("black", "azure3", "black", "red", "blue")) +
  labs(x=expression(log[2](FC)), y=expression(-log[10] ~ "(" ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
png(filename ="figures/MvsF_SS_HI_VOLCANOplotFC_brain.png")
p
dev.off()
dev.off()
# Getting a subset of genes to label
subset_data_Y <- subset(df.t, adj.P.Val<0.05 & (chr=="Y") & (Color=="5"))
subset_data_X <- subset(df.t, adj.P.Val<0.05 & (chr == "X") & (Color=="4"))
subset_data <- rbind(subset_data_Y,subset_data_X)

# Volcanoplot with gene labels. Change log10(p-value) significance line to match 
# the p-value correspondin to FDR adjusted p-value of 0.01
png(filename ="figures/MvsF_SS_HI_VOLCANOplotFC_brain.png",width=6,height=6,units="in",res=1200)
p+geom_text_repel(data=subset_data,aes(x = logFC, y = -log10(P.Value), label=gene, color=Color)
)+geom_hline(yintercept = 2, colour="#990000", linetype="dashed"
) + geom_vline(xintercept = 2, colour="#990000", linetype="dashed"
) + geom_vline(xintercept = -2, colour="#990000", linetype="dashed")
dev.off()
dev.off()

#--------------------------------
# log(CPM) male to female
#--------------------------------
# Getting the sex of each sample from the metada
males <- pheno[pheno$sex == 'male' & pheno$genome == 'SS' & pheno$aligner == 'HI',]
males <- as.vector(males["sampleid"])
males <- as.vector(unlist(males))
females <- pheno[pheno$sex == 'female' & pheno$genome == 'SS' & pheno$aligner == 'HI',]
females <- as.vector(females["sampleid"])
females <- as.vector(unlist(females))

# Subsetting male and female data
transformed_males <- as.data.frame(transformed)[males]
transformed_females <- as.data.frame(transformed)[females]

# Adding columns with mean expression values
transformed_males$mean = apply(transformed_males,1,mean,na.rm=TRUE)
transformed_females$mean = apply(transformed_females,1,mean,na.rm=TRUE)

# New dataframe with gene names and male & female mean expression values
mean_exp <- data.frame(matrix(nrow=nrow(transformed)))
mean_exp$gene <- v$genes$gene
mean_exp$chr <- v$genes$Chr
mean_exp$female <- transformed_females$mean
mean_exp$male <- transformed_males$mean
mean_exp_MvsF_SS_HI <- mean_exp

# Plotting with ggplot
# Assigning color to HIghlight X and Y linked genes

# Y chromsomes significant 
Y_sig_genes<-subset(df.Y, select=c(MvsF_SS_HI.gene))
names(Y_sig_genes)[1] <- "gene"
Y_sig_genes <- (Y_sig_genes$gene) 
mean_exp.y<- subset(mean_exp, gene %in% Y_sig_genes)
mean_exp.y <- subset(mean_exp.y, female >= 0.00001 | male >= 0.00001)
mean_exp.y <- cbind(mean_exp.y, rep(1, nrow(mean_exp.y)))
colnames(mean_exp.y)[6] <- "Color"

# X chromsomes significant 
X_sig_genes<-subset(df.X, select=c(MvsF_SS_HI.gene))
names(X_sig_genes)[1] <- "gene"
X_sig_genes <- (X_sig_genes$gene) 
mean_exp.x<- subset(mean_exp, gene %in% X_sig_genes)
mean_exp.x <- subset(mean_exp.x, female >= 0.00001 | male >= 0.00001)
mean_exp.x <- cbind(mean_exp.x, rep(2, nrow(mean_exp.x)))
colnames(mean_exp.x)[6] <- "Color"

# autosomes 
mean_exp.a <- subset(mean_exp, chr != "X" & chr != "Y" & chr !="X;X")
mean_exp.a<- cbind(mean_exp.a, rep(3, nrow(mean_exp.a)))
colnames(mean_exp.a)[6] <- "Color"

# female bias 
FB_sig_genes<-subset(df.FB, select=c(MvsF_SS_HI.gene))
names(FB_sig_genes)[1] <- "gene"
FB_sig_genes <- (FB_sig_genes$gene) 
mean_exp.FB <- subset(mean_exp, gene %in% FB_sig_genes)
mean_exp.FB <- subset(mean_exp.FB, female >= 0.00001 | male >= 0.00001)
mean_exp.FB <- cbind(mean_exp.FB, rep(4, nrow(mean_exp.FB)))
colnames(mean_exp.FB)[6] <- "Color"

# male bias 
MB_sig_genes<-subset(df.MB, select=c(MvsF_SS_HI.gene))
names(MB_sig_genes)[1] <- "gene"
MB_sig_genes <- (MB_sig_genes$gene) 
mean_exp.MB <- subset(mean_exp, gene %in% MB_sig_genes)
mean_exp.MB <- subset(mean_exp.MB, female >= 0.00001 | male >= 0.00001)
mean_exp.MB <- cbind(mean_exp.MB, rep(5, nrow(mean_exp.MB)))
colnames(mean_exp.MB)[6] <- "Color"

# Combine together 
mean_exp.c <- rbind(mean_exp.a, mean_exp.FB, mean_exp.MB, mean_exp.y, mean_exp.x)
mean_exp.c$Color <- as.factor(mean_exp.c$Color)
mean_exp.c_minCPM<-subset(mean_exp.c, female >= 0.00001 | male >= 0.00001 )
expressedIn_brain_SS_HISAT <- mean_exp.c_minCPM
#write.table(mean_exp.c_minCPM, "genelists/expressedIn_brain_SS_HISAT.txt", sep="\t", quote = FALSE,row.names=FALSE)

# if gene in
mean_exp.c.df<-data.frame(mean_exp.c_minCPM)
mean_exp.c.df_gene <- (mean_exp.c.df$gene)
DEG_CPM1_MvsF_SS_HI_brain <- subset(DEG_MvsF_SS_HI_brain, gene %in% mean_exp.c.df_gene)
DEG_CPM1_MvsF_SS_HI_brain <- DEG_MvsF_SS_HI_brain[-c(6)]
DEG_CPM1_MvsF_SS_HI_brain<-DEG_CPM1_MvsF_SS_HI_brain[duplicated(DEG_CPM1_MvsF_SS_HI_brain), ]
DEG_CPM1_MvsF_SS_HI_brain <- subset(DEG_MvsF_SS_HI_brain, gene %in% mean_exp.c.df_gene)
write.table(DEG_CPM1_MvsF_SS_HI_brain, "genelists/DEG_CPM1_MvsF_SS_HISAT_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)
DEG_CPM1_MvsF_SS_HI_maleBias_brain<-subset(DEG_CPM1_MvsF_SS_HI_brain, logFC>=0)
DEG_CPM1_MvsF_SS_HI_maleBias_brain<-DEG_CPM1_MvsF_SS_HI_maleBias_brain[c(4)]
write.table(DEG_CPM1_MvsF_SS_HI_maleBias_brain, "genelists/DEG_CPM1_MvsF_SS_HI_maleBias_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)
DEG_CPM1_MvsF_SS_HI_femaleBias_brain<-subset(DEG_CPM1_MvsF_SS_HI_brain, logFC<0)
DEG_CPM1_MvsF_SS_HI_femaleBias_brain<-DEG_CPM1_MvsF_SS_HI_femaleBias_brain[c(4)]
write.table(DEG_CPM1_MvsF_SS_HI_femaleBias_brain, "genelists/DEG_CPM1_MvsF_SS_HI_femaleBias_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)

# $E component of the output object contains a "numeric matrix of normalized expression values 
# on the log2 scale": negative values indicate low levels of (normalized) expression, <1 read.
mean_exp_plot <- ggplot(data = mean_exp.c, aes(x = female, y = male, color=Color )) +
  geom_point(alpha = 0.5, size = 1.75) +
  theme( legend.position = "none") +
  xlim(c(-7, 17)) + ylim(c(-7, 17)) +
  scale_color_manual(values = c("blue", "red", "azure3", "black", "black")) +
  labs(x="Female log(CPM)", y="Male log(CPM)") +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=12)) +
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=12))
#png(filename ="figures/MvsF_SS_HI_MEANEXPFC.png")
mean_exp_plot
dev.off()
dev.off()


# Getting a subset of genes to label
subset_data <- subset(mean_exp.c, abs(female - male) > 5 | abs(male - female) > 5)
subset_data <- mean_exp.c
subset_data<- subset(subset_data, (chr == "X" | chr=="Y") & (Color=="1" | Color=="2"))
# Plot with gene labels and regression
png(filename ="figures/MvsF_SS_HI_MEANEXP_XYFC_brain.png",width=6,height=6,units="in",res=1200)
mean_exp_plot+geom_text_repel(data=subset_data, aes(x = female, y = male, label=gene, color=Color))+geom_smooth(method="lm",se=FALSE, color="black")+theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+coord_cartesian(ylim=c(-6,15))+coord_cartesian(xlim=c(-6,15))+scale_y_continuous(breaks=seq(-7,15,2))+ scale_x_continuous(breaks=seq(-7,15,2))+ theme_bw()
mean_exp_plot+coord_cartesian(ylim=c(-6,15))+coord_cartesian(xlim=c(-6,15))
mean_exp_plot+geom_text_repel(data=subset_data, aes(x = female, y = male, label=gene, color=Color))+
  geom_abline(intercept = 0, slope = 1, color="gray", linetype="dashed", size=1.5)+
  theme(legend.position="none")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_y_continuous(breaks=seq(-7,15,2))+ 
  scale_x_continuous(breaks=seq(-7,15,2)) 
dev.off()
dev.off()

exp_FemalesMales_SS_HISAT<-cbind(transformed_females, transformed_males, mean_exp_MvsF_SS_HI)
geneticSEXgenes_SS_HISAT<- subset(exp_FemalesMales_SS_HISAT, gene=="XIST" | gene=="ZFY" |gene=="ZFX" | gene=="SRY" | gene=="PCDH11X" |gene=="PCDH11Y" | gene=="USP9X" |gene=="USP9Y" | gene=="UTY" |gene=="DDX3X" |gene=="DDX3Y"| gene=="KDM6A")
geneticSEXgenes_SS_HISAT<-subset(geneticSEXgenes_SS_HISAT, select = -c(mean,mean.1,matrix.nrow...nrow.transformed..,female,male))
write.table(geneticSEXgenes_SS_HISAT, "genelists/geneticSEXgenes_SS_HISAT_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)
exp_FemalesMales_SS_HISAT$CPM_female <-exp(exp_FemalesMales_SS_HISAT$female)
exp_FemalesMales_SS_HISAT$CPM_male <-exp(exp_FemalesMales_SS_HISAT$male)
write.table(exp_FemalesMales_SS_HISAT, "genelists/exp_FemalesMales_SS_HISAT_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)


#-----------
# combine genetic sex gene lists together 
#-----------
geneticSEX<-cbind(geneticSEXgenes_SS_HISAT,geneticSEXgenes_def_HISAT,geneticSEXgenes_SS_STAR,geneticSEXgenes_def_STAR)
write.table(geneticSEX, "genelists/geneticSEXgenes_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)

# expressed in brain
ExpressedIn_brain_ALL <- rbind(expressedIn_brain_SS_HISAT,expressedIn_brain_SS_STAR,expressedIn_brain_def_HISAT,expressedIn_brain_def_STAR) 
ExpressedIn_brain_ALLgenes <- ExpressedIn_brain_ALL[c(2)]
ExpressedIn_brain_ALLuniquegenes<-unique(ExpressedIn_brain_ALLgenes[1])
write.table(ExpressedIn_brain_ALLuniquegenes, "genelists/ExpressedIn_brain_ALLgenes.txt", sep="\t", quote = FALSE,row.names=FALSE)
#---------------------------------#---------------------------------#---------------------------------#---------------------------------
#---------------------------------#---------------------------------#---------------------------------#---------------------------------
#---------------------------------#---------------------------------#---------------------------------#---------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
#
#
#
#         X chromosome expression differences between default and sex-specific mapping 
#
#
#
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
# Female SS and female def STAR
mean_exp_MvsF_def_STAR$female_def_STAR_cpm<-((exp(mean_exp_MvsF_def_STAR$female)))
mean_exp_MvsF_SS_STAR$female_SS_STAR_cpm<-((exp(mean_exp_MvsF_SS_STAR$female)))
F_def_STAR <- mean_exp_MvsF_def_STAR[ c("gene","chr","female_def_STAR_cpm") ]
F_SS_STAR <- mean_exp_MvsF_SS_STAR[ c("female_SS_STAR_cpm") ]
FvsF_def_SS_STAR_cpm<-cbind(F_def_STAR,F_SS_STAR)
FvsF_def_SS_STAR_cpm$SS_def_diff<-((FvsF_def_SS_STAR_cpm$female_SS_STAR_cpm)-(FvsF_def_SS_STAR_cpm$female_def_STAR_cpm))
FvsF_def_SS_STAR_cpm$SS_def_ratio<-((FvsF_def_SS_STAR_cpm$female_SS_STAR_cpm)/(FvsF_def_SS_STAR_cpm$female_def_STAR_cpm))
FvsF_def_SS_STAR_cpm$SS_def_logratio<-(log((FvsF_def_SS_STAR_cpm$female_SS_STAR_cpm)/(FvsF_def_SS_STAR_cpm$female_def_STAR_cpm)))
FvsF_def_SS_STAR_cpm$SS_def_log2FC<-(log2((FvsF_def_SS_STAR_cpm$female_SS_STAR_cpm)/(FvsF_def_SS_STAR_cpm$female_def_STAR_cpm)))
FvsF_STAR_log2FC<-subset(FvsF_def_SS_STAR_cpm, SS_def_log2FC >=1 | SS_def_log2FC <= -1)
write.table(FvsF_STAR_log2FC, "genelists/FvsF_STAR_log2FC_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)
write.table(FvsF_def_SS_STAR_cpm,"genelists/FvsF_def_SS_STAR_cpm_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)

# order X chromosome genes
#chrX_FvsF_def_SS_STAR_cpm <-subset(FvsF_def_SS_STAR_cpm, gene %in% X_gene)
chrX_FvsF_def_SS_STAR_cpm<-subset(FvsF_def_SS_STAR_cpm, grepl("^X",chr))
chrX_FvsF_def_SS_STAR_cpm_order <- chrX_FvsF_def_SS_STAR_cpm[order(chrX_FvsF_def_SS_STAR_cpm$gene),]
chrX_FvsF_def_SS_STAR_cpm_order$gene <- factor(chrX_FvsF_def_SS_STAR_cpm_order$gene, levels = chrX_genes$gene )
write.table(chrX_FvsF_def_SS_STAR_cpm, "genelists/chrX_FvsF_def_SS_STAR_cpm_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)

PAR_FvsF<-subset(chrX_FvsF_def_SS_STAR_cpm_order, gene %in% PAR_gene)
PAR_FvsF_order <- PAR_FvsF[order(PAR_FvsF$gene),]
PAR_FvsF_order$gene <- factor(PAR_FvsF_order$gene, levels = chrX_genes$gene)
write.table(PAR_FvsF_order, "genelists/PAR_XTR_genes_FvsF_STAR_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)

# X_5Mb
X_5Mb_FvsF<-subset(chrX_FvsF_def_SS_STAR_cpm_order, gene %in% X_5Mb_gene)
X_5Mb_FvsF_order <- X_5Mb_FvsF[order(X_5Mb_FvsF$gene),]
X_5Mb_FvsF_order$gene <- factor(X_5Mb_FvsF_order$gene, levels = chrX_genes$gene)
X_5Mb_FvsF_order$gene <- droplevels(X_5Mb_FvsF_order$gene)

# PAR1 
PAR1_FvsF<-subset(chrX_FvsF_def_SS_STAR_cpm_order, gene %in% PAR1_gene)
PAR1_FvsF_order <- PAR1_FvsF[order(PAR1_FvsF$gene),]
PAR1_FvsF_order$gene <- factor(PAR1_FvsF_order$gene, levels = chrX_genes$gene)

# FvsF_def_SS_STAR_cpm
XYHOM_FvsF<-subset(FvsF_def_SS_STAR_cpm, gene %in% XYHOM_gene)
XYHOM_FvsF_order <- XYHOM_FvsF[order(XYHOM_FvsF$gene),]
#XYHOM_FvsF_order$gene <- factor(XYHOM_FvsF_order$gene, levels = chrX_genes$gene)
write.table(XYHOM_FvsF_order, "genelists/XYHOM_genes_FvsF_STAR_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)

# male SS and male def STAR
mean_exp_MvsF_def_STAR$male_def_STAR_cpm<-((exp(mean_exp_MvsF_def_STAR$male)))
mean_exp_MvsF_SS_STAR$male_SS_STAR_cpm<-((exp(mean_exp_MvsF_SS_STAR$male)))
M_def_STAR <- mean_exp_MvsF_def_STAR[ c("gene","chr","male_def_STAR_cpm") ]
M_SS_STAR <- mean_exp_MvsF_SS_STAR[ c("male_SS_STAR_cpm") ]
MvsM_def_SS_STAR_cpm<-cbind(M_def_STAR,M_SS_STAR)

MvsM_def_SS_STAR_cpm$SS_def_diff<-((MvsM_def_SS_STAR_cpm$male_SS_STAR_cpm)-(MvsM_def_SS_STAR_cpm$male_def_STAR_cpm))
MvsM_def_SS_STAR_cpm$SS_def_ratio<-((MvsM_def_SS_STAR_cpm$male_SS_STAR_cpm)/(MvsM_def_SS_STAR_cpm$male_def_STAR_cpm))
MvsM_def_SS_STAR_cpm$SS_def_logratio<-(log((MvsM_def_SS_STAR_cpm$male_SS_STAR_cpm)/(MvsM_def_SS_STAR_cpm$male_def_STAR_cpm)))
MvsM_def_SS_STAR_cpm$SS_def_log2FC<-(log2((MvsM_def_SS_STAR_cpm$male_SS_STAR_cpm)/(MvsM_def_SS_STAR_cpm$male_def_STAR_cpm)))
MvsM_STAR_log2FC<-subset(MvsM_def_SS_STAR_cpm, SS_def_log2FC >=1 | SS_def_log2FC <= -1)
write.table(MvsM_STAR_log2FC, "genelists/MvsM_STAR_log2FC_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)
write.table(MvsM_def_SS_STAR_cpm,"genelists/MvsM_def_SS_STAR_cpm_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)

# order X chromosome genes
#chrX_MvsM_def_SS_STAR_cpm <-subset(MvsM_def_SS_STAR_cpm, gene %in% X_gene)
chrX_MvsM_def_SS_STAR_cpm<-subset(MvsM_def_SS_STAR_cpm, grepl("^X",chr))
chrX_MvsM_def_SS_STAR_cpm_order <- chrX_MvsM_def_SS_STAR_cpm[order(chrX_MvsM_def_SS_STAR_cpm$gene),]
chrX_MvsM_def_SS_STAR_cpm_order$gene <- factor(chrX_MvsM_def_SS_STAR_cpm_order$gene, levels = chrX_genes$gene )
write.table(chrX_MvsM_def_SS_STAR_cpm, "genelists/chrX_MvsM_def_SS_STAR_cpm_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)

PAR_MvsM<-subset(chrX_MvsM_def_SS_STAR_cpm_order, gene %in% PAR_gene)
PAR_MvsM_order <- PAR_MvsM[order(PAR_MvsM$gene),]
PAR_MvsM_order$gene <- factor(PAR_MvsM_order$gene, levels = chrX_genes$gene)
write.table(PAR_MvsM_order, "genelists/PAR_XTR_genes_MvsM_STAR_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)

# X_5Mb
X_5Mb_MvsM<-subset(chrX_MvsM_def_SS_STAR_cpm_order, gene %in% X_5Mb_gene)
X_5Mb_MvsM_order <- X_5Mb_MvsM[order(X_5Mb_MvsM$gene),]
X_5Mb_MvsM_order$gene <- factor(X_5Mb_MvsM_order$gene, levels = chrX_genes$gene)
X_5Mb_MvsM_order$gene <- droplevels(X_5Mb_MvsM_order$gene)

# PAR1 
PAR1_MvsM<-subset(chrX_MvsM_def_SS_STAR_cpm_order, gene %in% PAR1_gene)
PAR1_MvsM_order <- PAR1_MvsM[order(PAR1_MvsM$gene),]
PAR1_MvsM_order$gene <- factor(PAR1_MvsM_order$gene, levels = chrX_genes$gene)

# MvsM_def_SS_STAR_cpm
XYHOM_MvsM<-subset(MvsM_def_SS_STAR_cpm, gene %in% XYHOM_gene)
XYHOM_MvsM_order <- XYHOM_MvsM[order(XYHOM_MvsM$gene),]
#XYHOM_FvsF_order$gene <- factor(XYHOM_FvsF_order$gene, levels = chrX_genes$gene)
write.table(XYHOM_MvsM_order, "genelists/XYHOM_genes_MvsM_STAR_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)


# brain - X chromosome fold change of the CPMs 
png(filename ="figures/Xchr_foldChangeOfCPM_STAR_brain.png",width=7,height=2,units="in",res=1200)
par(mar=c(0,2,0,0))
plot(chrX_FvsF_def_SS_STAR_cpm_order$gene, chrX_FvsF_def_SS_STAR_cpm_order$SS_def_log2FC, ylim=c(-1,10), xaxt='n', ann=FALSE, yaxt='n', ann=FALSE)
axis(2,cex.axis=.5)
points(chrX_FvsF_def_SS_STAR_cpm_order$gene, chrX_FvsF_def_SS_STAR_cpm_order$SS_def_log2FC, cex=1, pch=20, col=c("white"))
points(chrX_FvsF_def_SS_STAR_cpm_order$gene, chrX_FvsF_def_SS_STAR_cpm_order$SS_def_log2FC, cex=.64, pch=1, col=c("darkorange"))
points(PAR_FvsF_order$gene,PAR_FvsF_order$SS_def_log2FC, cex=.65,pch=1, col=c("white"))
points(PAR_FvsF_order$gene,PAR_FvsF_order$SS_def_log2FC, cex=.65,pch=1, col=c("white"))
points(PAR_FvsF_order$gene,PAR_FvsF_order$SS_def_log2FC, cex=.83,pch=21, col=c("darkorange2"))
#points(chrX_MvsM_def_SS_STAR_cpm_order$gene, chrX_MvsM_def_SS_STAR_cpm_order$SS_def_log2FC, cex=1, pch=20, col=c("white"))
points(chrX_MvsM_def_SS_STAR_cpm_order$gene, chrX_MvsM_def_SS_STAR_cpm_order$SS_def_log2FC, cex=.32, pch=0, col=c("light blue"))
points(PAR_MvsM_order$gene,PAR_MvsM_order$SS_def_log2FC, cex=.32,pch=0, col=c("white"))
points(PAR_MvsM_order$gene,PAR_MvsM_order$SS_def_log2FC, cex=.6,pch=0, col=c("blue"))
dev.off()
dev.off()


png(filename ="figures/Xchr_5MB_foldChangeOfCPM_STAR_brain.png",width=3,height=2,units="in",res=1200)
par(mar=c(0,2,0,0))
plot(X_5Mb_FvsF_order$gene, X_5Mb_FvsF_order$SS_def_log2FC, ylim=c(-1,10), xaxt='n', ann=FALSE, yaxt='n', ann=FALSE)
axis(2,cex.axis=.4)
points(X_5Mb_FvsF_order$gene, X_5Mb_FvsF_order$SS_def_log2FC, cex=4, pch="-", col=c("white"))
points(X_5Mb_FvsF_order$gene, X_5Mb_FvsF_order$SS_def_log2FC, cex=.65, pch=1, col=c("darkorange"))
points(PAR1_FvsF_order$gene, PAR1_FvsF_order$SS_def_log2FC, cex=.65,pch=1, col=c("white"))
points(PAR1_FvsF_order$gene, PAR1_FvsF_order$SS_def_log2FC, cex=.65,pch=1, col=c("white"))
points(PAR1_FvsF_order$gene, PAR1_FvsF_order$SS_def_log2FC, cex=.83,pch=21, col=c("darkorange2"))
points(X_5Mb_MvsM_order$gene, X_5Mb_MvsM_order$SS_def_log2FC, cex=.32, pch=0, col=c("light blue"))
points(PAR1_MvsM_order$gene,PAR1_MvsM_order$SS_def_log2FC, cex=.32,pch=0, col=c("white"))
points(PAR1_MvsM_order$gene,PAR1_MvsM_order$SS_def_log2FC, cex=.6,pch=0, col=c("blue"))
dev.off()
dev.off()

# Female SS and female def HISAT
mean_exp_MvsF_def_HI$female_def_HI_cpm<-((exp(mean_exp_MvsF_def_HI$female)))
mean_exp_MvsF_SS_HI$female_SS_HI_cpm<-((exp(mean_exp_MvsF_SS_HI$female)))
F_def_HI <- mean_exp_MvsF_def_HI[ c("gene","chr","female_def_HI_cpm") ]
F_SS_HI <- mean_exp_MvsF_SS_HI[ c("female_SS_HI_cpm") ]
FvsF_def_SS_HISAT_cpm<-cbind(F_def_HI,F_SS_HI)
FvsF_def_SS_HISAT_cpm$SS_def_diff<-((FvsF_def_SS_HISAT_cpm$female_SS_HI_cpm)-(FvsF_def_SS_HISAT_cpm$female_def_HI_cpm))
FvsF_def_SS_HISAT_cpm$SS_def_ratio<-((FvsF_def_SS_HISAT_cpm$female_SS_HI_cpm)/(FvsF_def_SS_HISAT_cpm$female_def_HI_cpm))
FvsF_def_SS_HISAT_cpm$SS_def_logratio<-(log((FvsF_def_SS_HISAT_cpm$female_SS_HI_cpm)/(FvsF_def_SS_HISAT_cpm$female_def_HI_cpm)))
FvsF_def_SS_HISAT_cpm$SS_def_log2FC<-(log2((FvsF_def_SS_HISAT_cpm$female_SS_HI_cpm)/(FvsF_def_SS_HISAT_cpm$female_def_HI_cpm)))
FvsF_HISAT_log2FC<-subset(FvsF_def_SS_HISAT_cpm, SS_def_log2FC >=1 | SS_def_log2FC <= -1)
write.table(FvsF_HISAT_log2FC, "genelists/FvsF_HISAT_log2FC_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)
write.table(FvsF_def_SS_HISAT_cpm,"genelists/FvsF_def_SS_HISAT_cpm_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)


# order X chromosome genes
#chrX_FvsF_def_SS_HISAT_cpm <-subset(FvsF_def_SS_HISAT_cpm, gene %in% X_gene)
chrX_FvsF_def_SS_HISAT_cpm<-subset(FvsF_def_SS_HISAT_cpm, grepl("^X",chr))
chrX_FvsF_def_SS_HISAT_cpm_order <- chrX_FvsF_def_SS_HISAT_cpm[order(chrX_FvsF_def_SS_HISAT_cpm$gene),]
chrX_FvsF_def_SS_HISAT_cpm_order$gene <- factor(chrX_FvsF_def_SS_HISAT_cpm_order$gene, levels = chrX_genes$gene )
write.table(chrX_FvsF_def_SS_HISAT_cpm, "genelists/chrX_FvsF_def_SS_HISAT_cpm_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)

PAR_FvsF<-subset(chrX_FvsF_def_SS_HISAT_cpm_order, gene %in% PAR_gene)
PAR_FvsF_order <- PAR_FvsF[order(PAR_FvsF$gene),]
PAR_FvsF_order$gene <- factor(PAR_FvsF_order$gene, levels = chrX_genes$gene)
write.table(PAR_FvsF_order, "genelists/PAR_XTR_genes_FvsF_HISAT_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)

# FvsF_def_SS_HISAT_cpm
XYHOM_FvsF<-subset(FvsF_def_SS_HISAT_cpm, gene %in% XYHOM_gene)
XYHOM_FvsF_order <- XYHOM_FvsF[order(XYHOM_FvsF$gene),]
#XYHOM_FvsF_order$gene <- factor(XYHOM_FvsF_order$gene, levels = chrX_genes$gene)
write.table(XYHOM_FvsF_order, "genelists/XYHOM_genes_FvsF_HISAT_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)

# X_5Mb
X_5Mb_FvsF<-subset(chrX_FvsF_def_SS_HISAT_cpm_order, gene %in% X_5Mb_gene)
X_5Mb_FvsF_order <- X_5Mb_FvsF[order(X_5Mb_FvsF$gene),]
X_5Mb_FvsF_order$gene <- factor(X_5Mb_FvsF_order$gene, levels = chrX_genes$gene)
X_5Mb_FvsF_order$gene <- droplevels(X_5Mb_FvsF_order$gene)

# PAR1 
PAR1_FvsF<-subset(chrX_FvsF_def_SS_HISAT_cpm_order, gene %in% PAR1_gene)
PAR1_FvsF_order <- PAR1_FvsF[order(PAR1_FvsF$gene),]
PAR1_FvsF_order$gene <- factor(PAR1_FvsF_order$gene, levels = chrX_genes$gene)

# male SS and male def HISAT
mean_exp_MvsF_def_HI$male_def_HI_cpm<-((exp(mean_exp_MvsF_def_HI$male)))
mean_exp_MvsF_SS_HI$male_SS_HI_cpm<-((exp(mean_exp_MvsF_SS_HI$male)))
M_def_HI <- mean_exp_MvsF_def_HI[ c("gene","chr","male_def_HI_cpm") ]
M_SS_HI <- mean_exp_MvsF_SS_HI[ c("male_SS_HI_cpm") ]
MvsM_def_SS_HISAT_cpm<-cbind(M_def_HI,M_SS_HI)

MvsM_def_SS_HISAT_cpm$SS_def_diff<-((MvsM_def_SS_HISAT_cpm$male_SS_HI_cpm)-(MvsM_def_SS_HISAT_cpm$male_def_HI_cpm))
MvsM_def_SS_HISAT_cpm$SS_def_ratio<-((MvsM_def_SS_HISAT_cpm$male_SS_HI_cpm)/(MvsM_def_SS_HISAT_cpm$male_def_HI_cpm))
MvsM_def_SS_HISAT_cpm$SS_def_logratio<-(log((MvsM_def_SS_HISAT_cpm$male_SS_HI_cpm)/(MvsM_def_SS_HISAT_cpm$male_def_HI_cpm)))
MvsM_def_SS_HISAT_cpm$SS_def_log2FC<-(log2((MvsM_def_SS_HISAT_cpm$male_SS_HI_cpm)/(MvsM_def_SS_HISAT_cpm$male_def_HI_cpm)))
MvsM_HISAT_log2FC<-subset(MvsM_def_SS_HISAT_cpm, SS_def_log2FC >=1 | SS_def_log2FC <= -1)
write.table(MvsM_HISAT_log2FC, "genelists/MvsM_HISAT_log2FC_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)
write.table(MvsM_def_SS_HISAT_cpm,"genelists/MvsM_def_SS_HISAT_cpm_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)

# order X chromosome genes
#chrX_MvsM_def_SS_HISAT_cpm <-subset(MvsM_def_SS_HISAT_cpm, gene %in% X_gene)
chrX_MvsM_def_SS_HISAT_cpm<-subset(MvsM_def_SS_HISAT_cpm, grepl("^X",chr))
chrX_MvsM_def_SS_HISAT_cpm_order <- chrX_MvsM_def_SS_HISAT_cpm[order(chrX_MvsM_def_SS_HISAT_cpm$gene),]
chrX_MvsM_def_SS_HISAT_cpm_order$gene <- factor(chrX_MvsM_def_SS_HISAT_cpm_order$gene, levels = chrX_genes$gene )
write.table(chrX_MvsM_def_SS_HISAT_cpm, "genelists/chrX_MvsM_def_SS_HISAT_cpm_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)

PAR_MvsM<-subset(chrX_MvsM_def_SS_HISAT_cpm_order, gene %in% PAR_gene)
PAR_MvsM_order <- PAR_MvsM[order(PAR_MvsM$gene),]
PAR_MvsM_order$gene <- factor(PAR_MvsM_order$gene, levels = chrX_genes$gene)
write.table(PAR_MvsM_order, "genelists/PAR_XTR_genes_MvsM_HISAT_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)

# X_5Mb
X_5Mb_MvsM<-subset(chrX_MvsM_def_SS_HISAT_cpm_order, gene %in% X_5Mb_gene)
X_5Mb_MvsM_order <- X_5Mb_MvsM[order(X_5Mb_MvsM$gene),]
X_5Mb_MvsM_order$gene <- factor(X_5Mb_MvsM_order$gene, levels = chrX_genes$gene)
X_5Mb_MvsM_order$gene <- droplevels(X_5Mb_MvsM_order$gene)

# PAR1 
PAR1_MvsM<-subset(chrX_MvsM_def_SS_HISAT_cpm_order, gene %in% PAR1_gene)
PAR1_MvsM_order <- PAR1_MvsM[order(PAR1_MvsM$gene),]
PAR1_MvsM_order$gene <- factor(PAR1_MvsM_order$gene, levels = chrX_genes$gene)

# MvsM_def_SS_HISAT_cpm
XYHOM_MvsM <-subset(MvsM_def_SS_HISAT_cpm, gene %in% XYHOM_gene)
XYHOM_MvsM_order <- XYHOM_MvsM[order(XYHOM_MvsM$gene),]
#XYHOM_FvsF_order$gene <- factor(XYHOM_FvsF_order$gene, levels = chrX_genes$gene)
write.table(XYHOM_MvsM_order, "genelists/XYHOM_genes_MvsM_HISAT_brain.txt", sep="\t", quote = FALSE,row.names=FALSE)


# brain - X chromosome fold change of the CPMs 

png(filename ="figures/Xchr_foldChangeOfCPM_HISAT_brain.png",width=7,height=2,units="in",res=1200)
par(mar=c(0,2,0,0))
plot(chrX_FvsF_def_SS_HISAT_cpm_order$gene, chrX_FvsF_def_SS_HISAT_cpm_order$SS_def_log2FC, ylim=c(-1,4), xaxt='n', ann=FALSE, yaxt='n', ann=FALSE)
axis(2,cex.axis=.5)
points(chrX_FvsF_def_SS_HISAT_cpm_order$gene, chrX_FvsF_def_SS_HISAT_cpm_order$SS_def_log2FC, cex=1, pch=20, col=c("white"))
points(chrX_FvsF_def_SS_HISAT_cpm_order$gene, chrX_FvsF_def_SS_HISAT_cpm_order$SS_def_log2FC, cex=.64, pch=1, col=c("darkorange"))
points(PAR_FvsF_order$gene,PAR_FvsF_order$SS_def_log2FC, cex=.65,pch=1, col=c("white"))
points(PAR_FvsF_order$gene,PAR_FvsF_order$SS_def_log2FC, cex=.65,pch=1, col=c("white"))
points(PAR_FvsF_order$gene,PAR_FvsF_order$SS_def_log2FC, cex=.83,pch=21, col=c("darkorange2"))
#points(chrX_MvsM_def_SS_HISAT_cpm_order$gene, chrX_MvsM_def_SS_HISAT_cpm_order$SS_def_log2FC, cex=1, pch=20, col=c("white"))
points(chrX_MvsM_def_SS_HISAT_cpm_order$gene, chrX_MvsM_def_SS_HISAT_cpm_order$SS_def_log2FC, cex=.32, pch=0, col=c("light blue"))
points(PAR_MvsM_order$gene,PAR_MvsM_order$SS_def_log2FC, cex=.32,pch=0, col=c("white"))
points(PAR_MvsM_order$gene,PAR_MvsM_order$SS_def_log2FC, cex=.6,pch=0, col=c("blue"))
dev.off()
dev.off()


png(filename ="figures/Xchr_5MB_foldChangeOfCPM_HISAT_brain.png",width=3,height=2,units="in",res=1200)
par(mar=c(0,2,0,0))
plot(X_5Mb_FvsF_order$gene, X_5Mb_FvsF_order$SS_def_log2FC, ylim=c(-1,4), xaxt='n', ann=FALSE, yaxt='n', ann=FALSE)
axis(2,cex.axis=.4)
points(X_5Mb_FvsF_order$gene, X_5Mb_FvsF_order$SS_def_log2FC, cex=4, pch="-", col=c("white"))
points(X_5Mb_FvsF_order$gene, X_5Mb_FvsF_order$SS_def_log2FC, cex=.65, pch=1, col=c("darkorange"))
points(PAR1_FvsF_order$gene, PAR1_FvsF_order$SS_def_log2FC, cex=.65,pch=1, col=c("white"))
points(PAR1_FvsF_order$gene, PAR1_FvsF_order$SS_def_log2FC, cex=.65,pch=1, col=c("white"))
points(PAR1_FvsF_order$gene, PAR1_FvsF_order$SS_def_log2FC, cex=.83,pch=21, col=c("darkorange2"))
points(X_5Mb_MvsM_order$gene, X_5Mb_MvsM_order$SS_def_log2FC, cex=.32, pch=0, col=c("light blue"))
points(PAR1_MvsM_order$gene,PAR1_MvsM_order$SS_def_log2FC, cex=.32,pch=0, col=c("white"))
points(PAR1_MvsM_order$gene,PAR1_MvsM_order$SS_def_log2FC, cex=.6,pch=0, col=c("blue"))
dev.off()
dev.off()

#---------------------------
# log2FC Female SS to female def in STAR and HISAT
#--------------------------

#FvsF_STAR_log2FC<-as.vector(FvsF_STAR_log2FC$gene)
#FvsF_HISAT_log2FC<-as.vector(FvsF_HISAT_log2FC$gene)

#venn.plot <- venn.diagram(list(FvsF_HISAT_log2FC,FvsF_STAR_log2FC), NULL, fill=c("pink","pink"), cex = 2, cat.fontface=1, main="", category.names=c("FvsF_HISAT","FvsF_STAR"))
#grid.draw(venn.plot)
#dev.copy(jpeg,filename="figures/Venn_FvsF_log2FC_STARtoHISAT_brain.jpg",width=6,height=4,units="in",res=1200);
#dev.off()
#dev.off()

#MvsM_STAR_log2FC<-as.vector(MvsM_STAR_log2FC$gene)
#MvsM_HISAT_log2FC<-as.vector(MvsM_HISAT_log2FC$gene)

#venn.plot <- venn.diagram(list(MvsM_HISAT_log2FC,MvsM_STAR_log2FC), NULL, fill=c("light blue","light blue"), cex = 2, cat.fontface=1, main="", category.names=c("MvsM_HISAT","MvsM_STAR"))
#grid.draw(venn.plot)
#dev.copy(jpeg,filename="figures/Venn_MvsM_log2FC_STARtoHISAT_brain.jpg",width=6,height=4,units="in",res=1200);
#dev.off()
#dev.off()

#-------------
# venn def to SS for male bias and female bias genes 
#-------------
DEG_CPM1_MvsF_SS_STAR_maleBias_brain<-as.vector(DEG_CPM1_MvsF_SS_STAR_maleBias_brain$gene)
DEG_CPM1_MvsF_SS_HI_maleBias_brain<-as.vector(DEG_CPM1_MvsF_SS_HI_maleBias_brain$gene)

venn.plot <- venn.diagram(list(DEG_CPM1_MvsF_SS_HI_maleBias_brain,DEG_CPM1_MvsF_SS_STAR_maleBias_brain), NULL, fill=c("light blue","light blue"), main="", category.names=c("SS_HISAT","SS_STAR"),cat.cex = .75, cex = 2,cat.fontface=1,cat.just=list(c(1,.5) , c(0,0)))
grid.draw(venn.plot)
dev.copy(jpeg,filename="figures/Venn_MaleBias_SS_STARtoHISAT_brain.jpg",width=6,height=4,units="in",res=1200);
dev.off()
dev.off()

DEG_CPM1_MvsF_SS_STAR_femaleBias_brain<-as.vector(DEG_CPM1_MvsF_SS_STAR_femaleBias_brain$gene)
DEG_CPM1_MvsF_SS_HI_femaleBias_brain<-as.vector(DEG_CPM1_MvsF_SS_HI_femaleBias_brain$gene)

venn.plot <- venn.diagram(list(DEG_CPM1_MvsF_SS_HI_femaleBias_brain,DEG_CPM1_MvsF_SS_STAR_femaleBias_brain), NULL, fill=c("pink","pink"), cex = 2, cat.fontface=1, main="", category.names=c("SS_HISAT","SS_STAR"),cat.cex = .75,cat.just=list(c(1,.5) , c(0,0)))
grid.draw(venn.plot)
dev.copy(jpeg,filename="figures/Venn_FemaleBias_SS_STARtoHISAT_brain.jpg",width=6,height=4,units="in",res=1200);
dev.off()
dev.off()

DEG_CPM1_MvsF_def_STAR_maleBias_brain<-as.vector(DEG_CPM1_MvsF_def_STAR_maleBias_brain$gene)
DEG_CPM1_MvsF_def_HI_maleBias_brain<-as.vector(DEG_CPM1_MvsF_def_HI_maleBias_brain$gene)

venn.plot <- venn.diagram(list(DEG_CPM1_MvsF_def_HI_maleBias_brain,DEG_CPM1_MvsF_def_STAR_maleBias_brain), NULL, fill=c("light blue","light blue"), cex = 2, cat.fontface=1, main="", category.names=c("def_HISAT","def_STAR"),cat.cex = .75)#,cat.just=list(c(1,.5) , c(0,0)))
grid.draw(venn.plot)
dev.copy(jpeg,filename="figures/Venn_MaleBias_def_STARtoHISAT_brain.jpg",width=6,height=4,units="in",res=1200);
dev.off()
dev.off()

DEG_CPM1_MvsF_def_STAR_femaleBias_brain<-as.vector(DEG_CPM1_MvsF_def_STAR_femaleBias_brain$gene)
DEG_CPM1_MvsF_def_HI_femaleBias_brain<-as.vector(DEG_CPM1_MvsF_def_HI_femaleBias_brain$gene)

venn.plot <- venn.diagram(list(DEG_CPM1_MvsF_def_HI_femaleBias_brain,DEG_CPM1_MvsF_def_STAR_femaleBias_brain), NULL, fill=c("pink","pink"), cex = 2, cat.fontface=1, main="", category.names=c("def_HISAT","def_STAR"),cat.cex = .75)#,cat.just=list(c(1,.5) , c(0,0)))
grid.draw(venn.plot)
dev.copy(jpeg,filename="figures/Venn_FemaleBias_def_STARtoHISAT_brain.jpg",width=6,height=4,units="in",res=1200);
dev.off()
dev.off()

venn.plot <- venn.diagram(list(DEG_CPM1_MvsF_def_STAR_femaleBias_brain,DEG_CPM1_MvsF_SS_STAR_femaleBias_brain), NULL, fill=c("pink","pink"), cex = 2, cat.fontface=1, main="", category.names=c("def_STAR","SS_STAR"),cat.cex = .75,cat.just=list(c(1,.5) , c(0,0)))
grid.draw(venn.plot)
dev.copy(jpeg,filename="figures/Venn_FemaleBias_defToSS_STAR_brain.jpg",width=6,height=4,units="in",res=1200);
dev.off()
dev.off()

venn.plot <- venn.diagram(list(DEG_CPM1_MvsF_def_STAR_maleBias_brain,DEG_CPM1_MvsF_SS_STAR_maleBias_brain), NULL, fill=c("light blue","light blue"), cex = 2, cat.fontface=1, main="", category.names=c("def_STAR","SS_STAR"),cat.cex = .75,cat.just=list(c(1,.5) , c(0,0)))
grid.draw(venn.plot)
dev.copy(jpeg,filename="figures/Venn_MaleBias_defToSS_STAR_brain.jpg",width=6,height=4,units="in",res=1200);
dev.off()
dev.off()

venn.plot <- venn.diagram(list(DEG_CPM1_MvsF_def_HI_femaleBias_brain,DEG_CPM1_MvsF_SS_HI_femaleBias_brain), NULL, fill=c("pink","pink"), cex = 2, cat.fontface=1, main="", category.names=c("def_HISAT","SS_HISAT"),cat.cex = .75)#,cat.just=list(c(1,.5) , c(0,0)))
grid.draw(venn.plot)
dev.copy(jpeg,filename="figures/Venn_FemaleBias_defToSS_HISAT_brain.jpg",width=6,height=4,units="in",res=1200);
dev.off()
dev.off()

venn.plot <- venn.diagram(list(DEG_CPM1_MvsF_def_HI_maleBias_brain,DEG_CPM1_MvsF_SS_HI_maleBias_brain), NULL, fill=c("light blue","light blue"), cex = 2, cat.fontface=1, main="", category.names=c("def_HISAT","SS_HISAT"),cat.cex = .75)#,cat.just=list(c(1,.5) , c(0,0)))
grid.draw(venn.plot)
dev.copy(jpeg,filename="figures/Venn_MaleBias_defToSS_HISAT_brain.jpg",width=6,height=4,units="in",res=1200);
dev.off()
dev.off()

# venn gene lists
## HISAT
malebias_defSS_HI<-intersect(DEG_CPM1_MvsF_def_HI_maleBias_brain,DEG_CPM1_MvsF_SS_HI_maleBias_brain)
malebias_def_HI<-setdiff(DEG_CPM1_MvsF_def_HI_maleBias_brain,DEG_CPM1_MvsF_SS_HI_maleBias_brain)
malebias_SS_HI<-setdiff(DEG_CPM1_MvsF_SS_HI_maleBias_brain,DEG_CPM1_MvsF_def_HI_maleBias_brain)

df_malebias_defSS_HI <- data.frame(malebias_defSS_HI)
df_malebias_defSS_HI <- cbind(df_malebias_defSS_HI, rep("m_sharedbt_defAndSS", nrow(df_malebias_defSS_HI)))
colnames(df_malebias_defSS_HI)[2] <- "sex_biased"
colnames(df_malebias_defSS_HI)[1] <- "gene"

df_malebias_def_HI <- data.frame(malebias_def_HI)
df_malebias_def_HI <- cbind(df_malebias_def_HI, rep("m_def_only", nrow(df_malebias_def_HI)))
colnames(df_malebias_def_HI)[2] <- "sex_biased"
colnames(df_malebias_def_HI)[1] <- "gene"

df_malebias_SS_HI <- data.frame(malebias_SS_HI)
df_malebias_SS_HI <- cbind(df_malebias_SS_HI, rep("m_SS_only", nrow(df_malebias_SS_HI)))
colnames(df_malebias_SS_HI)[2] <- "sex_biased"
colnames(df_malebias_SS_HI)[1] <- "gene"

malebias_HISAT<-rbind(df_malebias_defSS_HI,df_malebias_def_HI,df_malebias_SS_HI)
malebias_HISAT_chr<-cbind(M_def_HI,malebias_HISAT[match(M_def_HI$gene, malebias_HISAT$gene),])
malebias_HISAT_chr<-malebias_HISAT_chr[c("gene","chr","sex_biased")]
malebias_HISAT_df<-subset(malebias_HISAT_chr, gene %in% malebias_HISAT$gene)

femalebias_defSS_HI<-intersect(DEG_CPM1_MvsF_def_HI_femaleBias_brain,DEG_CPM1_MvsF_SS_HI_femaleBias_brain)
femalebias_def_HI<-setdiff(DEG_CPM1_MvsF_def_HI_femaleBias_brain,DEG_CPM1_MvsF_SS_HI_femaleBias_brain)
femalebias_SS_HI<-setdiff(DEG_CPM1_MvsF_SS_HI_femaleBias_brain,DEG_CPM1_MvsF_def_HI_femaleBias_brain)

df_femalebias_defSS_HI <- data.frame(femalebias_defSS_HI)
df_femalebias_defSS_HI <- cbind(df_femalebias_defSS_HI, rep("f_sharedbt_defAndSS", nrow(df_femalebias_defSS_HI)))
colnames(df_femalebias_defSS_HI)[2] <- "sex_biased"
colnames(df_femalebias_defSS_HI)[1] <- "gene"

df_femalebias_def_HI <- data.frame(femalebias_def_HI)
df_femalebias_def_HI <- cbind(df_femalebias_def_HI, rep("f_def_only", nrow(df_femalebias_def_HI)))
colnames(df_femalebias_def_HI)[2] <- "sex_biased"
colnames(df_femalebias_def_HI)[1] <- "gene"

df_femalebias_SS_HI <- data.frame(femalebias_SS_HI)
df_femalebias_SS_HI <- cbind(df_femalebias_SS_HI, rep("f_SS_only", nrow(df_femalebias_SS_HI)))
colnames(df_femalebias_SS_HI)[2] <- "sex_biased"
colnames(df_femalebias_SS_HI)[1] <- "gene"

femalebias_HISAT<-rbind(df_femalebias_defSS_HI,df_femalebias_def_HI,df_femalebias_SS_HI)
femalebias_HISAT_chr<-cbind(M_def_HI,femalebias_HISAT[match(M_def_HI$gene, femalebias_HISAT$gene),])
femalebias_HISAT_chr<-femalebias_HISAT_chr[c("gene","chr","sex_biased")]
femalebias_HISAT_df<-subset(femalebias_HISAT_chr, gene %in% femalebias_HISAT$gene)

SexBiasedGeneExression_HISAT<-rbind(femalebias_HISAT_df,malebias_HISAT_df)
SexBiasedGeneExression_HISATorder<-SexBiasedGeneExression_HISAT[order(SexBiasedGeneExression_HISAT$sex_biased),]
write.table(SexBiasedGeneExression_HISATorder, "genelists/SexBiasedGeneExression_HISAT_brain.txt",sep="\t", quote = FALSE,row.names=FALSE)


## STAR
malebias_defSS_STAR<-intersect(DEG_CPM1_MvsF_def_STAR_maleBias_brain,DEG_CPM1_MvsF_SS_STAR_maleBias_brain)
malebias_def_STAR<-setdiff(DEG_CPM1_MvsF_def_STAR_maleBias_brain,DEG_CPM1_MvsF_SS_STAR_maleBias_brain)
malebias_SS_STAR<-setdiff(DEG_CPM1_MvsF_SS_STAR_maleBias_brain,DEG_CPM1_MvsF_def_STAR_maleBias_brain)

df_malebias_defSS_STAR <- data.frame(malebias_defSS_STAR)
df_malebias_defSS_STAR <- cbind(df_malebias_defSS_STAR, rep("m_sharedbt_defAndSS", nrow(df_malebias_defSS_STAR)))
colnames(df_malebias_defSS_STAR)[2] <- "sex_biased"
colnames(df_malebias_defSS_STAR)[1] <- "gene"

df_malebias_def_STAR <- data.frame(malebias_def_STAR)
df_malebias_def_STAR <- cbind(df_malebias_def_STAR, rep("m_def_only", nrow(df_malebias_def_STAR)))
colnames(df_malebias_def_STAR)[2] <- "sex_biased"
colnames(df_malebias_def_STAR)[1] <- "gene"

df_malebias_SS_STAR <- data.frame(malebias_SS_STAR)
df_malebias_SS_STAR <- cbind(df_malebias_SS_STAR, rep("m_SS_only", nrow(df_malebias_SS_STAR)))
colnames(df_malebias_SS_STAR)[2] <- "sex_biased"
colnames(df_malebias_SS_STAR)[1] <- "gene"

malebias_STAR<-rbind(df_malebias_defSS_STAR,df_malebias_def_STAR,df_malebias_SS_STAR)
malebias_STAR_chr<-cbind(M_def_STAR,malebias_STAR[match(M_def_STAR$gene, malebias_STAR$gene),])
malebias_STAR_chr<-malebias_STAR_chr[c("gene","chr","sex_biased")]
malebias_STAR_df<-subset(malebias_STAR_chr, gene %in% malebias_STAR$gene)

femalebias_defSS_STAR<-intersect(DEG_CPM1_MvsF_def_STAR_femaleBias_brain,DEG_CPM1_MvsF_SS_STAR_femaleBias_brain)
femalebias_def_STAR<-setdiff(DEG_CPM1_MvsF_def_STAR_femaleBias_brain,DEG_CPM1_MvsF_SS_STAR_femaleBias_brain)
femalebias_SS_STAR<-setdiff(DEG_CPM1_MvsF_SS_STAR_femaleBias_brain,DEG_CPM1_MvsF_def_STAR_femaleBias_brain)

df_femalebias_defSS_STAR <- data.frame(femalebias_defSS_STAR)
df_femalebias_defSS_STAR <- cbind(df_femalebias_defSS_STAR, rep("f_sharedbt_defAndSS", nrow(df_femalebias_defSS_STAR)))
colnames(df_femalebias_defSS_STAR)[2] <- "sex_biased"
colnames(df_femalebias_defSS_STAR)[1] <- "gene"

df_femalebias_def_STAR <- data.frame(femalebias_def_STAR)
df_femalebias_def_STAR <- cbind(df_femalebias_def_STAR, rep("f_def_only", nrow(df_femalebias_def_STAR)))
colnames(df_femalebias_def_STAR)[2] <- "sex_biased"
colnames(df_femalebias_def_STAR)[1] <- "gene"

df_femalebias_SS_STAR <- data.frame(femalebias_SS_STAR)
df_femalebias_SS_STAR <- cbind(df_femalebias_SS_STAR, rep("f_SS_only", nrow(df_femalebias_SS_STAR)))
colnames(df_femalebias_SS_STAR)[2] <- "sex_biased"
colnames(df_femalebias_SS_STAR)[1] <- "gene"

femalebias_STAR<-rbind(df_femalebias_defSS_STAR,df_femalebias_def_STAR,df_femalebias_SS_STAR)
femalebias_STAR_chr<-cbind(M_def_STAR,femalebias_STAR[match(M_def_STAR$gene, femalebias_STAR$gene),])
femalebias_STAR_chr<-femalebias_STAR_chr[c("gene","chr","sex_biased")]
femalebias_STAR_df<-subset(femalebias_STAR_chr, gene %in% femalebias_STAR$gene)

SexBiasedGeneExression_STAR<-rbind(femalebias_STAR_df,malebias_STAR_df)
SexBiasedGeneExression_STARorder<-SexBiasedGeneExression_STAR[order(SexBiasedGeneExression_STAR$sex_biased),]
write.table(SexBiasedGeneExression_STARorder, "genelists/SexBiasedGeneExression_STAR_brain.txt",sep="\t",quote = FALSE,row.names=FALSE)

#--------------
# FEMALE
# STAR X chromosome regions mean and median 
#--------------
# call in PAR1 genes 
PAR1_FvsF<-subset(chrX_FvsF_def_SS_STAR_cpm_order, gene %in% PAR1_gene)
PAR1_FvsF_order <- PAR1_FvsF[order(PAR1_FvsF$gene),]
PAR1_FvsF_order$gene <- factor(PAR1_FvsF_order$gene, levels = chrX_genes$gene)
mean(PAR1_FvsF_order$female_def_STAR_cpm)
median(PAR1_FvsF_order$female_def_STAR_cpm)
mean(PAR1_FvsF_order$female_SS_STAR_cpm)
median(PAR1_FvsF_order$female_SS_STAR_cpm)
((mean(PAR1_FvsF_order$female_SS_STAR_cpm))/(mean(PAR1_FvsF_order$female_def_STAR_cpm)))

# call in PAR2 genes 
PAR2_FvsF<-subset(chrX_FvsF_def_SS_STAR_cpm_order, gene %in% PAR2_gene)
PAR2_FvsF_order <- PAR2_FvsF[order(PAR2_FvsF$gene),]
PAR2_FvsF_order$gene <- factor(PAR2_FvsF_order$gene, levels = chrX_genes$gene)
mean(PAR2_FvsF_order$female_def_STAR_cpm)
median(PAR2_FvsF_order$female_def_STAR_cpm)
mean(PAR2_FvsF_order$female_SS_STAR_cpm)
median(PAR2_FvsF_order$female_SS_STAR_cpm)
((mean(PAR2_FvsF_order$female_SS_STAR_cpm))/(mean(PAR2_FvsF_order$female_def_STAR_cpm)))

# call in XTR genes 
XTR_FvsF<-subset(chrX_FvsF_def_SS_STAR_cpm_order, gene %in% XTR_gene)
XTR_FvsF_order <- XTR_FvsF[order(XTR_FvsF$gene),]
XTR_FvsF_order$gene <- factor(XTR_FvsF_order$gene, levels = chrX_genes$gene)
mean(XTR_FvsF_order$female_def_STAR_cpm)
median(XTR_FvsF_order$female_def_STAR_cpm)
mean(XTR_FvsF_order$female_SS_STAR_cpm)
median(XTR_FvsF_order$female_SS_STAR_cpm)
((mean(XTR_FvsF_order$female_SS_STAR_cpm))/(mean(XTR_FvsF_order$female_def_STAR_cpm)))


# call in XDR genes 
XDR_FvsF<-subset(chrX_FvsF_def_SS_STAR_cpm_order, gene %in% XDR_gene)
XDR_FvsF_order <- XDR_FvsF[order(XDR_FvsF$gene),]
XDR_FvsF_order$gene <- factor(XDR_FvsF_order$gene, levels = chrX_genes$gene)
mean(XDR_FvsF_order$female_def_STAR_cpm)
median(XDR_FvsF_order$female_def_STAR_cpm)
mean(XDR_FvsF_order$female_SS_STAR_cpm)
median(XDR_FvsF_order$female_SS_STAR_cpm)
((mean(XDR_FvsF_order$female_SS_STAR_cpm))/(mean(XDR_FvsF_order$female_def_STAR_cpm)))

# call in XAR genes 
XAR_FvsF<-subset(chrX_FvsF_def_SS_STAR_cpm_order, gene %in% XAR_gene)
XAR_FvsF_order <- XAR_FvsF[order(XAR_FvsF$gene),]
XAR_FvsF_order$gene <- factor(XAR_FvsF_order$gene, levels = chrX_genes$gene)
mean(XAR_FvsF_order$female_def_STAR_cpm)
median(XAR_FvsF_order$female_def_STAR_cpm)
mean(XAR_FvsF_order$female_SS_STAR_cpm)
median(XAR_FvsF_order$female_SS_STAR_cpm)
((mean(XAR_FvsF_order$female_SS_STAR_cpm))/(mean(XAR_FvsF_order$female_def_STAR_cpm)))


# call in XCR genes 
XCR_FvsF<-subset(chrX_FvsF_def_SS_STAR_cpm_order, gene %in% XCR_gene)
XCR_FvsF_order <- XCR_FvsF[order(XCR_FvsF$gene),]
XCR_FvsF_order$gene <- factor(XCR_FvsF_order$gene, levels = chrX_genes$gene)
mean(XCR_FvsF_order$female_def_STAR_cpm)
median(XCR_FvsF_order$female_def_STAR_cpm)
mean(XCR_FvsF_order$female_SS_STAR_cpm)
median(XCR_FvsF_order$female_SS_STAR_cpm)
((mean(XCR_FvsF_order$female_SS_STAR_cpm))/(mean(XCR_FvsF_order$female_def_STAR_cpm)))

#--------------
# FEMALE
# HISAT X chromosome regions mean and median 
#--------------
# call in PAR1 genes 
PAR1_FvsF<-subset(chrX_FvsF_def_SS_HISAT_cpm_order, gene %in% PAR1_gene)
PAR1_FvsF_order <- PAR1_FvsF[order(PAR1_FvsF$gene),]
PAR1_FvsF_order$gene <- factor(PAR1_FvsF_order$gene, levels = chrX_genes$gene)
mean(PAR1_FvsF_order$female_def_HI_cpm)
median(PAR1_FvsF_order$female_def_HI_cpm)
mean(PAR1_FvsF_order$female_SS_HI_cpm)
median(PAR1_FvsF_order$female_SS_HI_cpm)
((mean(PAR1_FvsF_order$female_SS_HI_cpm))/(mean(PAR1_FvsF_order$female_def_HI_cpm)))

# call in PAR2 genes 
PAR2_FvsF<-subset(chrX_FvsF_def_SS_HISAT_cpm_order, gene %in% PAR2_gene)
PAR2_FvsF_order <- PAR2_FvsF[order(PAR2_FvsF$gene),]
PAR2_FvsF_order$gene <- factor(PAR2_FvsF_order$gene, levels = chrX_genes$gene)
mean(PAR2_FvsF_order$female_def_HI_cpm)
median(PAR2_FvsF_order$female_def_HI_cpm)
mean(PAR2_FvsF_order$female_SS_HI_cpm)
median(PAR2_FvsF_order$female_SS_HI_cpm)
((mean(PAR2_FvsF_order$female_SS_HI_cpm))/(mean(PAR2_FvsF_order$female_def_HI_cpm)))

# call in XTR genes 
XTR_FvsF<-subset(chrX_FvsF_def_SS_HISAT_cpm_order, gene %in% XTR_gene)
XTR_FvsF_order <- XTR_FvsF[order(XTR_FvsF$gene),]
XTR_FvsF_order$gene <- factor(XTR_FvsF_order$gene, levels = chrX_genes$gene)
mean(XTR_FvsF_order$female_def_HI_cpm)
median(XTR_FvsF_order$female_def_HI_cpm)
mean(XTR_FvsF_order$female_SS_HI_cpm)
median(XTR_FvsF_order$female_SS_HI_cpm)
((mean(XTR_FvsF_order$female_SS_HI_cpm))/(mean(XTR_FvsF_order$female_def_HI_cpm)))

# call in XDR genes 
XDR_FvsF<-subset(chrX_FvsF_def_SS_HISAT_cpm_order, gene %in% XDR_gene)
XDR_FvsF_order <- XDR_FvsF[order(XDR_FvsF$gene),]
XDR_FvsF_order$gene <- factor(XDR_FvsF_order$gene, levels = chrX_genes$gene)
mean(XDR_FvsF_order$female_def_HI_cpm)
median(XDR_FvsF_order$female_def_HI_cpm)
mean(XDR_FvsF_order$female_SS_HI_cpm)
median(XDR_FvsF_order$female_SS_HI_cpm)
((mean(XDR_FvsF_order$female_SS_HI_cpm))/(mean(XDR_FvsF_order$female_def_HI_cpm)))

# call in XAR genes 
XAR_FvsF<-subset(chrX_FvsF_def_SS_HISAT_cpm_order, gene %in% XAR_gene)
XAR_FvsF_order <- XAR_FvsF[order(XAR_FvsF$gene),]
XAR_FvsF_order$gene <- factor(XAR_FvsF_order$gene, levels = chrX_genes$gene)
mean(XAR_FvsF_order$female_def_HI_cpm)
median(XAR_FvsF_order$female_def_HI_cpm)
mean(XAR_FvsF_order$female_SS_HI_cpm)
median(XAR_FvsF_order$female_SS_HI_cpm)
((mean(XAR_FvsF_order$female_SS_HI_cpm))/(mean(XAR_FvsF_order$female_def_HI_cpm)))

# call in XCR genes 
XCR_FvsF<-subset(chrX_FvsF_def_SS_HISAT_cpm_order, gene %in% XCR_gene)
XCR_FvsF_order <- XCR_FvsF[order(XCR_FvsF$gene),]
XCR_FvsF_order$gene <- factor(XCR_FvsF_order$gene, levels = chrX_genes$gene)
mean(XCR_FvsF_order$female_def_HI_cpm)
median(XCR_FvsF_order$female_def_HI_cpm)
mean(XCR_FvsF_order$female_SS_HI_cpm)
median(XCR_FvsF_order$female_SS_HI_cpm)
((mean(XCR_FvsF_order$female_SS_HI_cpm))/(mean(XCR_FvsF_order$female_def_HI_cpm)))

#--------------
# MALE
# STAR X chromosome regions mean and median 
#--------------
# call in PAR1 genes 
PAR1_MvsM<-subset(chrX_MvsM_def_SS_STAR_cpm_order, gene %in% PAR1_gene)
PAR1_MvsM_order <- PAR1_MvsM[order(PAR1_MvsM$gene),]
PAR1_MvsM_order$gene <- factor(PAR1_MvsM_order$gene, levels = chrX_genes$gene)
mean(PAR1_MvsM_order$male_def_STAR_cpm)
median(PAR1_MvsM_order$male_def_STAR_cpm)
mean(PAR1_MvsM_order$male_SS_STAR_cpm)
median(PAR1_MvsM_order$male_SS_STAR_cpm)
((mean(PAR1_MvsM_order$male_SS_STAR_cpm))/(mean(PAR1_MvsM_order$male_def_STAR_cpm)))

# call in PAR2 genes 
PAR2_MvsM<-subset(chrX_MvsM_def_SS_STAR_cpm_order, gene %in% PAR2_gene)
PAR2_MvsM_order <- PAR2_MvsM[order(PAR2_MvsM$gene),]
PAR2_MvsM_order$gene <- factor(PAR2_MvsM_order$gene, levels = chrX_genes$gene)
mean(PAR2_MvsM_order$male_def_STAR_cpm)
median(PAR2_MvsM_order$male_def_STAR_cpm)
mean(PAR2_MvsM_order$male_SS_STAR_cpm)
median(PAR2_MvsM_order$male_SS_STAR_cpm)
((mean(PAR2_MvsM_order$male_SS_STAR_cpm))/(mean(PAR2_MvsM_order$male_def_STAR_cpm)))

# call in XTR genes 
XTR_MvsM<-subset(chrX_MvsM_def_SS_STAR_cpm_order, gene %in% XTR_gene)
XTR_MvsM_order <- XTR_MvsM[order(XTR_MvsM$gene),]
XTR_MvsM_order$gene <- factor(XTR_MvsM_order$gene, levels = chrX_genes$gene)
mean(XTR_MvsM_order$male_def_STAR_cpm)
median(XTR_MvsM_order$male_def_STAR_cpm)
mean(XTR_MvsM_order$male_SS_STAR_cpm)
median(XTR_MvsM_order$male_SS_STAR_cpm)
((mean(XTR_MvsM_order$male_SS_STAR_cpm))/(mean(XTR_MvsM_order$male_def_STAR_cpm)))


# call in XDR genes 
XDR_MvsM<-subset(chrX_MvsM_def_SS_STAR_cpm_order, gene %in% XDR_gene)
XDR_MvsM_order <- XDR_MvsM[order(XDR_MvsM$gene),]
XDR_MvsM_order$gene <- factor(XDR_MvsM_order$gene, levels = chrX_genes$gene)
mean(XDR_MvsM_order$male_def_STAR_cpm)
median(XDR_MvsM_order$male_def_STAR_cpm)
mean(XDR_MvsM_order$male_SS_STAR_cpm)
median(XDR_MvsM_order$male_SS_STAR_cpm)
((mean(XDR_MvsM_order$male_SS_STAR_cpm))/(mean(XDR_MvsM_order$male_def_STAR_cpm)))

# call in XAR genes 
XAR_MvsM<-subset(chrX_MvsM_def_SS_STAR_cpm_order, gene %in% XAR_gene)
XAR_MvsM_order <- XAR_MvsM[order(XAR_MvsM$gene),]
XAR_MvsM_order$gene <- factor(XAR_MvsM_order$gene, levels = chrX_genes$gene)
mean(XAR_MvsM_order$male_def_STAR_cpm)
median(XAR_MvsM_order$male_def_STAR_cpm)
mean(XAR_MvsM_order$male_SS_STAR_cpm)
median(XAR_MvsM_order$male_SS_STAR_cpm)
((mean(XAR_MvsM_order$male_SS_STAR_cpm))/(mean(XAR_MvsM_order$male_def_STAR_cpm)))


# call in XCR genes 
XCR_MvsM<-subset(chrX_MvsM_def_SS_STAR_cpm_order, gene %in% XCR_gene)
XCR_MvsM_order <- XCR_MvsM[order(XCR_MvsM$gene),]
XCR_MvsM_order$gene <- factor(XCR_MvsM_order$gene, levels = chrX_genes$gene)
mean(XCR_MvsM_order$male_def_STAR_cpm)
median(XCR_MvsM_order$male_def_STAR_cpm)
mean(XCR_MvsM_order$male_SS_STAR_cpm)
median(XCR_MvsM_order$male_SS_STAR_cpm)
((mean(XCR_MvsM_order$male_SS_STAR_cpm))/(mean(XCR_MvsM_order$male_def_STAR_cpm)))

#--------------
# MALE
# HISAT X chromosome regions mean and median 
#--------------
# call in PAR1 genes 
PAR1_MvsM<-subset(chrX_MvsM_def_SS_HISAT_cpm_order, gene %in% PAR1_gene)
PAR1_MvsM_order <- PAR1_MvsM[order(PAR1_MvsM$gene),]
PAR1_MvsM_order$gene <- factor(PAR1_MvsM_order$gene, levels = chrX_genes$gene)
mean(PAR1_MvsM_order$male_def_HI_cpm)
median(PAR1_MvsM_order$male_def_HI_cpm)
mean(PAR1_MvsM_order$male_SS_HI_cpm)
median(PAR1_MvsM_order$male_SS_HI_cpm)
((mean(PAR1_MvsM_order$male_SS_HI_cpm))/(mean(PAR1_MvsM_order$male_def_HI_cpm)))

# call in PAR2 genes 
PAR2_MvsM<-subset(chrX_MvsM_def_SS_HISAT_cpm_order, gene %in% PAR2_gene)
PAR2_MvsM_order <- PAR2_MvsM[order(PAR2_MvsM$gene),]
PAR2_MvsM_order$gene <- factor(PAR2_MvsM_order$gene, levels = chrX_genes$gene)
mean(PAR2_MvsM_order$male_def_HI_cpm)
median(PAR2_MvsM_order$male_def_HI_cpm)
mean(PAR2_MvsM_order$male_SS_HI_cpm)
median(PAR2_MvsM_order$male_SS_HI_cpm)
((mean(PAR2_MvsM_order$male_SS_HI_cpm))/(mean(PAR2_MvsM_order$male_def_HI_cpm)))

# call in XTR genes 
XTR_MvsM<-subset(chrX_MvsM_def_SS_HISAT_cpm_order, gene %in% XTR_gene)
XTR_MvsM_order <- XTR_MvsM[order(XTR_MvsM$gene),]
XTR_MvsM_order$gene <- factor(XTR_MvsM_order$gene, levels = chrX_genes$gene)
mean(XTR_MvsM_order$male_def_HI_cpm)
median(XTR_MvsM_order$male_def_HI_cpm)
mean(XTR_MvsM_order$male_SS_HI_cpm)
median(XTR_MvsM_order$male_SS_HI_cpm)
((mean(XTR_MvsM_order$male_SS_HI_cpm))/(mean(XTR_MvsM_order$male_def_HI_cpm)))

# call in XDR genes 
XDR_MvsM<-subset(chrX_MvsM_def_SS_HISAT_cpm_order, gene %in% XDR_gene)
XDR_MvsM_order <- XDR_MvsM[order(XDR_MvsM$gene),]
XDR_MvsM_order$gene <- factor(XDR_MvsM_order$gene, levels = chrX_genes$gene)
mean(XDR_MvsM_order$male_def_HI_cpm)
median(XDR_MvsM_order$male_def_HI_cpm)
mean(XDR_MvsM_order$male_SS_HI_cpm)
median(XDR_MvsM_order$male_SS_HI_cpm)
((mean(XDR_MvsM_order$male_SS_HI_cpm))/(mean(XDR_MvsM_order$male_def_HI_cpm)))

# call in XAR genes 
XAR_MvsM<-subset(chrX_MvsM_def_SS_HISAT_cpm_order, gene %in% XAR_gene)
XAR_MvsM_order <- XAR_MvsM[order(XAR_MvsM$gene),]
XAR_MvsM_order$gene <- factor(XAR_MvsM_order$gene, levels = chrX_genes$gene)
mean(XAR_MvsM_order$male_def_HI_cpm)
median(XAR_MvsM_order$male_def_HI_cpm)
mean(XAR_MvsM_order$male_SS_HI_cpm)
median(XAR_MvsM_order$male_SS_HI_cpm)
((mean(XAR_MvsM_order$male_SS_HI_cpm))/(mean(XAR_MvsM_order$male_def_HI_cpm)))

# call in XCR genes 
XCR_MvsM<-subset(chrX_MvsM_def_SS_HISAT_cpm_order, gene %in% XCR_gene)
XCR_MvsM_order <- XCR_MvsM[order(XCR_MvsM$gene),]
XCR_MvsM_order$gene <- factor(XCR_MvsM_order$gene, levels = chrX_genes$gene)
mean(XCR_MvsM_order$male_def_HI_cpm)
median(XCR_MvsM_order$male_def_HI_cpm)
mean(XCR_MvsM_order$male_SS_HI_cpm)
median(XCR_MvsM_order$male_SS_HI_cpm)
((mean(XCR_MvsM_order$male_SS_HI_cpm))/(mean(XCR_MvsM_order$male_def_HI_cpm)))

