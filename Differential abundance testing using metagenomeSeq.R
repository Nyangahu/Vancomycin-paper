#Load R packages (To avoid errors when sourcing the custom functions doc, load all packages)
library(phyloseq)
library(ape)
library(BAT)
library(cluster)
library(corrplot)#for cor.mtest
library(dplyr)
library(fifer)
library(ggplot2)
library(graphics)
library(matrixStats)#rowSds
library(NMF)#for heatmap function
library(metagenomeSeq)#differential abundance testing
library(pheatmap)
library(plyr)
library(psych)#corr.test
library(reshape2)
library(RColorBrewer)
library(scales)
library(vegan)
library(picante)
library(gridExtra)
library(knitr)
library(colorspace)
library(Biostrings)
library(permute)
library(lattice)
library(nlme)
library(gplots)
library(Heatplus)
library(BIOM.utils)
library(scales)
library(Rcpp)
library(MASS)
library(fifer)
library(ggtree)
library(randomForest)
#Importing data (OTU table generated in QIIME, now in the biom format)
phy <- import_biom("otus_table.tax.biom", verbose = TRUE)
#Data clean up
colnames(tax_table(phy))
#Replace "Rank1" with Kingdom etc,
colnames(tax_table(phy)) <-  c("Kingdom", "Phylum" , "Class" , "Order" , "Family" , "Genus", "Species")
# Clean taxonomic annotations
tax_table(phy)[,"Kingdom"] <- sub("k__","",tax_table(phy)[,"Kingdom"])
tax_table(phy)[,"Phylum"] <- sub("p__","",tax_table(phy)[,"Phylum"])
tax_table(phy)[,"Class"] <- sub("c__","",tax_table(phy)[,"Class"])
tax_table(phy)[,"Order"] <- sub("o__","",tax_table(phy)[,"Order"])
tax_table(phy)[,"Family"] <- sub("f__","",tax_table(phy)[,"Family"])
tax_table(phy)[,"Genus"] <- sub("g__","",tax_table(phy)[,"Genus"])
tax_table(phy)[,"Species"] <- sub("s__","",tax_table(phy)[,"Species"])
#IMPORT METADATA
meta <-  read.table ("mapping file combined.txt", sep = "\t", header =TRUE, row.names=1)
#Check if the same sample numbers in meta file as the phy object
length(rownames(meta))
#(check that the sample names match in all cases)
length(intersect(rownames(meta),sample_names(phy)))
# Add metadata to phyloseq object.
#assign the metadata to the phyloseq object 'phy' (phyloseq will put these in the right order)
sample_data(phy) <- meta

#Subsetting samples to include only pup samples
D <- subset_samples(phy, Sample.type=="pup")
M <-subset_samples(phy, Sample.type=="mother")
##Differential abundance testing by metagenomeseq package (Method by Katie Lennard)
#Load the custom functions"microbiome_custom_functions.R"
source( "/Users/dnyang/Desktop/R 2017/Vanco 1/microbiome_custom_functions.R")
outDir= paste0(getwd(),"/")
#Merge at lowest taxonomic annotation
phy.new <- tax_glom.kv(physeq = D)
#First identify sample names for subsets:
colnames(sample_data(phy.new))
head(sample_data(phy.new))
u=as.character(unlist(unique(sample_data(phy.new)[,"Description"])))#convert to character string for matching
u
# Description
# Pg     Pups born to mothers treated during gestation
# Pc                     Pups born to untreated mothers
# Pgn Pups born to moms treated during gest and nursing
# Pn       Pups born to mothers treated during nursing

#so we are looking for pups born to untreated mothers vs pups born to mothers treated during a) gestation b) nursing and c) gest and nursing
IDs.1 = rownames(sample_data(phy.new)[sample_data(phy.new)[,"Description"]==u[1] | sample_data(phy.new)[,"Description"]==u[2],]) #pg vs. pc
length(IDs.1)#N=23. Repeat for all comparisons
#now we subset the phyloseq obj according to these IDs.1
physeq.1 <- prune_samples(sample_names(phy.new)%in% IDs.1,phy.new)
head(sample_data(physeq.1))
nsamples(physeq.1)#19
(sample_data(physeq.1)[,"Description"])#Inspect the description column
#now we do the differential abundance testing using metagenomeSeq (The super.fitZig.kv function is included in the custom functions file)
super.fitZig.kv(physeq=physeq.1,factor = "Description",outDir = outDir,FileName =paste0("1_25FC_0.2_pc_vs_pn 2"),
                heatmap.descriptor=c("tax_annot"), main=paste0("pc vs pn "), subt=c("subt = FDR < 0.05,|coeff| >= 1.25, >20%+ in either group"), 
                ordered=TRUE, p=0.05,FC = 1.25,colours=NULL, perc=0.2)



