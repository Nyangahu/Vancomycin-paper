#Load in the phyloseq, metagenomeSeq & other packages
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

#Importing data (OTU table generated in QIIME, now in the biom format)
phy <- import_biom("vanco.biom", verbose = TRUE)
sample_names(phy)

# Number of OTUs in the phyloseq object
ntaxa(phy)
# Rename filenames as required ie pup.10.g to pupg10
sample_names(phy) <- sub("pup10.gn","pupgn10",sample_names(phy))
sample_names(phy) <- sub("pup1.gn","pupgn1",sample_names(phy))
sample_names(phy) <- sub("mom2.G.and.N","Mumgn1",sample_names(phy))
sample_names(phy) <- sub("mom1.G.and.N","Mumgn2",sample_names(phy))
sample_names(phy) <- sub("mom1.nurs","Mumn1",sample_names(phy))
sample_names(phy) <- sub("mom1.ctrl","Mumc1",sample_names(phy))
sample_names(phy) <- sub("Mom1.gest","Mumg1",sample_names(phy))
sample_names(phy) <- sub("pup2.n","pupn2",sample_names(phy))
sample_names(phy) <- sub("pup1.g","pupg1",sample_names(phy))
sample_names(phy) <- sub("pup1.n","pupn1",sample_names(phy))
sample_names(phy) <- sub("pup4.n","pupn4",sample_names(phy))
sample_names(phy) <- sub("pup3.n","pupn3",sample_names(phy))
sample_names(phy) <- sub("pup6.n","pupn6",sample_names(phy))
sample_names(phy) <- sub("pup5.n","pupn5",sample_names(phy))
sample_names(phy) <- sub("pupg7n","pupgn7",sample_names(phy))
sample_names(phy) <- sub("pup1.c","pupc1",sample_names(phy))
sample_names(phy) <- sub("pup.10.g","pupg10",sample_names(phy))
sample_names(phy) <- sub("pup2.c","pupc2",sample_names(phy))
sample_names(phy) <- sub("pup3.c","pupc3",sample_names(phy))
sample_names(phy) <- sub("pup5.c","pupc5",sample_names(phy))
sample_names(phy) <- sub("pup4.c","pupc4",sample_names(phy))
sample_names(phy) <- sub("pup6.c","pupc6",sample_names(phy))
sample_names(phy) <- sub("pup7.g","pupg7",sample_names(phy))
sample_names(phy) <- sub("pup6.g","pupg7",sample_names(phy))
sample_names(phy) <- sub("pup8.g","pup8g",sample_names(phy))
sample_names(phy) <- sub("pup5.g","pupg5",sample_names(phy))
sample_names(phy) <- sub("pup8.gn","pupgn8",sample_names(phy))
sample_names(phy) <- sub("pup4.gn","pupgn4",sample_names(phy))
sample_names(phy) <- sub("pup2.g","pupg2",sample_names(phy))
sample_names(phy) <- sub("pup3.gn","pupgn3",sample_names(phy))
sample_names(phy) <- sub("pup6.gn","pupgn6",sample_names(phy))
sample_names(phy) <- sub("pup3.g","pupg3",sample_names(phy))
sample_names(phy) <- sub("pup2.gn","pupgn2",sample_names(phy))
sample_names(phy) <- sub("pup9.g","pupg9",sample_names(phy))
sample_names(phy) <- sub("pup.9.gn","pupgn9",sample_names(phy))
sample_names(phy) <- sub("pup5.gn","pupgn5",sample_names(phy))
sample_names(phy) <- sub("pup4.g","pupg4",sample_names(phy))
sample_names(phy) <- sub("pup7.gn","pupgn7",sample_names(phy))
sample_names(phy) <- sub("pup6.g","pupg6",sample_names(phy))
sample_names(phy) <- sub("pup8g","pupg8",sample_names(phy))
sample_names(phy) <- sub("pupg8n","pupgn8",sample_names(phy))
sample_names(phy) <- sub("pupg2n","pupgn2",sample_names(phy))
sample_names(phy) <- sub("pupg5n","pupgn5",sample_names(phy))
sample_names(phy) <- sub("pupg7n","pupgn7",sample_names(phy))
sample_names(phy)
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
meta <-  read.table ("Vancomycin R metadata.txt", sep = "\t", header =TRUE, row.names=1)
head(meta)
tail(meta)
length(sample_names(phy))#37
#Check if the same sample numbers in meta file as the phy object
length(rownames(meta))
#(check that the sample names match in all cases)
length(intersect(rownames(meta),sample_names(phy)))
# Add metadata to phyloseq object.
#assign the metadata to the phyloseq object 'phy' (phyloseq will put these in the right order)
sample_data(phy) <- meta
phy
#Subsetting samples to include only pup samples
D <- subset_samples(phy, Sample.type=="Pup")
M <-subset_samples(phy, Sample.type=="Mother")


#Run the "microbiome_custom_functions.R". These has various functions for microbiome analysis including metagenomeSeq
#Provide the path to where the file is located in your computer
source("C:/Users/Donald/Desktop/Vancomycin paper Microbiome journal/Vanco 1/microbiome_custom_functions.R")
D

#Removing NA that you don't want to appear on bar plots
phy1 <- subset_taxa(D, !is.na(Phylum) & !Phylum %in% c("NA"))
phy2 <- subset_taxa(phy1, !is.na(Class) & !Phylum %in% c("NA"))
phy3 <- subset_taxa(phy2, !is.na(Order) & !Phylum %in% c("NA"))
phy4 <- subset_taxa(phy3, !is.na(Family) & !Phylum %in% c("NA"))
phy5 <- subset_taxa(phy4, !is.na(Genus) & !Phylum %in% c("NA"))
phy6 <- subset_taxa(phy5, !is.na(Species) & !Phylum %in% c("NA"))
phy6

phy6#43 taxa

#Now merge at the lowest taxonomic annotation and create a new phyloseq object
phy.new <- tax_glom.kv(physeq = phy6)
phy.new#26 taxa
phy.new#Phyloseq object that has been merged at the lowest taxonomic level.
#Only look at species level
phy.new <- tax_glom(phy.new, taxrank= "Species")
#Lets subset the phyloseq object for the two queries you want to do:
#Firt identify sample names for subsets:
colnames(sample_data(phy.new))
head(sample_data(phy.new))
u=as.character(unlist(unique(sample_data(phy.new)[,"Description"])))#convert to character string for matching
u
# Description
# pupg10     Pups born to mothers treated during gestation
# pupc2                     Pups born to untreated mothers
# pupgn1 Pups born to moms treated during gest and nursing
# pupn2        Pups born to mothers treated during nursing

#so we are looking for pups born to untreated mothers vs pups born to mothers treated during a) gestation and b) nursing
IDs.1 = rownames(sample_data(phy.new)[sample_data(phy.new)[,"Description"]==u[1] | sample_data(phy.new)[,"Description"]==u[2],]) #pg vs. pc
length(IDs.1)#N=16
#Repeat for other comparisons
IDs.1 = rownames(sample_data(phy.new)[sample_data(phy.new)[,"Description"]==u[2] | sample_data(phy.new)[,"Description"]==u[4],]) #pc vs. pn
length(IDs.1)#N=12
#Repeat for other comparisons
IDs.1 = rownames(sample_data(phy.new)[sample_data(phy.new)[,"Description"]==u[2] | sample_data(phy.new)[,"Description"]==u[3],]) #pc vs. pgn
length(IDs.1)#N=16

#now we subset the phyloseq obj according to these IDs.1
physeq.1 <- prune_samples(sample_names(phy.new)%in% IDs.1,phy.new)
head(sample_data(physeq.1))
nsamples(physeq.1)#16
#lets look at what the description column looks like (which we will use as 'factor' argument in metagenomeSeq)
(sample_data(physeq.1)[,"Description"])#Inspect the description column

#now run differential abundance testing by metagenomeSeq
source("C:/Users/Donald/Desktop/Vancomycin paper Microbiome journal/Vanco 1/microbiome_custom_functions.R")
outDir= paste0(getwd(),"/")

super.fitZig.kv(physeq=physeq.1,factor = "Description",outDir = outDir,FileName =paste0("1_25FC_0.2_pc_vs_pg"),
                heatmap.descriptor=c("tax_annot"), main=paste0("pc vs pg "), subt=c("subt = FDR < 0.05,|coeff| >= 1.25, >20%+ in either group"), 
                ordered=TRUE, p=0.05,FC = 1.25,colours=NULL, perc=0.2)



