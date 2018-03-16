#For Github, generating the stacked plot showing maternal microbiota (Fig S1C).
#Load phyloseq and ggplot2 packages

library(phyloseq)
library(ggplot2)

#Importing data (OTU table generated in QIIME, now in the biom format)
import_biom(BIOMfilename = "Vancomycin repeat.biom", verbose = TRUE)
phy <- import_biom("Vancomycin repeat.biom", verbose = TRUE)
#Check sample names in the phyloseq object
sample_names(phy)
# Number of OTUs in the phyloseq object
ntaxa(phy)
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
meta <-  read.table ("Mapping vancorepeat II.txt", sep = "\t", header =TRUE, row.names=1)
head(sample_names(phy))
#Check number of samples in the phyloseq object
length(sample_names(phy))#78
#Check if the same sample numbers in meta file match those in the phyloseq object
length(rownames(meta))
#(check that the sample names match in all cases)
length(intersect(rownames(meta),sample_names(phy)))
# Add metadata to phyloseq object.
#assign the metadata to the phyloseq object 'phy' (phyloseq will put these in the right order)
sample_data(phy) <- meta
phy
ntaxa(phy)
#Subsetting samples into various groups contained in the metadata
sample_variables(phy)
Genital <- subset_samples(phy, Sample.type=="Genital Tract")
BM4 <-subset_samples(phy, Sample.type=="Breast Milk D4")
BM14 <-subset_samples(phy, Sample.type=="Breast Milk D14")
FecalPup <- subset_samples(phy, Sample.type=="Feces Pups")
FecalMoms <- subset_samples(phy, Sample.type=="Feces Dam")
FecalMoms

#Now focus on "FecalMoms"
#Remove the NAs from FecalMoms
phy1 <- subset_taxa(FecalMoms, !is.na(Phylum) & !Phylum %in% c("NA"))
plot_bar(phy1, fill = "Phylum")
phy2 <- subset_taxa(phy1, !is.na(Class) & !Phylum %in% c("NA"))
plot_bar(phy2, fill = "Class")
phy3 <- subset_taxa(phy2, !is.na(Order) & !Phylum %in% c("NA"))
plot_bar(phy3, fill = "Order")
phy4 <- subset_taxa(phy3, !is.na(Family) & !Phylum %in% c("NA"))
plot_bar(phy4, fill = "Family")
phy5 <- subset_taxa(phy4, !is.na(Genus) & !Phylum %in% c("NA"))
phy6 <- subset_taxa(phy5, !is.na(Species) & !Phylum %in% c("NA"))

#Draw a stacked bar graph of maternal microbiota at phylum level
M1 <-psmelt(phy6)
#Initially, we used this script to create the figure in the paper but (labels = percent) is obsolete in the ggplot2 version am now using
ggplot(M1,aes(x = Description, y = Abundance,fill = Family)) + geom_bar(position = "fill",stat = "identity") + scale_y_continuous(labels = percent_format()) + theme(axis.text=element_text(size=12, face="bold"),axis.title=element_text(size=14,face="bold")) + theme(legend.title=element_text(size=14,face="bold")) + theme(legend.text=element_text(size=14))

#This has been replaced with this.

ggplot(M1,aes(x = Description, y = Abundance,fill = Phylum)) + geom_bar(position = "fill",stat = "identity") + scale_y_continuous(labels = scales::percent) + theme(axis.text=element_text(size=12, face="bold"),
                                                                                                                                                                     axis.title=element_text(size=14,face="bold")) + theme(legend.title=element_text(size=14,face="bold")) + theme(legend.text=element_text(size=14))


