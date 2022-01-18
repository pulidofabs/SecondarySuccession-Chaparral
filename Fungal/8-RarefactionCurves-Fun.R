#Nov 22, 2020

#Reset R's Brain
rm(list=ls())

#Set working directory-------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/3-HolyFire/1-HolyFire-Pub1/2-Fungi")


#Load librarires-------------------------------------------------------------------------------------------------------------------
library(phyloseq)# to load in qiime files and work w them based on .tza extension
library(qiime2R)
library(ape)#to build tree but I imported from QIIME
library(tidyverse) #required to load tax table
library(ggplot2)
library(ggpubr) 
library(ochRe)

#IMPORT PHYLOSEQ FILES FOR GRAPHS BELOW-------------**************

#LOAD QIIME DATA-NON-NORMALIZED BECUASE I WILL RAREFY IN PHYLOSEQ -------------------------------------------------------------------------------------------
Metadata<-read_tsv("Qiime/Metadata/Fungal-Metadata-Lib1-4E-MJ-Clean-NC.tsv")#Rare metadata
Table<-read_qza("Qiime/DADA2/Clean/NoControls/Fungal-Table-Merged-1-4E-MJ-NC.qza", tmp = "C:/tmp")#none rarefied table will rarefy here
Tree<-read_qza("Qiime/Phylo/Fungal-rooted-tree-NC.qza",tmp = "C:/tmp" )
Taxonomy<-read_qza("Qiime/Taxonomy/Fungal-taxonomy-paired-RC-Clean.qza",tmp = "C:/tmp")
Taxtable<-Taxonomy$data %>% as_tibble() %>% separate(Taxon, sep = ";",
                   c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

#CREATE PHYLOSEQ ARTIFACT................................................................
physeq<-phyloseq(
  otu_table(Table$data, taxa_are_rows = TRUE), 
  phy_tree(Tree$data), 
  tax_table(as.data.frame(Taxtable) %>% column_to_rownames("Feature.ID") %>% as.matrix()), 
  sample_data(Metadata%>% as.data.frame() %>% column_to_rownames("FungalSeqID")))


#----
#----
#********************************************************************************************************************----
#------------------------------------     QUALITY CONTROL        ---------------------------------------------
#********************************************************************************************************************----

rank_names(physeq)# Look at rank names

#Quality control: Remove the g__ from each rank number
colnames(tax_table(physeq))= c("Kingdom","Phylum","Class","Order","Family","Genus","Species", "Confidence")
tax_table(physeq)[, "Kingdom"] <- gsub("k__", "", tax_table(physeq)[, "Kingdom"])
tax_table(physeq)[, "Phylum"] <- gsub("p__", "", tax_table(physeq)[, "Phylum"])
tax_table(physeq)[, "Class"] <- gsub("c__", "", tax_table(physeq)[, "Class"])
tax_table(physeq)[, "Order"] <- gsub("o__", "", tax_table(physeq)[, "Order"])
tax_table(physeq)[, "Family"] <- gsub("f__", "", tax_table(physeq)[, "Family"])
tax_table(physeq)[, "Genus"] <- gsub("g__", "", tax_table(physeq)[, "Genus"])
tax_table(physeq)[, "Species"] <- gsub("s__", "", tax_table(physeq)[, "Species"])


#--SUBSET DATA BETWEEN BURNED AND UNBURNED SAMPLES ...........................................................
sample_data(physeq)$TSFdays<-factor(sample_data(physeq)$TSFdays) # make as factor

# * * Subset All to maintain only the burned samples of all data.................................
physeqAllB<-subset_samples(physeq,Treatment=="Burned");physeqAllB#14578x234
physeqAllUn<-subset_samples(physeq,Treatment=="Unburned");physeqAllUn#14578x116


# * * Subset data for soil only  ...................................................................
physeqSoil<-subset_samples(physeq,SampleType=="Soil");physeqSoil#14578x314
sample_names(physeqSoil); sample_data(physeqSoil)$TSFdays

PhySoilB<-subset_samples(physeqSoil,Treatment %in% "Burned");PhySoilB#14578x210
sample_names(PhySoilB);sample_data(PhySoilB)$TSFdays

PhySoilUn<-subset_samples(physeqSoil,Treatment=="Unburned");PhySoilUn#14578x104
sample_names(PhySoilUn); sample_data(PhySoilUn)$TSFdays

print(Metadata)
# * * Subset data for Rain only ........................................................----
physeqRain<-subset_samples(physeq,SampleType=="Rain");physeqRain##14578x36
sample_names(physeqRain) 

PhyRainB<-subset_samples(physeqRain,Treatment=="Burned");PhyRainB#14578x24
sample_names(PhyRainB);sample_data(PhyRainB)$TSFdays 

PhyRainUn<-subset_samples(physeqRain,Treatment=="Unburned");PhyRainUn#14578x12
sample_names(PhyRainUn);sample_data(PhyRainUn)$TSFdays


#---
#----
#********************************************************************************************************************----
#------------------------------   ESTIMATE SPECIES RICHNESS SOIL SAMPLES   ---------------------------------------------
#********************************************************************************************************************----

#RAREFY THE DATA TO SAME VALUES AS WE DID MANUALLY ...........................................
library(ranacapa)
RarefyPlotFun <- ggrare(physeqSoil, step = 11058, color = "Plot",  se = TRUE)#Rarefy data


TrtFun <- RarefyPlotFun + 
  facet_wrap(~Treatment, scales = "free")+
  theme()+
  theme(text = element_text(size=20), 
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.spacing = unit(1, "lines"),
        legend.position="bottom",
        strip.text = element_text(size = 18));TrtFun



SiteFun <- RarefyPlotFun + facet_wrap(~Plot, scales = "free")+
  theme()+
  theme(text = element_text(size=20), 
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 16),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.spacing = unit(1, "lines"),
        legend.position="bottom",strip.text = element_text(size = 18));SiteFun



TSFplotFun <- RarefyPlotFun + facet_wrap(~TSFdays, scales = "free")+
  theme()+
  theme(text = element_text(size=20), 
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.spacing = unit(1, "lines"),
        legend.position="bottom",strip.text = element_text(size = 18));TSFplotFun



#----
#----
#*********************************************************************************************----
#--EXPORT GRAPS ----------------------------------------------------------------------------------
#*********************************************************************************************----


pdf("1-Analysis/Diversity/Alpha/Fungi/Graphs/RarefactionCurves-Treatment.pdf", height=6, width=10.5)
TrtFun
dev.off()

pdf("1-Analysis/Diversity/Alpha/Fungi/Graphs/RarefactionCurves-Site.pdf", height=10, width=12)
SiteFun
dev.off()

pdf("1-Analysis/Diversity/Alpha/Fungi/Graphs/RarefactionCurves-TSF.pdf", height=10, width=12)
TSFplotFun
dev.off()





