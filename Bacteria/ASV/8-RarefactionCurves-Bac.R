#June 20, 2020

#Reset R's Brain
rm(list=ls())

#Set working directory-------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/3-HolyFire/1-HolyFire-Pub1/1-Bacteria")


#Load librarires-------------------------------------------------------------------------------------------------------------------
library(phyloseq)# to load in qiime files and work w them based on .tza extension
library(qiime2R)
library(ape)#to build tree but I imported from QIIME
library(tidyverse) #required to load tax table
library(ggplot2)
library(ggpubr) 
library(ochRe)

#LOAD QIIME DATA-NON-NORMALIZED BECUASE I WILL RAREFY IN PHYLOSEQ -------------------------------------------------------------------------------------------
Metadata<-read_tsv("1-Analysis/Metadata/Bacteria-Metadata-lib1-4E-MJ-Clean-NC.tsv")#Rare metadata
Table<-read_qza("Qiime/DADA2/Soil/Bacteria-Table-Soil-1-4EMJ-FilteredNC.qza", tmp = "C:/tmp")#unrarefied table
Tree<-read_qza("Qiime/Phylo/Bacteria-rooted-tree-NC.qza",tmp = "C:/tmp" )
Taxonomy<-read_qza("Qiime/Taxonomy/Bacteria-taxonomy-paired-Clean-2020-132.qza",tmp = "C:/tmp")
Taxtable<-Taxonomy$data %>% as_tibble() %>% separate(Taxon, sep = ";",
                   c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

#CREATE PHYLOSEQ ARTIFACT-----------------------------------------------------------------------------
physeq<-phyloseq(
  otu_table(Table$data, taxa_are_rows = TRUE), 
  phy_tree(Tree$data), 
  tax_table(as.data.frame(Taxtable) %>% column_to_rownames("Feature.ID") %>% as.matrix()), 
  sample_data(Metadata%>% as.data.frame() %>% column_to_rownames("#SampleID"))
)


#********************************************************************************************************************----
#------------------------------------     QUALITY CONTROL        ---------------------------------------------
#********************************************************************************************************************----

rank_names(physeq)# Look at rank names

#Quality control: Remove the g__ from each rank number
colnames(tax_table(physeq))= c("Kingdom","Phylum","Class","Order","Family","Genus","Species", "Confidence")
tax_table(physeq)[, "Kingdom"] <- gsub("D_0__", "", tax_table(physeq)[, "Kingdom"])
tax_table(physeq)[, "Phylum"] <- gsub("D_1__", "", tax_table(physeq)[, "Phylum"])
tax_table(physeq)[, "Class"] <- gsub("D_2__", "", tax_table(physeq)[, "Class"])
tax_table(physeq)[, "Order"] <- gsub("D_3__", "", tax_table(physeq)[, "Order"])
tax_table(physeq)[, "Family"] <- gsub("D_4__", "", tax_table(physeq)[, "Family"])
tax_table(physeq)[, "Genus"] <- gsub("D_5__", "", tax_table(physeq)[, "Genus"])
tax_table(physeq)[, "Species"] <- gsub("D_6__", "", tax_table(physeq)[, "Species"])



#--SUBSET DATA BETWEEN BURNED AND UNBURNED SAMPLES ...........................................................
sample_data(physeq)$TSFdays<-factor(sample_data(physeq)$TSFdays) # make as factor

physeqBurn<-subset_samples(physeq,Treatment=="Burned");physeqBurn#33078x215
sample_names(physeqBurn)

physeqUnburn<-subset_samples(physeq,Treatment=="Unburned");physeqUnburn#33078x101
sample_names(physeqUnburn)



#----
#----
#********************************************************************************************************************----
#------------------------------------     ESTIMATE SPECIES RICHNESS       ---------------------------------------------
#********************************************************************************************************************----

# using GlobalPatterns
library(ranacapa)

#Rarefy dataset.....................................................................
RarefyPlotFun <- ggrare(physeq, step = 7115, color = "Plot",  se = TRUE)

#Create graph.......................................................................
#* Treatment............................................................
TrtFun <- RarefyPlotFun + facet_wrap(~Treatment, scales = "free")+
  theme()+
  theme(text = element_text(size=20), 
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.spacing = unit(1, "lines"),
        legend.position="bottom",strip.text = element_text(size = 18));TrtFun


#* Plot............................................................
SiteFun <- RarefyPlotFun + facet_wrap(~Plot, scales = "free")+
  theme()+
  theme(text = element_text(size=20), 
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 16),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.spacing = unit(1, "lines"),
        legend.position="bottom",strip.text = element_text(size = 18));SiteFun


#* TSFdays............................................................
TSFplotFun <- RarefyPlotFun + facet_wrap(~TSFdays, scales = "free")+
  theme()+
  theme(text = element_text(size=20), 
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.spacing = unit(1, "lines"),
        legend.position="bottom",strip.text = element_text(size = 18));TSFplotFun




#*********************************************************************************************----
#--EXPORT GRAPS ----------------------------------------------------------------------------------
#*********************************************************************************************----


pdf("1-Analysis/Diversity/RarefactionCurves-Treatment.pdf", height=6, width=10.5)
TrtFun
dev.off()


pdf("1-Analysis/Diversity/Alpha/RarefactionCurves-Site.pdf", height=10, width=12)
SiteFun
dev.off()

pdf("1-Analysis/Diversity/Alpha/RarefactionCurves-TSF.pdf", height=10, width=12)
TSFplotFun
dev.off()





