#August 21,2021

#Reset R's Brain
rm(list=ls())

#Set working directory............................................................
setwd("C:/Users/fabipc/Dropbox/3-HolyFire/1-HolyFire-Pub1/1-Bacteria/")

#Load librarires..................................................................
library(phyloseq)#load qiime files and work w them based on .tza extension
library(qiime2R)
library(ape)#to build tree but I imported from QIIME
library(tidyverse) #required to load tax table
library(ggplot2)
library(ggpubr) 
library(ochRe)
library(MicEco)#tp create venndiagram

#Load rarefied data -------------------------------------------------------------------------------------------
Metadata<-read_tsv("1-Analysis/Metadata/MetaRare.tsv")#Rare metadata
Table<-read_qza("Qiime/core-metrics-results/rarefied_table.qza", tmp = "C:/tmp")#Rarefied table
Tree<-read_qza("Qiime/Phylo/Bacteria-rooted-tree-nc.qza",tmp = "C:/tmp" )
Taxonomy<-read_qza("Qiime/Taxonomy/Bacteria-taxonomy-paired-Clean-2020-132.qza",tmp = "C:/tmp")
Taxtable<-Taxonomy$data %>% as_tibble() %>% separate(Taxon, sep = ";",
                        c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

#CREATE PHYLOSEQ ARTIFACT-----------------------------------------------------------------------------
physeq<-phyloseq(
  otu_table(Table$data, taxa_are_rows = TRUE), 
  phy_tree(Tree$data), 
  tax_table(as.data.frame(Taxtable) %>% column_to_rownames("Feature.ID") %>% as.matrix()), 
  sample_data(Metadata%>% as.data.frame() %>% column_to_rownames("SampleID"))
)

#----
#----
#********************************************************************************************************************----
#------------------------------------     QUALITY CONTROL        ---------------------------------------------
#********************************************************************************************************************----
rank_names(physeq)# Look at rank names


#Quality control: Remove the g__ from each rank number..............................................................
colnames(tax_table(physeq))= c("Kingdom","Phylum","Class","Order","Family","Genus","Species", "Confidence")
tax_table(physeq)[, "Kingdom"] <- gsub("D_0__", "", tax_table(physeq)[, "Kingdom"])
tax_table(physeq)[, "Phylum"] <- gsub("D_1__", "", tax_table(physeq)[, "Phylum"])
tax_table(physeq)[, "Class"] <- gsub("D_2__", "", tax_table(physeq)[, "Class"])
tax_table(physeq)[, "Order"] <- gsub("D_3__", "", tax_table(physeq)[, "Order"])
tax_table(physeq)[, "Family"] <- gsub("D_4__", "", tax_table(physeq)[, "Family"])
tax_table(physeq)[, "Genus"] <- gsub("D_5__", "", tax_table(physeq)[, "Genus"])
tax_table(physeq)[, "Species"] <- gsub("D_6__", "", tax_table(physeq)[, "Species"])


#--SUBSET DATA BETWEEN BURNED AND UNBURNED SAMPLES ...........................................................
#First convert TSF to factor .....................................................
sample_data(physeq)$TSFdays<-factor(sample_data(physeq)$TSFdays) # make as factor

#Convert data to factors
df <- as.data.frame(lapply(sample_data(physeq),function (y) if(class(y)!="factor" ) as.factor(y) else y),stringsAsFactors=T)
row.names(df) <- sample_names(physeq)
sample_data(physeq) <- sample_data(df)
str(sample_data(physeq))



# * * Subset All to maintain only the burned samples of all data.................................
physeqBall<-subset_samples(physeq,Treatment=="Burned");physeqBall#29819x229
physeqUnAll<-subset_samples(physeq,Treatment=="Unburned");physeqUnAll#29819x110


# * * Subset data for soil only  ...................................................................
physeqSoil<-subset_samples(physeq,SampleType=="Soil");physeqSoil#29819x306
sample_names(physeqSoil); sample_data(physeqSoil)$TSFdays

PhySoilB<-subset_samples(physeqSoil,Treatment=="Burned");PhySoilB#29819x208
sample_names(PhySoilB);sample_data(PhySoilB)$TSFdays

PhySoilUn<-subset_samples(physeqSoil,Treatment=="Unburned");PhySoilUn#29819x98
sample_names(PhySoilUn); sample_data(PhySoilUn)$TSFdays


# * * Subset data for Rain only ........................................................----
physeqRain<-subset_samples(physeq,SampleType=="Rain");physeqRain#29819x33
sample_names(physeqRain) 

PhyRainB<-subset_samples(physeqRain,Treatment=="Burned");PhyRainB#29819x21
sample_names(PhyRainB);sample_data(PhyRainB)$TSFdays 

PhyRainUn<-subset_samples(physeqRain,Treatment=="Unburned");PhyRainUn#29819x12
sample_names(PhyRainUn);sample_data(PhyRainUn)$TSFday

#----
#----
#**************************************************************************************************************-----
# RELATIVE ABUNDANCE SOIL SAMPLES ONLY--------------------------------------------------------------------------
#**************************************************************************************************************-----

# * * * LOOK AT SOIL SAMPLES ONLY  -----------------------------------------------------------------------
#** Treatment analysis.................................................................................
GenusSoil <- tax_glom(physeqSoil, taxrank = 'Genus') # agglomerate taxa
(GenTrtSoil = merge_samples(GenusSoil, "Treatment")) # merge samples on sample variable of interest
RelGenTrtSoil <- transform_sample_counts(GenTrtSoil, function(x) x/sum(x)) #get abundance in %
RelGenTrtSoil1 <- psmelt(RelGenTrtSoil) # create dataframe from phyloseq object
RelGenTrtSoil1$Genus <- as.character(RelGenTrtSoil1$Genus) #convert to character
RelGenTrtSoil1$Genus[RelGenTrtSoil1$Abundance < 0.02] <- "< 2% abund." #rename genera with < 2% abundance


#Severity analysis.........................................................................................
GenusSoilSev <- tax_glom(physeqSoil, taxrank = 'Genus')
(GenSevSoil = merge_samples(GenusSoilSev, "AshSeverity"))
RelGenSevSoil <- transform_sample_counts(GenSevSoil, function(x) x/sum(x))
RelGenSevSoil1 <- psmelt(RelGenSevSoil)
RelGenSevSoil1$Genus <- as.character(RelGenSevSoil1$Genus)
RelGenSevSoil1$Genus[RelGenSevSoil1$Abundance < 0.02] <- "< 2% abund."


#TSF burned analysis.................................................................................
GenusSoilB <- tax_glom(PhySoilB, taxrank = 'Genus')
(GenTsfsSoilB = merge_samples(GenusSoilB, "TSFdays")) 
sample_data(GenTsfsSoilB)$TSFdays <- factor(sample_names(GenTsfsSoilB))
RelGenTsfSoilB <- transform_sample_counts(GenTsfsSoilB, function(x) x/sum(x)) 
RelGenTsfSoilB1 <- psmelt(RelGenTsfSoilB)
RelGenTsfSoilB1$Genus <- as.character(RelGenTsfSoilB1$Genus)
RelGenTsfSoilB1$Genus[RelGenTsfSoilB1$Abundance < 0.03] <- "< 3% abund."


#TSF Unburned analysis.................................................................................
GenusSoilUn <- tax_glom(PhySoilUn, taxrank = 'Genus')
(GenTsfSoilUn = merge_samples(GenusSoilUn, "TSFdays"))
sample_data(GenTsfsSoilUn)$TSFdays <- factor(sample_names(GenTsfsSoilUn))
RelGenTsfSoilUn <- transform_sample_counts(GenTsfSoilUn, function(x) x/sum(x))
RelGenTsfSoilUn1 <- psmelt(RelGenTsfSoilUn)
RelGenTsfSoilUn1$Genus <- as.character(RelGenTsfSoilUn1$Genus) 
RelGenTsfSoilUn1$Genus[RelGenTsfSoilUn1$Abundance < 0.03] <- "< 3% abund." 



#Phylum per treatment...................................................................................
PhylumSoil <- tax_glom(physeqSoil, taxrank = 'Phylum') # agglomerate taxa
(PhyTrtSoil = merge_samples(PhylumSoil, "Treatment")) # merge samples on sample variable of interest
RelPhyTrtSoil <- transform_sample_counts(PhyTrtSoil, function(x) x/sum(x)) #get abundance in %
RelPhyTrtSoil1 <- psmelt(RelPhyTrtSoil) # create dataframe from phyloseq object
RelPhyTrtSoil1$Phylum <- as.character(RelPhyTrtSoil1$Phylum) #convert to character
#RelPhyTrtSoil1$Phylum[RelPhyTrtSoil1$Abundance < 0.02] <- "< 2% abund." #rename Phyera with < 2% abundance



# TSF burned phylum analysis.................................................................................
PhylumSoilB <- tax_glom(PhySoilB, taxrank = 'Phylum')
(PhyTsfsSoilB = merge_samples(PhylumSoilB, "TSFdays")) 
sample_data(PhyTsfsSoilB)$TSFdays <- factor(sample_names(PhyTsfsSoilB))
RelPhyTsfSoilB <- transform_sample_counts(PhyTsfsSoilB, function(x) x/sum(x)) 
RelPhyTsfSoilB1 <- psmelt(RelPhyTsfSoilB)
RelPhyTsfSoilB1$Phylum <- as.character(RelPhyTsfSoilB1$Phylum)
#RelPhyTsfSoilB1$Phylum[RelPhyTsfSoilB1$Abundance < 0] <- "< 3% abund."



#TSF PHYLUM Unburned analysis.................................................................................
PhylumSoilUn <- tax_glom(PhySoilUn, taxrank = 'Phylum')
(PhyTsfSoilUn = merge_samples(PhylumSoilUn, "TSFdays"))
sample_data(PhyTsfSoilUn)$TSFdays <- factor(sample_names(PhyTsfSoilUn))
RelPhyTsfSoilUn <- transform_sample_counts(PhyTsfSoilUn, function(x) x/sum(x))
RelPhyTsfSoilUn1 <- psmelt(RelPhyTsfSoilUn)
RelPhyTsfSoilUn1$Phylum <- as.character(RelPhyTsfSoilUn1$Phylum) 
#RelPhyTsf1SoilUn$Phylum[RelPhyTsf1SoilUn$Abundance < 0.03] <- "< 3% abund." 



#----
#----
# **********************************************************************************************************----
#E---------EXPORT ALL RELATIVE ABUNDANCE FILES FROM ABOVE ------------------------------------------------------
# **********************************************************************************************************----

#Create directory for storage ...............................................................................
dir.create(file.path("1-Analysis/RelativeAbundance/Tables/Barplots/Soil"), recursive = TRUE)
dir.create(file.path("1-Analysis/RelativeAbundance/Tables/Barplots/Rain"), recursive = TRUE)#saved but not part of paper
dir.create(file.path("1-Analysis/RelativeAbundance/Tables/Barplots/SampleType"), recursive = TRUE)

dir.create(file.path("1-Analysis/RelativeAbundance/Graphs/Soil"), recursive = TRUE)
#dir.create(file.path("1-Analysis/RelativeAbundance/Graphs/Rain"), recursive = TRUE)..removed from analysis not part of paper

#Only soil samples ..........................................................................................................
#* Treatment...................
write.csv(RelGenTrtSoil1, "1-Analysis/RelativeAbundance/Tables/Barplots/Soil/RelativeAbund-Species-Trt-Soil.csv")
write.csv(RelGenSevSoil1, "1-Analysis/RelativeAbundance/Tables/Barplots/Soil/RelativeAbund-Genus-Sev-Soil-2Per.csv")

#TimeSinceFire,,,,,,,,,,,,,,,,,
write.csv(RelGenTsfSoilB1, "1-Analysis/RelativeAbundance/Tables/Barplots/Soil/RelativeAbund-Genus-TSF-Burned-Soil-3Per.csv")
write.csv(RelGenTsfSoilUn1, "1-Analysis/RelativeAbundance/Tables/Barplots/Soil/RelativeAbund-Genus-TSF-Unburned-Soil-3Per.csv")
write.csv(RelPhyTsfSoilB1, "1-Analysis/RelativeAbundance/Tables/Barplots/Soil/RelativeAbund-Phyla-TSF-Burned-Soil.csv")
write.csv(RelPhyTsfSoilUn1, "1-Analysis/RelativeAbundance/Tables/Barplots/Soil/RelativeAbund-Phyla-TSF-Unburned-Soil.csv")


#----
#-----
#********************************************************************************************************************************----
# ---------------------------------------           GRAPHS              ------------------------------------------------------
#********************************************************************************************************************************----

#READ IN THE FILES

#Only soil samples ..........................................................................................................
RelGenTrtSoil1<-read.csv("1-Analysis/RelativeAbundance/Tables/Barplots/Soil/RelativeAbund-Species-Trt-Soil.csv", row.names = 1)
RelGenSevSoil1<-read.csv("1-Analysis/RelativeAbundance/Tables/Barplots/Soil/RelativeAbund-Genus-Sev-Soil-2Per.csv", row.names = 1)
RelGenTsfSoilB1<-read.csv("1-Analysis/RelativeAbundance/Tables/Barplots/Soil/RelativeAbund-Genus-TSF-Burned-Soil-3Per.csv", row.names = 1)
RelGenTsfSoilUn1<-read.csv("1-Analysis/RelativeAbundance/Tables/Barplots/Soil/RelativeAbund-Genus-TSF-Unburned-Soil-3Per.csv", row.names = 1)
RelPhyTsfSoilB1<-read.csv("1-Analysis/RelativeAbundance/Tables/Barplots/Soil/RelativeAbund-Phyla-TSF-Burned-Soil.csv", row.names = 1)
RelPhyTsfSoilUn1<-read.csv("1-Analysis/RelativeAbundance/Tables/Barplots/Soil/RelativeAbund-Phyla-TSF-Unburned-Soil.csv", row.names = 1)



#Only Soil samples..........................
length(unique(RelGenTrtSoil1$Genus))#13
length(unique(RelGenSevSoil1$Genus))#13
length(unique(RelGenTsfSoilB1$Genus))#18
length(unique(RelGenTsfSoilUn1$Genus))#6
length(unique(RelPhyTsfSoilB1$Phylum))#35
length(unique(RelPhyTsfSoilUn1$Phylum))#35



#ORDER DATA ...........................................................................

# RENAMING TAXA FOR GRAPHING AND REORDERING --------------------------------------------------------------------------------------------
#* * * SOIL SAMPLES.............................................................................................
str(RelGenTsfSoilUn1$Genus)

RelGenTrtSoil1$Genus[RelGenTrtSoil1$Genus == "Candidatus Udaeobacter"] <- "Can. Udaeobacter" 
RelGenSevSoil1$Genus[RelGenSevSoil1$Genus == "Candidatus Udaeobacter"] <- "Can. Udaeobacter" 
RelGenTsfSoilB1$Genus[RelGenTsfSoilB1$Genus == "Candidatus Udaeobacter"] <- "Can. Udaeobacter" 


RelGenTsfSoilUn1$Genus[RelGenTsfSoilUn1$Genus =="Burkholderia-Caballeronia-Paraburkholderia" ]<-"Burkholderia" 
RelGenTsfSoilUn1$Genus[RelGenTsfSoilUn1$Genus == "Candidatus Udaeobacter"] <- "Can. Udaeobacter" 

#Reorder soil values .......................................................................................
RelGenSevSoil1$Sample<-ordered(RelGenSevSoil1$Sample, c("Unburned","Low","Moderate","High"))

#Genus TSF
RelGenTsfSoilB1$Sample <- factor(RelGenTsfSoilB1$Sample , levels = c("17","25","34","67","95","131","187","286","376"))
RelGenTsfSoilUn1$Sample <- factor(RelGenTsfSoilUn1$Sample , levels = c("17","25","34","67","95","131","187","286","376"))

#Phylum TSF
RelPhyTsfSoilB1$Sample <- factor(RelPhyTsfSoilB1$Sample , levels = c("17","25","34","67","95","131","187","286","376"))
RelPhyTsfSoilUn1$Sample <- factor(RelPhyTsfSoilUn1$Sample , levels = c("17","25","34","67","95","131","187","286","376"))


#CREATE COLOR SCHEME FOR PLOTS................................................................................
Trt<-c("< 2% abund."="#dee0df","Bacillus"="#34113d","Bryobacter"="#9C5850",
        "Can. Udaeobacter"="#E3BDA5","Cohnella"="#608a70","Domibacillus"="#DAE7DF",
        "Gemmatimonas"="#4d4b2d","Massilia"="#706f41","Noviherbaspirillum"="#d2dba8",
        "Paenibacillus"="#ebda96","RB41"="#cd5725","Sphingomonas"="#f0d1e1",
        "uncultured"="#3d3d40")
        

Burned<-c("< 3% abund."="#dee0df","Adhaeribacter"="#928a94","Bacillus"="#34113d",
        "Blastococcus"="#632a2a","Can. Udaeobacter"="#E3BDA5","Cohnella"="#608a70",
        "Conexibacter"="#9CBFA9","Domibacillus"="#DAE7DF","Gemmatimonas"="#4d4b2d",
        "Massilia"="#706f41","Mucilaginibacter"="#9a9858","Noviherbaspirillum"="#d2dba8",
        "Paenibacillus"="#ebda96","Pedobacter"="#7b4f2c","RB41"="#cd5725",
        "Solirubrobacter"="#cfcfcf","Tumebacillus"="#949494","uncultured"="#3d3d40")


Unburned<-c("< 3% abund."="#dee0df","Bryobacter"="#9C5850","Can. Udaeobacter"="#E3BDA5",
            "RB41"="#cd5725","Sphingomonas"="#f0d1e1","uncultured"="#3d3d40")



#----
#----
#*****************************************************************************************************----
# ------------------------- GRAPHS PER SOIL SAMPLES ONLY -------------------------------------------------
#*****************************************************************************************************----

#TREATMENT...............................................................................
AbunTrtSoil<-ggplot(data=RelGenTrtSoil1, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = Trt)+
  theme_bw()+ 
  theme(plot.margin = margin(1.5, 1.5, 0.8, 0.8, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 24, colour = "black"), 
        axis.title.x = element_text(size = 24, colour = "black"), 
        axis.text.y = element_text(size=18, hjust=0.5, colour = "black"),
        axis.text.x = element_text(size=18, colour = "black"),
        legend.position="right",
        legend.text = element_text(face = "italic", size=18),
        legend.title=element_text(size=21, face="bold"))+ 
  xlab("Time Since Fire (days)")+
  ylab("Relative Abundance")+ 
  guides(fill=guide_legend(ncol = 1));AbunTrtSoil



#BURNED PLOTS ...............................................................................
AbunBsoil<-ggplot(data=RelGenTsfSoilB1, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = Burned)+
  theme_bw()+ 
  theme(plot.margin = margin(1.5, 1.5, 0.8, 0.8, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 24, colour = "black"), 
        axis.title.x = element_text(size = 24, colour = "black"), 
        axis.text.y = element_text(size=18,angle=90,hjust=0.60,colour = "black"),
        axis.text.x = element_text(size=18, colour = "black"),
        legend.position="bottom",
        legend.text = element_text(face = "italic", size=18),
        #legend.title=element_text(size=21, face="bold"))+ 
        legend.title=element_blank())+ #size=26, face="bold"
  xlab("Time Since Fire (days)")+
  ylab("Relative Abundance")+ 
  guides(fill=guide_legend(ncol = 1));AbunBsoil


#UNBURNED PLOTS ...............................................................................
AbunUnSoil<-ggplot(data=RelGenTsfSoilUn1, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = Unburned)+
  theme_bw()+ 
  theme(plot.margin = margin(1.5, 1.5, 0.8, 0.8, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 24, colour = "black"), 
        axis.title.x = element_text(size = 24, colour = "black"), 
        axis.text.y = element_text(size=18,angle=90,hjust=0.60,colour = "black"),
        axis.text.x = element_text(size=18,colour = "black"),
        legend.position="bottom",
        legend.text = element_text(face = "italic", size=18),
        legend.title=element_text(size=21, face="bold"))+ 
  xlab("Time Since Fire (days)")+
  ylab("Relative Abundance")+ 
  guides(fill=guide_legend(ncol = 1));AbunUnSoil


#Export individual plots before making panels ...........................................
pdf("1-Analysis/RelativeAbundance/Graphs/Soil/TRT-Genus.pdf", height=7, width=8)
AbunTrtSoil
dev.off()

pdf("1-Analysis/RelativeAbundance/Graphs/Soil/TSF-Genus-Unburned.pdf", height=8, width=8)
AbunUnSoil
dev.off()

pdf("1-Analysis/RelativeAbundance/Graphs/Soil/TSF-Genus-Burned.pdf", height=8, width=8)
AbunBsoil
dev.off()



#PLOT DESIGN ........................................................................................
GenTSF<-ggarrange(AbunUnSoil, AbunBsoil, ncol=2, nrow=1, labels= c("A Unburned","B Burned"), 
              font.label = list(size = 22, color ="black"),common.legend = FALSE,
              legend = "bottom", align = "hv");GenTSF


#Rerun the graphs, but make the guide_legend(ncol-1)
GenTSF2<-ggarrange(AbunUnSoil, AbunBsoil, ncol=1, nrow=2,labels= c("A Unburned","B Burned"), 
               font.label = list(size = 22, color ="black"),
                   common.legend = FALSE,legend = "right", align = "hv");GenTSF2



#Export Panels....................................................................................................
pdf("1-Analysis/RelativeAbundance/Graphs/Soil/TSF-Genus-Unburned.vs.Burned-Hori.pdf", height=8, width=12)
GenTSF
dev.off()

pdf("1-Analysis/RelativeAbundance/Graphs/Soil/TSF-Genus-Unburned.vs.Burned-vertical.pdf", height=12, width=10)
GenTSF2
dev.off()





#PHYLUM PLOTS ...............................................................................------
PhyCol<-c("Acidobacteria"= "","Actinobacteria","Armatimonadetes", "Bacteroidetes",
          "BRC1","Chlamydiae","Chloroflexi","Cyanobacteria","Deinococcus-Thermus",
          "Dependentiae","Elusimicrobia","Entotheonellaeota","Euryachaeota","FBP",
          "fcpu426","Fibrobacteres","Firmicutes","Fusobacteria","Gemmatimonadetes",
          "Hydrogenedentes","Kiritimatiellaeota","Latescibacteria","Margulisbacteria")
          "Nitrospirae")

TrtSoilPhy<-ggplot(data=RelPhyTrtSoil1, aes(x=Sample, y=Abundance, fill=Phylum))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = PhyCol)+
  theme_bw()+ 
  theme(plot.margin = margin(1.5, 1.5, 0.8, 0.8, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 24, colour = "black"), 
        axis.title.x = element_text(size = 24, colour = "black"), 
        axis.text.y = element_text(size=24,angle=90,hjust=0.60,colour = "black"),
        axis.text.x = element_text(size=24, colour = "black"),
        legend.position="right",
        legend.text = element_text(face = "italic", size=26),
        legend.title=element_text(size=24, face="bold"))+ 
  xlab("Treatment")+
  ylab("Relative Abundance")+ 
  guides(fill=guide_legend(ncol = 1));TrtSoilPhy



AbunBsoil<-ggplot(data=RelPhyTsfSoilB1, aes(x=Sample, y=Abundance, fill=Phylum))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = PhyCol)+
  theme_bw()+ 
  theme(plot.margin = margin(1.5, 1.5, 0.8, 0.8, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 18, colour = "black"), 
        axis.title.x = element_text(size = 18, colour = "black"), 
        axis.text.y = element_text(size=20,angle=90,hjust=0.60,colour = "black"),
        axis.text.x = element_text(size=20,colour = "black"),
        legend.position="bottom",
        legend.text = element_text(face = "italic", size=20),
        legend.title=element_blank())+ #size=26, face="bold"
  xlab("Time Since Fire (days)")+
  ylab("Relative Abundance")+ 
  guides(fill=guide_legend(ncol = 3));AbunBsoil


#UNBURNED Phyla PLOTS ...............................................................................
AbunUnSoil<-ggplot(data=RelPhyTsfSoilUn1, aes(x=Sample, y=Abundance, fill=Phylum))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = PhyCol)+
  theme_bw()+ 
  theme(plot.margin = margin(1.5, 1.5, 0.8, 0.8, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 30, colour = "black"), 
        axis.title.x = element_text(size = 30, colour = "black"), 
        axis.text.y = element_text(size=20,angle=90,hjust=0.60,colour = "black"),
        axis.text.x = element_text(size=20,colour = "black"),
        legend.position="bottom",
        legend.text = element_text(face = "italic", size=20),
        legend.title=element_text(size=26, face="bold"))+ 
  xlab("Time Since Fire (days)")+
  ylab("Relative Abundance")+ 
  guides(fill=guide_legend(ncol = 3));AbunUnSoil








##----
####
#*****************************************************************************************************----
# ------------------------- GRAPHS PER RAIN SAMPLES ONLY -------------------------------------------------
#*****************************************************************************************************----

#TREATMENT-------------------------------------------------------- # NOT SURE Y THEY ARE AT 200-------------------
#Look only at soil samples
AbunTrtRain<-ggplot(data=RelGenTrtRain1, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = GenusCol)+
  theme_bw()+ 
  theme(plot.margin = margin(1.5, 1.5, 0.8, 0.8, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 38, colour = "black",
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(size = 36, colour = "black"), 
        axis.text.y = element_text(size=30, hjust=0.5, colour = "black"),
        axis.text.x = element_text(size=30, colour = "black",
                                   margin = margin(t = 0, r =0, b = 30, l = 0)),
        legend.position="right",
        legend.text = element_text(face = "italic", size=28),
        legend.title=element_text(size=30, face="bold"))+ 
  xlab("Time Since Fire (days)")+
  ylab("Relative Abundance")+ 
  guides(fill=guide_legend(ncol = 1));AbunTrtRain


#BURNED RAIN PLOTS-----------------------------------------------------------------------------------
AbunBrain<-ggplot(data=RelGenTsf1RainB, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = GenusCol)+
  theme_bw()+ 
  theme(plot.margin = margin(1.5, 1.5, 0.8, 0.8, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 30, colour = "black",
              margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.title.x = element_text(size = 30, colour = "black"), 
        axis.text.y = element_text(size=26,angle=90, hjust=0.5, colour = "black"),
        axis.text.x = element_text(size=28, angle=45, hjust=1,colour = "black",
               margin = margin(t = 0, r =0, b = 30, l = 0)),
        legend.position="right",
        legend.text = element_text(face = "italic", size=18),
        legend.title=element_text(size=19, face="bold"))+ 
  xlab("Time Since Fire (days)")+
  ylab("Relative Abundance")+ 
  guides(fill=guide_legend(ncol = 1));AbunBrain

AbunBrain$data$Sample<- factor(AbunBrain$data$Sample, 
      levels = c("17","25","34","67")); AbunBrain

#UNBURNED RAIN PLOTS-----------------------------------------------------------------------------------
AbunUnRain<-ggplot(data=RelGenTsf1RainUn, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = GenusCol)+
  theme_bw()+ 
  theme(plot.margin = margin(1.5, 1.5, 0.8, 0.8, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 30, colour = "black",
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.title.x = element_text(size = 30, colour = "black"), 
        axis.text.y = element_text(size=26,angle=90, hjust=0.5,colour = "black"),
        axis.text.x = element_text(size=28, angle=45,hjust=1,colour = "black",
                                   margin = margin(t = 0, r =0, b = 30, l = 0)),
        legend.position="right",
        legend.text = element_text(face = "italic", size=18),
        legend.title=element_text(size=19, face="bold"))+ 
  xlab("Time Since Fire (days)")+
  ylab("Relative Abundance")+ 
  guides(fill=guide_legend(ncol = 1));AbunUnRain

AbunUnRain$data$Sample <- factor(AbunUnRain$data$Sample, 
            levels = c("17","25","34","67")); AbunUnRain


#PLOT DESIGN---------------------------------------------------------------------------------------
GenTSF<-ggarrange(AbunBrain,AbunUnRain,  ncol=2, nrow=1, labels= c("a Burned","b Unburned"), 
                  font.label = list(size = 22, color ="black"),
                  common.legend = FALSE,legend = "bottom", align = "hv");GenTSF

GenTSF2<-ggarrange(AbunBrain,AbunUnRain,  ncol=1, nrow=2, labels= c("a Burned","b Unburned"), 
                   font.label = list(size = 22, color ="black"),
                   common.legend = FALSE,legend = "bottom", align = "hv");GenTSF2



#**************************************************************************************************************----
# -----------------------------EXPORT GRAPHS-----------------------------------------------------------------------
#**************************************************************************************************************----

pdf("1-Analysis/RelativeAbundance/Fungi/Graphs/Rain/TRT-Genus.pdf", height=10, width=12)
AbunTrtRain
dev.off()

pdf("1-Analysis/RelativeAbundance/Fungi/Graphs/Rain/TSF-Genus-Unburned.pdf", height=10, width=12)
AbunUnRain
dev.off()

pdf("1-Analysis/RelativeAbundance/Fungi/Graphs/Rain/TSF-Genus-Burned.pdf", height=10, width=11)
AbunBrain
dev.off()



#Export Panels...............................................................................................
#Legend size were 28 header, 24 for text
pdf("1-Analysis/RelativeAbundance/Fungi/Graphs/Rain/TSF-Genus-Unburned.vs.Burned-Hori.pdf", height=8, width=18)#28 x 24
GenTSF
dev.off()

pdf("1-Analysis/RelativeAbundance/Fungi/Graphs/Rain/TSF-Genus-Unburned.vs.Burned-vertical.pdf", height=18, width=10)#26x22
GenTSF2
dev.off()





