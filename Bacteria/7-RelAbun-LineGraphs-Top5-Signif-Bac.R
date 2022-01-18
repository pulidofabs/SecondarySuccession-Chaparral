#Nov 18, 2020

#Reset R's Brain
rm(list=ls())

setwd("C:/Users/fabipc/Dropbox/3-HolyFire/1-HolyFire-Pub1/1-Bacteria/")

#Load librarires-------------------------------------------------------------------------------------------------------------------
library(phyloseq)# to load in qiime files and work w them based on .tza extension
library(qiime2R)
library(ape)#to build tree but I imported from QIIME
library(tidyverse) #required to load tax table
library(ggplot2)
library(ggpubr) 


#To create the graphs by just importing the previously quantified phyloseq data scroll down and skip calculations--------


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
#--------
#********************************************************************************************************************----
#------------------------------------     QUALITY CONTROL        ---------------------------------------------
#********************************************************************************************************************----

rank_names(physeq)# Look at rank names

#Quality control: Remove the g__ from each rank number
colnames(tax_table(physeq))= c("Kingdom","Phylum","Class","Order","Family","Genus","Species",    "Confidence")
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

# * * Subset All to maintain only the burned samples of all data..........................
physeqAllB<-subset_samples(physeq,Treatment=="Burned");physeqAllB#29819x229
physeqAllUn<-subset_samples(physeq,Treatment=="Unburned");physeqAllUn#29819x110

# * * Subset data for soil only  ........................................................
physeqSoil<-subset_samples(physeq,SampleType=="Soil");physeqSoil#29819x306
sample_names(physeqSoil); sample_data(physeqSoil)$TSFdays

PhySoilB<-subset_samples(physeqSoil,Treatment=="Burned");PhySoilB#29819x208
sample_names(PhySoilB);sample_data(PhySoilB)$TSFdays

PhySoilUn<-subset_samples(physeqSoil,Treatment=="Unburned");PhySoilUn#298195x98
sample_names(PhySoilUn); sample_data(PhySoilUn)$TSFdays




#.----
#.----
#********************************************************************************************************************----
# RELATIVE ABUNDANCE SOIL  --------------------------------------------------------------------------------------
#********************************************************************************************************************----
# Treatment.............................................................................
GenAglom <- tax_glom(physeqSoil, "Genus")#Agglomerate by taxa of interest 
GenTop5<-names(sort(taxa_sums(GenAglom), TRUE)[1:5])#Top taxa
Gen5prune<-prune_taxa(GenTop5, physeqSoil) #Prune 
Gen5RA <- transform_sample_counts(Gen5prune, function(x) x / sum(x))#Rel Abund
Gen5Trt <- merge_samples(Gen5RA, "TSFdays") #Merge by cat of interest
RAGen5Trt <- psmelt(Gen5Trt)##Melt samples to create table to work in base R.............
names(RAGen5Trt)[2]<-"TSF"

# Burned.............................................................................
GenB <- tax_glom(PhySoilB, "Genus")#Agglomerate by taxa of interest 
GenTop5B<-names(sort(taxa_sums(GenB), TRUE)[1:5])#get top 5 most abundant genus
Gen5PruneB<-prune_taxa(GenTop5B, PhySoilB) #Prune 
Gen5RAB <- transform_sample_counts(Gen5PruneB, function(x) x / sum(x))#Calculate relative abundance
Gen5TSFB <- merge_samples(Gen5RAB, "TSFdays") #Mersge samples based on category of interest
Gen5RAtsfB<- transform_sample_counts(Gen5TSFB, function(x) x / sum(x))#calculate raltive abundance per TSF
Gen5Burn<-psmelt(Gen5RAtsfB)#Crate table to work in base R
names(Gen5Burn)[2]<-"TSF" #add label to sample column containing TSF created by phyloseq


#unburned.............................................................................
GenUn <- tax_glom(PhySoilUn, "Genus")#Agglomerate by taxa of interest 
GenTop5Un<-names(sort(taxa_sums(GenUn), TRUE)[1:5])#get top 5 most abundant genus
Gen5PruneUn<-prune_taxa(GenTop5Un, PhySoilUn) #Prune 
Gen5RAUn <- transform_sample_counts(Gen5PruneUn, function(x) x / sum(x))#Calculate relative abundance
Gen5TSFUn <- merge_samples(Gen5RAUn, "TSFdays") #Mersge samples based on category of interest
Gen5RAtsfUn<- transform_sample_counts(Gen5TSFUn, function(x) x / sum(x))#calculate raltive abundance per TSF
Gen5Unburn<-psmelt(Gen5RAtsfUn)#Create table to work in base R
names(Gen5Unburn)[2]<-"TSF" #add label to sample column containing TSF created by phyloseq



# EXPORT DATA FILES --------------------------------------------------------------------------------------------------------------
dir.create(file.path("1-Analysis/RelativeAbundance/Tables/LineGraph/Soil"), recursive=TRUE)


#Export soil samples ......................................................................................................
write.csv(RAGen5Trt,"1-Analysis/RelativeAbundance/Tables/LineGraph/Soil/Data-Trt-Top5Genus-Soil.csv")
write.csv(Gen5Burn , "1-Analysis/RelativeAbundance/Tables/LineGraph/Soil/Data-TSF-Top5Genus-Burned-Soil.csv")
write.csv(Gen5Unburn, "1-Analysis/RelativeAbundance/Tables/LineGraph/Soil/Data-TSF-Top5Genus-Unburned-Soil.csv")




#----
#----
#*************************************************************************************************************----
# CREATE PLOTS SOIL ONLY  --------------------------------------------------------------------------------------
#*************************************************************************************************************----
library(dplyr)


#** * IF NEED BE REIMPORT THE FILES FOR ANALYSIS----
RAGen5Trt<-read.csv("1-Analysis/RelativeAbundance/Tables/LineGraph/Soil/Data-Trt-Top5Genus-Soil.csv")
Gen5Burn<-read.csv("1-Analysis/RelativeAbundance/Tables/LineGraph/Soil/Data-TSF-Top5Genus-Burned-Soil.csv")
Gen5Unburn<-read.csv("1-Analysis/RelativeAbundance/Tables/LineGraph/Soil/Data-TSF-Top5Genus-Unburned-Soil.csv")


# *** Reorder labels for plotting......................................................................... 
Gen5Burn$TSF<-factor(Gen5Burn$TSF, levels = c("17","25","34","67","95","131","187","286","376"))
Gen5Unburn$TSF<-factor(Gen5Unburn$TSF, levels = c("17","25","34","67","95","131","187","286","376"))


# * * Make a table containing the relative abundance of each of the 5 taxa of interest ........................
Gen5BurnTable <-Gen5Burn%>%group_by(TSF, Genus)
Gen5BurnTable$TSF<-factor(Gen5BurnTable$TSF);Gen5BurnTable

Gen5UnburnTable<- Gen5Unburn%>% group_by(TSF, Genus)
Gen5UnburnTable$TSF<-factor(Gen5UnburnTable$TSF);Gen5UnburnTable


#Create manual color scale for graphs......................................................................
GenColB<- c("#8a3500","#ffbb1c","#c7b797","#8fc3ff","#1f5498")
GenColUn<-c("#594a2d","#557A95","#83677B","#a34126","#a7c4ba")

GenColors<-c("Bacillus"="#7f6e99","Massilia"="#45a340",
            "Noviherbaspirillum"="#acd1e6","uncultured"="#757473", 
            "Paenibacillus"="#734c3d")

# * PLOTS -------------------------------------------------------------------------------------------------------------------
Top5GenBsoil<-ggplot(Gen5Burn, aes(x=TSF, y=Abundance, group=Genus))+
  geom_line(aes(linetype=Phylum,color=Genus), size=1.2)+
  geom_point(aes(color=Genus),size=3)+ 
  scale_linetype_manual(values=c("solid","twodash"))+
  scale_color_manual(values = GenColors)+
  theme_bw()+
  theme(panel.grid= element_blank(),
        text = element_text(size=28), 
        axis.title.y = element_text(size = 28), 
        axis.text.y = element_text(size=26, angle = 90, hjust=0.5, colour = "black"),
        axis.text.x = element_text(size=26,colour = "black"),
        legend.text=element_text(size=26,face = "italic"), 
        legend.position = "right")+
  labs(x = "Time-Since-Fire (days)", y="Relative Abundance");Top5GenBsoil




#Does not quite work, various changes...........................................
Top5GenUsoil<-ggplot(Gen5Unburn, aes(x=TSF, y=Abundance, group=Genus))+
  geom_line(aes(linetype=Phylum,color=Genus), size=1.2)+
  geom_point(aes(color=Genus),size=3)+ 
  scale_linetype_manual(values=c("solid","twodash"))+
  scale_color_manual(values = GenColUn)+
  theme_bw()+
  theme(panel.grid= element_blank(),
        text = element_text(size=28), 
        axis.title.y = element_text(size = 28), 
        axis.text.y = element_text(size=26, angle = 90, hjust=0.5, colour = "black"),
        axis.text.x = element_text(size=26,colour = "black"),
        legend.text=element_text(size=26,face = "italic"), 
        legend.position = "right")+
  labs(x = "Time-Since-Fire (days)", y="Relative Abundance");Top5GenUsoil



#EXPORT FILES AS PDF------------------------------------------------------------------------------------------
dir.create(file.path("1-Analysis/RelativeAbundance/Graphs/Soil/Lineplots"), recursive=TRUE)


pdf("1-Analysis/RelativeAbundance/Graphs/Soil/Lineplots/Top5-Genus-Burned-linearSoil.pdf", height=6.5, width=10)
Top5GenBsoil
dev.off()

pdf("1-Analysis/RelativeAbundance/Graphs/Soil/Lineplots/Top5-Genus-Unburned-linearSoil.pdf", height=8, width=12)
Top5GenUsoil
dev.off()



#----
#------
#***********************************************************************************************************************-----
#----------------------TEST FOR SIGNIFICANCE --------------------------------------------------------------------------
#**********************************************************************************************************************----
#Load data---------------------------------------------------------------------------------------------
#Asco<-read.csv("Analysis/RelativeAbundance/TSF-Line/Tables/Ascos-Top3-Genus.csv")
#Asco$TSF<- factor(Asco$TSF, levels = c("17","25","34","67","95","131","187","286","376"))


Gen5Burn<-read.csv("1-Analysis/RelativeAbundance/Tables/LineGraph/Soil/Data-TSF-Top5Genus-Burned-Soil.csv")
Gen5Burn$TSF<- factor(Gen5Burn$TSF, levels = c("17","25","34","67","95","131","187","286","376"))



# * * Subset .............................................................................
attach(Gen5Burn)
Pyronema<-Gen5Burn[which(Genus == "Pyronema"), ]; head(Pyronema)
Geminibasidium<-Gen5Burn[which(Genus == "Geminibasidium"), ]; head(Geminibasidium)
Inocybe<-Gen5Burn[which(Genus == "Inocybe"), ]; head(Inocybe)
Aspergillus<-Gen5Burn[which(Genus == "Aspergillus"), ]; head(Aspergillus)
Penicillium<-Gen5Burn[which(Genus == "Penicillium"), ]; head(Penicillium)


#Quality control--check for normality.......................................................
par(mfrow=c(3,2))
plot(Pyronema$TSF, Pyronema$Abundance)
plot(Geminibasidium$TSF, Geminibasidium$Abundance)
plot(Inocybe$TSF, Inocybe$Abundance)
plot(Aspergillus$TSF, Aspergillus$Abundance)
plot(Penicillium$TSF, Penicillium$Abundance)





#.----
#.----
#**********************************************************************************************----
# NEGATIVE BINOMIAL REGRESSION --------------------------------------------------------------------
#**********************************************************************************************----

#RUN REGRESSION -------------------------------------------------------------------
library(MuMIn) #to get pseudo R
library(lme4)
library(lmerTest)#Gives us P values 
#Using lmer to use site as the random effect

#The totals have been consolidated per timepoint, no need
#to control for sampling

#https://www.jaredknowles.com/journal/2013/11/25/getting-started-with-mixed-effect-models-in-r
#Random intercep for nested: fit nested group effect terms: Our samples are nested within Site and Cardinal.Direction


PyroLM<-glm(Abundance ~ TSFdays,data = Pyronema);summary(PyroLM)
LMPyroSum<-summary(PyroLM);LMPyroSum#0.8956

InoLM<-glm(Inocybe$Abundance~TSFdays, data = Inocybe)
LMInoSum<-summary(InoLM);LMInoSum#0.02767 *

GeminiLM<-glm(Geminibasidium$Abundance~TSFdays, data = Geminibasidium)
LMGemiSum<-summary(GeminiLM);LMGemiSum# 0.0706 .

AsperLM<-glm(Aspergillus$Abundance~TSFdays, data = Aspergillus)
LMAsperSum<-summary(AsperLM);LMAsperSum#0.000636 ***

PeniLM<-glm(Penicillium$Abundance~TSFdays, data = Penicillium)
LMPeniSum<-summary(PeniLM);LMPeniSum #0.013003 * 



# * *  Export results of LME4 ----------------------------------------------------------------------------------
dir.create(file.path("1-Analysis/RelativeAbundance/Tables/LineGraph/Soil"), recursive=TRUE)


capture.output(LMPyroSum, file="1-Analysis/RelativeAbundance/Tables/LineGraph/Soil/Top5Pyronema-TSF.csv")
capture.output(LMInoSum, file="1-Analysis/RelativeAbundance/Tables/LineGraph/Soil/Top5Inocybe-TSF.csv")
capture.output(LMGemiSum, file = "1-Analysis/RelativeAbundance/Tables/LineGraph/Soil/Top5Geminibasidium-TSF.csv")
capture.output(LMAsperSum, file="1-Analysis/RelativeAbundance/Tables/LineGraph/Soil/Top5Aspergillus-TSF.csv")
capture.output(LMPeniSum, file = "1-Analysis/RelativeAbundance/Tables/LineGraph/Soil/Top5Penicillium-TSF.csv")















