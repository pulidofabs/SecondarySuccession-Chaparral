#Nov 18, 2020

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


##Re-Import soil samples below, do not need to recalculate values


#LOAD QIIME DATA- -------------------------------------------------------------------------------------------
Metadata<-read_tsv("Metadata/Fungi/MetaRare.tsv")#Rare metadata
Table<-read_qza("Qiime/core-metrics-results/rarefied_table.qza", tmp = "C:/tmp")#Rarefied table
Tree<-read_qza("Qiime/Phylo/Fungal-rooted-tree-NC.qza",tmp = "C:/tmp" )
Taxonomy<-read_qza("Qiime/Taxonomy/Fungal-taxonomy-paired-RC-Clean.qza",tmp = "C:/tmp")
Taxtable<-Taxonomy$data %>% as_tibble() %>% separate(Taxon, sep = ";",
           c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

#CREATE PHYLOSEQ ARTIFACT-..........................................................................
physeq<-phyloseq(
  otu_table(Table$data, taxa_are_rows = TRUE), 
  phy_tree(Tree$data), 
  tax_table(as.data.frame(Taxtable) %>% column_to_rownames("Feature.ID") %>% as.matrix()), 
  sample_data(Metadata%>% as.data.frame() %>% column_to_rownames("SampleID")))



#----
#--------
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


#--
#----
#***************************************************************************************************************----
#-- SUBSET DATA (sample type and treatment) ------------------------------------------------------------------------
#***************************************************************************************************************----sample_data(physeq)$TSFdays<-factor(sample_data(physeq)$TSFdays) # make as factor

# * * Subset data for soil only  .................................................
physeqSoil<-subset_samples(physeq,SampleType=="Soil");physeqSoil#11595x302
sample_names(physeqSoil); sample_data(physeqSoil)$TSFdays

PhySoilB<-subset_samples(physeqSoil,Treatment=="Burned");PhySoilB#11595x199
sample_names(PhySoilB);sample_data(PhySoilB)$TSFdays

PhySoilUn<-subset_samples(physeqSoil,Treatment=="Unburned");PhySoilUn#11595x103
sample_names(PhySoilUn); sample_data(PhySoilUn)$TSFdays


#First convert TSF to factor ............
sample_data(physeqSoil)$TSFdays<-factor(sample_data(physeqSoil)$TSFdays)#as factor


#.----
#.----
#***********************************************************************************************************-----
# RELATIVE ABUNDANCE CALCULATIONS -------------------------------------------------------------------------------
#***********************************************************************************************************-----

#* * Treatment......................................................................
GenSoil <- tax_glom(physeqSoil, "Genus")
RelGenSoilTrt <- transform_sample_counts(GenSoil, function(x) x / sum(x))
RelGenSoilTrt1<-psmelt(RelGenSoilTrt)

#* * Burned ........................................................................
GenSoilB<- tax_glom(PhySoilB, "Genus")
RelGenSoilB <- transform_sample_counts(GenSoilB, function(x) x / sum(x))
RelGenSoilB1<-psmelt(RelGenSoilB)

#* * Unburned ......................................................................
GenSoilUn<- tax_glom(PhySoilUn, "Genus")
RelGenSoilUn<- transform_sample_counts(GenSoilUn, function(x) x / sum(x))
RelGenSoilUn1<-psmelt(RelGenSoilUn)



#.----
#.----
#**************************************************************************************************************************-----
# EXPORT DATA FILES --------------------------------------------------------------------------------------------------------------
#**************************************************************************************************************************-----
dir.create(file.path("1-Analysis/RelativeAbundance/Fungi/Tables/DeStats/Soil"), recursive=TRUE)

#Export soil samples ......................................................................................................
write.csv(RelGenSoilTrt1, "1-Analysis/RelativeAbundance/Fungi/Tables/DeStats/Soil/Trt-PerChange-Data-Soil.csv")
write.csv(RelGenSoilB1, "1-Analysis/RelativeAbundance/Fungi/Tables/DeStats/Soil/TSF-PerChange-Data-Burned-Soil.csv")
write.csv(RelGenSoilUn1, "1-Analysis/RelativeAbundance/Fungi/Tables/DeStats/Soil/TSF-PerChange-Data-Unburned-Soil.csv")


#Re-Import soil samples ......................................................................................................
RelGenSoilTrt1<-read.csv("1-Analysis/RelativeAbundance/Fungi/Tables/DeStats/Soil/Trt-PerChange-Data-Soil.csv")
RelGenSoilB1<-read.csv("1-Analysis/RelativeAbundance/Fungi/Tables/DeStats/Soil/TSF-PerChange-Data-Burned-Soil.csv")
RelGenSoilUn1<-read.csv("1-Analysis/RelativeAbundance/Fungi/Tables/DeStats/Soil/TSF-PerChange-Data-Unburned-Soil.csv")



#* Quality control ....................................................................................
#Set reference level to Unburn..................................................
RelGenSoilTrt1$Treatment<-as.factor(RelGenSoilTrt1$Treatment)
RelGenSoilTrt1$Treatment <- try(relevel(RelGenSoilTrt1$Treatment , "Unburned"))
levels(RelGenSoilTrt1$Treatment)




#---
#-----
#********************************************************************************************************************----
#---DESCRIPTIVE STATISTICS -------------------------------------------------------------------------------------------------
#********************************************************************************************************************----
#Treatment..............................................................
TrtStatSoil<-RelGenSoilTrt1 %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus,Treatment) %>%
  summarise(n_obs = n(),
            TotalASV=sum(Abundance,na.rm = TRUE),
            Average = mean(Abundance, na.rm = TRUE), 
            sd = sd(Abundance, na.rm=TRUE),
            min=min(Abundance, na.rm = TRUE),
            max =max(Abundance, na.rm = TRUE));TrtStatSoil

# *  * TIME SINCE FIRE SOIL ONLY........................................
TsfTrtStatSoilB<-RelGenSoilB1 %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus, Treatment, TSFdays) %>%
  summarise(n_obs = n(),
            TotalASV=sum(Abundance,na.rm = TRUE),
            Average = mean(Abundance, na.rm = TRUE), 
            sd = sd(Abundance, na.rm=TRUE),
            min=min(Abundance, na.rm = TRUE),
            max =max(Abundance, na.rm = TRUE));TsfTrtStatSoilB

TsfStatSoilB<-RelGenSoilB1 %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus, TSFdays) %>%
  summarise(n_obs = n(),
            TotalASV=sum(Abundance,na.rm = TRUE),
            Average = mean(Abundance, na.rm = TRUE), 
            sd = sd(Abundance, na.rm=TRUE),
            min=min(Abundance, na.rm = TRUE),
            max =max(Abundance, na.rm = TRUE));TsfStatSoilB


TsfStatSoilUn<-RelGenSoilUn1 %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus, TSFdays) %>%
  summarise(n_obs = n(),
            TotalASV=sum(Abundance,na.rm = TRUE),
            Average = mean(Abundance, na.rm = TRUE), 
            sd = sd(Abundance, na.rm=TRUE),
            min=min(Abundance, na.rm = TRUE),
            max =max(Abundance, na.rm = TRUE));TsfStatSoilUn




#----
#
#***********************************************************************************************************************----
#--- PERCENT CHANGE CALCULATIONS -------------------------------------------------------------------------------------------
#***********************************************************************************************************************----
library(dplyr)

#* * * SOIL SAMPLES ONLY ....................................................
TSFTrtPerChgSoil<-RelGenSoilTrt%>%
 # filter(!is.na(Genus)) %>%
  group_by(Genus,TSFdays,Treatment) %>%
  summarise(TotalASV=sum(Abundance,na.rm = TRUE),
            mean = mean(Abundance,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 2))%>%
  as.data.frame();TSFTrtPerChgSoil


TrtPerChgSoil<-RelGenSoilTrt%>%
  filter(!is.na(Genus)) %>%
  group_by(Genus, Treatment) %>%
  summarise(TotalASV=sum(Abundance,na.rm = TRUE),
            mean = mean(Abundance,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 2))%>%
  as.data.frame();TrtPerChgSoil


TsfPerChgSoilB<-RelGenSoilB1%>%
  filter(!is.na(Genus)) %>%
  group_by(Genus, TSFdays) %>%
  summarise(TotalASV=sum(Abundance,na.rm = TRUE),
  mean = mean(Abundance,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 2))%>%
  as.data.frame();TsfPerChgSoilB

TsfPerChgSoilUn<-RelGenSoilUn1%>%
  filter(!is.na(Genus)) %>%
  group_by(Genus, TSFdays) %>%
  summarise(TotalASV=sum(Abundance,na.rm = TRUE),
            mean = mean(Abundance,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 2))%>%
  as.data.frame();TsfPerChgSoilUn



#.----
#.----
#**************************************************************************************************************************-----
# EXPORT RESULTS OF DESCRIPTIVE STATISTICS--------------------------------------------------------------------------------------------------------------
#**************************************************************************************************************************-----
dir.create(file.path("1-Analysis/RelativeAbundance/Fungi/Tables/PerChange/Soil"), recursive=TRUE)

#Export descriptive statistics ......................................................................................................
write.csv(TrtStatSoil, "1-Analysis/RelativeAbundance/Fungi/Tables/DeStats/Soil/Trt-Soil-DeStats-Results.csv")
write.csv(TsfStatSoilB, "1-Analysis/RelativeAbundance/Fungi/Tables/DeStats/Soil/TSF-Soil-DeStats-Results-Burned.csv")
write.csv(TsfStatSoilUn , "1-Analysis/RelativeAbundance/Fungi/Tables/DeStats/Soil/TSF-Soil-DeStats-Results-Unburned.csv")
write.csv(TsfTrtStatSoilB , "1-Analysis/RelativeAbundance/Fungi/Tables/DeStats/Soil/Trt-TSF-Soil-DeStats-Results-Burned.csv")


#Export soil samples ......................................................................................................
write.csv(TrtPerChgSoil, "1-Analysis/RelativeAbundance/Fungi/Tables/PerChange/Soil/Trt-PerChange-Soil.csv")
write.csv(TSFTrtPerChgSoil, "1-Analysis/RelativeAbundance/Fungi/Tables/PerChange/Soil/TRT-TSF-PerChange-Soil.csv")
write.csv(TsfPerChgSoilB, "1-Analysis/RelativeAbundance/Fungi/Tables/PerChange/Soil/TSF-PerChange-Burned-Soil.csv")
write.csv(TsfPerChgSoilUn, "1-Analysis/RelativeAbundance/Fungi/Tables/PerChange/Soil/TSF-PerChange-Unburned-Soil.csv")




#----
#-----
#***********************************************************************************************************************-----
# CALCULATE SIGNIFICANCE----------------------------------------------------------------------------------------------------
#***********************************************************************************************************************-----
#Re-Import soil samples ......................................................................................................
RelGenSoilTrt<-read.csv("1-Analysis/RelativeAbundance/Fungi/Tables/DeStats/Soil/Trt-PerChange-Data-Soil.csv")
RelGenSoilB1<-read.csv("1-Analysis/RelativeAbundance/Fungi/Tables/DeStats/Soil/TSF-PerChange-Data-Burned-Soil.csv")
RelGenSoilUn1<-read.csv("1-Analysis/RelativeAbundance/Fungi/Tables/DeStats/Soil/TSF-PerChange-Data-Unburned-Soil.csv")





#SUBSET THE DATA FOR TREATMENT ........................................................
attach(RelGenSoilTrt)

Aspergillus<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "Aspergillus"),]; head(Aspergillus$Genus)
Balsamia<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "Balsamia"),]; head(Balsamia$Genus)
Cenococcum<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "Cenococcum"),]; head(Cenococcum$Genus)
Cladophialophora<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "Cladophialophora"),]; head(Cladophialophora$Genus)
Cortinarius<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "Cortinarius"),]; head(Cortinarius$Genus)
Geminibasidium<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "Geminibasidium"),]; head(Geminibasidium$Genus)
Geopora<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "Geopora"),]; head(Geopora$Genus)
Hyaloscypha<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "Hyaloscypha"),]; head(Hyaloscypha$Genus)
Inocybe<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "Inocybe"),]; head(Inocybe$Genus)
Penicillium<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "Penicillium"),]; head(Penicillium$Genus)
Pyronema<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "Pyronema"),]; head(Pyronema$Genus)
Tephrocybe<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "Tephrocybe"),]; head(Tephrocybe$Genus)
Thelephora<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "Thelephora"),]; head(Thelephora$Genus)
Tomentella<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "Tomentella"),]; head(Tomentella$Genus)
Trechispora<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "Trechispora"),]; head(Trechispora$Genus)
Venturia<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "Venturia"),]; head(Venturia$Genus)



#CHECK DATA FOR NORMALITY..........................................................................
par(mfrow=c(4,4))

#Check histogram and run shapiro test ...........................................
hist(Aspergillus$Abundance); shapiro.test(Aspergillus$Abundance)
hist(Balsamia$Abundance); shapiro.test(Balsamia$Abundance)
hist(Cenococcum$Abundance); shapiro.test(Cenococcum$Abundance)
hist(Cladophialophora$Abundance); shapiro.test(Cladophialophora$Abundance)
hist(Cortinarius$Abundance); shapiro.test(Cortinarius$Abundance)
hist(Geminibasidium$Abundance); shapiro.test(Geminibasidium$Abundance)
hist(Geopora$Abundance); shapiro.test(Geopora$Abundance)
hist(Hyaloscypha$Abundance); shapiro.test(Hyaloscypha$Abundance)
hist(Inocybe$Abundance); shapiro.test(Inocybe$Abundance)
hist(Penicillium$Abundance); shapiro.test(Penicillium$Abundance)
hist(Pyronema$Abundance); shapiro.test(Pyronema$Abundance)
hist(Tephrocybe$Abundance); shapiro.test(Tephrocybe$Abundance)
hist(Tomentella$Abundance); shapiro.test(Tomentella$Abundance)
hist(Trechispora$Abundance); shapiro.test(Trechispora$Abundance)
hist(Venturia$Abundance); shapiro.test(Venturia$Abundance)




#-----
#----
#*****************************************************************************************************************************************----
#SIGNIFICANCE PER TAXA PER TREATMENT----------------------------------------------------------------------------------------------
#*****************************************************************************************************************************************----
library(MASS)
library(lmerTest)
library(lme4)
library(MuMIn)


#Explore Temporal correlations to ensure that we are capturing as much of the variance............... 
ASVneg1 <- glmer(Aspergillus$Abundance ~ Aspergillus$Treatment + (1|Plot)+(1|Subplot)+(1|TSFdays), data = Aspergillus);Anova(ASVneg1)
ASVneg2 <- glmer(Aspergillus$Abundance ~ Aspergillus$Treatment +  (1|Plot)+(1|TSFdays), data = Aspergillus);Anova(ASVneg2)
ASVneg3 <- glmer(Aspergillus$Abundance ~ Aspergillus$Treatment + (1|Subplot)+(1|TSFdays), data = Aspergillus);Anova(ASVneg3)


AICc(ASVneg1,ASVneg2,ASVneg3)#ASVneg3 better


#Test for significance
library(car)#anova for glmer
Aspergillusglmer<- glmer(Aspergillus$Abundance ~ Aspergillus$Treatment + (1|Subplot)+(1|TSFdays), data = Aspergillus)
Anova(ASVneg1, type=3)#3.289e-07 ***

Balsamiaglmer<-glmer(Balsamia$Abundance ~ Balsamia$Treatment + (1|Subplot)+(1|TSFdays), data = Balsamia)
Anova(Balsamiaglmer, type=3)#0.01454 *

Cenococcumglmer<-glmer(Cenococcum$Abundance ~ Cenococcum$Treatment + (1|Subplot)+(1|TSFdays), data = Cenococcum)
Anova(Cenococcumglmer, type=3)#2.176e-05 ***

Cladophialophoraglmer<-glmer(Cladophialophora$Abundance ~ Cladophialophora$Treatment + (1|Subplot)+(1|TSFdays), data = Cladophialophora)
Anova(Cladophialophoraglmer)#2.9e-12 ***

Cortinariusglmer<-glmer(Cortinarius$Abundance ~ Cortinarius$Treatment + (1|Subplot)+(1|TSFdays), data = Cortinarius)
Anova(Cortinariusglmer, type=3)#2.621e-05 ***

Geminibasidiumglmer<-glmer(Geminibasidium$Abundance ~ Geminibasidium$Treatment + (1|Subplot)+(1|TSFdays), data = Geminibasidium)
Anova(Geminibasidiumglmer, type=3)#0.028511 *

Geoporaglmer<-glmer(Geopora$Abundance ~ Geopora$Treatment + (1|Subplot)+(1|TSFdays), data = Geopora)
Anova(Geoporaglmer, type=3)#1.038e-06 ***

Hyaloscyphaglmer<-glmer(Hyaloscypha$Abundance ~ Hyaloscypha$Treatment + (1|Subplot)+(1|TSFdays), data = Hyaloscypha)
Anova(Hyaloscyphaglmer, type=3)#1.246e-07 ***

Inocybeglmer<-glmer(Inocybe$Abundance ~ Inocybe$Treatment + (1|Subplot)+(1|TSFdays), data = Inocybe)
Anova(Inocybeglmer, type=3)#0.0009642 **

Penicilliumglmer<-glmer(Penicillium$Abundance ~ Penicillium$Treatment + (1|Subplot)+(1|TSFdays), data = Penicillium)
Anova(Penicilliumglmer, type=3)#8.942e-06 ***

Pyronemaglmer<-glmer(Pyronema$Abundance ~ Pyronema$Treatment + (1|Subplot)+(1|TSFdays), data = Pyronema)
Anova(Pyronemaglmer, type=3)# 2.130e-08 

Tephrocybeglmer<-glmer(Tephrocybe$Abundance ~ Tephrocybe$Treatment + (1|Subplot)+(1|TSFdays), data = Tephrocybe)
Anova(Tephrocybeglmer, type=3)#0.12641 

Tomentellaglmer<-glmer(Tomentella$Abundance ~ Tomentella$Treatment + (1|Subplot)+(1|TSFdays), data = Tomentella)
Anova(Tomentellaglmer, type=3)#0.05293 .

Trechisporaglmer<-glmer(Trechispora$Abundance ~ Trechispora$Treatment + (1|Subplot)+(1|TSFdays), data = Trechispora)
Anova(Trechisporaglmer, type=3)#5.86e-07 ***

Venturiaglmer<-glmer(Venturia$Abundance ~ Venturia$Treatment + (1|Subplot)+(1|TSFdays), data = Venturia)
Anova(Venturiaglmer, type=3)#4.055e-06 ***



# EXPORT RESULTS --------------------------------------------------------------------------------------------------------------------
dir.create(file.path("1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignifnif/Soil/Trt"), recursive=TRUE)

#Export ANOVA values ...............................................................................................................
capture.output(Anova(Balsamiaglmer), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignifnif/Soil/Trt/Balsamiaglmer.csv")
capture.output(Anova(Cenococcumglmer), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignifnif/Soil/Trt/Cenococcumglmer.csv")
capture.output(Anova(Cladophialophoraglmer), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignifnif/Soil/Trt/Cladophialophoraglmer.csv")
capture.output(Anova(Cortinariusglmer), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignifnif/Soil/Trt/Cortinariusglmer.csv")
capture.output(Anova(Geminibasidiumglmer), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignifnif/Soil/Trt/Geminibasidiumglmer.csv")
capture.output(Anova(Geoporaglmer), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignifnif/Soil/Trt/Geoporaglmer.csv")
capture.output(Anova(Hyaloscyphaglmer), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignifnif/Soil/Trt/Hyaloscyphaglmer.csv")
capture.output(Anova(Inocybeglmer), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignifnif/Soil/Trt/Inocybeglmer.csv")
capture.output(Anova(Penicilliumglmer), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignifnif/Soil/Trt/Penicilliumglmer.csv")
capture.output(Anova(Pyronemaglmer), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignifnif/Soil/Trt/Pyronemaglmer.csv")
capture.output(Anova(Tephrocybeglmer), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignifnif/Soil/Trt/Tephrocybeglmer.csv")
capture.output(Anova(Tomentellaglmer), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignifnif/Soil/Trt/Tomentellaglmer.csv")
capture.output(Anova(Trechisporaglmer), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignifnif/Soil/Trt/Trechisporaglmer.csv")
capture.output(Anova(Venturiaglmer), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignifnif/Soil/Trt/Venturiaglmer.csv")

detach(RelGenSoilTrt)

#-----
#-------
#************************************************************************************************************************----
# -- CHANGES OF TAXA WITH TIME SINCE FIRE -BURNED  ----------------------------------------------------------
##***********************************************************************************************************************----
attach(RelGenSoilB1)

#Subset the taxa of interest...................................................................................
AspergillusB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Aspergillus"),]; head(AspergillusB$Genus)
BalsamiaB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Balsamia"),]; head(BalsamiaB$Genus)
CatenuliferaB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Catenulifera"),]; head(CatenuliferaB$Genus)
ConiochaetaB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Coniochaeta"),]; head(ConiochaetaB$Genus)
CortinariusB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Cortinarius"),]; head(CortinariusB$Genus)
CoprinellusB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Coprinellus"),]; head(CoprinellusB$Genus)
GeminibasidiumB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Geminibasidium"),]; head(GeminibasidiumB$Genus)
InocybeB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Inocybe"),]; head(InocybeB$Genus)
LachnumB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Mallocybe"),]; head(LachnumB$Genus)
MycenaB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Mycena"),]; head(MycenaB$Genus)
PenicilliumB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Penicillium"),]; head(PenicilliumB$Genus)
PyronemaB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Pyronema"),]; head(PyronemaB$Genus)
RasamsoniaB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Rasamsonia"),]; head(RasamsoniaB$Genus)
TephrocybeB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Tephrocybe"),]; head(TephrocybeB$Genus)
TomentellaB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Tomentella"),]; head(TomentellaB$Genus)


#Explore Temporal correlations to ensure that we are capturing as much of the variance................................................................ 
ASVneg1 <- glmer(AspergillusB$Abundance ~ AspergillusB$TSFdays + (1|Plot)+(1|Subplot)+(1|TSFdays), data = AspergillusB);Anova(ASVneg1, Type=3)
ASVneg2 <- glmer(AspergillusB$Abundance ~ AspergillusB$TSFdays +  (1|Plot)+(1|TSFdays), data = AspergillusB);Anova(ASVneg2, Type=3)
ASVneg3 <- glmer(AspergillusB$Abundance ~ AspergillusB$TSFdays + (1|Subplot)+(1|TSFdays), data = AspergillusB);Anova(ASVneg3, Type=3)

AICc(ASVneg1,ASVneg2,ASVneg3)#ASVneg3 better


#TEST EACH TAXA INDEPENDENTLY...............................................................................................
library(lme4)
AspergillusLM<-glmer(AspergillusB$Abundance ~ AspergillusB$TSFdays + (1|Subplot)+(1|TSFdays), data = AspergillusB)
Anova(AspergillusLM)#0.954

BalsamiaLM<-glmer(BalsamiaB$Abundance ~ BalsamiaB$TSFdaysk + (1|Subplot)+(1|TSFdays), data = BalsamiaB)
Anova(BalsamiaLM)#0.2883

CatenuliferaLM<-glmer(CatenuliferaB$Abundance ~ CatenuliferaB$TSFdays + (1|Subplot)+(1|TSFdays), data = CatenuliferaB)
Anova(CatenuliferaLM)#0.963

ConiochaetaLM<-glmer(ConiochaetaB$Abundance ~ ConiochaetaB$TSFdays + (1|Subplot)+(1|TSFdays), data = ConiochaetaB)
Anova(ConiochaetaLM)#0.7219

CortinariusLM<-glmer(CortinariusB$Abundance ~ CortinariusB$TSFdays + (1|Subplot)+(1|TSFdays), data = CortinariusB)
Anova(CortinariusLM)#0.2104

CoprinellusLM<-glmer(CoprinellusB$Abundance ~ CoprinellusB$TSFdays + (1|Subplot)+(1|TSFdays), data = CoprinellusB)
Anova(CoprinellusLM)#0.9671 

GeminibasidiumLM<-glmer(GeminibasidiumB$Abundance ~ GeminibasidiumB$TSFdays + (1|Subplot)+(1|TSFdays), data = GeminibasidiumB)
Anova(GeminibasidiumLM)#0.5636

InocybeLM<-glmer(InocybeB$Abundance ~ InocybeB$TSFdays + (1|Subplot)+(1|TSFdays), data = InocybeB)
Anova(InocybeLM)#0.7269

MycenaLM<-glmer(MycenaUn$Abundance ~ MycenaUn$TSFdays + (1|Subplot)+(1|TSFdays), data = MycenaUn)
Anova(MycenaLM)#0.4622 

PenicilliumLM<-glmer(PenicilliumB$Abundance ~ PenicilliumB$TSFdays + (1|Subplot)+(1|TSFdays), data = PenicilliumB)
Anova(PenicilliumLM)#1.116e-06 ***

PyronemaLM<-glm(PyronemaB$Abundance ~ PyronemaB$TSFdays + (1/Plot)+ (1/Subplot)+(1/TSFdays), data = PyronemaB)
Anova(PyronemaLM)#0.9684

RasamsoniaLM<-glmer(RasamsoniaB$Abundance ~ RasamsoniaB$TSFdays + (1|Subplot)+(1|TSFdays), data = RasamsoniaB)
Anova(RasamsoniaLM)#0.8238

TephrocybeLM<-glmer(TephrocybeB$Abundance ~ TephrocybeB$TSFdays + (1|Subplot)+(1|TSFdays), data = TephrocybeB)
Anova(TephrocybeLM)#0.8377  
 
TomentellaLM<-glmer(TomentellaB$Abundance ~ TomentellaB$TSFdays + (1|Subplot)+(1|TSFdays), data = TomentellaB)
Anova(TomentellaLM)#0.2064


#----
#-----
#----------
#********************************************************************************************************************************----
# EXPORT RESULTS --------------------------------------------------------------------------------------------------------------------
#********************************************************************************************************************************----
dir.create(file.path("1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFBurned/"), recursive = TRUE)

#Generalized Regression ....................................................................................................................
capture.output(Anova(AspergillusLM), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFBurned/Aspergillus-LmBurn.csv")
capture.output(Anova(BalsamiaLM), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFBurned/Balsamia-LmBurn.csv")
capture.output(Anova(CatenuliferaLM), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFBurned/Catenulifera-LmBurn.csv")
capture.output(Anova(ConiochaetaLM), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFBurned/Coniochaeta-LmBurn.csv")
capture.output(Anova(CortinariusLM), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFBurned/Cortinarius-LmBurn.csv")
capture.output(Anova(CoprinellusLM), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFBurned/Coprinellus-LmBurn.csv")
capture.output(Anova(GeminibasidiumLM), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFBurned/Geminibasidium-LmBurn.csv")
capture.output(Anova(InocybeLM), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFBurned/Inocybe-LmBurn.csv")
capture.output(Anova(MycenaLM), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFBurned/Mycena-LmBurn.csv")
capture.output(Anova(PenicilliumLM), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFBurned/Penicillium-LmBurn.csv")
capture.output(Anova(PyronemaLM), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFBurned/Pyronema-LmBurn.csv")
capture.output(Anova(RasamsoniaLM), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFBurned/Rasamsonia-LmBurn.csv")
capture.output(Anova(TephrocybeLM), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFBurned/Tephrocybe-LmBurn.csv")
capture.output(Anova(TomentellaLM), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFBurned/Tomentella-LmBurn.csv")


detach(RelGenSoilB1)

#----
#-----
#************************************************************************************************************************----
# -- CHANGES OF TAXA WITH TIME SINCE FIRE UNBURNED DATA ----------------------------------------------------------
##***********************************************************************************************************************----
#
attach(RelGenSoilUn1)

BalsamiaUn<-RelGenSoilUn1[which(RelGenSoilUn1$Genus == "Balsamia"),]
CenococcumUn<-RelGenSoilUn1[which(RelGenSoilUn1$Genus == "Cenococcum"),]
CortinariusUn<-RelGenSoilUn1[which(RelGenSoilUn1$Genus == "Cortinarius"),]
CladophialophoraUn<-RelGenSoilUn1[which(RelGenSoilUn1$Genus == "Cladophialophora"),]
GeoporaUn<-RelGenSoilUn1[which(RelGenSoilUn1$Genus == "Geopora"),]
HyaloscyphaUn<-RelGenSoilUn1[which(RelGenSoilUn1$Genus == "Hyaloscypha"),]
HygrocybeUn<-RelGenSoilUn1[which(RelGenSoilUn1$Genus == "Hygrocybe"),]
InocybeUn<-RelGenSoilUn1[which(RelGenSoilUn1$Genus == "Inocybe"),]
LachnumUn<-RelGenSoilUn1[which(RelGenSoilUn1$Genus == "Lachnum"),]
MycenaUn<-RelGenSoilUn1[which(RelGenSoilUn1$Genus == "Mycena"),]
OchroconisUn<-RelGenSoilUn1[which(RelGenSoilUn1$Genus == "Ochroconis"),]
PoculumUn<-RelGenSoilUn1[which(RelGenSoilUn1$Genus == "Poculum"),]
ThelephoraUn<-RelGenSoilUn1[which(RelGenSoilUn1$Genus == "Thelephora"),]
TomentellaUn<-RelGenSoilUn1[which(RelGenSoilUn1$Genus == "Tomentella"),]
TrechisporaUn<-RelGenSoilUn1[which(RelGenSoilUn1$Genus == "Trechispora"),]
VenturiaUn<-RelGenSoilUn1[which(RelGenSoilUn1$Genus == "Venturia"),]


#Explore Temporal correlations to ensure that we are capturing as much of the variance............... 
ASVneg1 <- glmer(BalsamiaB$Abundance ~ BalsamiaB$TSFdays + (1|Subplot)+(1|Subplot)+(1|TSFdays), data = BalsamiaB);Anova(ASVneg1)
ASVneg2 <- glmer(BalsamiaB$Abundance ~ BalsamiaB$TSFdays +  (1|Subplot)+(1|TSFdays), data = BalsamiaB);Anova(ASVneg2)
ASVneg3 <- glmer(BalsamiaB$Abundance ~ BalsamiaB$TSFdays + (1|Subplot)+(1|TSFdays), data = BalsamiaB);Anova(ASVneg3)

AICc(ASVneg1,ASVneg2,ASVneg3)#ASVneg3 better



#TEST SIGNIFICANCE ....................................................................................................
BalsamiaLMU<-glmer(BalsamiaUn$Abundance ~ BalsamiaUn$TSFdays + (1|Subplot)+(1|TSFdays), data = BalsamiaUn)
Anova(BalsamiaLMU)#0.5278 

CenococcumLMU<-glmer(CenococcumUn$Abundance ~ CenococcumUn$TSFdays + (1|Subplot)+(1|TSFdays), data = CenococcumUn)
Anova(CenococcumLMU)#0.06747 

CortinariusLMU<-glmer(CortinariusUn$Abundance ~ CortinariusUn$TSFdays + (1|Subplot)+(1|TSFdays), data = CortinariusUn)
Anova(CortinariusLMU)# 0.2574

CladophialophoraLMU<-glmer(CladophialophoraUn$Abundance ~ CladophialophoraUn$TSFdays + (1|Subplot)+(1|TSFdays), data = CladophialophoraUn)
Anova(CladophialophoraLMU)#3.971e-05 ***

GeoporaLMU<-glmer(GeoporaUn$Abundance ~ GeoporaUn$TSFdays + (1|Subplot)+(1|TSFdays), data = GeoporaUn)
Anova(GeoporaLMU)#.6861

HyaloscyphaLMU<-glmer(HyaloscyphaUn$Abundance ~ HyaloscyphaUn$TSFdays + (1|Subplot)+(1|TSFdays), data = HyaloscyphaUn)
Anova(HyaloscyphaLMU)#0.7032

HygrocybeLMU<-glmer(HygrocybeUn$Abundance ~ HygrocybeUn$TSFdays + (1|Subplot)+(1|TSFdays), data = HygrocybeUn)
Anova(HygrocybeLMU)#0.2359

InocybeLMU<-glmer(InocybeUn$Abundance ~ InocybeUn$TSFdays + (1|Subplot)+(1|TSFdays), data = InocybeUn)
Anova(InocybeLMU)#0.001826 **

LachnumLMU<-glmer(LachnumUn$Abundance ~ LachnumUn$TSFdays + (1|Subplot)+(1|TSFdays), data = LachnumUn)
Anova(LachnumLMU)#0.7587

MycenaLMU<-glmer(MycenaUn$Abundance ~ MycenaUn$TSFdays + (1|Subplot)+(1|TSFdays), data = MycenaUn)
Anova(MycenaLMU)#0.02746 *

OchroconisLMU<-glmer(OchroconisUn$Abundance ~ OchroconisUn$TSFdays + (1|Subplot)+(1|TSFdays), data = OchroconisUn)
Anova(OchroconisLMU)#0.9873

PoculumLMU<-glmer(PoculumUn$Abundance ~ PoculumUn$TSFdays + (1|Subplot)+(1|TSFdays), data = PoculumUn)
Anova(PoculumLMU)#0.3762

ThelephoraLMU<-glmer(ThelephoraUn$Abundance ~ ThelephoraUn$TSFdays + (1|Subplot)+(1|TSFdays), data = ThelephoraUn)
Anova(ThelephoraLMU)#0.6737 

TomentellaLMU<-glmer(TomentellaUn$Abundance ~ TomentellaUn$TSFdays + (1|Subplot)+(1|TSFdays), data = TomentellaUn)
Anova(TomentellaLMU)#0..2084 

TrechisporaLMU<-glmer(TrechisporaUn$Abundance ~ TrechisporaUn$TSFdays + (1|Subplot)+(1|TSFdays), data = TrechisporaUn)
Anova(TrechisporaLMU)# 0.917  

VenturiaLMU<-glmer(VenturiaUn$Abundance ~ VenturiaUn$TSFdays + (1|Subplot)+(1|TSFdays), data = VenturiaUn)
Anova(VenturiaLMU)#0.5351




#********************************************************************************************************************************----
# EXPORT RESULTS --------------------------------------------------------------------------------------------------------------------
#********************************************************************************************************************************----
dir.create(file.path("1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFUn/"), recursive = TRUE)

#Generalized Regression ....................................................................................................................
capture.output(Anova(BalsamiaLMU), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFUn/Balsamia-LmUn.csv")
capture.output(Anova(CenococcumLMU), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFUn/Cenoccocum-LmUn.csv")
capture.output(Anova(CortinariusLMU), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFUn/Cortinarius-LmUn.csv")
capture.output(Anova(CladophialophoraLMU), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFUn/Cladophialophora-LmUn.csv")
capture.output(Anova(GeoporaLMU), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFUn/Geopora-LmUn.csv")
capture.output(Anova(HyaloscyphaLMU), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFUn/Hyaloscypha-LmUn.csv")

capture.output(Anova(HygrocybeLMU), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFUn/Hygrocybe-LmUn.csv")
capture.output(Anova(InocybeLMU), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFUn/Inocybe-LmUn.csv")
capture.output(Anova(LachnumLMU), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFUn/Lachnum-LmUn.csv")
capture.output(Anova(MycenaLMU), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFUn/Mycena-LmUn.csv")

capture.output(Anova(OchroconisLMU), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFUn/Ochroconis-LmUn.csv")
capture.output(Anova(PoculumLMU), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFUn/Poculum-LmUn.csv")
capture.output(Anova(ThelephoraLMU), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFUn/Thelephora-LmUn.csv")
capture.output(Anova(TomentellaLMU), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFUn/Tomentella-LmUn.csv")
capture.output(Anova(TrechisporaLMU), file="1-Analysis/RelativeAbundance/Fungi/Tables/TaxaSignif/TSFUn/Trechispora-LmUn.csv")



