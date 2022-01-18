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


#--
#----
#***************************************************************************************************************----
#-- SUBSET DATA (sample type and treatment) ------------------------------------------------------------------------
#***************************************************************************************************************----sample_data(physeq)$TSFdays<-factor(sample_data(physeq)$TSFdays) # make as factor

# * * Subset data for soil only  .................................................
physeqSoil<-subset_samples(physeq,SampleType=="Soil");physeqSoil#29819x306
sample_names(physeqSoil); sample_data(physeqSoil)$TSFdays

PhySoilB<-subset_samples(physeqSoil,Treatment=="Burned");PhySoilB#29819x208
sample_names(PhySoilB);sample_data(PhySoilB)$TSFdays

PhySoilUn<-subset_samples(physeqSoil,Treatment=="Unburned");PhySoilUn#29819x98
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
dir.create(file.path("1-Analysis/RelativeAbundance/Tables/DeStats/Soil"), recursive=TRUE)

#Export soil samples ......................................................................................................
write.csv(RelGenSoilTrt1, "1-Analysis/RelativeAbundance/Tables/DeStats/Soil/Trt-PerChange-Data-Soil.csv")
write.csv(RelGenSoilB1, "1-Analysis/RelativeAbundance/Tables/DeStats/Soil/TSF-PerChange-Data-Burned-Soil.csv")
write.csv(RelGenSoilUn1, "1-Analysis/RelativeAbundance/Tables/DeStats/Soil/TSF-PerChange-Data-Unburned-Soil.csv")


#Re-Import soil samples ......................................................................................................
RelGenSoilTrt1<-read.csv("1-Analysis/RelativeAbundance/Tables/DeStats/Soil/Trt-PerChange-Data-Soil.csv")
RelGenSoilB1<-read.csv("1-Analysis/RelativeAbundance/Tables/DeStats/Soil/TSF-PerChange-Data-Burned-Soil.csv")
RelGenSoilUn1<-read.csv("1-Analysis/RelativeAbundance/Tables/DeStats/Soil/TSF-PerChange-Data-Unburned-Soil.csv")



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
TSFTrtPerChgSoil<-RelGenSoilTrt1%>%
  filter(!is.na(Genus)) %>%
  group_by(Genus,TSFdays,Treatment) %>%
  summarise(TotalASV=sum(Abundance,na.rm = TRUE),
            mean = mean(Abundance,na.rm = T),mean = round(mean, 6)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>% 
  mutate(percent = round (percent, 2))%>%
  as.data.frame();TSFTrtPerChgSoil


TrtPerChgSoil<-RelGenSoilTrt1%>%
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
dir.create(file.path("1-Analysis/RelativeAbundance/Tables/PerChange/Soil"), recursive=TRUE)

#Export descriptive statistics ......................................................................................................
write.csv(TrtStatSoil, "1-Analysis/RelativeAbundance/Tables/DeStats/Soil/Trt-Soil-DeStats-Results.csv")
write.csv(TsfStatSoilB, "1-Analysis/RelativeAbundance/Tables/DeStats/Soil/TSF-Soil-DeStats-Results-Burned.csv")
write.csv(TsfStatSoilUn , "1-Analysis/RelativeAbundance/Tables/DeStats/Soil/TSF-Soil-DeStats-Results-Unburned.csv")
write.csv(TsfTrtStatSoilB , "1-Analysis/RelativeAbundance/Tables/DeStats/Soil/Trt-TSF-Soil-DeStats-Results-Burned.csv")


#Export soil samples ......................................................................................................
write.csv(TrtPerChgSoil, "1-Analysis/RelativeAbundance/Tables/PerChange/Soil/Trt-PerChange-Soil.csv")
write.csv(TSFTrtPerChgSoil, "1-Analysis/RelativeAbundance/Tables/PerChange/Soil/TRT-TSF-PerChange-Soil.csv")
write.csv(TsfPerChgSoilB, "1-Analysis/RelativeAbundance/Tables/PerChange/Soil/TSF-PerChange-Burned-Soil.csv")
write.csv(TsfPerChgSoilUn, "1-Analysis/RelativeAbundance/Tables/PerChange/Soil/TSF-PerChange-Unburned-Soil.csv")




#----
#-----
#***********************************************************************************************************************-----
# CALCULATE SIGNIFICANCE----------------------------------------------------------------------------------------------------
#***********************************************************************************************************************-----
#Re-Import soil samples ......................................................................................................
RelGenSoilTrt<-read.csv("1-Analysis/RelativeAbundance/Tables/DeStats/Soil/Trt-PerChange-Data-Soil.csv")
RelGenSoilB1<-read.csv("1-Analysis/RelativeAbundance/Tables/DeStats/Soil/TSF-PerChange-Data-Burned-Soil.csv")
RelGenSoilUn1<-read.csv("1-Analysis/RelativeAbundance/Tables/DeStats/Soil/TSF-PerChange-Data-Unburned-Soil.csv")





#SUBSET THE DATA FOR TREATMENT ........................................................
attach(RelGenSoilTrt)

Bryobacter<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "Bryobacter"),]; head(Bryobacter$Genus)
Bryobacter<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "Bryobacter"),]; head(Bryobacter$Genus)
Udaeobacter<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "Candidatus Udaeobacter"),]; head(Udaeobacter$Genus)
RB41<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "RB41"),]; head(RB41$Genus)
DomiBryobacter<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "DomiBryobacter"),]; head(DomiBryobacter$Genus)
Gemmatimonas<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "Gemmatimonas"),]; head(Gemmatimonas$Genus)
Massilia<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "Massilia"),]; head(Massilia$Genus)
Noviherbaspirillum<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "Noviherbaspirillum"),]; head(Noviherbaspirillum$Genus)
PaeniBryobacter<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "PaeniBryobacter"),]; head(PaeniBryobacter$Genus)
RB41<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "RB41"),]; head(RB41$Genus)
Sphingomonas<-RelGenSoilTrt[which(RelGenSoilTrt$Genus == "Sphingomonas"),]; head(Sphingomonas$Genus)

#CHECK DATA FOR NORMALITY..........................................................................
par(mfrow=c(4,4))

#Check histogram and run shapiro test ...........................................
hist(Bryobacter$Abundance); shapiro.test(Bryobacter$Abundance)
hist(Bryobacter$Abundance); shapiro.test(Bryobacter$Abundance)
hist(Udaeobacter$Abundance); shapiro.test(Udaeobacter$Abundance)
hist(RB41$Abundance); shapiro.test(RB41$Abundance)
hist(DomiBryobacter$Abundance); shapiro.test(DomiBryobacter$Abundance)
hist(Gemmatimonas$Abundance); shapiro.test(Gemmatimonas$Abundance)
hist(Massilia$Abundance); shapiro.test(Massilia$Abundance)
hist(Noviherbaspirillum$Abundance); shapiro.test(Noviherbaspirillum$Abundance)
hist(PaeniBryobacter$Abundance); shapiro.test(PaeniBryobacter$Abundance)
hist(RB41$Abundance); shapiro.test(RB41$Abundance)
hist(Sphingomonas$Abundance); shapiro.test(Sphingomonas$Abundance)




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
ASVneg1 <- glmer(Bryobacter$Abundance ~ Bryobacter$Treatment + (1|Plot)+(1|Subplot)+(1|TSFdays), data = Bryobacter);Anova(ASVneg1)
ASVneg2 <- glmer(Bryobacter$Abundance ~ Bryobacter$Treatment +  (1|Plot)+(1|TSFdays), data = Bryobacter);Anova(ASVneg2)
ASVneg3 <- glmer(Bryobacter$Abundance ~ Bryobacter$Treatment + (1|Subplot)+(1|TSFdays), data = Bryobacter);Anova(ASVneg3)


AICc(ASVneg1,ASVneg2,ASVneg3)#ASVneg1 better


#Test for significance
library(car)#anova for glmer
Bryobacterglmer<- glmer(Bryobacter$Abundance ~ Bryobacter$Treatment + (1|Plot) +(1|Subplot)+(1|TSFdays), data = Bryobacter)
Anova(Bryobacterglmer, type=3)#1.627e-06 ***

Bryobacterglmer<- glmer(Bryobacter$Abundance ~ Bryobacter$Treatmen + (1|Plot) + (1|Subplot)+(1|TSFdays), data = Bryobacter)
Anova(Bryobacterglmer, type=3)#2.174e-05 ***

Udaeobacterglmer<- glmer(Udaeobacter$Abundance ~ Udaeobacter$Treatment + (1|Plot) +(1|Subplot)+(1|TSFdays), data = Udaeobacter)
Anova(Udaeobacterglmer, type=3)#0.0003398 ***

RB41glmer<- glmer(RB41$Abundance ~ RB41$Treatment + (1|Plot) +(1|Subplot)+(1|TSFdays), data = RB41)
Anova(RB41glmer, type=3)#1.333e-05 ***

DomiBryobacterglmer<- glmer(DomiBryobacter$Abundance ~ DomiBryobacter$Treatment + (1|Plot) +(1|Subplot)+(1|TSFdays), data = DomiBryobacter)
Anova(DomiBryobacterglmer, type=3)#2.135e-06 ***

Gemmatimonasglmer<- glmer(Gemmatimonas$Abundance ~ Gemmatimonas$Treatment + (1|Plot) +(1|Subplot)+(1|TSFdays), data = Gemmatimonas)
Anova(Gemmatimonasglmer, type=3)#0.0364841 *  

Massiliaglmer<- glmer(Massilia$Abundance ~ Massilia$Treatment + (1|Plot) +(1|Subplot)+(1|TSFdays), data = Massilia)
Anova(Massiliaglmer, type=3)#1.526e-08 ***

Noviherbaspirillumglmer<- glmer(Noviherbaspirillum$Abundance ~ Noviherbaspirillum$Treatment + (1|Plot) +(1|Subplot)+(1|TSFdays), data = Noviherbaspirillum)
Anova(Noviherbaspirillumglmer, type=3)#5.603e-06 ***

PaeniBryobacterglmer<- glmer(PaeniBryobacter$Abundance ~ PaeniBryobacter$Treatment + (1|Plot) +(1|Subplot)+(1|TSFdays), data = PaeniBryobacter)
Anova(PaeniBryobacterglmer, type=3)#3.667e-07 ***

RB41glmer<- glmer(RB41$Abundance ~ RB41$Treatment + (1|Plot) +(1|Subplot)+(1|TSFdays), data = RB41)
Anova(RB41glmer, type=3)# 0.0014454 ** 

Sphingomonasglmer<- glmer(Sphingomonas$Abundance ~ Sphingomonas$Treatment + (1|Plot) +(1|Subplot)+(1|TSFdays), data = Sphingomonas)
Anova(Sphingomonasglmer, type=3)#3.799e-06 ***



# EXPORT RESULTS --------------------------------------------------------------------------------------------------------------------
dir.create(file.path("1-Analysis/RelativeAbundance/Tables/TaxaSignifnif/Soil/Trt"), recursive=TRUE)

#Export ANOVA values ...............................................................................................................
capture.output(Anova(Bryobacterglmer), file="1-Analysis/RelativeAbundance/Tables/TaxaSignifnif/Soil/Trt/Bryobacterglmer.csv")
capture.output(Anova(Bryobacterglmer), file="1-Analysis/RelativeAbundance/Tables/TaxaSignifnif/Soil/Trt/Can.Udaeobacterglmer.csv")
capture.output(Anova(Udaeobacterglmer), file="1-Analysis/RelativeAbundance/Tables/TaxaSignifnif/Soil/Trt/Sphingomonasglmer.csv")
capture.output(Anova(RB41glmer), file="1-Analysis/RelativeAbundance/Tables/TaxaSignifnif/Soil/Trt/RB41glmer.csv")
capture.output(Anova(DomiBryobacterglmer), file="1-Analysis/RelativeAbundance/Tables/TaxaSignifnif/Soil/Trt/DomiBryobacterglmer.csv")
capture.output(Anova(Gemmatimonasglmer), file="1-Analysis/RelativeAbundance/Tables/TaxaSignifnif/Soil/Trt/unculturedglmer.csv")
capture.output(Anova(Massiliaglmer), file="1-Analysis/RelativeAbundance/Tables/TaxaSignifnif/Soil/Trt/Hyaloscyphaglmer.csv")
capture.output(Anova(Noviherbaspirillumglmer), file="1-Analysis/RelativeAbundance/Tables/TaxaSignifnif/Soil/Trt/Gemmatimonasglmer.csv")
capture.output(Anova(PaeniBryobacterglmer), file="1-Analysis/RelativeAbundance/Tables/TaxaSignifnif/Soil/Trt/Noviherbaspirillumglmer.csv")
capture.output(Anova(RB41glmer), file="1-Analysis/RelativeAbundance/Tables/TaxaSignifnif/Soil/Trt/PaeniBryobacterglmer.csv")
capture.output(Anova(Sphingomonasglmer), file="1-Analysis/RelativeAbundance/Tables/TaxaSignifnif/Soil/Trt/RB41glmer.csv")


#-----
#-------
#************************************************************************************************************************----
# -- CHANGES OF TAXA WITH TIME SINCE FIRE -BURNED  ----------------------------------------------------------
##***********************************************************************************************************************----
attach(RelGenSoilB1)

#Subset the taxa of interest...................................................................................
AdhaeribacterB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Adhaeribacter"),]; head(AdhaeribacterB$Genus)
BryobacterB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Bryobacter"),]; head(BryobacterB$Genus)
BlastococcusB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Blastococcus"),]; head(BlastococcusB$Genus)
Can.UdaeobacterB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Candidatus Udaeobacter"),]; head(Can.UdaeobacterB$Genus)
RB41B<-RelGenSoilB1[which(RelGenSoilB1$Genus == "RB41"),]; head(RB41B$Genus)
ConexibacterB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Conexibacter"),]; head(ConexibacterB$Genus)
DomiBryobacterB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "DomiBryobacter"),]; head(DomiBryobacterB$Genus)
GemmatimonasB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Gemmatimonas"),]; head(GemmatimonasB$Genus)
MassiliaB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Massilia"),]; head(MassiliaB$Genus)
MucilaginibacterB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Mucilaginibacter"),]; head(MucilaginibacterB$Genus)
NoviherbaspirillumB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Noviherbaspirillum"),]; head(NoviherbaspirillumB$Genus)
PaeniBryobacterB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "PaeniBryobacter"),]; head(PaeniBryobacterB$Genus)
PedobacterB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Pedobacter"),]; head(PedobacterB$Genus)
RB41B<-RelGenSoilB1[which(RelGenSoilB1$Genus == "RB41"),]; head(RB41B$Genus)
SolirubrobacterB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "Solirubrobacter"),]; head(SolirubrobacterB$Genus)
TumeBryobacterB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "TumeBryobacter"),]; head(TumeBryobacterB$Genus)
unculturedB<-RelGenSoilB1[which(RelGenSoilB1$Genus == "uncultured"),]; head(unculturedB$Genus)


#Explore Temporal correlations to ensure that we are capturing as much of the variance................................................................ 
ASVneg1 <- glmer(AdhaeribacterB$Abundance ~ AdhaeribacterB$TSFdays + (1|Plot)+(1|Subplot)+(1|TSFdays), data = AdhaeribacterB);Anova(ASVneg1, Type=3)
ASVneg2 <- glmer(AdhaeribacterB$Abundance ~ AdhaeribacterB$TSFdays +  (1|Plot)+(1|TSFdays), data = AdhaeribacterB);Anova(ASVneg2, Type=3)
ASVneg3 <- glmer(AdhaeribacterB$Abundance ~ AdhaeribacterB$TSFdays + (1|Subplot)+(1|TSFdays), data = AdhaeribacterB);Anova(ASVneg3, Type=3)

AICc(ASVneg1,ASVneg2,ASVneg3)#ASVneg1 better


#TEST EACH TAXA INDEPENDENTLY...............................................................................................
library(lme4)
AdhaeribacterLM<-glmer(AdhaeribacterB$Abundance ~ AdhaeribacterB$TSFdays + (1|Plot)+ (1|Subplot)+(1|TSFdays), data = AdhaeribacterB)
Anova(AdhaeribacterLM)#0.002596 **

BryobacterLM<-glmer(BryobacterB$Abundance ~ BryobacterB$TSFdays + (1|Plot)+(1|Subplot)+(1|TSFdays), data = BryobacterB)
Anova(BryobacterLM)#0.05579 

BlastococcusLM<-glmer(BlastococcusB$Abundance ~ BlastococcusB$TSFdays + (1|Plot)+(1|Subplot)+(1|TSFdays), data = BlastococcusB)
Anova(BlastococcusLM)#0.963

Can.UdaeobacterLM<-glmer(Can.UdaeobacterB$Abundance ~ Can.UdaeobacterB$TSFdays + (1|Plot)+(1|Subplot)+(1|TSFdays), data = Can.UdaeobacterB)
Anova(Can.UdaeobacterLM)#0.4191

RB41LM<-glmer(RB41B$Abundance ~ RB41B$TSFdays + (1|Plot)+(1|Subplot)+(1|TSFdays), data = RB41B)
Anova(RB41LM)#0.001428 **

ConexibacterLM<-glmer(ConexibacterB$Abundance ~ ConexibacterB$TSFdays + (1|Plot)+(1|Subplot)+(1|TSFdays), data = ConexibacterB)
Anova(ConexibacterLM)#0.5339

DomiBryobacterLM<-glmer(DomiBryobacterB$Abundance ~ DomiBryobacterB$TSFdays + (1|Plot)+(1|Subplot)+(1|TSFdays), data = DomiBryobacterB)
Anova(DomiBryobacterLM)#0.27

GemmatimonasLM<-glmer(GemmatimonasB$Abundance ~ GemmatimonasB$TSFdays + (1|Plot)+(1|Subplot)+(1|TSFdays), data = GemmatimonasB)
Anova(GemmatimonasLM)#4.27e-07 ***

MucilaginibacterLM<-glmer(MucilaginibacterB$Abundance ~ MucilaginibacterB$TSFdays + (1|Plot)+(1|Subplot)+(1|TSFdays), data = MucilaginibacterB)
Anova(MucilaginibacterLM)#0.0005327 ***

NoviherbaspirillumLM<-glmer(NoviherbaspirillumB$Abundance ~ NoviherbaspirillumB$TSFdays + (1|Plot)+(1|Subplot)+(1|TSFdays), data = NoviherbaspirillumB)
Anova(NoviherbaspirillumLM)#0.0004527 ***

PaeniBryobacterLM<-glmer(PaeniBryobacterB$Abundance ~ PaeniBryobacterB$TSFdays + (1|Plot)+(1|Subplot)+(1|TSFdays), data = PaeniBryobacterB)
Anova(PaeniBryobacterLM)#0.03667 *

PedobacterLM<-glmer(PedobacterB$Abundance ~ PedobacterB$TSFdays + (1|Plot)++(1|TSFdays), data = PedobacterB)
Anova(PedobacterLM)#0.1854

RB41LM<-glmer(RB41B$Abundance ~ RB41B$TSFdays + (1|Plot)+(1|Subplot)+(1|TSFdays), data = RB41B)
Anova(RB41LM)#0.9945 
 
SolirubrobacterLM<-glmer(SolirubrobacterB$Abundance ~ SolirubrobacterB$TSFdays + (1|Plot)+(1|Subplot)+(1|TSFdays), data = SolirubrobacterB)
Anova(SolirubrobacterLM)#0.3038

TumeBryobacterLM<-glmer(TumeBryobacterB$Abundance ~ TumeBryobacterB$TSFdays + (1|Plot)+(1|Subplot)+(1|TSFdays), data = TumeBryobacterB)
Anova(TumeBryobacterLM)#0.03372 *



#----
#-----
#----------
#********************************************************************************************************************************----
# EXPORT RESULTS --------------------------------------------------------------------------------------------------------------------
#********************************************************************************************************************************----
dir.create(file.path("1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/"), recursive = TRUE)

#Generalized Regression ....................................................................................................................
capture.output(Anova(AdhaeribacterLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/Adhaeribacter-GlmerBurn.csv")
capture.output(Anova(BryobacterLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/Bryobacter-GlmerBurn.csv")
capture.output(Anova(BlastococcusLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/Blastococcus-GlmerBurn.csv")
capture.output(Anova(Can.UdaeobacterLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/Can.Udaeobacter-GlmerBurn.csv")
capture.output(Anova(RB41LM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/RB41-GlmerBurn.csv")
capture.output(Anova(ConexibacterLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/Conexibacter-GlmerBurn.csv")
capture.output(Anova(DomiBryobacterLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/DomiBryobacter-GlmerBurn.csv")
capture.output(Anova(GemmatimonasLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/Gemmatimonas-GlmerBurn.csv")
capture.output(Anova(MucilaginibacterLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/Mucilaginibacter-GlmerBurn.csv")
capture.output(Anova(NoviherbaspirillumLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/Noviherbaspirillum-GlmerBurn.csv")
capture.output(Anova(PaeniBryobacterLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/PaeniBryobacter-GlmerBurn.csv")
capture.output(Anova(PedobacterLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/Pedobacter-GlmerBurn.csv")
capture.output(Anova(RB41LM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/RB41-GlmerBurn.csv")
capture.output(Anova(SolirubrobacterLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/Solirubrobacter-GlmerBurn.csv")
capture.output(Anova(TumeBryobacterLM), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFBurned/TumeBryobacter-GlmerBurn.csv")
detach(RelGenSoilB1)



#----
#-----
#************************************************************************************************************************----
# -- CHANGES OF TAXA WITH TIME SINCE FIRE UNBURNED DATA ----------------------------------------------------------
##***********************************************************************************************************************----
#
attach(RelGenSoilUn1)

BryobacterUn<-RelGenSoilUn1[which(RelGenSoilUn1$Genus == "Bryobacter"),]
Can.UdaeobacterUn<-RelGenSoilUn1[which(RelGenSoilUn1$Genus == "Candidatus Udaeobacter"),]
RB41Un<-RelGenSoilUn1[which(RelGenSoilUn1$Genus == "RB41"),]
SphingomonasUn<-RelGenSoilUn1[which(RelGenSoilUn1$Genus == "Sphingomonas"),]
unculturedUn<-RelGenSoilUn1[which(RelGenSoilUn1$Genus == "uncultured"),]

#Explore Temporal correlations to ensure that we are capturing as much of the variance............... 
ASVneg1 <- glmer(BryobacterUn$Abundance ~ BryobacterUn$TSFdays + (1|Subplot)+(1|Subplot)+(1|TSFdays), data = BryobacterUn);Anova(ASVneg1)
ASVneg2 <- glmer(BryobacterUn$Abundance ~ BryobacterUn$TSFdays +  (1|Subplot)+(1|TSFdays), data = BryobacterUn);Anova(ASVneg2)
ASVneg3 <- glmer(BryobacterUn$Abundance ~ BryobacterUn$TSFdays + (1|Subplot)+(1|TSFdays), data = BryobacterUn);Anova(ASVneg3)

AICc(ASVneg1,ASVneg2,ASVneg3)#ASVneg1 better



#TEST SIGNIFICANCE ....................................................................................................
BryobacterLMU<-glmer(BryobacterUn$Abundance ~ BryobacterUn$TSFdays + (1|Plot)+(1|Subplot)+(1|TSFdays), data = BryobacterUn)
Anova(BryobacterLMU)#0.41 

Can.UdaeobacterLMU<-glmer(Can.UdaeobacterUn$Abundance ~ Can.UdaeobacterUn$TSFdays + (1|Plot)+(1|Subplot)+(1|TSFdays), data = Can.UdaeobacterUn)
Anova(Can.UdaeobacterLMU)#0.1988 

RB41LMU<-glmer(RB41Un$Abundance ~ RB41Un$TSFdays + (1|Plot)+(1|Subplot)+(1|TSFdays), data = RB41Un)
Anova(RB41LMU)# 0.58056

SphingomonasLMU<-glmer(SphingomonasUn$Abundance ~ SphingomonasUn$TSFdays + (1|Plot)+(1|Subplot)+(1|TSFdays), data = SphingomonasUn)
Anova(SphingomonasLMU)#0.4314

unculturedLMU<-glmer(unculturedUn$Abundance ~ unculturedUn$TSFdays +(1|Plot)+ (1|Subplot)+(1|TSFdays), data = unculturedUn)
Anova(unculturedLMU)#0.4167





#********************************************************************************************************************************----
# EXPORT RESULTS --------------------------------------------------------------------------------------------------------------------
#********************************************************************************************************************************----
dir.create(file.path("1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFUn/"), recursive = TRUE)

#Generalized Regression ....................................................................................................................
capture.output(Anova(BryobacterLMU), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFUn/Bryobacter-LmUn.csv")
capture.output(Anova(Can.UdaeobacterLMU), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFUn/Cenoccocum-LmUn.csv")
capture.output(Anova(RB41LMU), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFUn/RB41-LmUn.csv")
capture.output(Anova(SphingomonasLMU), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFUn/Sphingomonas-LmUn.csv")
capture.output(Anova(unculturedLMU), file="1-Analysis/RelativeAbundance/Tables/TaxaSignif/TSFUn/uncultured-LmUn.csv")



