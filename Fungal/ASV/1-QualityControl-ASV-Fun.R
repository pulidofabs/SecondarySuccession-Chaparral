#Reset R's Brain
#8/27/2021

rm(list=ls())

#Set working directory....................................................................
setwd("C:/Users/fabipc/Dropbox/3-HolyFire/1-HolyFire-Year1/2-Fungi")

#Load packages...................................................................
library(dplyr)
library(EcolUtils) 
library(SPECIES)
library(BiodiversityR)
library(scales)
library(plyr)


#Export Dirty Metadata..............................................................................
metadata<-read.csv("Qiime/Metadata/Fungal-Metadata-Lib1-4E-MJ-Clean.csv",
                   na.strings = "N/A", header = TRUE, row.names = 1)  
dim(metadata)#379x59

# * * LOAD OTU TABLES (no singletons or contaminated samples).......................................
RawOtu<-read.csv("Qiime/ASVtables/Fungal-table-with-taxonomy.csv", row.names=1,check.names=FALSE) 
dim(RawOtu)#14651x380


#
###
#*********************************************************************************************************************************------
#-------------------CREATE CLEAN ASV TABLE AND ID/TAXA LABEL TABLES  ---------------------------------------------
#*********************************************************************************************************************************------

# * Remove FeatureID and Taxonomy Labels to use later-...........................................................
lastcolumn <- ncol(RawOtu); lastcolumn #--380---What is the last column in the dataset
taxononomy<-RawOtu[,380]; taxononomy #Extract taxonomy column from the dataset
featureID <- row.names(RawOtu); featureID #Extract feature ID, from the dataset
TaxonID <- data.frame(featureID,taxononomy); head(TaxonID) #Create dataframe of featureID and taxonomy

# * Clean Taxonomy labels by separating them by (;)-..................................................................................
   #divide last row into subsets #separate taxonomy w spaces 
SplitTaxonomy<-ldply(str_split(string = TaxonID$taxononomy, pattern=";"), rbind)#Divide columns using ";" & convert  to data frame
names(SplitTaxonomy)<-c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") #rename the column names
SplitTaxonomy2 <- as.data.frame(lapply(SplitTaxonomy, gsub, pattern=" ", replacement=""));SplitTaxonomy2
tail(SplitTaxonomy2);head(SplitTaxonomy2)

# * Change the format of the TaxonID dataframe created above so that is is split
head(TaxonID)
TaxonID <- cbind(TaxonID[,1:2 ], SplitTaxonomy2);head(TaxonID)
rownames(TaxonID) <- featureID;head(TaxonID)#add feaure ID to table


# * * Look at your dataset to make sure that they have the same length................................................................
dim(TaxonID);dim(RawOtu) # should have same row length 

# * ATTACH TAXON ID CHANGES TO RAWOTU TABLE TO CREATE A CLEAN TABLE
RawOtu2<- cbind(RawOtu, TaxonID);head(RawOtu2[1:2,])

# CREATE SV TABLE--NO TAXONOMY .......................................................................................................
RawOtuTable  <- RawOtu2[ ,1:379];head(RawOtuTable)#number as shown above-1

#Transpose raw otu table to create actual OTU table...................................................................................
OtuTrans <- t(RawOtuTable);dim(OtuTrans) #379x14651 


#EXPORT DATA ...................................................................
dir.create(file.path("1-Analysis/ASVtables/Fungi"), recursive = TRUE)
dir.create(file.path("1-Analysis/QualityControl/Fungi"), recursive=TRUE)

write.csv(RawOtuTable, "1-Analysis/ASVtables/Fungi/RawOtuTable.csv")
write.csv(OtuTrans, "1-Analysis/ASVtables/Fungi/OtuTrans.csv")
write.csv(TaxonID, "1-Analysis/QualityControl/Fungi/TaxononomyLabels.csv")


#
#
#************************************************************************************------
# -------------------  QUALITY CONTROL---------------------------------------------
#************************************************************************************------

# * look at Negative DNA controls ----------------------
NegControls <- which(colnames(RawOtuTable) %in% 
               c("F057HolyNE02","FN1HolyNE01","F122HolyNeg5T4",
              "F146HolyNeg6T3","F170HolyNeg7T2","F186HolyNeg8T2",
              "F206HolyNeg9T1","F256HolyNeg11T7","F277HolyNeg12T8",
              "F294HolyNeg13T8","F74HolyNeg3bT5","F98HolyNeg4bT4",
              "FNeg_HolyLib3N","FHolyNE14","FPCRNeg","FPCRNegCtrl1",
              "FNegExtCol1","FNegExtCol2","FNegExtCol3","FNegExtCol4",
              "FNeg_HolyPcrNegLin2-Neg","F237HolyNeg10T7","F074HolyNE03",
              "FPCRNegCtrlFabi","FPCRNegCtrl2")); NegControls

MockComm<-which(colnames(RawOtuTable) %in%
          c("FMC1Holy","FPCRMoc","FPCRMockPositiveCtrl1",
          "FPCRMockPositiveCtrl2","FMC1Holy","FNP1Holy","FMockHolyLib3M",
          "FMockHolyMockLib2-M","FCNF08ETP7","FMockHolyMockLib2-M")); MockComm


# make otu table of Negative DNA controls & Mock Community .........................
OtuNegControl<-RawOtuTable[ ,NegControls]; head(OtuNegControl)
OtuMockComm<-RawOtuTable[,MockComm]; head(OtuMockComm)

SumNegs<-sort(rowSums(OtuNegControl), decreasing = TRUE);SumNegs
SumMocks<-sort(rowSums(OtuMockComm), decreasing = TRUE);SumMocks

#Remove Negcontrols and Moc comunities from the RawOtu Table........................
RawOtuTable<-RawOtuTable[ ,-NegControls]; dim(RawOtuTable)#14671x357
RawOtuTable<-RawOtuTable[ ,-MockComm]; dim(RawOtuTable)#14671X350


#EXPORT tables............................................................................................
write.csv(RawOtuTable, "1-Analysis/ASVtables/Fungi/RawOtuTable-NC.csv")


#ReImport raw Otu Table because it did not remove all moc samples--removed manually......
RawOtuTableNC<-read.csv("1-Analysis/ASVtables/Fungi/RawOtuTable-NC.csv", row.names=1,na.strings = "N/A")
dim(RawOtuTableNC)#14651X344

# TRANSPOSE RAW OTU TABLE.................................................................................
OtuTransNC <- t(RawOtuTableNC);dim(OtuTransNC) #344X14651
write.csv(OtuTrans, "1-Analysis/ASVtables/Fungi/OtuTrans-NC.csv")

#Redo RepseQ to see what you get when you run in R---they match-- do n
TotSeq<-sort(rowSums(OtuTransNC), decreasing = TRUE);TotSeq


#Export file..........................................................................
write.csv(TotSeq, "1-Analysis/QualityControl/Fungi-TotalSequence-Samples.csv")



#.----
#.----
#***********************************************************************************************************************------------
# RAREFY TABLE- FULL FUNGAL--------------------------------------------------------------------------------------------------------
#***********************************************************************************************************************------------
# Chose rarefaction value using qiime feature table, removed bottom 10%
# * * Trial 1 Remove samples w not enough sequences.......................................
OtuTrans2<-OtuTransNC[rowSums(OtuTransNC)>11058,]#Maintain rows w totals> rarefaction
dim(OtuTrans2)#335X14651 ---Loose 7 samples

# * * Rarefy, remove low abund, normalize to 11058 seq/sample, normalize 1x...............
OtuRare1 <- rrarefy(OtuTrans2,  11058)
dim(OtuRare1)#335X14651

# *  * Normalize the data 100x and get the mean...........................................
OtuRare2<- rrarefy.perm(OtuTrans2, sample =11058, n = 250, round.out = T)
dim(OtuRare2) #335X14651

# * * * TEST IF RAREFACTION METHOD MATTERS................................................
mantel(OtuRare1, OtuRare2)#0.9993

# REMOVE COLUMNS WITH ZERO READS (OTURARE2)...............................................
# * Create an object containing zeros..........
zeroes <- which(colSums(OtuRare2)==0)
head(zeroes)

# * Remove the columns containing zero's from OTU table ..................................
#Import OTU Rare3, without the controls---
OtuRare3NC <- OtuRare2[ , -zeroes]
head(OtuRare3NC[,1:2]);dim(OtuRare3NC)#335X12344


# * * EXPORT RAREFACTION TABLES...........................................................
#Tables ........................................................................
write.csv(OtuTrans2, "1-Analysis/ASVtables/Fungi/OtuTrans2-NC.csv")
write.csv(OtuRare2, "1-Analysis/ASVtables/Fungi/OtuRare2-NC.csv")
write.csv(OtuRare1, "1-Analysis/ASVtables/Fungi/OtuRare1-NC.csv")
write.csv(OtuRare3NC, "1-Analysis/ASVtables/Fungi/OtuRare3-NC.csv")#table w/o 0


#-----
#------
#***************************************************************************************************************************----
#RAREFY TABLE TO MATCH THE RAREFIED TABLE L-------------------------------------------------------------------------------------
#***************************************************************************************************************************----
#Load OtuRare3 table and then remove samples to match btwn metadata and OtuRare3
metadata<-read.csv("Qiime/Metadata/Fungal-Metadata-Lib1-4E-MJ-Clean.csv", na.strings = "N/A", header = TRUE)
OtuRare3NC<-read.csv("1-Analysis/ASVTables/Fungi/OtuRare3-NC.csv", check.names = FALSE)

names(OtuRare3NC)[1]<- "SampleID"
names(metadata)[1]<- "SampleID"

#Verify that column1 has the same name in both tables...................................
names(metadata[,1:2]); names(OtuRare3NC[,1:2])
dim(metadata);dim(OtuRare3NC)

MetaRare<-metadata %>% semi_join(OtuRare3NC, by = "SampleID")#KeepRows w matchingIDs....
dim(MetaRare)#335X60


#Export rarefied metadata...............................................................
dir.create(file.path("Metadata/Fungi"), recursive=TRUE)
write.csv(MetaRare, "Metadata/Fungi/MetaRare.csv")


#******************************************************************************************----
# CALCULATE ALPHA DIVERSITY -------------------------------------------------------------------
#******************************************************************************************----

# RICHNESS METHOD 1: BIODIVERSITY CALCULATED..............................................
# * Load metadata that corresponds to R-rarefied table crea.ted above -------
MetaRare<-read.csv("Metadata/Fungi/MetaRare.csv", row.names = 1);dim(MetaRare)#335x60
OtuRare3NC<-read.csv("1-Analysis/ASVtables/Fungi/OtuRare3-NC.csv", row.names = 1,
                     check.names = FALSE);dim(OtuRare3NC)#335 12311

#Verify that rownames match
all(rownames(MetaRare)==rownames(OtuRare3NC))


# * Richness calculations using Biodivesity package......................
OtuRichness <- estimateR(OtuRare3NC)
OtuRichness<- t(estimateR(OtuRare3NC))#transformed table

# Create a function to estimate diversity ...............................
estimates_plot_function <- function(datatrans,name){
  #plot S.chao vs S.obs to see if they are correlated
  estimates2<- t(estimateR(datatrans))
  par(mfrow=c(1,2))
  plot(estimates2[,2] ~estimates2[,1], xlab="S.obs", 
       ylab="chao", col=alpha("red", 0.5),pch=16)
  mtext("A",side=3,adj=0)
  #plot S. Ace vs S.obs to see if they are correlated
  plot(estimates2[,4] ~estimates2[,1], xlab="S.obs",
       ylab="ACE",col=alpha("black", 0.5),pch=16)
  mtext("B",side=3,adj=0)
  #dev.off()
}

# DIVERSITY INDICES -...................................................
library(BiodiversityR)
shanEntro <- diversity(OtuRare3NC)#ShannonEntropy
shannon <- exp(OtuRare3NC)#ShannonDiversity
simpson <- diversity(OtuRare3NC, "simpson")#SimpsonDiversity
simpEven<- diversity(OtuRare3NC, "inv")#SimpsonInverse(Evenness)

# * Dataframe of shannon entropy, diversity &  simpson diversity.........
otu.richness <- data.frame(shanEntro, shannon, simpson, simpEven)

# * Add above diversity metrics to S obs, chao1, ACE
OtuSppRichness <- cbind(otu.richness, OtuRichness)
dim(OtuSppRichness)#335X12319

# * * * Export Species Richness....................................................
dir.create("1-Analysis/Diversity")
dir.create(file.path("1-Analysis/Diversity/Alpha/Fungi"),recursive = TRUE)

write.csv(OtuSppRichness, "1-Analysis/Diversity/Alpha/Fungi/CoreMetricsAll.csv") 


#Subset the species richness ......................................................
ncol(OtuSppRichness)#12319

#1 column and 6th column
CoreMetrics1<-OtuSppRichness[,12313:12319];head(CoreMetrics1)#12319-6
CoreMetrics2<-as.data.frame(OtuSppRichness[,1]);head(CoreMetrics2)
names(CoreMetrics2)[1]<- "ShannonH"

#Bind the data into one dataframe
CoreMetrics<-cbind(CoreMetrics1,CoreMetrics2);head(CoreMetrics)

#Export core metrics results......................................................
write.csv(CoreMetrics,"1-Analysis/Diversity/Alpha/Fungi/CoreMetrics.csv")


#----
#----
#***************************************************************************************************----
#Merge table containing species richness(CoreMetrics) to MetaRare -------------------------------------
#**************************************************************************************************----

#Reload data to make sure we can merge the datasets
CoreMetrics<-read.csv("1-Analysis/Diversity/Alpha/Fungi/CoreMetrics.csv", check.names = FALSE)
MetaRare<-read.csv("Metadata/Fungi/MetaRare.csv");dim(MetaRare)#335x60
names(CoreMetrics)[1]<- "SampleID"

#Add the coremetrics file to the metadata file
MetaRareSpp<-merge(MetaRare,CoreMetrics,by="SampleID")
dim(MetaRareSpp)#335x68

#Export rarefied metadata w core metrics attached---------------------
write.csv(MetaRareSpp, "Metadata/Fungi/MetaRareSpp.csv")



#----
#----
#************************************************************************************************----
# Subset Metadata to create Indv Files for Soil, Rain, etc-------------------------------------------
#************************************************************************************************----
attach(MetaRareSpp)


#SOIL SAMPLES ..................................................
MetaSoil<-MetaRareSpp[which(SampleType== "Soil"), ]
head(MetaSoil[,1:2]);dim(MetaSoil)#302x68

#MASONJAR SAMPLES ...............................................
MetaRain<-MetaRareSpp[which(SampleType== "Rain"), ]
head(MetaRain[,1:2]);dim(MetaRain)#33x68
detach(MetaRareSpp)

# SUBSET BY  SOIL TREATMENT ......................................
attach(MetaSoil)
# * * UNBURNED PLOTS ...................................
Unburned<-MetaSoil[which(Treatment== "Unburned"), ]
head(Unburned[,1:2]);dim(Unburned)#103x68

# * * BURNED PLOTS ...................................
Burned<-MetaSoil[which(Treatment== "Burned"), ]
head(Burned[,1:2]);dim(Burned)#199x68
detach(MetaSoil)

# SUBSET BY  RAIN WATER TREATMENT ...................................
# * * UNBURNED PLOTS MJ ...................................
attach(MetaRain)
UnburnedRain<-MetaRain[which(Treatment== "Unburned"), ]
head(UnburnedRain[,1:2]);dim(UnburnedRain)#12x6

# * * BURNED PLOTS MJ ...................................
BurnedRain<-MetaRain[which(Treatment== "Burned"), ]
head(BurnedRain[,1:2]);dim(BurnedRain)#21x68
detach(MetaRain)


#Export tables................................................................
#Analysis of manuscripts are performed solely on soil samples only 
write.csv(MetaSoil, file = "Metadata/Fungi/MetaSoilSpp.csv")
write.csv(MetaRain, file = "Metadata/Fungi/MetaRainSpp.csv")

write.csv(Burned, file = "Metadata/Fungi/MetaSoil-BurnedSpp.csv")
write.csv(Unburned, file = "Metadata/Fungi/MetaSoil-UnburnedSpp.csv")

write.csv(BurnedRain, file = "Metadata/Fungi/MetaRain-BurnedSpp.csv")
write.csv(UnburnedRain, file = "Metadata/Fungi/MetaRain-UnburnedSpp.csv")








