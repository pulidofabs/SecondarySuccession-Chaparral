#Reset R's Brain
#11-18-2020

rm(list=ls())

#Set working directory....................................................................
setwd("C:/Users/fabipc/Dropbox/3-HolyFire/1-HolyFire-Year1/1-Bacteria/")

#Load packages...................................................................
library(dplyr)
library(EcolUtils) 
library(SPECIES)
library(BiodiversityR)
library(scales)
library(plyr)
library(stringr)


#LOAD DATA ................................................................................
Metadata<-read.csv("Qiime/Metadata/Bacteria-Metadata-Lib1-4E-MJ-Clean.csv", 
           row.names = 1,na.strings = "N/A", header = TRUE);dim(Metadata)#382x59

# * * LOAD OTU TABLES (no singletons, contm'd samples removed..............................
RawOtu<-read.csv("Qiime/DADA2/FilteredTable/Bacteria-table-with-taxonomy.csv",
        row.names = 1, check.names = FALSE);dim(RawOtu)#34238x383


#*********************************************************************************************************************************------
#-------------CREATE CLEAN ASV TABLE AND ID/TAXA LABEL TABLES  ---------------------------------------------
#*********************************************************************************************************************************------

# * Remove FeatureID and Taxonomy Labels to use later-...........................................................
lastcolumn <- ncol(RawOtu); lastcolumn #--383---What is the last column in the dataset, need for removing taxa
taxononomy<-RawOtu[,383]; taxononomy #Extract taxonomy column from the dataset
featureID <- row.names(RawOtu); featureID #Extract feature ID, from the dataset
TaxonID <- data.frame(featureID,taxononomy); head(TaxonID) #Create a dataframe of the featureID and taxonomy

# * Clean Taxonomy labels by separating them by (;)-..................................................................................
#divide last row into subsets #separate taxonomy w spaces 
SplitTaxonomy<-ldply(str_split(string = TaxonID$taxononomy, pattern=";"), rbind) # Divide a column using ";"and convert list to data frame
names(SplitTaxonomy)<-c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") #rename the column names
SplitTaxonomy2 <- as.data.frame(lapply(SplitTaxonomy, gsub, pattern=" ", replacement=""));SplitTaxonomy2
tail(SplitTaxonomy2);head(SplitTaxonomy2)

# * Change the format of the TaxonID dataframe created above so that is is split
head(TaxonID)
TaxonID <- cbind(TaxonID[,1:2 ], SplitTaxonomy2);head(TaxonID)
rownames(TaxonID) <- featureID;head(TaxonID)#add feaure ID to table

# * * Look at your dataset to make sure that they have the same length................................................................
dim(TaxonID);dim(RawOtu) # should have same row length 

# * ATTACH TAXON ID CHANGES TO RAWOTU TABLE TO CREATE A CLEAN TABLE...................................................................
RawOtu2<- cbind(RawOtu, TaxonID);head(RawOtu2[1:2,])

# CREATE SV TABLE--NO TAXONOMY .......................................................................................................
RawOtuTable  <- RawOtu2[ ,1:382];head(RawOtuTable)#number as shown above-1

#Transpose raw otu table to create actual OTU table...................................................................................
OtuTrans <- t(RawOtuTable);dim(OtuTrans) #382x34238

#EXPORT DATA ...................................................................
dir.create(file.path("Output/Abundance/Tables"))
dir.create(file.path("1-Analysis/ASVtables"), recursive = TRUE)
dir.create(file.path("1-Analysis/QualityControl"), recursive=TRUE)

write.csv(RawOtuTable, "1-Analysis/ASVtables/RawOtuTable.csv")
write.csv(OtuTrans, "1-Analysis/ASVtables/OtuTrans.csv")
write.csv(TaxonID, "1-Analysis/QualityControl/TaxononomyLabels.csv")


#----
#----
#*********************************************************************************************************************************------
# --- -----------   QUALITY CONTROL  ---------------------------------------------
#*********************************************************************************************************************************------

# * look at Negative DNA controls .....................................................
NegControls <- which(colnames(RawOtuTable) %in% 
             c("B122HolyNeg5T4","B146HolyNeg6T3","B170HolyNeg7T2", "B186HolyNeg8T2",
             "B206HolyNeg9T1","B237HolyNeg10T7","B256HolyNeg11T7","BHolyNE14",
             "BNeg_HolyPcrNegLin2-Neg", "BNeg_HolyLib3N","B277HolyNeg12T8",
             "B294HolyNeg13T8","B98HolyNeg4bT4", "BNegExtCol4","BNegExtCol1",
             "BPCRNeg","BPCRNegCtrl1","BPCRNegCtrl2", "BPCRNegCtrlFabi","BNegExtCol2",
             "B057HolyNE02","BN1HolyNE01","BPCRNeg"));NegControls

# * look at Mock DNA controls ..........................................................
MockComm<-which(colnames(RawOtuTable) %in% 
           c("BMockHolyMockLib2.M","BMockHolyMockLib2-M","BMockHolyLib3M","BPCRMoc",
             "BPCRMockPositiveCtrl1","BPCRMockPositiveCtrl2","BMC1Holy",
             "BPCRMoc","BNP1Holy")); MockComm



#Remove Negcontrols and Moc comunities from the RawOtu Table.............................................
OtuNegControl<-RawOtuTable[, NegControls]; head(OtuNegControl)
RawOtuTable<-RawOtuTable[,-NegControls]; dim(RawOtuTable)#34238x360


# make otu table of Negative DNA controls & Mock Community...............................................
OtuMockComm<-RawOtuTable[, MockComm]; head(OtuMockComm)
RawOtuTable<-RawOtuTable[,-MockComm]; dim(RawOtuTable)#34238x353


SumNegs<-sort(rowSums(OtuNegControl), decreasing = TRUE);SumNegs
SumMocks<-sort(rowSums(OtuMockComm), decreasing = TRUE);SumMocks



#EXPORT .................................................................................................
#For some reason it is not removing everything so removed the manually and reimported
write.csv(RawOtuTable, "1-Analysis/ASVtables/RawOtuTable-NC.csv")


#ReImport raw Otu Table because it did not remove all moc samples--removed manually......................
RawOtuTableNC<-read.csv("1-Analysis/ASVtables/RawOtuTable-NC.csv", row.names=1,na.strings = "N/A")
dim(RawOtuTableNC)#34238x347


# TRANSPOSE RAW OTU TABLE---------------------------------------------------------------------------------------------------------
OtuTransNC <- t(RawOtuTableNC);dim(OtuTransNC) #347 x 34238
write.csv(OtuTrans, "1-Analysis/ASVtables/OtuTrans-NC.csv")

#Redo RepseQ to see what you get when you run in R---they match-- do n
TotSeq<-sort(rowSums(OtuTransNC), decreasing = TRUE);TotSeq
write.csv(TotSeq, "1-Analysis/QualityControl/TotalSequence-Sample.csv")
#Rarefaction will be 7115 loosing 6 samples



#.----
#.----
#***********************************************************************************************************************------------
# RAREFY TABLE- FULL Bacteria--------------------------------------------------------------------------------------------------------
#***********************************************************************************************************************------------
# * Chose rarefaction value using qiime feature table, removed botto 10%

OtuTrans2<-OtuTransNC[rowSums(OtuTransNC)>7115,] #Subset data to maintain only rows with totals above rarefaction depth
dim(OtuTrans2)#340x34238 ---

# * * Rarefy table to remove low abundance, normalize to 7115 seq/sample, normalize 1x...............................................
OtuRare1 <- rrarefy(OtuTrans2, 7115)
dim(OtuRare1)#340x34238

# *  * Normalize the data 100x and get the mean......................................................................................
OtuRare2<- rrarefy.perm(OtuTrans2, sample = 7115, n = 250, round.out = T)
dim(OtuRare2)#340x34238

# * * * TEST IF RAREFACTION METHOD MATTERS...........................................................................................
mantel(OtuRare1, OtuRare2)#0.9985

# REMOVE COLUMNS WITH ZERO READS (OTURARE2)..........................................................................................
# * Create an object containing zeros -------------------
zeroes <- which(colSums(OtuRare2)==0)
head(zeroes)

# * Remove the columns containing zero's from OTU table..............................................................................
OtuRare3NC <- OtuRare2[ , -zeroes]
head(OtuRare3NC[,1:2]);dim(OtuRare3NC)#340 x 31858


# * * EXPORT RAREFACTION TABLES ......................................................................................................
#Tables ...............................................................................
write.csv(OtuTrans2, "1-Analysis/ASVtables/OtuTrans2-NC.csv")
write.csv(OtuRare2, "1-Analysis/ASVtables/OtuRare2-NC.csv")
write.csv(OtuRare1, "1-Analysis/ASVtables/OtuRare1-NC.csv")
write.csv(OtuRare3NC, "1-Analysis/ASVtables/OtuRare3-NC.csv")#table w/o zeroes



#
#
#***************************************************************************************************************************************----
#RAREFY TABLE TO MATCH THE RAREFIED TABLE L-------------------------------------------------------------------------------------------------
#***************************************************************************************************************************************----
#Load OtuRare3 table and then remove samples to match btwn metadata and OtuRare3
metadata<-read.csv("Qiime/Metadata/Bacteria-Metadata-Lib1-4E-MJ-Clean.csv",na.strings = "N/A", header = TRUE);dim(metadata)#382X59
OtuRare3<-read.csv("1-Analysis/ASVtables/OtuRare3-NC.csv",check.names = FALSE,na.strings = "N/A");dim(OtuRare3)#343x33107

names(OtuRare3)[1]<- "SampleID"
names(metadata)[1]<- "SampleID"

#Verify that column1 has the same name in both tables----------------------------
names(metadata[,1:2]); names(OtuRare3[,1:2])
dim(metadata);dim(OtuRare3)

# keep rows with matching ID---------------------------------------------------
MetaRare<-metadata %>% semi_join(OtuRare3, by = "SampleID") 
dim(MetaRare)#340x60


#EXPORT RAREFIED TABLES ---------------------------------------
dir.create(file.path("1-Analysis/Metadata"))
write.csv(MetaRare, "1-Analysis/Metadata/MetaRare.csv")




#***************************************************************************************************************************************----
# CALCULATE ALPHA DIVERSITY -------------------------------------------------------------------------------------------------
#***************************************************************************************************************************************----
# RICHNESS METHOD 1: BIODIVERSITY CALCULATED------------------------------------------------------------------------------
# * Load metadata that correspands to R-rarefied table created above -------
#row names were verified to match in excel
library(BiodiversityR)

MetaRare<-read.csv("1-Analysis/Metadata/MetaRare.csv", na.strings = "N/A", header = TRUE, row.names = 1);dim(MetaRare)#338x67
OtuRare3NC<-read.csv("1-Analysis/ASVtables/OtuRare3-NC.csv", check.names = FALSE, row.names = 1);dim(OtuRare3NC)#338x 34405

#Verify names match.............................................................................................
names(OtuRare3)[1]==names(metadata)[1]

names(metadata)[1]<- "SampleID"
names(metadata)[1]<- "SampleID"

# * Richness calculations using Biodivesity package.............................................................
OtuRichness <- estimateR(OtuRare3NC)
OtuRichness<- t(estimateR(OtuRare3NC))#run function on transformed table to switch rows and columns


# CREATE FUNCTION TO ESTIMATE DIVERSITY METRICS................................................................
estimates_plot_function <- function(datatrans,name){
  #plot S.chao vs S.obs to see if they are correlated
  estimates2<- t(estimateR(datatrans))
  # pdf(paste("figures/richnesscores_",name,".pdf"))
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

# DIVERSITY INDICES ..........................................................................................
shanEntro <- diversity(OtuRare3NC) # Shannon entropy
shannon <- exp(OtuRare3NC) ## Shannon number of diversity
simpson <- diversity(OtuRare3NC, "simpson")#Simpson diversity
simpEven<- diversity(OtuRare3NC, "inv") ## Simpson Inverse (Eveness)

# * Dataframe of shannon entropy, diversity &  simpson diversity................
otu.richness <- data.frame(shanEntro, shannon, simpson, simpEven)

# * Add above diversity metrics to S obs, chao1, ACE............................
OtuSppRichness <- cbind(otu.richness, OtuRichness)
dim(OtuSppRichness)#340x31854


# * * * Export Biodiversity package spp richness-............................................................
dir.create(file.path("1-Analysis/Diversity/Alpha"), recursive = TRUE)
write.csv(OtuSppRichness, "1-Analysis/Diversity/Alpha/CoreMetricsAll.csv") #all above species calculations

#Subset the richness metrics..................................................................................
ncol(OtuSppRichness)#31854

#1 column and 6th column........................................................
CoreMetrics1<-OtuSppRichness[,31848:31854];head(CoreMetrics1)#31854-6
CoreMetrics2<-as.data.frame(OtuSppRichness[,1]);head(CoreMetrics2)
names(CoreMetrics2)[1]<- "ShannonH"


CoreMetrics<-cbind(CoreMetrics1,CoreMetrics2);head(CoreMetrics)

#Export core metrics results....................................................
write.csv(CoreMetrics, "1-Analysis/Diversity/Alpha/CoreMetrics.csv")



#----
#----
#*******************************************************************************************************----
#Merge table containing species richness(CoreMetrics) to MetaRare ------------------------------------------
#*******************************************************************************************************----

#Reload data to make sure we can merge the dataset..............................................
CoreMetrics<-read.csv("1-Analysis/Diversity/Alpha/CoreMetrics.csv", check.names = FALSE)
MetaRare<-read.csv("1-Analysis/Metadata/MetaRare.csv");dim(MetaRare)#340x60
names(CoreMetrics)[1]<- "SampleID"

#Confirm that names match.......................................................................
names(CoreMetrics)[1]==names(MetaRare)[1]

#Add the coremetrics file to the metadata file..................................................
MetaRareSpp<-merge(MetaRare,CoreMetrics,by="SampleID")
dim(MetaRareSpp)#340x68

#Export rarefied metadata w core metrics attached...............................................
write.csv(MetaRareSpp, "1-Analysis/Metadata/MetaRareSpp.csv")



#----
#----
#****************************************************************************************************************----
# Subset Metadata to create Indv Files for Soil, Rain, etc---------------------------------------------------
#****************************************************************************************************************----
attach(MetaRareSpp)


#SOIL SAMPLES ..................................................
MetaSoil<-MetaRareSpp[which(SampleType== "Soil"), ]
head(MetaSoil[,1:2]);dim(MetaSoil)#307x68

#MASONJAR SAMPLES ...............................................
MetaRain<-MetaRareSpp[which(SampleType== "Rain"), ]
head(MetaRain[,1:2]);dim(MetaRain)#33x68
detach(MetaRareSpp)

# SUBSET BY  SOIL TREATMENT ......................................
attach(MetaSoil)
# * * UNBURNED PLOTS ...................................
Unburned<-MetaSoil[which(Treatment== "Unburned"), ]
head(Unburned[,1:2]);dim(Unburned)#99x68

# * * BURNED PLOTS ...................................
Burned<-MetaSoil[which(Treatment== "Burned"), ]
head(Burned[,1:2]);dim(Burned)#208x68
detach(MetaSoil)

# SUBSET BY  RAIN WATER TREATMENT ...................................
# * * UNBURNED PLOTS MJ ...................................
attach(MetaRain) summary(Metadata)
UnburnedRain<-MetaRain[which(Treatment== "Unburned"), ]
head(UnburnedRain[,1:2]);dim(UnburnedRain)#12x68

# * * BURNED PLOTS MJ ...................................
BurnedRain<-MetaRain[which(Treatment== "Burned"), ]
head(BurnedRain[,1:2]);dim(BurnedRain)#21x68
detach(MetaRain)


#Export tables...................................................................
write.csv(MetaSoil, file = "1-Analysis/Metadata/MetaSoilSpp.csv")
write.csv(MetaRain, file = "1-Analysis/Metadata/MetaRainSpp.csv")

write.csv(Burned, file = "1-Analysis/Metadata/MetaSoil-BurnedSpp.csv")
write.csv(Unburned, file = "1-Analysis/Metadata/MetaSoil-UnburnedSpp.csv")

write.csv(BurnedRain, file = "1-Analysis/Metadata/MetaRain-BurnedSpp.csv")
write.csv(UnburnedRain, file = "1-Analysis/Metadata/MetaRain-UnburnedSpp.csv")
