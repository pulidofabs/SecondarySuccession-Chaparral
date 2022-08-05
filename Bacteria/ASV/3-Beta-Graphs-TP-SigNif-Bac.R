#8-25-2021

#Set working directory...........................................................
setwd("~/Dropbox/1-Dissertation/1-Research/HolyFire/1-HolyFire-Pub1/1-Bacteria/")


#Load libary.....................................................................
library(vegan)#required for adonis2
library(tidyverse)# to change column to rownames
library(permute)# to be able to use the strata option in adonis2


#Load data.......................................................................
MetaSoil<- read.csv("1-Analysis/Metadata/MetaSoilSpp.csv", header = TRUE)

OtuRare1Soil<-read.csv("1-Analysis/ASVtables/OtuRare1-NC-Soil.csv")
OtuTrans2soil<-read.csv("1-Analysis/ASVtables/OtuTrans2-NC-Soil.csv")

#For future analysis import this files............................





#Verify rownames match...................................
all(rownames(MetaSoil)== rownames(OtuTrans2soil))
all(rownames(MetaSoil)== rownames(OtuRare1soil))


#************************************************************************************************************************----
#-------------- QUALITY CONTROL -------------------------------------------------------------------------------------------
#************************************************************************************************************************----
#Important,run this code the first time to subset data, afterwhich you can just 
#reimport the files, see below

attach(MetaSoil)

#Verify rownames match..........................................................
all(rownames(MetaSoil)== rownames(OtuTrans2soil));nrow(MetaSoil)


#Subset Metadata by Timepoints..................................................
TP1<-MetaSoil[which(Timepoint== "T1"), ];head(TP1[,1:2]);dim(TP1)#29x69
TP2<-MetaSoil[which(Timepoint== "T2"), ];head(TP1[,1:2]);dim(TP2)#36x69
TP3<-MetaSoil[which(Timepoint== "T3"), ];head(TP1[,1:2]);dim(TP3)#36x69
TP4<-MetaSoil[which(Timepoint== "T4"), ];head(TP1[,1:2]);dim(TP4)#35x69
TP5<-MetaSoil[which(Timepoint== "T5"), ];head(TP1[,1:2]);dim(TP5)#35x69
TP6<-MetaSoil[which(Timepoint== "T6"), ];head(TP1[,1:2]);dim(TP6)#36x69
TP7<-MetaSoil[which(Timepoint== "T7"), ];head(TP1[,1:2]);dim(TP7)#33x69
TP8<-MetaSoil[which(Timepoint== "T8"), ];head(TP1[,1:2]);dim(TP8)#32x69
TP9<-MetaSoil[which(Timepoint== "T9"), ];head(TP1[,1:2]);dim(TP9)#36x69

#Name change of the data if they do not match the rest of the da................
#Rename cells as they don't match the rest......................
OtuRare1Soil$X <- gsub("BCNF06WTP4", "BCNF06WT4",OtuRare1Soil$X)
OtuRare1Soil$X <- gsub("BCNF01WTP2", "BCNF01WT2",OtuRare1Soil$X)


#Subset ASV tables, by timepoint................................................
#columns = 35742-subset...................................................... 
ASVTP1<-OtuRare1Soil[grep("T1", OtuRare1Soil$X), ];dim(ASVTP1)#29
ASVTP2<-OtuRare1Soil[grep("T2", OtuRare1Soil$X), ];dim(ASVTP2)#36
ASVTP3<-OtuRare1Soil[grep("T3", OtuRare1Soil$X), ];dim(ASVTP3)#36
ASVTP4<-OtuRare1Soil[grep("T4", OtuRare1Soil$X), ];dim(ASVTP4)#35
ASVTP5<-OtuRare1Soil[grep("T5", OtuRare1Soil$X), ];dim(ASVTP5)#35
ASVTP6<-OtuRare1Soil[grep("T6", OtuRare1Soil$X), ];dim(ASVTP6)#36
ASVTP7<-OtuRare1Soil[grep("T7", OtuRare1Soil$X), ];dim(ASVTP7)#33
ASVTP8<-OtuRare1Soil[grep("T8", OtuRare1Soil$X), ];dim(ASVTP8)#32
ASVTP9<-OtuRare1Soil[grep("T9", OtuRare1Soil$X), ];dim(ASVTP9)#35


#Make the site ID as row names before exporting.............................
ASVTP1<-ASVTP1 %>% remove_rownames %>% column_to_rownames(var="X")
ASVTP2<-ASVTP2 %>% remove_rownames %>% column_to_rownames(var="X")
ASVTP3<-ASVTP3 %>% remove_rownames %>% column_to_rownames(var="X")
ASVTP4<-ASVTP4 %>% remove_rownames %>% column_to_rownames(var="X")
ASVTP5<-ASVTP5 %>% remove_rownames %>% column_to_rownames(var="X")
ASVTP6<-ASVTP6 %>% remove_rownames %>% column_to_rownames(var="X")
ASVTP7<-ASVTP7 %>% remove_rownames %>% column_to_rownames(var="X")
ASVTP8<-ASVTP8 %>% remove_rownames %>% column_to_rownames(var="X")
ASVTP9<-ASVTP9 %>% remove_rownames %>% column_to_rownames(var="X")




#----
#----

attach(MetaSoil)

#Verify rownames match.......................................
all(rownames(MetaSoil)== rownames(OtuTrans2soil))
nrow(MetaSoil)

#Subset Metadata by Timepoints.............................................
TP1<-MetaSoil[which(Timepoint== "T1"), ];head(TP1[,1:2]);dim(TP1)#29x68
TP2<-MetaSoil[which(Timepoint== "T2"), ];head(TP1[,1:2]);dim(TP2)#36x68
TP3<-MetaSoil[which(Timepoint== "T3"), ];head(TP1[,1:2]);dim(TP3)#36x68
TP4<-MetaSoil[which(Timepoint== "T4"), ];head(TP1[,1:2]);dim(TP4)#35x68
TP5<-MetaSoil[which(Timepoint== "T5"), ];head(TP1[,1:2]);dim(TP5)#35x68
TP6<-MetaSoil[which(Timepoint== "T6"), ];head(TP1[,1:2]);dim(TP6)#36x68
TP7<-MetaSoil[which(Timepoint== "T7"), ];head(TP1[,1:2]);dim(TP7)#33x68
TP8<-MetaSoil[which(Timepoint== "T8"), ];head(TP1[,1:2]);dim(TP8)#32x68
TP9<-MetaSoil[which(Timepoint== "T9"), ];head(TP1[,1:2]);dim(TP9)#35x68


#SUBSET THE TABLES BY TIMEPOINT .............................................
#Rename cells as they don't match the rest......................
OtuRare1Soil$X <- gsub("FCNF06WTP4", "FCNF06WT4",OtuRare1Soil$X)
OtuRare1Soil$X <- gsub("FCNF01WTP2", "FCNF01WT2",OtuRare1Soil$X)
OtuRare1Soil$X <- gsub("FCNF09ETP8", "FCNF09ET8",OtuRare1Soil$X)
OtuRare1Soil$X <- gsub("FCNF09STP8", "FCNF09ST8",OtuRare1Soil$X)
OtuRare1Soil$X <- gsub("FCNF09WTP8", "FCNF09WT8",OtuRare1Soil$X)


#columns = 35742-subset...................................................... 
ASVTP1<-OtuRare1Soil[grep("T1", OtuRare1Soil$X), ];dim(ASVTP1)#22
ASVTP2<-OtuRare1Soil[grep("T2", OtuRare1Soil$X), ];dim(ASVTP2)#31
ASVTP3<-OtuRare1Soil[grep("T3", OtuRare1Soil$X), ];dim(ASVTP3)#36
ASVTP4<-OtuRare1Soil[grep("T4", OtuRare1Soil$X), ];dim(ASVTP4)#35
ASVTP5<-OtuRare1Soil[grep("T5", OtuRare1Soil$X), ];dim(ASVTP5)#36
ASVTP6<-OtuRare1Soil[grep("T6", OtuRare1Soil$X), ];dim(ASVTP6)#36
ASVTP7<-OtuRare1Soil[grep("T7", OtuRare1Soil$X), ];dim(ASVTP7)#35
ASVTP8<-OtuRare1Soil[grep("T8", OtuRare1Soil$X), ];dim(ASVTP8)#36
ASVTP9<-OtuRare1Soil[grep("T9", OtuRare1Soil$X), ];dim(ASVTP9)#35


#Make the site ID as row names before exporting.............................
ASVTP1<-ASVTP1 %>% remove_rownames %>% column_to_rownames(var="X")
ASVTP2<-ASVTP2 %>% remove_rownames %>% column_to_rownames(var="X")
ASVTP3<-ASVTP3 %>% remove_rownames %>% column_to_rownames(var="X")
ASVTP4<-ASVTP4 %>% remove_rownames %>% column_to_rownames(var="X")
ASVTP5<-ASVTP5 %>% remove_rownames %>% column_to_rownames(var="X")
ASVTP6<-ASVTP6 %>% remove_rownames %>% column_to_rownames(var="X")
ASVTP7<-ASVTP7 %>% remove_rownames %>% column_to_rownames(var="X")
ASVTP8<-ASVTP8 %>% remove_rownames %>% column_to_rownames(var="X")
ASVTP9<-ASVTP9 %>% remove_rownames %>% column_to_rownames(var="X")



#Subset unrarefied OTUtrans table for NMDS--------------------------------------
#SUBSET THE TABLES BY TIMEPOINT .............................................
#Rename cells as they don't match the rest......................
OtuTrans2soil$X <- gsub("FCNF06WTP4", "FCNF06WT4",OtuTrans2soil$X)
OtuTrans2soil$X <- gsub("FCNF01WTP2", "FCNF01WT2",OtuTrans2soil$X)
OtuTrans2soil$X <- gsub("FCNF09ETP8", "FCNF09ET8",OtuTrans2soil$X)
OtuTrans2soil$X <- gsub("FCNF09STP8", "FCNF09ST8",OtuTrans2soil$X)
OtuTrans2soil$X <- gsub("FCNF09WTP8", "FCNF09WT8",OtuTrans2soil$X)


#columns = 35742-subset...................................................... 
TransTP1<-OtuTrans2soil[grep("T1", OtuTrans2soil$X), ];dim(TransTP1)#29
TransTP2<-OtuTrans2soil[grep("T2", OtuTrans2soil$X), ];dim(TransTP2)#35
TransTP3<-OtuTrans2soil[grep("T3", OtuTrans2soil$X), ];dim(TransTP3)#36
TransTP4<-OtuTrans2soil[grep("T4", OtuTrans2soil$X), ];dim(TransTP4)#34
TransTP5<-OtuTrans2soil[grep("T5", OtuTrans2soil$X), ];dim(TransTP5)#35
TransTP6<-OtuTrans2soil[grep("T6", OtuTrans2soil$X), ];dim(TransTP6)#36
TransTP7<-OtuTrans2soil[grep("T7", OtuTrans2soil$X), ];dim(TransTP7)#33
TransTP8<-OtuTrans2soil[grep("T8", OtuTrans2soil$X), ];dim(TransTP8)#32
TransTP9<-OtuTrans2soil[grep("T9", OtuTrans2soil$X), ];dim(TransTP9)#35


#Make the site ID as row names before exporting.............................
TransTP1<-TransTP1 %>% remove_rownames %>% column_to_rownames(var="X")
TransTP2<-TransTP2 %>% remove_rownames %>% column_to_rownames(var="X")
TransTP3<-TransTP3 %>% remove_rownames %>% column_to_rownames(var="X")
TransTP4<-TransTP4 %>% remove_rownames %>% column_to_rownames(var="X")
TransTP5<-TransTP5 %>% remove_rownames %>% column_to_rownames(var="X")
TransTP6<-TransTP6 %>% remove_rownames %>% column_to_rownames(var="X")
TransTP7<-TransTP7 %>% remove_rownames %>% column_to_rownames(var="X")
TransTP8<-TransTP8 %>% remove_rownames %>% column_to_rownames(var="X")
TransTP9<-TransTP9 %>% remove_rownames %>% column_to_rownames(var="X")



#Export tables .................................................................
#and reimport with row.name = 1 in order to run distance matrix

write.csv(ASVTP1, "1-Analysis/ASVtables/ASVTP1.csv")
write.csv(ASVTP2, "1-Analysis/ASVtables/ASVTP2.csv")
write.csv(ASVTP3, "1-Analysis/ASVtables/ASVTP3.csv")
write.csv(ASVTP4, "1-Analysis/ASVtables/ASVTP4.csv")
write.csv(ASVTP5, "1-Analysis/ASVtables/ASVTP5.csv")
write.csv(ASVTP6, "1-Analysis/ASVtables/ASVTP6.csv")
write.csv(ASVTP7, "1-Analysis/ASVtables/ASVTP7.csv")
write.csv(ASVTP8, "1-Analysis/ASVtables/ASVTP8.csv")
write.csv(ASVTP9, "1-Analysis/ASVtables/ASVTP9.csv")


#Tranposed, unrarefied table.....................................
write.csv(TransTP1, "1-Analysis/ASVtables/TransTP1.csv")
write.csv(TransTP2, "1-Analysis/ASVtables/TransTP2.csv")
write.csv(TransTP3, "1-Analysis/ASVtables/TransTP3.csv")
write.csv(TransTP4, "1-Analysis/ASVtables/TransTP4.csv")
write.csv(TransTP5, "1-Analysis/ASVtables/TransTP5.csv")
write.csv(TransTP6, "1-Analysis/ASVtables/TransTP6.csv")
write.csv(TransTP7, "1-Analysis/ASVtables/TransTP7.csv")
write.csv(TransTP8, "1-Analysis/ASVtables/TransTP8.csv")
write.csv(TransTP9, "1-Analysis/ASVtables/TransTP9.csv")


#----
#----
#****************************************************************************************----
#-------Reimport ASV tables and run analysys---------------------------------------------
#****************************************************************************************----

#Reimport OtuTrans, w/o names
OtuTrans2soil<-read.csv("1-Analysis/ASVtables/OtuTrans2-NC-Soil.csv",  row.names = 1, check.names = FALSE)

#Reimport oturate, subsetted tables, per timepoint, needed to test significance..........
ASVTP1<-read.csv("1-Analysis/ASVtables/ASVTP1.csv", row.names = 1)
ASVTP2<-read.csv("1-Analysis/ASVtables/ASVTP2.csv", row.names = 1)
ASVTP3<-read.csv("1-Analysis/ASVtables/ASVTP3.csv", row.names = 1)
ASVTP4<-read.csv("1-Analysis/ASVtables/ASVTP4.csv", row.names = 1)
ASVTP5<-read.csv("1-Analysis/ASVtables/ASVTP5.csv", row.names = 1)
ASVTP6<-read.csv("1-Analysis/ASVtables/ASVTP6.csv", row.names = 1)
ASVTP7<-read.csv("1-Analysis/ASVtables/ASVTP7.csv", row.names = 1)
ASVTP8<-read.csv("1-Analysis/ASVtables/ASVTP8.csv", row.names = 1)
ASVTP9<-read.csv("1-Analysis/ASVtables/ASVTP9.csv", row.names = 1)

#Reimport Tranposed, unrarefied table....................................................
#Needed for distance matrix and NMDS (we rarefy withtin dist matrix)

TransTP1<-read.csv("1-Analysis/ASVtables/TransTP1.csv", row.names = 1)
TransTP2<-read.csv("1-Analysis/ASVtables/TransTP2.csv", row.names = 1)
TransTP3<-read.csv("1-Analysis/ASVtables/TransTP3.csv", row.names = 1)
TransTP4<-read.csv("1-Analysis/ASVtables/TransTP4.csv", row.names = 1)
TransTP5<-read.csv("1-Analysis/ASVtables/TransTP5.csv", row.names = 1)
TransTP6<-read.csv("1-Analysis/ASVtables/TransTP6.csv", row.names = 1)
TransTP7<-read.csv("1-Analysis/ASVtables/TransTP7.csv", row.names = 1)
TransTP8<-read.csv("1-Analysis/ASVtables/TransTP8.csv", row.names = 1)
TransTP9<-read.csv("1-Analysis/ASVtables/TransTP9.csv", row.names = 1)





#.----
#.----
#**************************************************************************************************************----
# BETA DIVERSITY (NMDS) ALL SOIL SAMPLES --------------------------------------------------------------------------
#**************************************************************************************************************----
# CALCULATE DISTANCE .....................................................................................


#Soil ALL distance ....................................................................................
DistSoil<-avgdist(OtuTrans2soil,7115, iterations=100, meanfun=median, transf= sqrt, dmethod="bray")
NMDSsoil<-metaMDS(DistSoil, autotransform = FALSE, engine = "monoMDS",k=3, 
                  weakties =TRUE, model="global",maxit = 400, try=80, trymax=100)
print(NMDSsoil)#Stress=0.11



# * * NMDS PER TIMEPOINT-------------------------------------------------------------------------------
#TimePoint1
DisTT1<-avgdist(TransTP1, 7115, iterations=100, meanfun=median, transf= sqrt, dmethod="bray")
NMDS1<-metaMDS(DisTT1, autotransform = FALSE, engine = "monoMDS",k=3, 
               weakties =TRUE, model="global",maxit = 400, try=80, trymax=100)
print(NMDS1)#Stress=0.07


#TimePoint2
DisTT2<-avgdist(TransTP2,7115, iterations=100, meanfun=median, transf= sqrt, dmethod="bray")
NMDS2<-metaMDS(DisTT2, autotransform = FALSE, engine = "monoMDS",k=3, 
               weakties =TRUE, model="global",maxit = 400, try=80, trymax=100)
print(NMDS2)#Stress=0.06


#TimePoint3
DisTT3<-avgdist(TransTP3, 7115, iterations=100, meanfun=median, transf= sqrt, dmethod="bray")
NMDS3<-metaMDS(DisTT3, autotransform = FALSE, engine = "monoMDS",k=3, 
               weakties =TRUE, model="global",maxit = 400, try=80, trymax=100)
print(NMDS3)#Stress=0.06


#TimePoint4
DisTT4<-avgdist(TransTP4, 7115, iterations=100, meanfun=median, transf= sqrt, dmethod="bray")
NMDST4<-metaMDS(DisTT4, autotransform = FALSE, engine = "monoMDS",k=3, 
                weakties =TRUE, model="global",maxit = 400, try=80, trymax=100)
print(NMDST4)#Stress 0.07


#TimePoint5
DisTT5<-avgdist(TransTP5, 7115, iterations=100, meanfun=median, transf= sqrt, dmethod="bray")
NMDST5<-metaMDS(DisTT5, autotransform = FALSE, engine = "monoMDS",k=3, 
                weakties =TRUE, model="global",maxit = 400, try=80, trymax=100)
print(NMDST5)#Stress=0.07



#TimePoint6
DisTT6<-avgdist(TransTP6, 7115, iterations=100, meanfun=median, transf= sqrt, dmethod="bray")
NMDST6<-metaMDS(DisTT6, autotransform = FALSE, engine = "monoMDS",k=3, 
                weakties =TRUE, model="global",maxit = 400, try=80, trymax=100)
print(NMDST6)#Stress=0.07



#TimePoint7
DisTT7<-avgdist(TransTP7, 7115, iterations=100, meanfun=median, transf= sqrt, dmethod="bray")
NMDST7<-metaMDS(DisTT7, autotransform = FALSE, engine = "monoMDS",k=3, 
                weakties =TRUE, model="global",maxit = 400, try=80, trymax=100)
print(NMDST7)#Stress=0.06



#TimePoint8
DisTT8<-avgdist(TransTP8, 7115, iterations=100, meanfun=median, transf= sqrt, dmethod="bray")
NMDST8<-metaMDS(DisTT8, autotransform = FALSE, engine = "monoMDS",k=3, 
                weakties =TRUE, model="global",maxit = 400, try=80, trymax=100)
print(NMDST8)#Stress=0.07


#TimePoint9
DisTT9<-avgdist(TransTP9, 7115, iterations=100, meanfun=median, transf= sqrt, dmethod="bray")
NMDST9<-metaMDS(DisTT9, autotransform = FALSE, engine = "monoMDS",k=3, 
                weakties =TRUE, model="global",maxit = 400, try=80, trymax=100)
print(NMDST9)#Stress=0.06




# * * Export all distance matrix for easy analysis in future, edit grapsh ....................................
dir.create("Analysis/Diversity/Beta/Tables", recursive = TRUE)#Create folder

DistanceMatrixSoil<-as.data.frame(as.matrix(DistSoil))

DistMatrixTP1<-as.data.frame(as.matrix(DisTT1))
DistMatrixTP2<-as.data.frame(as.matrix(DisTT2))
DistMatrixTP3<-as.data.frame(as.matrix(DisTT3))
DistMatrixTP4<-as.data.frame(as.matrix(DisTT4))
DistMatrixTP5<-as.data.frame(as.matrix(DisTT5))
DistMatrixTP6<-as.data.frame(as.matrix(DisTT6))
DistMatrixTP7<-as.data.frame(as.matrix(DisTT7))
DistMatrixTP8<-as.data.frame(as.matrix(DisTT8))
DistMatrixTP9<-as.data.frame(as.matrix(DisTT9))


# * * Export distance matrix....................................................................................
write.csv(DistanceMatrixSoil, "1-Analysis/Diversity/Beta/Tables/DistMatrix-Soil.csv")

write.csv(DistMatrixTP1, "1-Analysis/Diversity/Beta/Tables/DistMatrix-TP1-SoilBac.csv")
write.csv(DistMatrixTP2, "1-Analysis/Diversity/Beta/Tables/DistMatrix-TP2-SoilBac.csv")
write.csv(DistMatrixTP3, "1-Analysis/Diversity/Beta/Tables/DistMatrix-TP3-SoilBac.csv")
write.csv(DistMatrixTP4, "1-Analysis/Diversity/Beta/Tables/DistMatrix-TP4-SoilBac.csv")
write.csv(DistMatrixTP5, "1-Analysis/Diversity/Beta/Tables/DistMatrix-TP5-SoilBac.csv")
write.csv(DistMatrixTP6, "1-Analysis/Diversity/Beta/Tables/DistMatrix-TP6-SoilBac.csv")
write.csv(DistMatrixTP7, "1-Analysis/Diversity/Beta/Tables/DistMatrix-TP7-SoilBac.csv")
write.csv(DistMatrixTP8, "1-Analysis/Diversity/Beta/Tables/DistMatrix-TP8-SoilBac.csv")
write.csv(DistMatrixTP9, "1-Analysis/Diversity/Beta/Tables/DistMatrix-TP9-SoilBac.csv")




#----
#-----
#********************************************************************************************************----
#-------------------------------- SIGNIFICANCE ..........................................................
#********************************************************************************************************----        
attach(MetaSoil)
library(vegan)


# * * PERMANOVA ..............................................................................................
# * * * FULL MODEL PERMANOVA w INTERATIONS ........................................................
#Permanova for full model to see variance partitioned, controlling for nestedness (random effects)
OtuRare1Soil<-read.csv("1-Analysis/ASVtables/OtuRare1-NC-Soil.csv", row.names = 1)


#Define the strata (control for nestedness)......................
perm <- how(nperm = 9999)
`setBlocks<-`(perm,Plot)

#Run Full model ........................................................................................................................
MetaSoil$TSFdays<-as.numeric(MetaSoil$TSFdays)


AdonFull<-adonis2(OtuRare1Soil~ Treatment*TSFdays*InitialAshDepth*TotPrecip, method= "bray",permutations=perm,data=MetaSoil);AdonFull


#Pairwise adonis using using a distance matrix as an input..............................................................................................................
Distance<-vegdist(OtuTrans2soil, "bray", upper=TRUE, diag=TRUE)
FullPW<-pairwise.adonis2(OtuRare1soil ~Treatment*as.factor(TSFdays), strata = 'Plot',
                data = MetaSoil,permutation=9999);FullPW$Burned_vs_Unburned
PW<-pairwise.adonis(Distance,factors=TSFdays,p.adjust.m='bonferroni')




# * * * EXPORTED RESULTS ........................................................................................
dir.create(file.path("1-Analysis/Diversity/Beta/Tables/Signif"),recursive = TRUE)

capture.output(AdonFull,file="1-Analysis/Diversity/Beta/Tables/Signif/Adonis-FullModel-BacALL5var.csv")



# * * TEST FOR HOMOGENEITY --------------------------------------------------------------------------------------
# permanova significant so we have to test for homogeneity
#_______HOMOGENEITY

#If PERMANOVA is significant but PERMDISP IS NOT, then you can infer that there is only a location effect. 
#If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a 
#location effect.


# * BETADISPER ------------------------------------------------------------------------------------------------------------
#  * * TREATMENT ..............................................................................................

Trt.dips<-betadisper(d=vegdist(OtuRare1soil, "bray"), group=Treatment, type="centroid")
Trt.dips$distances# show the distances

#Test for differences in distance and do pairwise comparisons using permutation................................
Trt.disp.Sign<-permutest(Trt.dips, pairwise = TRUE, permutations = 9999);Trt.disp.Sign #3e-04 ***
#Perm and Disp significant= Treatment and maybe location significant



# * * TSF .....................................................................................................
TSF.dips<-betadisper(d=vegdist(OtuRare1soil, "bray"), group= TSFdays, type="centroid"); TSF.dips
TSF.dips$distances# show the distances
permutest(TSF.dips)
plot(TSF.dips, hull=FALSE, ellipse=TRUE, label = FALSE)##sd ellipse, hull=FALSE, ellipse=TRUE)#sd ellipse

#Test for differences in distance and do pairwise comparisons using permutation ...............................
TSF.dips.Sign<-permutest(TSF.dips, pairwise = TRUE, permutations = 9999);TSF.dips.Sign 
#Pemanova & permdist NotSign = Only location effect



# * *  EXPORT THE RESULTS  ..................................................................
capture.output(Trt.disp.Sign,file="1-Analysis/Diversity/Beta/Tables/Signif/BetaDisp-TRT-Bac-sig.csv")
capture.output(TSF.dips.Sign,file="1-Analysis/Diversity/Beta/Tables/Signif/BetaDisp-TSF-Bac-sig.csv")
capture.output(Sev.dips.Sign,file="1-Analysis/Diversity/Beta/Tables/Signif/BetaDisp-Sev-Bac-sig.csv")



#----
#------
#******************************************************************************************************************----
#---DETERMINE SIGNIFICANCE PER TIMEPOINT ..........................................................
#*******************************************************************************************************************----        


# DETERMINE SIGNIFICANCE PER TIMEPOINT.................................................................................
#need to figure out why strata is not working

AdonTP1<-adonis2(ASVTP1 ~ Treatment,method = "bray",block="Plot",permutations=9999, data=TP1);AdonTP1
AdonTP2<-adonis2(ASVTP2 ~ Treatment,method = "bray",block="Plot",permutations=9999, data=TP2);AdonTP2
AdonTP3<-adonis2(ASVTP3 ~ Treatment,method = "bray",block="Plot",permutations=9999, data=TP3);AdonTP3
AdonTP4<-adonis2(ASVTP4 ~ Treatment,method = "bray",block="Plot",permutations=9999, data=TP4);AdonTP4
AdonTP5<-adonis2(ASVTP5 ~ Treatment,method = "bray",block="Plot",permutations=9999, data=TP5);AdonTP5
AdonTP6<-adonis2(ASVTP6 ~ Treatment,method = "bray",block="Plot",permutations=9999, data=TP6);AdonTP6
AdonTP7<-adonis2(ASVTP7 ~ Treatment,method = "bray",block="Plot",permutations=9999, data=TP7);AdonTP7
AdonTP8<-adonis2(ASVTP8 ~ Treatment,method = "bray",block="Plot",permutations=9999, data=TP8);AdonTP8
AdonTP9<-adonis2(ASVTP9 ~ Treatment,method = "bray",block="Plot",permutations=9999, data=TP9);AdonTP9



# * * Goodnes of fit-test .................................................
par(mfrow=c(1,1))
gof2<-goodness(NMDSsoil);gof2



#EXPORT RESULTS ..........................................................................................
write.csv(AdonTP1, "1-Analysis/Diversity/Beta/Tables/Signif/AdonTP1-Bac.csv")
write.csv(AdonTP2, "1-Analysis/Diversity/Beta/Tables/Signif/AdonTP2-Bac.csv")
write.csv(AdonTP3, "1-Analysis/Diversity/Beta/Tables/Signif/AdonTP3-Bac.csv")
write.csv(AdonTP4, "1-Analysis/Diversity/Beta/Tables/Signif/AdonTP4-Bac.csv")
write.csv(AdonTP5, "1-Analysis/Diversity/Beta/Tables/Signif/AdonTP5-Bac.csv")
write.csv(AdonTP6, "1-Analysis/Diversity/Beta/Tables/Signif/AdonTP6-Bac.csv")
write.csv(AdonTP7, "1-Analysis/Diversity/Beta/Tables/Signif/AdonTP7-Bac.csv")
write.csv(AdonTP8, "1-Analysis/Diversity/Beta/Tables/Signif/AdonTP8-Bac.csv")
write.csv(AdonTP9, "1-Analysis/Diversity/Beta/Tables/Signif/AdonTP9-Bac.csv")



#----
#--------
#******************************************************************************************************----
# NMDS PLOTS ----------------------------------------------------------------------------------------------
#***************************************************************************************************----
# **** Check that names match in NMDS and metadata
#Extract NMDS scores (x and y coordinates)...................................................
DataScoresSoil <- as.data.frame(scores(NMDSsoil))

# ** EXPORT SCORES FOR PLOTTING DATA ........................................................
dir.create(file.path("1-Analysis/Diversity/Beta/Tables"), recursive = TRUE)

write.csv(DataScoresSoil, "1-Analysis/Diversity/Beta/Tables/Bacteria-DataScoresSoils.csv")


#To run graphs only without running distance again ...........................................
#Reimport datascore files for "all" and "soil only" ................................
DataScoresSoil<-read.csv("1-Analysis/Diversity/Beta/Tables/Bacteria-DataScoresSoils.csv")


#Add columns to data frame to plot ................................
DataScoresSoil$Treatment<- MetaSoil$Treatment
DataScoresSoil$TimePoint<- MetaSoil$TimePoint
DataScoresSoil$TSFdays<- MetaSoil$TSFdays
DataScoresSoil$Plot<- MetaSoil$Plot



#----
#-----
#****************************************************************************************************----
#-------------------------------- CREATE GRAPHS ---------------------------------------------------------
#****************************************************************************************************----        

NMDSTP<-ggplot(DataScoresSoil, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(shape = as.factor(TSFdays), colour = Treatment))+ 
  scale_shape_manual(values=c(15,24,19,1,18,25,8,17,13))+
  theme_bw()+ theme(panel.grid  = element_blank(),
      plot.margin = margin(0.9, 1.5, 0.3, 0.3, "cm"),
      axis.text.y = element_text(size =17,hjust=.5, angle=90,colour = "black"), 
      axis.text.x = element_text(size = 17,colour = "black"), 
      axis.title = element_text(size = 20, colour = "black"), 
      legend.position = "bottom",
      legend.title = element_text(size = 20), 
      legend.text = element_text(size = 18)) +
  labs(x = "NMDS1", colour = "Treatment", 
       y = "NMDS2", shape = "Time Since Fire (days)")+ ylim(-0.5, 0.5)+
  guides(color=guide_legend(title.position="top", nrow=2, byrow=TRUE),
         shape = guide_legend(title.position="top",
        override.aes = list(size = 6),nrow=2, byrow=TRUE))+
  scale_colour_manual(values = c("#a2673f","#45877f"));NMDSTP

#NMDS at each specific timepoint..................................................
NMDSTP2<-NMDSTP+facet_wrap(~TSFdays)+
  theme(strip.text = element_text(face="bold", size=16),
        strip.background = element_rect(fill="white", colour="black"));NMDSTP2




NMDSTSFsoil<-ggplot(DataScoresSoil, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3.5, aes(colour = as.factor(TSFdays), shape=Treatment))+ 
  scale_shape_manual(values=c(15,24,19,1,18,25,8,17,13))+ 
  scale_color_manual(values = c("#70aeb8","#2e56a6","#a99cb8","#564666",
      "#abc4b5", "#487359","#bf9363","#7a4d04","red"))+
  theme_bw()+ theme(panel.grid  = element_blank(),
        plot.margin = margin(0.3, 1.5, 0.3, 0.3, "cm"),
        axis.text.y = element_text(size =20,hjust=.5, angle=90,colour = "black"), 
        axis.text.x = element_text(size = 20,colour = "black"), 
        axis.title = element_text(size = 18, colour = "black"), 
        legend.position = "bottom",
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18)) +
  labs(x = "NMDS1", colour = "Treatment", 
       y = "NMDS2", shape = "Time Since Fire (days)")+ ylim(-0.5, 0.5)+
  guides(color=guide_legend(title.position="top", nrow=2, byrow=TRUE),
         shape = guide_legend(title.position="top",
         override.aes = list(size = 6),nrow=2, byrow=TRUE));NMDSTSFsoil




#--EXPORT GRAPHS .........................................................................
dir.create(file.path("1-Analysis/Diversity/Beta/Graphs"), recursive = TRUE)

pdf("1-Analysis/Diversity/Beta/Graphs/Soil-NMDS-Trt-TSF.pdf", width=8.5, height=7.5)
NMDSTP
dev.off()


pdf("1-Analysis/Diversity/Beta/Graphs/NMDS-TimeSinceFire-Site-All.pdf", width=10, height=12)
NMDSTP2
dev.off()


pdf("1-Analysis/Diversity/Beta/Graphs/Soil-NMDS-TSF.pdf", width=8.5, height=7.5)
NMDSTSFsoil
dev.off()

#Script last ran, cleaned and saved on 1/7/2022

