#12-7-2002

rm(list=ls())#reset working directory


#Load libraries -------------------------------------------------------------------
library(vegan)
library(knitr)
library(tidyverse)

#Set working directory ------------------------------------------------------------
setwd("~/Dropbox/3-HolyFire/1-HolyFire-Year1/1-Bacteria")
setwd("C:/Users/juchu/Dropbox/3-HolyFire/1-HolyFire-Year1/1-Bacteria/")

#----
#----
#******************************************************************************************************************----
# (1) RAREFY TABLE TO MATCH THE RAREFIED TABLE ONLY DO ONCE AFTERWARDS SKIP TO THE STEP (2-----------------------------
#******************************************************************************************************************----
#Load data for this step only ................................................................................
MetaSoil<- read.csv("1-Analysis/Metadata/MetaRareSpp-Soil.csv", na.strings = "N/A",header = TRUE)
OtuTrans2Soil<-read.csv("1-Analysis/ASVtables/OtuTrans2-NC-Soil.csv",check.names = FALSE)

#Subset metadata by treatment....................................
attach(MetaSoil)
MetaBurn<-MetaSoil[which(Treatment== "Burned"), ];dim(MetaBurn)#210x68
MetaUn<-MetaSoil[which(Treatment== "Unburned"), ];dim(MetaUn)#99x68
detach(MetaSoil)

#Check names of datafiles.............................................
names(MetaBurn)[1];names(OtuTrans2Soil)[1];names(MetaUn)[1]

#Add name to file to make sure they match.............................
names(OtuTrans2Soil)[1]<- "SampleID"#rename ASVtable

#Verify that column1 has the same name in both tables.................
names(MetaBurn)[1];names(OtuTrans2Soil)[1];names(MetaUn)[1]
dim(MetaBurn);dim(OtuTrans2Soil);dim(MetaUn)
#210X68,            309X35742,       99x68

# keep rows with matching ID.........................................................................
#Row collumns should match the Metadata
OtuTrans2SoilB<-OtuTrans2Soil %>% semi_join(MetaBurn, by = "SampleID");dim(OtuTrans2SoilB)#210x35742
OtuTrans2SoilUn<-OtuTrans2Soil %>% semi_join(MetaUn, by = "SampleID");dim(OtuTrans2SoilUn)#99x25742


#Move column one into rownames since when you export it add an extra column....................
OtuTrans2SoilB1<-OtuTrans2SoilB %>% remove_rownames %>% column_to_rownames(var="SampleID")
        rownames(OtuTrans2SoilB1)[1];dim(OtuTrans2SoilB1)
    
OtuTrans2SoilUn1<-OtuTrans2SoilUn %>% remove_rownames %>% column_to_rownames(var="SampleID")
         rownames(OtuTrans2SoilUn1)[1];dim(OtuTrans2SoilUn1)

#Do the same to the metadata
MetaBurn1<-MetaBurn %>% remove_rownames %>% column_to_rownames(var="SampleID")
         rownames(MetaBurn1)[1];dim(MetaBurn1)
         
MetaUn1<-MetaUn %>% remove_rownames %>% column_to_rownames(var="SampleID")
         rownames(MetaUn1)[1];dim(MetaUn1)
         
         
#EXPORT RAREFIED TABLES...........................................................
write.csv(OtuTrans2SoilB1, "1-Analysis/ASVtables/OtuTrans2-NC-Soil-Burned.csv")
write.csv(OtuTrans2SoilUn1, "1-Analysis/ASVtables/OtuTrans2-NC-Soil-Unburned.csv")
write.csv(MetaBurn1, "1-Analysis/Metadata/MetaRareSppSoil-Burned.csv")
write.csv(MetaUn1, "1-Analysis/Metadata/MetaRareSppSoil-Unburned.csv")



#----
#----
#*****************************************************************************************************************(----
#IMPORT CORRECT FILES FOR ANALYSIS-------------------------------------------------------------------------------------
#******************************************************************************************************************----

MetaSoilUn<- read.csv("1-Analysis/Metadata/MetaRareSppSoil-Unburned.csv", na.strings = "N/A",header = TRUE, row.name=1)
OtuTrans2SoilUn<-read.csv("1-Analysis/ASVtables/OtuTrans2-NC-Soil-Unburned.csv",check.names = FALSE, row.name=1)

#Check dimensions of dataframes, they should match...................
dim(MetaSoilUn); dim(OtuTrans2SoilUn)#
all(rownames(MetaSoilUn)== rownames(OtuTrans2SoilUn))

#----
#----
#***********************************************************************************************-----
#*#--------- CALCULATE DISTANCE MATRIX BOTH BAC DATA AND TEMPORAL DATA ------------------------------
#*#*********************************************************************************************-----
library(ecodist)
#Import rarefied Burned data ...................................................................................

#Burned dataset...................................................................
MetaSoilB<- read.csv("1-Analysis/Metadata/MetaRareSppSoil-Burned.csv",na.strings = "N/A",row.name=1)
OtuTrans2SoilB<-read.csv("1-Analysis/ASVtables/OtuTrans2-NC-Soil-Burned.csv",check.names = FALSE, row.name=1)

#Unburned dataset ................................................................
MetaSoilUn<- read.csv("1-Analysis/Metadata/MetaRareSppSoil-Unburned.csv",na.strings = "N/A",row.name=1)
OtuTrans2SoilUn<-read.csv("1-Analysis/ASVtables/OtuTrans2-NC-Soil-Unburned.csv",check.names = FALSE, row.name=1)


#Verify dimensions and verify that rows match.............
dim(MetaSoilB); dim(OtuTrans2SoilB)#210 rows
all(rownames(MetaSoilB)== rownames(OtuTrans2SoilB))

dim(MetaSoilUn); dim(OtuTrans2SoilUn)#99 rows
all(rownames(MetaSoilUn)== rownames(OtuTrans2SoilUn))


#Create the distance matrix for the burned and unburned data individually................................
par(mfrow=c(1,1),mar=c(2, 2, 0.5, 0.5))

#Taxonomic distance.............................................................
BacDistB<-avgdist(OtuTrans2SoilB,7115,iterations=250,meanfun=median,transf=sqrt,dmethod="bray")
BacDistUn<-avgdist(OtuTrans2SoilUn, 7115, iterations=250,meanfun=median,transf= sqrt,dmethod="bray")


#Temporal distance (Euclidean)...................................................
TempDistB<-distance(MetaSoilB$TSFdays)
TempDistUn<-distance(MetaSoilUn$TSFdays)


#Export data.......................................................................
dir.create(file.path("1-Analysis/Diversity/Beta/Mantel/Table"), recursive = TRUE)
write.csv(as.matrix(BacDistB), "1-Analysis/Diversity/Beta/Mantel/Table/BrayCurtis-Burned-Bac.csv")
write.csv(as.matrix(BacDistUn), "1-Analysis/Diversity/Beta/Mantel/Table/BrayCurtis-Unburned-Bac.csv")



#Create a simple graph just to see what data looks like------------------------- 
color=rgb(0,0,0,alpha=0.5) 

#Burned Dataset.................................................................    
plot(jitter(TempDistB),BacDistB,xlab=" ",ylab=" ",col=color, cex=0.5)
    #text(9, 0.95, sprintf("%s-%s", Sp, treat), col="red", cex =2.5, font =2)
    MantelB<-mantel(BacDistB ~ TempDistB, nperm = 999)
    text(350, 0.4, sprintf("R=%.3f, P=%.3f", MantelB[1], MantelB[2]),
         col="black", cex =1, font =2)


#Unburned Dataset...............................................................
plot(jitter(TempDistUn),BacDistUn, xlab=" ",ylab=" ",col=color, cex=0.5)
    #text(9, 0.95, sprintf("%s-%s", Sp, treat), col="red", cex =2.5, font =2)
    MantelUn<-mantel(BacDistUn ~ TempDistUn, nperm = 999)
    text(350, 0.5, sprintf("R=%.3f, P=%.3f", MantelUn[1], MantelUn[2]),
         col="black", cex =1, font =2)
    
    
   
#Export mantel results................................................................
write.csv(MantelUn, "1-Analysis/Diversity/Beta/Mantel/Table/MantelResults-Unburned.csv")
write.csv(MantelB, "1-Analysis/Diversity/Beta/Mantel/Table/MantelResults-Burned.csv")
    
    
     
#----
#----
#*******************************************************************************----
#---- -------------- BURNED PLOTS  -------------------------------------------------
#*******************************************************************************----
#abundance vs temperature
   library(ggplot2)    
BDB = as.vector(BacDistB)
TDB = as.vector(TempDistB)
    
#new data frame with vectorized distance matrices
DistMatB<-data.frame(BDB,TDB)#distance matrix of BacDistB and TempDistB
  
#Export distance matrix in case we need to remake graphs
write.csv(DistMatB, "1-Analysis/Diversity/Beta/Mantel/Table/DistMatBurned-Bac-Graph.csv")

#Create graph using distance matrix and keep X-axis as matrix    
MantelTSFBacB<-ggplot(DistMatB, aes(y = BDB, x = TDB)) + 
  geom_point(size=1, aes(color = TDB), show.legend = F) + 
  geom_smooth(method = "lm", colour = "black") + 
  labs(x = "", #Difference in Time Since Fire (days)
       y = "")+ #Bray-Curtis Dissimilarity
  theme(plot.margin = margin(1, 0.8, 0, 0, "cm"),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA,colour = "black"), 
        axis.title= element_text(face = "bold", size = 22, colour = "black"),
        axis.text.x = element_text(size = 19,colour = "black"), 
        axis.text.y = element_text(size = 19,colour = "black",hjust=.6,angle=90))+ 
  scale_color_continuous(high = "#b8a17b", low = "#5c3216")+
  annotate("text", x=300,  y=0.20, 
           label ="R = 0.333, P = 0.001", size=6);MantelTSFBacB


#Secondary graph, but X axis shows the actual sampling timepoints  
MantelTSFBacB1<-ggplot(DistMatB, aes(y = BDB, x = TDB)) + 
  geom_point(size=1, aes(color = TDB),show.legend = F) + 
  geom_smooth(method = "lm", colour = "black") + 
  labs(x = "",#Difference in Time Since Fire (days) 
       y = "")+ #Bray-Curtis Dissimilarity)  
  theme(plot.margin = margin(1, 0.8, 0, 0, "cm"),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA,colour = "black"), 
        axis.title= element_text(face = "bold", size = 22, colour = "black"),
        axis.text.x = element_text(size = 18,colour = "black",angle=90,hjust=1), 
        axis.text.y = element_text(size = 19,colour = "black",hjust=.6,angle=90))+ 
  scale_color_continuous(high = "#b8a17b", low = "#5c3216")+
  annotate("text", x=300,  y=0.20,label ="R = 0.333, P = 0.001", size=6)+
  scale_x_continuous(limits=c(0, 376), 
       breaks=c(0,17,25,34,67,95,131,187,286,376));MantelTSFBacB1
    

#----
#----
#*******************************************************************************----
#---- -------------- UNBURNED PLOTS  -------------------------------------------------
#*******************************************************************************----
#abundance vs temperature
BDU = as.vector(BacDistUn)
TDU = as.vector(TempDistUn)

#new data frame with vectorized distance matrices
DistMatUn<-data.frame(BDU,TDU)#distance matrix of BacDistB and TempDistB
write.csv(DistMatUn, "1-Analysis/Diversity/Beta/Mantel/DistMatUnburned-Bac-Graph.csv")


#Create graph using distance matrix and keep X-axis as matrix    
MantelTSFBacUn<-ggplot(DistMatUn, aes(y = BDU, x = TDU)) + 
  geom_point(size=1, aes(color = TDU), show.legend=F) + 
  geom_smooth(method = "lm", colour = "black") + 
  labs(x = "Difference in Time Since Fire (days)", 
       y = "Bray-Curtis Dissimilarity") + 
  theme(plot.margin = margin(1, 0.8, 0, 0, "cm"),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA,colour = "black"), 
        axis.title= element_text(face = "bold", size = 22, colour = "black"),
        axis.text.x = element_text(size = 18,colour = "black"), 
        axis.text.y = element_text(size = 19,colour = "black",hjust=.6,angle=90))+ 
  scale_color_continuous(low = "#004a41", high = "#82a8a1")+
  annotate("text", x=300,  y=0.40, 
           label ="R = 0.046, P = 0.122", size=6);MantelTSFBacUn


#Secondary graph, but X axis shows the actual sampling timepoints  
MantelTSFBacUn1<-ggplot(DistMatUn, aes(y = BDU, x = TDU)) + 
  geom_point(size=1, aes(color = TDU), show.legend = F) + 
  geom_smooth(method = "lm", colour = "black") + 
  labs(x = "Difference in Time Since Fire (days)", 
       y = "Bray-Curtis Dissimilarity", 
       color = "") + 
  theme(plot.margin = margin(1, 0.8, 0, 0, "cm"),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA,colour = "black"), 
        axis.title= element_text(face = "bold", size = 22, colour = "black"),
        axis.text.x = element_text(size = 18,colour = "black",angle=90,hjust=1), 
        axis.text.y = element_text(size = 19,colour = "black",hjust=.6,angle=90))+ 
  scale_color_continuous(low = "#004a41", high = "#82a8a1")+
  annotate("text",x=300,y=0.40,label ="R = 0.046, P = 0.122",size=6)+
  scale_x_continuous(limits=c(0, 376), 
     breaks=c(0,17,25,34,67,95,131,187,286,376));MantelTSFBacUn1




#Create panels of the graphs between unburned and unburned plots................
library(ggpubr)
# TSF MATRIX AS IS............................................................
ManTSFa<-ggarrange(MantelTSFFunUn,MantelTSFFunB, ncol=1, nrow=2,align="hv", 
                   common.legend = TRUE, font.label = list(size=18,face="bold"),
                   legend="bottom", labels = c("a Unburned", "b Burned"),
                   hjust = c(-0.2, -0.3));ManTSFa

ManTSFa1<-ggarrange(MantelTSFFunUn,MantelTSFFunB, ncol=2, nrow=1,align="hv", 
                    common.legend = TRUE, font.label = list(size=18,face="bold"),
                    legend="bottom", labels = c("a Unburned", "b Burned"),
                    hjust = c(-0.2, -0.3));ManTSFa1


#W actual sampling timepoint on  x asix.......................................
ManTSF<-ggarrange(MantelTSFBacUn1,MantelTSFBacB1, ncol=1, nrow=2,align="hv", 
                  common.legend = TRUE, font.label = list(size=18,face="bold"),
                  legend="bottom", labels = c("a Unburned", "b Burned"),
                  hjust = c(-.8,-1));ManTSF

ManTSF1<-ggarrange(MantelTSFBacUn1,MantelTSFBacB1, ncol=2, nrow=1,align="hv", 
                   common.legend = TRUE, font.label = list(size=18,face="bold"),
                   legend="bottom", labels = c("a Unburned", "b Burned"),
                   hjust = c(-.8,-1));ManTSF1




#*******************************************************************************----
#Export graphs----------------------------------------------------------------------    
#*******************************************************************************----

#Re run graphs but remove y anx x labels from burned plots
pdf("1-Analysis/Diversity/Beta/Mantel/Graphs/Mantel-BrayEucFun-BvsUn-Vert.pdf",height=12,width=8)
ManTSFa
dev.off()

pdf("1-Analysis/Diversity/Beta/Mantel/Graphs/Mantel-BrayEucFun-BvsUn-Hor.pdf",height=8,width=16)
ManTSFa1
dev.off()

pdf("1-Analysis/Diversity/Beta/Mantel/Graphs/Mantel-BrayEucFun-BvsUn-TSFhor.pdf",height=12,width=8)
ManTSF
dev.off()

pdf("1-Analysis/Diversity/Beta/Mantel/Graphs/Mantel-BrayEucFun-BvsUn-TSFvert.pdf",height=8,width=16)
ManTSF1
dev.off()









