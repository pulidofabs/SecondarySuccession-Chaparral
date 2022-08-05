rm(list=ls())#reset working directory

#Set working directory
setwd("C:/Users/juchu/Dropbox/3-HolyFire/1-HolyFire-Year1/1-Bacteria")


#load libraries ................................................................
library(vegan)
library(ggplot2)
library(plyr)

#Load data ................................................................................
MetaBurnB<-read.csv("1-Analysis/Metadata/MetaRareSppSoil-Burned.csv",row.names = 1)
MetaUnB<-read.csv("1-Analysis/Metadata/MetaRareSppSoil-Unburned.csv",row.names = 1)

OtuTrans2BacB<-read.csv("1-Analysis/ASVtables/OtuTrans2-NC-Soil-Burned.csv",row.names = 1)
OtuTrans2BacUn<-read.csv("1-Analysis/ASVtables/OtuTrans2-NC-Soil-Unburned.csv",row.names = 1)



#----
#----
#********************************************************************************************************----
# ------------  ANALYSIS ------------------------------------------------------------------------
#********************************************************************************************************----

#Calculate distance--will export one and reimport from here on out to just work with this instead........
DistBacB<-avgdist(OtuTrans2BacB, 7115, iterations=250, meanBac=median,transf= sqrt, dmethod="bray")
DistBacUn<-avgdist(OtuTrans2BacUn, 7115, iterations=250,  ppomeanBac=median,transf= sqrt, dmethod="bray")

dir.create(file.path("1-Analysis/Diversity/Beta/BetaDisp/Tables"),recursive = TRUE)
write.csv(as.matrix(DistBacB),"1-Analysis/Diversity/Beta/BetaDisp/Tables/DistMatrixBurn.csv")
write.csv(as.matrix(DistBacUn),"1-Analysis/Diversity/Beta/BetaDisp/Tables/DistMatrixUnbubrn.csv")


#Calculate Betadisper...................................................................................
BetaBacB<-betadisper(d=DistBacB, group= MetaBurnB$TSFdays, type="centroid");BetaBacB
BetaBacUn<-betadisper(d=DistBacUn, group= MetaUnB$TSFdays, type="centroid");BetaBacUn

AovBetaBB<-anova(BetaBacB);AovBetaBB # Pr(>F)=0.001718 **
AovBetaBU<-anova(BetaBacUn);AovBetaBU # Pr(>F)=0.2253

#Permutational test to get the F-value for analysis ...............................................
BetaBBsig<-permutest(BetaBacB, pairwise = T, permutations = 9999,strata="Plot");BetaBBsig#p=0.002 
BetaBUsig<-permutest(BetaBacB, pairwise = T, permutations = 9999,strata="Plot");BetaBUsig#p=0.001 


#Create boxplor to look at the data..............................................
BoxCol<-c("#CA0000","#fc7b03","#fcbe03","#fcf45b",
        "#f8ff99","#90C9D6","#0296bf","#CD96CD","#800080")

BacB<-boxplot(BetaBacB, xlab="Time Since Fire (days)", notch=FALSE, col=BoxCol)
BacUn<-boxplot(BetaBacUn, xlab="Time Since Fire (days)", notch=FALSE, col=BoxCol)


#Tukeys test.........................................................
BetaHSDBB <- TukeyHSD(BetaBacB);plot(BetaHSDBB)
BetaHSDBU<-TukeyHSD(BetaBacUn);plot(BetaHSDBU)

#Quick look at plot.. will create a nicer one using ggplot...........
plot(BetaBacB, hull=FALSE, ellipse=TRUE, conf=0.68)
plot(BetaBacUn, hull=FALSE, ellipse=TRUE, conf=0.68)

#Export all results from the above calculations.......
capture.output(BetaBacB, file="1-Analysis/Diversity/Beta/BetaDisp/Tables/BetaBacBurn.csv")
capture.output(BetaBacUn, file="1-Analysis/Diversity/Beta/BetaDisp/Tables/BetaBacUnburn.csv")
capture.output(BetaBBsig, file="1-Analysis/Diversity/Beta/BetaDisp/Tables/BetaBacBSignif.csv")
capture.output(BetaBUsig, file="1-Analysis/Diversity/Beta/BetaDisp/Tables/BetaBacUnSignif.csv")
capture.output(BetaHSDBB, file="1-Analysis/Diversity/Beta/BetaDisp/Tables/BetaBac-TukeyUn.csv")
capture.output(BetaHSDBU, file="1-Analysis/Diversity/Beta/BetaDisp/Tables/BetaBac-TukeyB.csv")



dir.create(file.path("1-Analysis/Diversity/Beta/BetaDisp/Graphs"),recursive = TRUE)
pdf( "1-Analysis/Diversity/Beta/BetaDisp/Graphs/TukeysBetaHSDFburn.pdf", width=5, height=5)
plot(BetaHSDFB)
dev.off()

pdf("1-Analysis/Diversity/Beta/BetaDisp/Graphs/TukeysBetaHSDBacburn.pdf", height = 5, width = 5)
plot(BetaHSDFU)
dev.off()


#----
#----
#************************************************************************************************----
#---------------- CREATE GRAPHS BURNED PLOTS ONLY ---------------------------------------------------------------------
#************************************************************************************************----

#Burn........................................................
TSFburn<-MetaBurnB[,21];TSFburn #Extract TSF column
TSBac<-MetaUnB[,21];TSBac #Extract TSF column

#Extractscore values from the Betadisp.......................
ScoresBB<-scores(BetaBacB)
ScoresBU<-scores(BetaBacUn)

#Extract the PCOA1 and PCOA2 points (Sites)..................
SitesBurn<-as.data.frame(ScoresBB$sites)
SitesBurn1<-cbind(SitesBurn,TSFburn); head(SitesBurn1)


#Extract centroids from Betadisp
CentroidB<-as.data.frame(ScoresBB$centroids)

#Add row names and specify the number of replicatesp per...........................
#CentroidB$SiteID<-rownames(CentroidB)
CentroidB$TSFdays<-rep(c("17","25","34","67","95","131","187","286","376"),each=1)

# Calculate segments...............................................................
names(SitesBurn1)[3]<-"TSFdays"#remane the third column
SegBac<-SitesBurn1[, c("PCoA1", "PCoA2", "TSFdays")];head(SegBac[1:2,])

#Make a copy of the centroids and rename columns..................................
#rename the columns
tmp1<-CentroidB
names(tmp1)[1]<-"PCoA1_ctr"
names(tmp1)[2]<-"PCoA2_ctr"

SegBac1<-join(SegBac, tmp1, "TSFdays");SegBac1[1:2,]


#Create a dataframe of the centroids and PCOA points...............................
#Dataframe needs to have same # of columns, and column names, if they dont
#have one, add a column and values should be NA, can remove later
#Example............CentroidB$SiteID = NA

#Add column to specify what the data is............................................
CentroidB$data.type = "centroid"
SitesBurn1$data.type = "site"

#verify that dataframes have the same column in same order
all(colnames(CentroidB)==(colnames(SitesBurn1))) #true = proceed to bind

#Bind centroids and sites..........................................................
SiteCentroid<-rbind(CentroidB, SitesBurn1); head(SiteCentroid[1:2,])
SiteCentroid$data.type


#Export files
write.csv(SegBac1,"1-Analysis/Diversity/Beta/BetaDisp/SegmentsBac.csv")
write.csv(SiteCentroid,"1-Analysis/Diversity/Beta/BetaDisp/SiteCentroids.csv")

#----
#----
#***************************************************************************************----
#-------------------- QUALITY CONTROL ------------------------------------------------------
#***************************************************************************************----

#Turn variables into factor but it does not help with getting number in legend in order
SegBac1$TSFdays<-as.factor(SegBac1$TSFdays)
SiteCentroid$TSFdays<-as.factor(SiteCentroid$TSFdays)


#change the cols vector definitions at the beginning of code to this
cols.colour <- c("#540203","#CA0000","#fc7b03","#fcbe03","#86b387",
                 "#90C9D6","#0296bf","#CD96CD","#800080")

shapeCol<-c("#540203","#CA0000","#fc7b03","#fcbe03","#86b387",
            "#90C9D6","#0296bf","#CD96CD","#800080")
#segment must go before point so points are in front of lines
BetaBacB<-ggplot() + 
  geom_segment(data = SegBac1,aes(x = PCoA1, y = PCoA2, xend = PCoA1_ctr,
            yend = PCoA2_ctr, colour = TSFdays)) +
  geom_point(data = SiteCentroid,aes(x = PCoA1, y = PCoA2, colour = TSFdays,
          fill = TSFdays, shape = TSFdays), size = 3) +
  scale_colour_manual(values = cols.colour) +
  scale_fill_manual(values = shapeCol) +
  scale_shape_manual(values = c(15,24,19,1,18,25,8,17,13)) + 
  coord_equal() +
  theme_bw()+
  theme(plot.margin = margin(0.8,0.8, 0.8, 0.8, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size =20, angle=90, hjust=0.5, colour = "black"), 
        axis.text.x = element_text(size = 20,colour = "black"), 
        axis.title.y = element_text(size = 24, hjust=0.5,colour = "black"), 
        axis.title.x = element_text(size = 24, colour = "black"),
        legend.title = element_text(size=24),
        legend.text=element_text(size = 20));BetaBacB



#----
#----
#************************************************************************************----
#*------------- MEAN AND STD. ERR BETADISPER -------------------------------------------
#*#**********************************************************************************----
SegBac1<-read.csv("1-Analysis/Diversity/Beta/BetaDisp/SegmentsBac.csv")


#Create mean and sd OF GRAPHS
PCoA1B<-SegBac1$PCoA1;PCoA1B
PCoA2B<-SegBac1$PCoA2;PCoA2B
TSF<-SegBac1$TSFdays
DataBacB <- data.frame(TSF, PCoA1B, PCoA2B)
DataBacB2<- merge(DataBacB,aggregate(cbind(mean.x=PCoA1B, mean.y=PCoA2B)~TSF,DataBacB,mean),by="TSF")


CentBacB <- aggregate(cbind(PCoA1B,PCoA2B) ~ TSF, DataBacB2, mean)
FB        <- function(z)sd(z)/sqrt(length(z)) # Function to calculate std.err
seB        <- aggregate(cbind(se.x=PCoA1B, se.y=PCoA2B)~TSF, DataBacB2, FB)
CentBacSEB <- merge(CentBacB,seB, by="TSF")# add std.err column to centroids


#make TSF a factor
CentBacSEB$TSF<-as.factor(CentBacSEB$TSF)

TSFBacB<-ggplot(CentBacSEB, aes(PCoA1B, PCoA2B, color=TSF))+
  geom_point(data=CentBacSEB, size=4)+
  geom_errorbar(data=CentBacSEB, aes(ymin= PCoA2B-se.y, ymax= PCoA2B+se.y), width=0.01)+
  geom_errorbarh(data=CentBacSEB, aes(xmin= PCoA1B-se.x, xmax= PCoA1B+se.x), height=0.01)+
  scale_colour_manual(values = cols.colour) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size=24),
        legend.position = "bottom",
        legend.text=element_text(size = 20),
        axis.text.y = element_text(size =20, angle=90, hjust=0.5, colour = "black"), 
        axis.text.x = element_text(size = 20,colour = "black"), 
        axis.title.y = element_text(size = 24, hjust=0.5,colour = "black"), 
        axis.title.x = element_text(size = 24, colour = "black"));TSFBacB




#----
#----
#************************************************************************************************----
#---------------- CREATE GRAPHS UNBURNED PLOTS ONLY ---------------------------------------------------------------------
#************************************************************************************************----

#Un........................................................
TSBacUn<-MetaUnB[,21];TSBacUn #Extract TSF column
TSBacUn<-MetaUnB[,21];TSBacUn #Extract TSF column

#Extractscore values from the Betadisp.......................
ScoresBU<-scores(BetaBacUn)

#Extract the PCOA1 and PCOA2 points (Sites)..................
SitesUn<-as.data.frame(ScoresBU$sites)
SitesUn1<-cbind(SitesUn,TSBacUn); head(SitesUn1)


#Extract centroids from Betadisp
CentroidUn<-as.data.frame(ScoresBU$centroids)

#Add row names and specify the number of replicatesp per...........................
#CentroidB$SiteID<-rownames(CentroidB)
CentroidUn$TSFdays<-rep(c("17","25","34","67","95","131","187","286","376"),each=1)

# Calculate segments...............................................................
names(SitesUn1)[3]<-"TSFdays"#remane the third column
SegBacUn<-SitesUn1[, c("PCoA1", "PCoA2", "TSFdays")];head(SegBacUn[1:2,])

#Make a copy of the centroids and rename columns..................................
#rename the columns
tmp1un<-CentroidUn
names(tmp1un)[1]<-"PCoA1_ctr"
names(tmp1un)[2]<-"PCoA2_ctr"

SegBacUn1<-join(SegBacUn, tmp1, "TSFdays");SegBacUn1[1:2,]


#Create a dataframe of the centroids and PCOA points...............................
#Dataframe needs to have same # of columns, and column names, if they dont
#have one, add a column and values should be NA, can remove later
#Example............CentroidB$SiteID = NA

#Add column to specify what the data is............................................
CentroidUn$data.type = "centroid"
SitesUn1$data.type = "site"

#verify that dataframes have the same column in same order
all(colnames(CentroidUn)==(colnames(SitesUn1))) #true = proceed to bind

#Bind centroids and sites..........................................................
SiteCentroidUn<-rbind(CentroidUn, SitesUn1); head(SiteCentroidUn[1:2,])
SiteCentroidUn$data.type


#Export files
write.csv(SegBacUn1,"1-Analysis/Diversity/Beta/BetaDisp/Tables/SegmentsBac-Unburn.csv")
write.csv(SiteCentroidUn,"1-Analysis/Diversity/Beta/BetaDisp/Tables/SiteCentroids-Unburn.csv")

#----
#----
#***************************************************************************************----
#-------------------- QUALITY CONTROL ------------------------------------------------------
#***************************************************************************************----

#Turn variables into factor but it does not help with getting number in legend in order
SegBacUn1$TSFdays<-as.factor(SegBacUn1$TSFdays)
SiteCentroidUn$TSFdays<-as.factor(SiteCentroidUn$TSFdays)


#change the cols vector definitions at the beginning of code to this
cols.colour <- c("#540203","#CA0000","#fc7b03","#fcbe03","#86b387",
                 "#90C9D6","#0296bf","#CD96CD","#800080")

shapeCol<-c("#540203","#CA0000","#fc7b03","#fcbe03","#86b387",
            "#90C9D6","#0296bf","#CD96CD","#800080")

#segment must go before point so points are in front of lines
BetaBacUn<-ggplot() + 
  geom_segment(data = SegBacUn1,aes(x = PCoA1, y = PCoA2, xend = PCoA1_ctr,
          yend = PCoA2_ctr, colour = TSFdays)) +
  geom_point(data = SiteCentroidUn,aes(x = PCoA1, y = PCoA2, colour = TSFdays,
             fill = TSFdays, shape = TSFdays), size = 3) +
  scale_colour_manual(values = cols.colour) +
  scale_fill_manual(values = shapeCol) +
  scale_shape_manual(values = c(15,24,19,1,18,25,8,17,13)) + 
  coord_equal() +
  theme_bw()+
  theme(plot.margin = margin(0.8,0.8, 0.8, 0.8, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size =20, angle=90, hjust=0.5, colour = "black"), 
        axis.text.x = element_text(size = 20,colour = "black"), 
        axis.title.y = element_text(size = 24, hjust=0.5,colour = "black"), 
        axis.title.x = element_text(size = 24, colour = "black"),
        legend.title = element_text(size=24),
        legend.text=element_text(size = 20));BetaBacUn



#----
#----
#************************************************************************************----
#*------------- MEAN AND STD. ERR BETADISPER -------------------------------------------
#*#**********************************************************************************----
SegBacUn1<-read.csv("1-Analysis/Diversity/Beta/BetaDisp/Tables/SegmentsBac-Unburn.csv")


#Create mean and sd OF GRAPHS
PCoA1Un<-SegBacUn1$PCoA1;PCoA1Un
PCoA2Un<-SegBacUn1$PCoA2;PCoA2Un
TSF<-SegBacUn1$TSFdays
DataBacU <- data.frame(TSF, PCoA1Un, PCoA2Un)
DataBacU2<- merge(DataBacU,aggregate(cbind(mean.x=PCoA1Un, mean.y=PCoA2Un)~TSF,DataBacU,mean),by="TSF")


CentBacUn <- aggregate(cbind(PCoA1Un,PCoA2Un) ~ TSF, DataBacU2, mean)
FU       <- function(z)sd(z)/sqrt(length(z)) # Bacction to calculate std.err
seU        <- aggregate(cbind(se.x=PCoA1Un, se.y=PCoA2Un)~TSF, DataBacU2, FU)
CentBacSEU <- merge(CentBacUn,seU, by="TSF")# add std.err column to centroids



#make TSF a factor
CentBacSEU$TSF<-as.factor(CentBacSEU$TSF)

TSFBacUn<-ggplot(CentBacSEU, aes(PCoA1Un, PCoA2Un, color=TSF))+
  geom_point(data=CentBacSEU, size=4)+
  geom_errorbar(data=CentBacSEU, aes(ymin= PCoA2Un-se.y, ymax= PCoA2Un+se.y), width=0.01)+
  geom_errorbarh(data=CentBacSEU, aes(xmin= PCoA1Un-se.x, xmax= PCoA1Un+se.x), height=0.01)+
  scale_colour_manual(values = cols.colour) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size=24),
        legend.position = "bottom",
        legend.text=element_text(size = 20),
        axis.text.y = element_text(size =20, angle=90, hjust=0.5, colour = "black"), 
        axis.text.x = element_text(size = 20,colour = "black"), 
        axis.title.y = element_text(size = 24, hjust=0.5,colour = "black"), 
        axis.title.x = element_text(size = 24, colour = "black"));TSFBacUn


#----
#----
#**********************************************************************************************----
#------------------------- EXPORT GRAPHS -------------------------------------------------
#**********************************************************************************************----
#Export graphs........................................................................................
pdf("1-Analysis/Diversity/Beta/BetaDisp/Graphs/BetaDisp-Unburn.pdf", height = 8, width = 8)
BetaBacUn
dev.off()

pdf("1-Analysis/Diversity/Beta/BetaDisp/Graphs/BetaDisp-Burn.pdf", height = 8, width = 8)
BetaBacB
dev.off()


#Create panels for exporting...................................................
BetaBU<-ggarrange(TSFBacUn,TSFBacB, ncol=2, nrow=1, common.legend = TRUE, 
          legend="bottom", labels = c("a", "b"), hjust = c(-.8,-1),
          align="hv",font.label = list(size=16,face="bold"));BetaBU

pdf("1-Analysis/Diversity/Beta/BetaDisp/Graphs/BetaDispBac-Burn-Unburn.pdf", height = 8, width = 13)
BetaBU
dev.off()
