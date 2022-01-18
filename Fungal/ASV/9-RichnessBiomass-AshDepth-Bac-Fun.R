#3-27-2021

rm(list=ls())#reset working directory

#Set working directory=---------------------------------------------------------
setwd("C:/Users/fabipc/Dropbox/3-HolyFire/1-HolyFire-Pub1")

#LOAD PACKAGES ------------------------------------------------------------------------------------
library(tidyverse) #required to load tax table
library(ggpubr) 
library(ggplot2)

# LOAD DATA (rarefied table and metadata)-------------------------------------------------------------------
MetaRareFun<- read.csv("2-Fungi/Metadata/Fungi/MetaSoilSpp.csv", na.strings = "NA", header = TRUE)
MetaRareBac<- read.csv("1-Bacteria/1-Analysis/Metadata/MetaSoilSpp.csv", na.strings = "N/A", header = TRUE)
MetaBio<- read.csv("4-Qpcr/HolyfireT1toT9/Calculating_Copy_Number_HolyFire.csv",na.strings = "NA", header = TRUE)



#set as factor first to set base level..............................................
MetaBio$Treatment<-as.factor(MetaBio$Treatment)
MetaRareBac$Treatment<-as.factor(MetaRareBac$Treatment)
MetaRareFun$Treatment<-as.factor(MetaRareFun$Treatment)

MetaBio$Treatment <- try(relevel(MetaBio$Treatment , "Unburned"))
MetaRareBac$Treatment <- try(relevel(MetaRareBac$Treatment , "Unburned"))
MetaRareFun$Treatment <- try(relevel(MetaRareFun$Treatment , "Unburned"))


#.----
#.----
#**********************************************************************************************************----
# QUALITY CONTROL )-----------------------------------------------------------------------------------------
#**********************************************************************************************************----
#Filter out the extreme outliers....why did we choose this number?????????????????
MetaBio<- filter(MetaBio,CopyNumber < 2000000)


attach(MetaBio)
MetaBioFun<-MetaBio[which(Kingdom== "Fungi"), ];dim(MetaBioFun)#306x25
MetaBioBac<-MetaBio[which(Kingdom== "Bacteria"), ];dim(MetaBioBac)#303x25
detach(MetaBio)

#Separate by treatment between kingdoms..........................................
BioFunB<-MetaBioFun[which(Treatment== "Burned"), ];dim(BioFunB)#205x25
BioFunUn<-MetaBioFun[which(Treatment== "Unburned"), ];dim(BioFunUn)#101x25

BioBacB<-MetaBioBac[which(Treatment== "Burned"), ];dim(BioBacB)#205x25
BioBacUn<-MetaBioBac[which(Treatment== "Unburned"), ];dim(BioBacUn)#101x25





#Quality control.............................................................................
#Need data as numerical, need to convert interval to number
str(MetaRareFun$CopyNumber); str(MetaRareBac$CopyNumber)
str(MetaBioFun$CopyNumber);str(MetaBioBac$CopyNumber)


#Conversion.............Species richness............................
MetaRareFun$S.obs<-as.numeric(MetaRareFun$S.obs)
MetaRareBac$S.obs<-as.numeric(MetaRareBac$S.obs)
MetaRareFun$TSFdays<-as.factor(MetaRareFun$TSFdays)
MetaRareBac$TSFdays<-as.factor(MetaRareBac$TSFdays)

#Conversion.............Biomass...................................
MetaBioFun$CopyNumber<-as.numeric(MetaBioFun$CopyNumber)
MetaBioBac$CopyNumber<-as.numeric(MetaBioBac$CopyNumber)
MetaBioFun$TSFdays<-as.factor(MetaBioFun$TSFdays)
MetaBioBac$TSFdays<-as.factor(MetaBioBac$TSFdays)


#----
#----
#*****************************************************************************************************-----
# ----------------------------- CREATE PLOTS --------------------------------------------------------------
#*****************************************************************************************************-----
# * * CREATE Color scheme ...............................................
#Scheme for timepoitns ..................................................
colorTP<-c("#de0025","#b05800","#ff6a00","#ffd000","#fac2e3","black","#0a2bff","#00b5f7","#007d2a")

#Maintain same color for all graphs .....................................
colorPlot<-c("black","#5e4839","#947f70","#c49f68","#bfb3aa","#ebddc7","#235978","#45877f","#aed4ca")


#Create graphs ..............................................................................
attach(MetaRareFun)
SppAshFun<-ggplot(MetaRareFun, aes(x=InitialAshDepth, y=S.obs))  +
 geom_smooth(method = MASS::glm.nb, se = TRUE)+ 
  scale_shape_manual(values=c(15,24,19,1,18,25,8,17,13))+
  geom_point(size=2,aes(shape = as.factor(TSFdays), colour = Treatment))+
  stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label..,
      sep = "~~~~")), formula = y~x, label.x = 5, label.y =600, size=6)+
  stat_cor(label.y = 600)+
  scale_color_manual(values=  c("#45877f","#a2673f")) +
  scale_fill_manual(values=colorTP) +
  theme_bw()+ 
  theme(plot.margin = margin(.6, 0.4, 0.4, 0, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),  
        legend.key=element_blank(),
        axis.text.y = element_text(size =22, angle=90,hjust=0.5,colour = "black"), 
        axis.text.x = element_text(size = 22,colour = "black"), 
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 25, colour = "black"),
        legend.title = element_text(size=25),
        legend.text=element_text(size = 22),
        legend.position = "bottom")+ 
 labs(colour="Treatment", shape="Time Since Fire (days)")+
  labs(x = "Ash Depth (cm)",y="Species Richness");SppAshFun
detach(MetaRareFun)



attach(MetaRareBac)
SppAshBac<-ggplot(MetaRareBac, aes(x=InitialAshDepth, y=S.obs))  +
  geom_smooth(method = MASS::glm.nb, se = TRUE)+ 
  scale_shape_manual(values=c(15,24,19,1,18,25,8,17,13))+
  geom_point(size=2,aes(shape = as.factor(TSFdays), colour = Treatment))+
  stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label..,
        sep = "~~~~")), formula = y~x, label.x = 5, label.y =1500, size=6)+
  scale_color_manual(values=  c("#45877f","#a2673f")) +
  scale_fill_manual(values=colorTP) +
  theme_bw()+ 
  theme(plot.margin = margin(.6, 0.4, 0.4, 0, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),  
        legend.key=element_blank(),
        axis.text.y = element_text(size =22, angle=90,hjust=0.5,colour = "black"), 
        axis.text.x = element_text(size = 22,colour = "black"), 
        axis.title.y = element_text(size = 25, colour = "black"),
        axis.title.x = element_text(size = 25, colour = "black"),
        legend.title = element_text(size=25),
        legend.text=element_text(size = 22),
        legend.position = "bottom")+ 
  labs(colour="Treatment", shape="Time Since Fire (days)")+
  labs(x = "Ash Depth (cm)",y="Species Richness");SppAshBac
detach(MetaRareBac)



#BIOMASS
attach(MetaBioBac)
SppAshBioBac<-ggplot(MetaBioBac, aes(x=InitialAshDepth, y=CopyNumber))  +
  geom_smooth(method = MASS::glm.nb, se = TRUE)+ 
  scale_shape_manual(values=c(15,24,19,1,18,25,8,17,13))+
  geom_point(size=2,aes(shape = as.factor(TSFdays), colour = Treatment))+
  stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label..,
           sep = "~~~~")), formula = y~x, label.x = 3, label.y =1550000, size=6)+
  scale_color_manual(values=  c("#45877f","#a2673f")) +
  scale_fill_manual(values=colorTP) +
  theme_bw()+ 
  theme(plot.margin = margin(1.5, 0.4, 0, 0, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),  
        legend.key=element_blank(),
        axis.text.y = element_text(size =18, angle=90,hjust=0.5,colour = "black"), 
        axis.text.x = element_text(size = 22,colour = "black"), 
        axis.title.y = element_text(size = 25, colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_text(size=25),
        legend.text=element_text(size = 22),
        legend.position = "bottom")+ 
  labs(colour="Treatment", shape="Time Since Fire (days)")+
  labs(x = "Ash Depth (cm)",y="Gene Copy Number");SppAshBioBac
detach(MetaBioBac)



attach(MetaBioFun)
SppAshBioFun<-ggplot(MetaBioFun, aes(x=InitialAshDepth, y=CopyNumber))  +
  geom_smooth(method = MASS::glm.nb, se = TRUE)+ 
  scale_shape_manual(values=c(15,24,19,1,18,25,8,17,13))+
  geom_point(size=2,aes(shape = as.factor(TSFdays), colour = Treatment))+
  stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label..,
       sep = "~~~~")), formula = y~x, label.x = 3, label.y =600000, size=6)+
  scale_color_manual(values=  c("#45877f","#a2673f"))+
  scale_fill_manual(values=colorTP) +
  theme_bw()+ 
  theme(plot.margin = margin(1.5, 0.4, 0, 0, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),  
        legend.key=element_blank(),
        axis.text.y = element_text(size =22, angle=90,hjust=0.5,colour = "black"), 
        axis.text.x = element_text(size = 22,colour = "black"), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_text(size=25),
        legend.text=element_text(size = 22),
        legend.position = "bottom")+ 
  labs(colour="Treatment", shape="Time Since Fire (days)")+
  labs(x = "Ash Depth (cm)",y="Gene Copy Number");SppAshBioFun
detach(MetaBioFun)




#Create panel
SppAshPlots<-ggarrange(SppAshBioBac,SppAshBioFun,SppAshBac, SppAshFun, ncol=2, nrow=2, 
              common.legend = TRUE, legend="bottom", hjust = c(-.8,-1),align="hv",
              labels = c("A Bacteria","B Fungi","C Bacteria","D Fungi"), 
              font.label = list(size=24,face="bold"));SppAshPlots

pdf("1a-PanelsBac-Fun/Graphs/SppBio-AhDepth-Bac-FunNEW-Legend.pdf", height=12, width=12)
SppAshPlots
dev.off()

#----
#----
#*****************************************************************************************************-----
#------------------------- CALCULATE SIGNIFICANCE --------------------------------------------------------
#*****************************************************************************************************-----
#cannot normalize, log, sqrt, cubed

hist(log(MetaBioBac$CopyNumber))
shapiro.test(sqrt(MetaBioBac$CopyNumber))

hist(log(MetaBioFun$CopyNumber))
shapiro.test(log(MetaBioFun$CopyNumber))

hist((BioBacB$CopyNumber)^3)
shapiro.test((BioBacB$CopyNumber)^3)

hist((BioFunB$CopyNumber)^3)
shapiro.test((BioFunB$CopyNumber)^3)




#Check for data normality and over disperssion..................................................
#Checking mean and variance of the data 

#Bacerial data..............................................................................
BacBioVar <-MetaBioBac %>%#metadata or your file
  group_by(Treatment) %>% #the variable you are interested in (what you are grouping by)
  summarise(means = mean(CopyNumber),
            variance = var(CopyNumber));BacBioVar 


BacBioVar2 <-MetaBioBac %>%#metadata or your file
  group_by(TSFdays) %>% #the variable you are interested in (what you are grouping by)
  summarise(means = mean(CopyNumber),
            variance = var(CopyNumber));BacBioVar2

#Burned data
BioBacBvar<-BioBacB %>%#metadata or your file
  group_by(TSFdays) %>% #the variable you are interested in (what you are grouping by)
  summarise(means = mean(CopyNumber),
            variance = var(CopyNumber));BioBacBvar 




#Fungal data...............................................................................
FunBioVar <-MetaBioFun %>%#metadata or your file
  group_by(Treatment) %>% #the variable you are interested in (what you are grouping by)
  summarise(means = mean(CopyNumber),
            variance = var(CopyNumber));FunBioVar


FunBioVar2 <-MetaBioFun %>%#metadata or your file
  group_by(TSFdays) %>% #the variable you are interested in (what you are grouping by)
  summarise(means = mean(CopyNumber),
            variance = var(CopyNumber));FunBioVar2


#Burned data
BioFunBvar<-BioFunB %>%#metadata or your file
  group_by(TSFdays) %>% #the variable you are interested in (what you are grouping by)
  summarise(means = mean(CopyNumber),
            variance = var(CopyNumber));BioFunBvar 


  
