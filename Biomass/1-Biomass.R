#T1-T9
#2021

#Reset R's Brain
rm(list=ls())

#set working directory
setwd("~/Dropbox/1-Dissertation/1-Research/HolyFire/1-HolyFire-Pub1/HolyfireT1toT9")

#load libraries
library("tidyverse")
library("scales")
library("gridExtra")
library("ggpubr")
library("ggplot2")

#set par bac to normal
par(mfrow=c(1,1))

#load in data
Metadata <- read.csv("Metadata/Calculating_Copy_Number_HolyFire.csv",na.strings = "NA", header = TRUE)

str(Metadata)

#****************************************************************************************************************-----
#---------------------- QUALITY CONTROL ------------------------------------------------------------------------------
#****************************************************************************************************************-----
Metadata$Plot <- as.factor(Metadata$Plot);levels(Metadata$Plot)
Metadata$TSFdays <- as.factor(Metadata$TSFdays)
Metadata$Treatment <- as.factor(Metadata$Treatment)
Metadata$Timepoint <- as.factor(Metadata$Timepoint)
Metadata$AshSeverity<-ordered(Metadata$AshSeverity, c("Unburned","Low","Moderate", "High"))
Metadata$Date<-ordered(Metadata$Date, c("9/30/2018","10/8/2018","10/17/2018", "11/19/2018",
               "12/17/2018","1/22/2019","3/19/2019","6/26/2019","9/24/2019"))
Metadata$Timepoint<-ordered(Metadata$Timepoint, c("TP1","TP2","TP3", "TP4",
                      "TP5","TP6","TP7", "TP8","TP9"))


#CHANGE BASE LEVEL ...........................................................
Metadata$Treatment <- try(relevel(Metadata$Treatment , "Unburned"))



#..................................................................................................
#-----------SUBSET METADATA -----------------------------------------------------------------------
attach(Metadata)

#Filtering out the extremem outliers---------------------------------------------------
Metadata<- filter(Metadata,CopyNumber < 2000000)#We had 2

MetaFun <- filter(Metadata, Kingdom == "Fungi");levels(MetaFun$Treatment) 
MetaBac<- filter(Metadata, Kingdom == "Bacteria"); levels(MetaBac$Treatment)

#Subset metadata ............,...................................................
MetaBurn<-filter(Metadata, Treatment == "Burned"); head(MetaBurn)

BacBurn<-filter(MetaBac, Treatment == "Burned"); head(BacBurn)
FunBurn<-filter(MetaFun, Treatment == "Burned"); head(FunBurn)

BacUn<-filter(MetaBac, Treatment == "Unburned"); head(BacUn)
FunUn<-filter(MetaFun, Treatment == "Unburned"); head(FunUn)



#Export tables for graphing purposes in future, clean graphs....
write.csv(MetaFun,"Metadata/MetaFun-Biomass.csv")
write.csv(MetaBac,"Metadata/MetaBac-Biomass.csv")

#Burn & Un
write.csv(BacBurn,"Metadata/MetaBac-Biomass-Burn.csv")
write.csv(BacUn,"Metadata/MetaBac-Biomass-Unburn.csv")
write.csv(FunBurn,"Metadata/MetaFun-Biomass-Burn.csv")
write.csv(FunUn,"Metadata/MetaFun-Biomass-Unburn.csv")



#----
#----
#*********************************************************************************************----
# NORMALITY TEST ,...........................................................................
#*********************************************************************************************----

# FULL DATASET ................................................
shapiro.test(Metadata$CopyNumber)#data not normal
hist(Metadata$CopyNumber)

shapiro.test(MetaFun$CopyNumber)
hist(MetaFun$CopyNumber)

shapiro.test(MetaBac$CopyNumber)
hist(MetaBac$CopyNumber)


#Bacteria dataset........................#cannot normalize.
shapiro.test(sqrt(BacBurn$CopyNumber))
hist(log(BacBurn$CopyNumber))

shapiro.test(sqrt(BacUn$CopyNumber))
hist(sqrt(BacUn$CopyNumber))


#Fungi dataset.........................#cannot normalize
shapiro.test(log(FunBurn$CopyNumber))
hist(FunBurn$CopyNumber)

shapiro.test(sqrt(FunUn$CopyNumber))
hist(sqrt(FunUn$CopyNumber))


#----
#----
#*****************************************************************************************************************************----
#DESRIPTIVE STATISTICS----------------------------------------------------------------------------------------------------------
#*****************************************************************************************************************************----

TrtStats<-Metadata %>%
  filter(!is.na(CopyNumber)) %>%
  group_by(Treatment, Kingdom) %>%
  summarise(n_obs = n(),
            TotalGenes =sum(CopyNumber, na.rm = TRUE),
            Avg_Genes = mean(CopyNumber, na.rm = TRUE), 
            Variance = var(CopyNumber, na.rm = TRUE), 
            sd = sd(CopyNumber), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 0) %>% as.data.frame();TrtStats 

TSFstats<-Metadata %>%
  filter(!is.na(CopyNumber)) %>%
  group_by(Treatment, Kingdom, TSFdays) %>%
  summarise(n_obs = n(),
            TotalGenes =sum(CopyNumber, na.rm = TRUE),
            Avg_ASV = mean(CopyNumber, na.rm = TRUE), 
            Variance = var(CopyNumber, na.rm = TRUE), 
            sd = sd(CopyNumber), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 0) %>% as.data.frame();TSFstats


Sevstats<-Metadata %>%
  filter(!is.na(CopyNumber)) %>%
  group_by(Severity, Kingdom) %>%
  summarise(n_obs = n(),
            TotalGenes =sum(CopyNumber, na.rm = TRUE),
            Avg_ASV = mean(CopyNumber, na.rm = TRUE), 
            Variance = var(CopyNumber, na.rm = TRUE), 
            sd = sd(CopyNumber), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 0) %>% as.data.frame();Sevstats 

#FUNGAL ................................................................
TrtStatsFun<-MetaFun%>%
  filter(!is.na(CopyNumber)) %>%
  group_by(Treatment) %>%
  summarise(n_obs = n(),
            TotalGenes =sum(CopyNumber, na.rm = TRUE),
            Avg_Genes = mean(CopyNumber, na.rm = TRUE), 
            Variance = var(CopyNumber, na.rm = TRUE), 
            sd = sd(CopyNumber), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 0) %>% as.data.frame();TrtStatsFun

TSFstatsFun<-MetaFun %>%
  filter(!is.na(CopyNumber)) %>%
  group_by(TSFdays,Treatment) %>%
  summarise(n_obs = n(),
            TotalGenes =sum(CopyNumber, na.rm = TRUE),
            Avg_ASV = mean(CopyNumber, na.rm = TRUE), 
            Variance = var(CopyNumber, na.rm = TRUE), 
            sd = sd(CopyNumber), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 0) %>% as.data.frame();TSFstatsFun

TSFstatsFunB<-FunBurn %>%
  filter(!is.na(CopyNumber)) %>%
  group_by(TSFdays) %>%
  summarise(n_obs = n(),
            TotalGenes =sum(CopyNumber, na.rm = TRUE),
            Avg_ASV = mean(CopyNumber, na.rm = TRUE), 
            Variance = var(CopyNumber, na.rm = TRUE), 
            sd = sd(CopyNumber), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 0) %>% as.data.frame();TSFstatsFunB 

SevstatsFunB<-FunBurn%>%
  filter(!is.na(CopyNumber)) %>%
  group_by(Severity) %>%
  summarise(n_obs = n(),
            TotalGenes =sum(CopyNumber, na.rm = TRUE),
            Avg_ASV = mean(CopyNumber, na.rm = TRUE), 
            Variance = var(CopyNumber, na.rm = TRUE), 
            sd = sd(CopyNumber), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 0) %>% as.data.frame();SevstatsFunB 


#Bacteria .................................................................
TrtStatsBac<-MetaBac%>%
  filter(!is.na(CopyNumber)) %>%
  group_by(Treatment) %>%
  summarise(n_obs = n(),
            TotalGenes =sum(CopyNumber, na.rm = TRUE),
            Avg_Genes = mean(CopyNumber, na.rm = TRUE), 
            Variance = var(CopyNumber, na.rm = TRUE), 
            sd = sd(CopyNumber), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 0) %>% as.data.frame();TrtStatsBac

TSFstatsBac<-MetaBac %>%
  filter(!is.na(CopyNumber)) %>%
  group_by(TSFdays,Treatment) %>%
  summarise(n_obs = n(),
            TotalGenes =sum(CopyNumber, na.rm = TRUE),
            Avg_ASV = mean(CopyNumber, na.rm = TRUE), 
            Variance = var(CopyNumber, na.rm = TRUE), 
            sd = sd(CopyNumber), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 0) %>% as.data.frame();TSFstatsBac

TSFstatsBacB<-BacBurn %>%
  filter(!is.na(CopyNumber)) %>%
  group_by(TSFdays) %>%
  summarise(n_obs = n(),
            TotalGenes =sum(CopyNumber, na.rm = TRUE),
            Avg_ASV = mean(CopyNumber, na.rm = TRUE), 
            Variance = var(CopyNumber, na.rm = TRUE), 
            sd = sd(CopyNumber), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 0) %>% as.data.frame();TSFstatsBacB

SevstatsBacB<-BacBurn%>%
  filter(!is.na(CopyNumber)) %>%
  group_by(Severity) %>%
  summarise(n_obs = n(),
            TotalGenes =sum(CopyNumber, na.rm = TRUE),
            Avg_ASV = mean(CopyNumber, na.rm = TRUE), 
            Variance = var(CopyNumber, na.rm = TRUE), 
            sd = sd(CopyNumber), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 0) %>% as.data.frame();SevstatsBacB


# EXPORT FILES ......................................................................................
dir.create(file.path("Analysis/Tables/DescriptiveStats/Bac"), recursive = TRUE)
dir.create(file.path("Analysis/Tables/DescriptiveStats/Fun"), recursive = TRUE)
dir.create(file.path("Analysis/Graphs/Fun"), recursive = TRUE)
dir.create(file.path("4Analysis/Graphs/Bac"), recursive = TRUE)

#DESCRIPTIVE STATS.......................................................................................................
#Full Metadata
write.csv(TrtStats, file="Analysis/Tables/DescriptiveStats/Results-Biomass-DesStats-TrtAll.csv") 
write.csv(TSFstats, file="Analysis/Tables/DescriptiveStats/Results-Biomass-DesStats-TSFall.csv") 
write.csv(Sevstats, file="Analysis/Tables/DescriptiveStats/Results-Biomass-DesStats-SEVall.csv") 

#Fungal only
write.csv(TrtStatsFun, file="Analysis/Tables/DescriptiveStats/Fun/Results-Biomass-DesStats-Trt.csv") 
write.csv(TSFstatsFun, file="Analysis/Tables/DescriptiveStats/Fun/Results-Biomass-DesStats-TSF-Trt.csv") 
write.csv(TSFstatsFunB, file="Analysis/Tables/DescriptiveStats/Fun/Results-Biomass-DesStats-TSF.csv") 
write.csv(SevstatsFunB, file="Analysis/Tables/DescriptiveStats/Fun/Results-Biomass-DesStats-Sev.csv") 

#Bacteria only
write.csv(TrtStatsBac, file="Analysis/Tables/DescriptiveStats/Bac/Results-Biomass-DesStats-Trt.csv") 
write.csv(TSFstatsBac, file="Analysis/Tables/DescriptiveStats/Bac/Results-Biomass-DesStats-TSF-Trt.csv") 
write.csv(TSFstatsBacB, file="Analysis/Tables/DescriptiveStats/Bac/Results-Biomass-DesStats-TSF.csv") 
write.csv(SevstatsBacB, file="Analysis/Tables/DescriptiveStats/Bac/Results-Biomass-DesStats-Sev.csv")






#----
#----
#----
#******************************************************************************************************************************----
#CALCULATE PERCENT CHANGE BETWEEN TREATMENTS AND TIMEPOINTS------------------------------------------------------------------------
#******************************************************************************************************************************----
Trt <-  Metadata%>% 
  group_by( Kingdom,Treatment) %>% 
  summarise(mean = mean(CopyNumber,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();Trt 

Sev<-Metadata%>% 
  group_by(Kingdom,AshSeverity) %>% 
  summarise(mean = mean(CopyNumber,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();Sev


TrtTSF <-  Metadata%>% 
  group_by(Kingdom,TSFdays,Treatment) %>% 
  summarise(mean = mean(CopyNumber,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();TrtTSF

#Fungi.....................................
TrtTSFfun<-  MetaFun%>% 
  group_by(TSFdays,Treatment) %>% 
  summarise(mean = mean(CopyNumber,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();TrtTSFfun

TsfBFun <- FunBurn%>% 
  group_by(TSFdays) %>% 
  summarise(mean = mean(CopyNumber,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();TsfBFun

TsfUnFun <- FunUn%>% 
  group_by(TSFdays) %>% 
  summarise(mean = mean(CopyNumber,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();TsfUnFun


#Bacteria.....................................
TrtTSFbac <-  MetaBac%>% 
  group_by(TSFdays,Treatment) %>% 
  summarise(mean = mean(CopyNumber,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();TrtTSFbac

TsfBbac <- BacBurn%>% 
  group_by(TSFdays) %>% 
  summarise(mean = mean(CopyNumber,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();TsfBbac

TsfUnBac <- BacUn%>% 
  group_by(TSFdays) %>% 
  summarise(mean = mean(CopyNumber,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();TsfUnBac



#Create directory........................................................................
dir.create(file.path("Analysis/Tables/PercentChange/Bac"), recursive = TRUE)
dir.create(file.path("Analysis/Tables/PercentChange/Fun"), recursive = TRUE)

#Percent change All (bac and fun together).....................................................................................
write.csv(Trt,"Analysis/Tables/PercentChange/PerCh-Trt-All.csv")
write.csv(Sev, "Analysis/Tables/PercentChange/PerCh-TSF-All.csv")
write.csv(TrtTSF,"Analysis/Tables/PercentChange/PerCh-TRT-TST-All.csv")

#Percent Bacteria ....................................................
write.csv(TrtTSFbac,"Analysis/Tables/PercentChange/Bac/PerCh-Trt-Tsf-Bac.csv")
write.csv(TsfBbac,"Analysis/Tables/PercentChange/Bac/PerCh-Tsf-Burned-Bac.csv")
write.csv(TsfUnBac,"Analysis/Tables/PercentChange/Bac/PerCh-Tsf-Unburned-Bac.csv")


#Percent Fungi....................................................
write.csv(TrtTSFfun,"Analysis/Tables/PercentChange/Fun/PerCh-Trt-Tsf-Fun.csv")
write.csv(TsfBFun,"Analysis/Tables/PercentChange/Fun/PerCh-Tsf-Burned-Fun.csv")
write.csv(TsfUnFun,"Analysis/Tables/PercentChange/Fun/PerCh-Tsf-Unburned-Fun.csv")






#----
#----
#----
#***********************************************************************************************************----
#-------------------------- CREATE GRAPHS ---------------------------------------------------------------------
#***********************************************************************************************************----
#RED COLORS "#a2673f","#45877f"
#Figure 1 Burned v Unburned Plots

attach(Metadata)
mycol<-c("#40271f","#5e4839","#947f70","#c49f68","#bfb3aa","#ebddc7","#235978","#45877f","#aed4ca")

TSFtrtPlot<-ggplot(Metadata, aes(x=Timepoint, y=CopyNumber, group=Treatment, col=Treatment, shape=Treatment)) +
  stat_summary(fun=mean,geom="line", size=1, aes(linetype = Treatment),show.legend=TRUE)+
  stat_summary(fun.data = mean_se,geom = "errorbar",size=.85,alpha=0.7,position = position_dodge(0.01))+
  scale_linetype_manual(values=c("solid","twodash"))+
  scale_color_manual(values=c("#45877f","#a2673f"))+ 
  theme_bw()+ 
  theme(panel.grid = element_blank(),
        legend.key=element_blank(),
        axis.text.y = element_text(size =14,angle=90,hjust=0.5,colour = "black"), 
        axis.text.x = element_text(size = 14,colour = "black"), 
        axis.title.y = element_text(size = 16, colour = "black"), 
        axis.title.x = element_text(size = 16, colour = "black"),
        legend.title = element_text(size=15),
        legend.text=element_text(size = 14),legend.position = "bottom")+ 
  labs(x = "Time Since Fire (days)",y="Gene Copy Number");TSFtrtPlot


(TSFtrtPlot2<- TSFtrtPlot + facet_wrap( ~ Kingdom, scale = "free", ncol = 1)+
    theme(strip.text = element_text(size=24)))



#SINGLE LINE TREATMENT AND TIMESINCEFIRE....................................................................
TSFtrt<-ggplot(Metadata, aes(x=Timepoint, y=CopyNumber,
              group=Treatment, col=Treatment, label="Sample")) +
stat_summary(fun=mean,geom="line", size=1, aes(linetype = Treatment),show.legend=TRUE)+
  stat_summary(fun.data = mean_se,geom = "errorbar",size=.85, 
               alpha=0.7,position = position_dodge(0.01))+
  scale_linetype_manual(values=c("solid","twodash"))+
  scale_color_manual(values=c("#45877f","#a2673f"))+
  theme_bw()+ 
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),  
        legend.key=element_blank(),
        axis.text.y = element_text(size =14, angle=90, hjust=0.5, colour = "black"), 
        axis.text.x = element_text(size = 14, colour = "black"), 
        axis.title.y = element_text(size = 16, colour = "black"), 
        axis.title.x = element_text(size = 16, colour = "black"),
        legend.title = element_text(size=18),
        legend.text=element_text(size = 18),
        legend.position = "bottom")+ 
  labs(x = "Time Since Fire (days)",y="Gene Copy Number") +
  scale_y_continuous(labels = scientific);TSFtrt


#Add  kruskal.wallis signfiicance at each time point
TSFtrtStat<- TSFtrt + stat_compare_means(method = "kruskal.test", 
             label="p.signif", label.y=650000, 
             show.legend = FALSE );TSFtrtStat 


(TSFtrtStat2<- TSFtrtStat + facet_wrap( ~ Kingdom, scale = "free")+
    theme(strip.text = element_text(size=24)));TSFtrtStat2
        
pdf("Analysis/Graphs/TSF-Trt-BacFun-TP1-9.pdf", height=7, width=12)
TSFtrtStat2
dev.off()

#TSF SOIL BURN SEVERITY FOR BURN DATA ONLY .........................................
Sub<-ggboxplot(FunUn, x="Subplot", y="CopyNumber") +   
  geom_boxplot(aes(fill=Treatment))+ 
  theme_bw() +  ylab("Species Richness") + 
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        plot.margin = margin(0.3, 1.5, 0.3, 0.3, "cm"),
        legend.position = "bottom",
        text = element_text(size=20),
        axis.title.x=element_blank(),
        axis.text.y = element_text(size=16,angle=90, hjust=.5,colour = "black"),
        axis.text.x = element_text(size=16, colour = "black"))+ 
  scale_fill_manual(values=c("#45877f","#a2673f"));Sub


(Sub2<- Sub  + facet_wrap( ~ Timepoint, scale = "free")+
    theme(strip.text = element_text(size=20)));Sub2

pdf("Analysis/Graphs/SubplotCopyNum.pdf", height=8, width=10)
Sub2
dev.off()









#SOIL BURN SEVERITY ....................................................................
TSFSevAll<-ggplot(Metadata, aes(x=TSFdays, y=CopyNumber,
           group=AshSeverity, col=AshSeverity, label="Sample")) +
  stat_summary(fun=mean,geom="line", size=1,show.legend=TRUE)+
  stat_summary(fun.data = mean_se,geom = "errorbar",size=.85, 
               alpha=0.7,position = position_dodge(0.01))+
  #scale_linetype_manual(values=c("solid","twodash"))+
  scale_color_manual(values=c("#027300","#93b396","#d9cd2c","#8B1919"))+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),  
        legend.key=element_blank(),
        axis.text.y = element_text(size =18, angle=90, hjust=0.5, colour = "black"), 
        axis.text.x = element_text(size = 18, colour = "black"), 
        axis.title.y = element_text(size = 26, colour = "black"), 
        axis.title.x = element_text(size = 26, colour = "black"),
        legend.title = element_text(size=22),
        legend.text=element_text(size = 18),
        legend.position = "bottom")+ 
  labs(x = "Time Since Fire (days)",y="Gene Copy Number") +
  scale_y_continuous(labels = scientific);TSFSevAll

(TSFSevAll2<- TSFSevAll  + facet_wrap( ~ Kingdom, scale = "free")+
    theme(strip.text = element_text(size=24)));TSFSevAll2


  Sevbp<-ggboxplot(Metadata, x="AshSeverity", y="CopyNumber") +   
  geom_boxplot(aes(fill=AshSeverity))+ 
  theme_bw() +  
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.position = "bottom",
        text = element_text(size=30),
        axis.title.x=element_blank(),
        axis.text.y = element_text(size=25, angle=90, hjust=.5, colour = "black"),
        axis.text.x = element_text(size=25, colour = "black"))+ 
  scale_fill_manual(values=c("#027300","#93b396","#d9cd2c","#8B1919"))+
  ylab("Gene Copy Number")+
    scale_y_continuous(labels = scientific);Sevbp







#----
#----
#******************************************************************************************----
#Export graphs ..................................................................................
#******************************************************************************************----
pdf("Analysis/Graphs/Biomass-BacFun_TrtTsf-Plot.pdf", height=8, width=14)
TSFtrtPlot2
dev.off()


pdf("Analysis/Graphs/Biomass-BacFun_KW-TrtTsf-TP13.pdf", height=8, width=14)
TSFtrtStat2
dev.off()

pdf("Analysis/Graphs/Biomass-BacFun_SevTSF.pdf", height=8, width=14)
TSFSev2
dev.off()

pdf("Analysis/Graphs/Biomass-BacFun_SevTSF-All.pdf", height=8, width=14)
TSFSevAll2
dev.off()



#***********************************************************************************************************************
#EXTRA GRAPHS WITH ADDED DATA POINTS JUST TO LOOK AT -----------------------------------------------------------
#Fig2, similar plot to Fig1 except adding in the data points with transparency.........................................
Fig2 <- ggplot(Metadata, aes(x=TSFdays, y=CopyNumber, group=Treatment, col=Treatment, label="Sample")) +
  stat_summary(fun=mean,geom="line", size=1, aes(linetype = Treatment),show.legend=TRUE)+
  stat_summary(fun.data = mean_se,geom = "errorbar",size=.85, 
               alpha=0.7,position = position_dodge(0.01))+
  geom_point(size=1, alpha=0.5)+ #Adds in the data points
  scale_linetype_manual(values=c("solid","twodash"))+
  scale_color_manual(values=c("#45877f","#a2673f"))+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),  
        legend.key=element_blank(),
        plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"),
        axis.text.y = element_text(size =18, angle=90, hjust=0.5, colour = "black"), 
        axis.text.x = element_text(size = 18, colour = "black",
              margin = margin(t = 0, r = 0, b = 20, l = 0)), 
        axis.title.y = element_text(size = 26, colour = "black",
              margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(size = 26, colour = "black"),
        legend.title = element_text(size=22),
        legend.text=element_text(size = 18),
        legend.position = "bottom")+ 
  labs(x = "Time Since Fire (days)", y="Copy Number");Fig2

#add in the significance
Fig2_Stat <- Fig2 + stat_compare_means(method = "kruskal.test", label="p.signif");Fig2_Stat

# This shows the variation between the burned and unburnded plots, free scales
Fig2_Facet <- Fig2_Stat +facet_wrap( ~ Kingdom, scale = "free")+
  theme(strip.text = element_text(size=24));Fig2_Facet



#Export graphs ..................................................................................................
dir.create("Analysis/Graphs")

pdf("Analysis/Graphs/Biomass-BacFun_KW-points.pdf", height=6, width=8)
Fig2_Stat
dev.off()

pdf("Analysis/Graphs/HolyFireQPCR-BacFun-Panels-points.pdf", height=8, width=12)
Fig2_Facet
dev.off()

  
#----  
#----
#*********************************************************************************************----
#MetaFun=-------------------------------------------------------------------------------------
#*********************************************************************************************----
TSFtrtFun<-ggplot(MetaFun, aes(x=TSFdays, y=CopyNumber, group=Plot, col=Plot, shape=Treatment, label="Sample")) +
  stat_summary(fun=mean,geom="line", size=1, aes(linetype = Treatment),show.legend=TRUE)+
  stat_summary(fun.data = mean_se,geom = "errorbar",size=.85,alpha=0.7,position = position_dodge(0.01))+
  scale_linetype_manual(values=c("solid","twodash"))+
  scale_color_manual(values=c("#40271f","#5e4839","#947f70",
       "#c49f68","#bfb3aa","#ebddc7","#235978","#45877f","#aed4ca"))+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),  
        legend.key=element_blank(),
        axis.text.y = element_text(size =17, angle=90, hjust=0.9, colour = "black"), 
        axis.text.x = element_text(size = 20, colour = "black"), 
        axis.title.y = element_text(size = 22, colour = "black"), 
        axis.title.x = element_text(size = 22, colour = "black"),
        legend.title = element_text(size=22),
        legend.text=element_text(size = 18),
        legend.position = "right")+ 
  labs(x = "",y="");TSFtrtFun

#MetaBac
TSFtrtBac<-ggplot(MetaBac, aes(x=TSFdays, y=CopyNumber, group=Plot, col=Plot, shape=Treatment, label="Sample")) +
  stat_summary(fun=mean,geom="line", size=1, aes(linetype = Treatment),show.legend=TRUE)+
  stat_summary(fun.data = mean_se,geom = "errorbar",size=.85,alpha=0.7,position = position_dodge(0.01))+
  scale_linetype_manual(values=c("solid","twodash"))+
  scale_color_manual(values=c("#40271f","#5e4839","#947f70",
         "#c49f68","#bfb3aa","#ebddc7","#235978","#45877f","#aed4ca"))+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),  
        legend.key=element_blank(),
        axis.text.y = element_text(size =17, angle=90, hjust=0.9, colour = "black"), 
        axis.text.x = element_text(size = 20, colour = "black"), 
        axis.title.y = element_text(size = 22, colour = "black"), 
        axis.title.x = element_text(size = 22, colour = "black"),
        legend.title = element_text(size=22),
        legend.text=element_text(size = 18),
        legend.position = "right")+ 
  labs(x = "",y="Gene Copy Number");TSFtrtBac


TSFtrt<-ggarrange(TSFtrtBac,TSFtrtFun, ncol=2, nrow=1, common.legend = TRUE, 
                  legend="right", labels = c("a Bacteria", "b Fungi"), hjust = c(-.8,-1),
                  align="hv",font.label = list(size=16,face="bold"));TSFtrt



pdf("Analysis/Graphs/Biomass-BacFun-Panels-TSF-trt.pdf", height=6, width=13)
TSFtrt
dev.off()


#PANELS WITH SPECIES RICHNESS
TSFtrt<-ggarrange(TSFtrtBac,TSFtrtFun,SppTSFBac,SppTSFfun,ncol=2,nrow = 2, 
                  common.legend = TRUE, legend="right", align="hv",
                  labels = c("a Bacteria", "b Fungi","c Bacteria", "c Fungi"),
                  hjust = c(-.8,-1),font.label = list(size=16,face="bold"));TSFtrt



TSFtrt2<-ggarrange(TSFtrtBac,TSFtrtFun,SppTSFBac,SppTSFfun,ncol=2,nrow = 2, 
                  common.legend = TRUE, legend="bottom", align="hv",
                  labels = c("a Bacteria", "b Fungi","c Bacteria", "c Fungi"),
                  hjust = c(-.8,-1),font.label = list(size=16,face="bold"));TSFtrt2



pdf("1a-PanelsBac-Fun/Graphs/Diversity/BiomassRichness-BacFun.pdf", height=12, width=14)
TSFtrt
dev.off()


pdf("1a-PanelsBac-Fun/Graphs/Diversity/BiomassRichness-BacFun2.pdf", height=12, width=12)
TSFtrt2
dev.off()


#PANELS WITH SPECIES RICHNESS--- 1 LINE PER TREATMENT
TSFtrt<-ggarrange(TSFtrtBac,TSFtrtFun,SppTSFBac,SppTSFfun,ncol=2,nrow = 2, 
                  common.legend = TRUE, legend="right", align="hv",
                  labels = c("a Bacteria", "b Fungi","c Bacteria", "c Fungi"),
                  hjust = c(-.8,-1),font.label = list(size=16,face="bold"));TSFtrt




TSFtrt2<-ggarrange(TSFtrtBac,TSFtrtFun,SppTSFBac,SppTSFfun,ncol=2,nrow = 2, 
                   common.legend = TRUE, legend="bottom", align="hv",
                   labels = c("a Bacteria", "b Fungi","c Bacteria", "c Fungi"),
                   hjust = c(-.8,-1),font.label = list(size=16,face="bold"));TSFtrt2



pdf("1a-PanelsBac-Fun/Graphs/Diversity/BiomassRichness-BacFun.pdf", height=12, width=14)
TSFtrt
dev.off()


pdf("1a-PanelsBac-Fun/Graphs/Diversity/BiomassRichness-BacFun2.pdf", height=12, width=12)
TSFtrt2
dev.off()






#----  
#----
#*********************************************************************************************----
# one line per graphs-------------------------------------------------------------------------------------
#*********************************************************************************************----
TSFFun<-ggplot(MetaFun, aes(x=TSFdays, y=CopyNumber,
        group=Treatment, col=Treatment, label="Sample")) +
  stat_summary(fun=mean,geom="line", size=1, aes(linetype = Treatment),show.legend=TRUE)+
  stat_summary(fun.data = mean_se,geom = "errorbar",size=.85, 
               alpha=0.7,position = position_dodge(0.01))+
  scale_linetype_manual(values=c("solid","twodash"))+
  scale_color_manual(values=c("#45877f","#a2673f"))+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),  
        legend.key=element_blank(),,
        axis.text.y = element_text(size =18, angle=90, hjust=0.5, colour = "black"), 
        axis.text.x = element_text(size = 18, colour = "black",vjust=0.3), 
        axis.title.y = element_text(size = 26, colour = "black"), 
        axis.title.x = element_text(size = 26, colour = "black"),
        legend.title = element_text(size=22),
        legend.text=element_text(size = 18),
        legend.position = "bottom")+ 
  labs(x = "",y="") +
  scale_y_continuous(labels = scientific);TSFFun


#Add  kruskal.wallis signfiicance at each time point
TSFFun2<- TSFFun + stat_compare_means(method = "kruskal.test", 
              label="p.signif", label.y=3e05, 
                show.legend = FALSE, size=8);TSFFun2

#MetaBac
TSFBac<-ggplot(MetaBac, aes(x=TSFdays, y=CopyNumber,
      group=Treatment, col=Treatment, label="Sample")) +
  stat_summary(fun=mean,geom="line", size=1, aes(linetype = Treatment),show.legend=TRUE)+
  stat_summary(fun.data = mean_se,geom = "errorbar",size=.85, 
               alpha=0.7,position = position_dodge(0.01))+
  scale_linetype_manual(values=c("solid","twodash"))+
  scale_color_manual(values=c("#45877f","#a2673f"))+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),  
        legend.key=element_blank(),
        axis.text.y = element_text(size =18, angle=90, hjust=0.5, colour = "black"), 
        axis.text.x = element_text(size = 18, colour = "black",vjust=0.3), 
        axis.title.y = element_text(size = 26, colour = "black"), 
        axis.title.x = element_text(size = 26, colour = "black"),
        legend.title = element_text(size=22),
        legend.text=element_text(size = 18),
        legend.position = "bottom")+ 
  labs(x = "",y="Gene Copy Number") +
  scale_y_continuous(labels = scientific);TSFBac


#Add  kruskal.wallis signfiicance at each time point
TSFBac2<- TSFBac + stat_compare_means(method = "kruskal.test", 
             label="p.signif", label.y=9e05, 
              show.legend = FALSE, size=8);TSFBac2


TSF<-ggarrange(TSFBac2,TSFFun2, ncol=2, nrow=1, common.legend = TRUE, 
        legend="right", labels = c("a Bacteria", "b Fungi"), 
        hjust = c(-.8,-1),align="hv",font.label = list(size=16,
         face="bold"));TSF



pdf("Analysis/Graphs/Biomass-BacFun-Panels-TSF.pdf", height=6, width=14)
TSF
dev.off()





#PANELS WITH SPECIES RICHNESS
TSFBacFun<-ggarrange(TSFBac2,TSFFun2,SppTSF2Bac2 ,SppTSF2Fun2, ncol=2,nrow = 2, 
            common.legend = TRUE, legend="right", align="hv",
            labels = c("a Bacteria", "b Fungi","c Bacteria", "c Fungi"),
            hjust = c(-.8,-1),font.label = list(size=16,face="bold"));TSFBacFun



TSFBacFun2<-ggarrange(TSFBac2,TSFFun2,SppTSF2Bac2 ,SppTSF2Fun2,ncol=2,nrow = 2, 
            common.legend = TRUE, legend="bottom", align="hv",
            labels = c("a Bacteria", "b Fungi","c Bacteria", "c Fungi"),
             hjust = c(-.8,-1),font.label = list(size=16,face="bold"));TSFBacFun2



pdf("1a-PanelsBac-Fun/Graphs/Diversity/BiomassRichness-BacFun-1Line.pdf", height=14, width=16)
TSFBacFun
dev.off()


pdf("1a-PanelsBac-Fun/Graphs/Diversity/BiomassRichness-BacFun2-1line.pdf", height=14, width=14)
TSFBacFun2
dev.off()


#PANELS WITH SPECIES RICHNESS--- 1 LINE PER TREATMENT
TSFtrt<-ggarrange(TSFtrtBac,TSFtrtFun,SppTSFBac,SppTSFfun,ncol=2,nrow = 2, 
                  common.legend = TRUE, legend="right", align="hv",
                  labels = c("a Bacteria", "b Fungi","c Bacteria", "c Fungi"),
                  hjust = c(-.8,-1),font.label = list(size=16,face="bold"));TSFtrt



TSFtrt2<-ggarrange(TSFtrtBac,TSFtrtFun,SppTSFBac,SppTSFfun,ncol=2,nrow = 2, 
                   common.legend = TRUE, legend="bottom", align="hv",
                   labels = c("a Bacteria", "b Fungi","c Bacteria", "c Fungi"),
                   hjust = c(-.8,-1),font.label = list(size=16,face="bold"));TSFtrt2



pdf("1a-PanelsBac-Fun/Graphs/Diversity/BiomassRichness-BacFun.pdf", height=12, width=14)
TSFtrt
dev.off()


pdf("1a-PanelsBac-Fun/Graphs/Diversity/BiomassRichness-BacFun2.pdf", height=12, width=12)
TSFtrt2
dev.off()