#T1-T9 #2021

#Reset R's Brain
rm(list=ls())

#set working directory
setwd("~/Dropbox/1-Dissertation/1-Research/HolyFire/1-HolyFire-Pub1/1c-Qpcr/Year1")

#load libraries
library("tidyverse")
library("scales")

#load in data
Metadata <- read.csv("Analysis/Metadata/CopyNumberHolyFireNormalized.csv",na.strings = "NA", header = TRUE)


#****************************************************************************************************************-----
#---------------------- QUALITY CONTROL ------------------------------------------------------------------------------
#****************************************************************************************************************-----
attach(Metadata)
Metadata$Plot <- as.factor(Metadata$Plot);levels(Metadata$Plot)
Metadata$TSFdays <- as.factor(Metadata$TSFdays)
Metadata$Treatment <- as.factor(Metadata$Treatment)
Metadata$Timepoint <- as.factor(Metadata$Timepoint)
Metadata$AshSeverity<-ordered(Metadata$AshSeverity, c("Unburned","Low","Moderate", "High"))
Metadata$Date<-ordered(Metadata$Date, c("9/30/2018","10/8/2018","10/17/2018", "11/19/2018",
               "12/17/2018","1/22/2019","3/19/2019","6/26/2019","9/24/2019"))
Metadata$Timepoint<-ordered(Metadata$Timepoint, c("TP1","TP2","TP3", "TP4",
                      "TP5","TP6","TP7", "TP8","TP9"))



#Filter out the extreme outliers..................................................we had 2 crazy huge numbers
Metadata<- filter(Metadata,CopyNumber.ul < 2000000)

#Subet by kingdom to separate bac vs. fun
MetaBac<-filter(Metadata, Kingdom == "Bacteria"); head(MetaBac);dim(MetaBac)#303x26
MetaFun<-filter(Metadata, Kingdom == "Fungi"); head(MetaFun);dim(MetaFun)#306x26

#Subset by treatment....................................................................
FunBurn<-filter(MetaFun, Treatment == "Burned"); head(FunBurn);dim(FunBurn)#205x26
FunUn<-filter(MetaFun, Treatment == "Unburned"); head(MetaFun);dim(FunUn)#98x26

BacBurn<-filter(MetaBac, Treatment == "Burned"); head(BacBurn);dim(BacBurn)#205x26
BacUn<-filter(MetaBac, Treatment == "Unburned"); head(BacUn);dim(BacUn)#101x26

#Export tables for graphing purposes in future, clean graphs.............................
write.csv(Metadata,"Analysis/Metadata/MetadataFilter.csv")
write.csv(MetaFun,"Analysis/Metadata/MetaFun-Biomass.csv")
write.csv(MetaBac,"Analysis/Metadata/MetaBac-Biomass.csv")
#.................................................................Burn & Un
write.csv(BacBurn,"Analysis/Metadata/MetaBac-Biomass-Burn.csv")
write.csv(BacUn,"Analysis/Metadata/MetaBac-Biomass-Unburn.csv")
write.csv(FunBurn,"Analysis/Metadata/MetaFun-Biomass-Burn.csv")
write.csv(FunUn,"Analysis/Metadata/MetaFun-Biomass-Unburn.csv")



#----
#----
#*********************************************************************************************----
# NORMALITY TEST ,...........................................................................
#*********************************************************************************************----
#Use normalize dataset (per gram soil) for all analysis
# Full dataset ......................................................
shapiro.test(Metadata$CopyNumberPerGram)#data not normal
hist(Metadata$CopyNumberPerGram)

shapiro.test(MetaFun$CopyNumberPerGram)
hist(MetaFun$CopyNumberPerGram)

shapiro.test(MetaBac$CopyNumberPerGram)
hist(MetaBac$CopyNumberPerGram)

#Bacteria dataset...............................#cannot normalize.
shapiro.test(sqrt(BacBurn$CopyNumberPerGram))
hist(log(BacBurn$CopyNumberPerGram))

shapiro.test(sqrt(BacUn$CopyNumberPerGram))
hist(sqrt(BacUn$CopyNumberPerGram))

#Fungi dataset..................................#cannot normalize
shapiro.test(log(FunBurn$CopyNumberPerGram))
hist(FunBurn$CopyNumberPerGram)

shapiro.test(sqrt(FunUn$CopyNumberPerGram))
hist(sqrt(FunUn$CopyNumberPerGram))


#----
#----
#*****************************************************************************************************************************----
#DESRIPTIVE STATISTICS----------------------------------------------------------------------------------------------------------
#*****************************************************************************************************************************----
#set as factor first to set base level............................................................-----.
Metadata$Treatment<-as.factor(Metadata$Treatment)
Metadata$Treatment <- try(relevel(Metadata$Treatment , "Unburned"))

MetaFun$Treatment<-as.factor(MetaFun$Treatment)
MetaFun$Treatment <- try(relevel(MetaFun$Treatment , "Unburned"))

MetaBac$Treatment<-as.factor(MetaBac$Treatment)
MetaBac$Treatment <- try(relevel(MetaBac$Treatment , "Unburned"))

Metadata$AshSeverity<-ordered(Metadata$AshSeverity, c("Unburned","Low","Moderate", "High"))
MetaFun$AshSeverity<-ordered(MetaFun$AshSeverity, c("Unburned","Low","Moderate", "High"))
MetaBac$AshSeverity<-ordered(MetaBac$AshSeverity, c("Unburned","Low","Moderate", "High"))

#----descriptive stats...........................................................................-----
TrtStats<-Metadata %>%
  filter(!is.na(CopyNumberPerGram)) %>%
  group_by(Treatment, Kingdom) %>%
  summarise(n_obs = n(),
            TotalGenes =sum(CopyNumberPerGram, na.rm = TRUE),
            Avg_Genes = mean(CopyNumberPerGram, na.rm = TRUE), 
            Variance = var(CopyNumberPerGram, na.rm = TRUE), 
            sd = sd(CopyNumberPerGram), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 0) %>% as.data.frame();TrtStats 

TSFstats<-Metadata %>%
  filter(!is.na(CopyNumberPerGram)) %>%
  group_by(Treatment, Kingdom, TSFdays) %>%
  summarise(n_obs = n(),
            TotalGenes =sum(CopyNumberPerGram, na.rm = TRUE),
            Avg_ASV = mean(CopyNumberPerGram, na.rm = TRUE), 
            Variance = var(CopyNumberPerGram, na.rm = TRUE), 
            sd = sd(CopyNumberPerGram), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 0) %>% as.data.frame();TSFstats

Sevstats<-Metadata %>%
  filter(!is.na(CopyNumberPerGram)) %>%
  group_by(Severity, Kingdom) %>%
  summarise(n_obs = n(),
            TotalGenes =sum(CopyNumberPerGram, na.rm = TRUE),
            Avg_ASV = mean(CopyNumberPerGram, na.rm = TRUE), 
            Variance = var(CopyNumberPerGram, na.rm = TRUE), 
            sd = sd(CopyNumberPerGram), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 0) %>% as.data.frame();Sevstats 

#FUNGAL ................................................................
TrtStatsFun<-MetaFun%>%
  filter(!is.na(CopyNumberPerGram)) %>%
  group_by(Treatment) %>%
  summarise(n_obs = n(),
            TotalGenes =sum(CopyNumberPerGram, na.rm = TRUE),
            Avg_Genes = mean(CopyNumberPerGram, na.rm = TRUE), 
            Variance = var(CopyNumberPerGram, na.rm = TRUE), 
            sd = sd(CopyNumberPerGram), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 0) %>% as.data.frame();TrtStatsFun

TSFstatsFun<-MetaFun %>%
  filter(!is.na(CopyNumberPerGram)) %>%
  group_by(TSFdays,Treatment) %>%
  summarise(n_obs = n(),
            TotalGenes =sum(CopyNumberPerGram, na.rm = TRUE),
            Avg_ASV = mean(CopyNumberPerGram, na.rm = TRUE), 
            Variance = var(CopyNumberPerGram, na.rm = TRUE), 
            sd = sd(CopyNumberPerGram), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 0) %>% as.data.frame();TSFstatsFun

TSFstatsFunB<-FunBurn %>%
  filter(!is.na(CopyNumberPerGram)) %>%
  group_by(TSFdays) %>%
  summarise(n_obs = n(),
            TotalGenes =sum(CopyNumberPerGram, na.rm = TRUE),
            Avg_ASV = mean(CopyNumberPerGram, na.rm = TRUE), 
            Variance = var(CopyNumberPerGram, na.rm = TRUE), 
            sd = sd(CopyNumberPerGram), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 0) %>% as.data.frame();TSFstatsFunB 

SevstatsFunB<-FunBurn%>%
  filter(!is.na(CopyNumberPerGram)) %>%
  group_by(Severity) %>%
  summarise(n_obs = n(),
            TotalGenes =sum(CopyNumberPerGram, na.rm = TRUE),
            Avg_ASV = mean(CopyNumberPerGram, na.rm = TRUE), 
            Variance = var(CopyNumberPerGram, na.rm = TRUE), 
            sd = sd(CopyNumberPerGram), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 0) %>% as.data.frame();SevstatsFunB 


#Bacteria .................................................................
TrtStatsBac<-MetaBac%>%
  filter(!is.na(CopyNumberPerGram)) %>%
  group_by(Treatment) %>%
  summarise(n_obs = n(),
            TotalGenes =sum(CopyNumberPerGram, na.rm = TRUE),
            Avg_Genes = mean(CopyNumberPerGram, na.rm = TRUE), 
            Variance = var(CopyNumberPerGram, na.rm = TRUE), 
            sd = sd(CopyNumberPerGram), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 0) %>% as.data.frame();TrtStatsBac

TSFstatsBac<-MetaBac %>%
  filter(!is.na(CopyNumberPerGram)) %>%
  group_by(TSFdays,Treatment) %>%
  summarise(n_obs = n(),
            TotalGenes =sum(CopyNumberPerGram, na.rm = TRUE),
            Avg_ASV = mean(CopyNumberPerGram, na.rm = TRUE), 
            Variance = var(CopyNumberPerGram, na.rm = TRUE), 
            sd = sd(CopyNumberPerGram), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 0) %>% as.data.frame();TSFstatsBac

TSFstatsBacB<-BacBurn %>%
  filter(!is.na(CopyNumberPerGram)) %>%
  group_by(TSFdays) %>%
  summarise(n_obs = n(),
            TotalGenes =sum(CopyNumberPerGram, na.rm = TRUE),
            Avg_ASV = mean(CopyNumberPerGram, na.rm = TRUE), 
            Variance = var(CopyNumberPerGram, na.rm = TRUE), 
            sd = sd(CopyNumberPerGram), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 0) %>% as.data.frame();TSFstatsBacB

SevstatsBacB<-BacBurn%>%
  filter(!is.na(CopyNumberPerGram)) %>%
  group_by(Severity) %>%
  summarise(n_obs = n(),
            TotalGenes =sum(CopyNumberPerGram, na.rm = TRUE),
            Avg_ASV = mean(CopyNumberPerGram, na.rm = TRUE), 
            Variance = var(CopyNumberPerGram, na.rm = TRUE), 
            sd = sd(CopyNumberPerGram), se =sd / sqrt(n_obs)) %>%
  filter(n_obs > 0) %>% as.data.frame();SevstatsBacB


# EXPORT FILES ......................................................................................
dir.create(file.path("Analysis/Tables/DescriptiveStats/Bac"), recursive = TRUE)
dir.create(file.path("Analysis/Tables/DescriptiveStats/Fun"), recursive = TRUE)
dir.create(file.path("Analysis/Tables/DescriptiveStats/Trt"), recursive = TRUE)
dir.create(file.path("Analysis/Graphs/Fun"), recursive = TRUE)
dir.create(file.path("Analysis/Graphs/Bac"), recursive = TRUE)

#DESCRIPTIVE STATS.......................................................................................................
#Full Metadata
write.csv(TrtStats, file="Analysis/Tables/DescriptiveStats/Trt/Results-Biomass-DesStats-Trt-Norm.csv") 
write.csv(TSFstats, file="Analysis/Tables/DescriptiveStats/Trt/Results-Biomass-DesStats-TSF-Norm.csv") 
write.csv(Sevstats, file="Analysis/Tables/DescriptiveStats/Trt/Results-Biomass-DesStats-SEV-Norm.csv") 

#Fungal only
write.csv(TrtStatsFun, file="Analysis/Tables/DescriptiveStats/Fun/Results-Biomass-DesStats-Trt-Norm.csv") 
write.csv(TSFstatsFun, file="Analysis/Tables/DescriptiveStats/Fun/Results-Biomass-DesStats-TSF-Trt-Norm.csv") 
write.csv(TSFstatsFunB, file="Analysis/Tables/DescriptiveStats/Fun/Results-Biomass-DesStats-TSF-Norm.csv") 
write.csv(SevstatsFunB, file="Analysis/Tables/DescriptiveStats/Fun/Results-Biomass-DesStats-Sev-Norm.csv") 

#Bacteria only
write.csv(TrtStatsBac, file="Analysis/Tables/DescriptiveStats/Bac/Results-Biomass-DesStats-Trt-Norm.csv") 
write.csv(TSFstatsBac, file="Analysis/Tables/DescriptiveStats/Bac/Results-Biomass-DesStats-TSF-Trt-Norm.csv") 
write.csv(TSFstatsBacB, file="Analysis/Tables/DescriptiveStats/Bac/Results-Biomass-DesStats-TSF-Norm.csv") 
write.csv(SevstatsBacB, file="Analysis/Tables/DescriptiveStats/Bac/Results-Biomass-DesStats-Sev-Norm.csv")


#----
#----
#----
#******************************************************************************************************************************----
#CALCULATE PERCENT CHANGE BETWEEN TREATMENTS AND TIMEPOINTS------------------------------------------------------------------------
#******************************************************************************************************************************----
Trt <-  Metadata%>% 
  group_by( Kingdom,Treatment) %>% 
  summarise(mean = mean(CopyNumberPerGram,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();Trt 

Sev<-Metadata%>% 
  group_by(Kingdom,AshSeverity) %>% 
  summarise(mean = mean(CopyNumberPerGram,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();Sev

TrtTSF <-  Metadata%>% 
  group_by(Kingdom,TSFdays,Treatment) %>% 
  summarise(mean = mean(CopyNumberPerGram,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();TrtTSF

#Fungi.....................................
TrtTSFfun<-  MetaFun%>% 
  group_by(TSFdays,Treatment) %>% 
  summarise(mean = mean(CopyNumberPerGram,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();TrtTSFfun

TsfBFun <- FunBurn%>% 
  group_by(TSFdays) %>% 
  summarise(mean = mean(CopyNumberPerGram,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();TsfBFun

TsfUnFun <- FunUn%>% 
  group_by(TSFdays) %>% 
  summarise(mean = mean(CopyNumberPerGram,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();TsfUnFun


#Bacteria.....................................
TrtTSFbac <-  MetaBac%>% 
  group_by(TSFdays,Treatment) %>% 
  summarise(mean = mean(CopyNumberPerGram,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();TrtTSFbac

TsfBbac <- BacBurn%>% 
  group_by(TSFdays) %>% 
  summarise(mean = mean(CopyNumberPerGram,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();TsfBbac

TsfUnBac <- BacUn%>% 
  group_by(TSFdays) %>% 
  summarise(mean = mean(CopyNumberPerGram,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();TsfUnBac



#Create directory........................................................................
dir.create(file.path("Analysis/Tables/PercentChange/Bac"), recursive = TRUE)
dir.create(file.path("Analysis/Tables/PercentChange/Fun"), recursive = TRUE)
dir.create(file.path("Analysis/Tables/PercentChange/Trt"), recursive = TRUE)


#Percent change All (bac and fun together).....................................................................................
write.csv(Trt,"Analysis/Tables/PercentChange/Trt/PerCh-Trt-All-Norm.csv")
write.csv(Sev, "Analysis/Tables/PercentChange/Trt/PerCh-TSF-All-Norm.csv")
write.csv(TrtTSF,"Analysis/Tables/PercentChange/Trt/PerCh-TRT-TST-All-Norm.csv")

#Percent Bacteria ....................................................
write.csv(TrtTSFbac,"Analysis/Tables/PercentChange/Bac/PerCh-Trt-Tsf-Bac-Norm.csv")
write.csv(TsfBbac,"Analysis/Tables/PercentChange/Bac/PerCh-Tsf-Burned-Bac-Norm.csv")
write.csv(TsfUnBac,"Analysis/Tables/PercentChange/Bac/PerCh-Tsf-Unburned-Bac-Norm.csv")

#Percent Fungi....................................................
write.csv(TrtTSFfun,"Analysis/Tables/PercentChange/Fun/PerCh-Trt-Tsf-Fun-Norm.csv")
write.csv(TsfBFun,"Analysis/Tables/PercentChange/Fun/PerCh-Tsf-Burned-Fun-Norm.csv")
write.csv(TsfUnFun,"Analysis/Tables/PercentChange/Fun/PerCh-Tsf-Unburned-Fun-Norm.csv")






#----
#----
#----
