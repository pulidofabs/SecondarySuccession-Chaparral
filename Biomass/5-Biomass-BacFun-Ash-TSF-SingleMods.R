#8-12-2022
rm(list=ls())#reset working directory

#Notes
#Significance testing for single TSF and Trt interactions 
#Significance testing for single TSF*Trt*InitialAshDept interactions
#Model selection has been donw in a different script


#Set working directory................................................................................
setwd("~/Dropbox/1-Dissertation/1-Research/HolyFire/1-HolyFire-Pub1/1c-Qpcr/Year1")

#Load packages 
library(phia)

# Load data ............................................................................................
MetaBioBac<-read.csv("Analysis/Metadata/MetaBac-Biomass.csv", header = TRUE)
MetaBioFun<-read.csv("Analysis/Metadata/MetaFun-Biomass.csv", header = TRUE)

#Quality control...............................................................................................................................
MetaBioBac$Treatment<-as.factor(MetaBioBac$Treatment)
MetaBioBac$TSF<-as.factor(MetaBioBac$TSFdays)
MetaBioBac$Treatment <- try(relevel(MetaBioBac$Treatment , "Unburned"))#relevel so that unburned is the base level for comparison

MetaBioFun$Treatment<-as.factor(MetaBioFun$Treatment)
MetaBioFun$TSF<-as.factor(MetaBioFun$TSFdays)
MetaBioFun$Treatment <- try(relevel(MetaBioFun$Treatment , "Unburned"))#relevel so that unburned is the base level for comparison



#.----
#.----
#******************************************************************************************************************----
# .................. SIGNIFICANCE TREATMENT * TSF only ...................................................
#******************************************************************************************************************----
#Wil use package phia since it can use glmer models
#adjust.p = “holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”

#Bacteria.............................................................................................................
BacTSF<-glmer.nb(CopyNumberPerGram ~ Treatment*TSF* (1|Plot) + (1|TSF), MetaBioBac, nAGQ=0);summary(BacTSF)
BTR2<-r.squaredGLMM(BacTSF);BTR2 #0.4332728 0.5173094
BacHoc<-testInteractions(BacTSF, pairwise ="Treatment",fixed = "TSF", adjustment = "BH");BacHoc

capture.output(summary(BacTSF), file="Analysis/Tables/Significance/Bac-NegBinom-TrtTSF.csv")
capture.output(BTR2, file= "Analysis/Tables/Significance/Bac-NegBinom-TrtTSF.csv")
capture.output(BacHoc, file= "Analysis/Tables/Significance/Bac-NegBinom-PostHoc-TrtTSF.csv")

#Fungi.............................................................................................................
FunTSF<-glmer.nb(CopyNumberPerGram ~ Treatment*TSF*TotPrecip + (1|Plot) + (1|TSF), MetaBioFun, nAGQ=0);summary(FunTSF)
BTR2<-r.squaredGLMM(FunTSF);BTR2 #0.5627252 0.6051276
FunHoc<-testInteractions(FunTSF, pairwise ="Treatment",fixed = "TSF", adjustment = "BH");FunHoc

capture.output(summary(FunTSF), file="Analysis/Tables/Significance/Fun-NegBinom-TrtTSF.csv")
capture.output(BTR2, file= "Analysis/Tables/Significance/Fun-NegBinom-TrtTSF.csv")
capture.output(FunHoc, file= "Analysis/Tables/Significance/Fun-NegBinom-PostHoc-TrtTSF.csv")




#.----
#.----
#***********************************************************************************************************************----
# .................. SIGNIFICANCE  TSF * ASH  ...................................................
#***********************************************************************************************************************----
#Wil use package phia since it can use glmer models
#adjust.p = “holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”

MetaBioBac$TSF<-scale(MetaBioBac$TSFdays)
MetaBioFun$TSF<-scale(MetaBioFun$TSFdays)

#Bacteria.............................................................................................................
#Time since fire-* ash depth- biomass
BioBAsh<-glmer.nb(CopyNumberPerGram ~ InitialAshDepth*TSF + (1|Plot)+ (1|TSF), nAGQ=0, MetaBioBac);summary(BioBAsh)
R2AshBB<-r.squaredGLMM(BioBAsh);R2AshBB#delta 0.09726821 0.4304912

#Fungi...............................................................................................................
BioFAsh<-glmer.nb(CopyNumberPerGram ~ InitialAshDepth*TSF + (1|Plot)+ (1|TSF), nAGQ=0, MetaBioFun);summary(BioFAsh)
R2AshBF<-r.squaredGLMM(BioFAsh);R2AshBF#delta 0.04900840 0.5460041

#Capture ash signif values......................................................................
capture.output(summary(BioBAsh), file="Analysis/Tables/Significance/Bac-NegBinom-AshTSF.csv")
capture.output(R2AshBB, file= "Analysis/Tables/Significance/Bac-NegBinom-AshTSF.csv")
capture.output(summary(BioFAsh), file="Analysis/Tables/Significance/Fun-NegBinom-AshTSF.csv")
capture.output(R2AshBF, file= "Analysis/Tables/Significance/Fun-NegBinom-AshTSF.csv")



