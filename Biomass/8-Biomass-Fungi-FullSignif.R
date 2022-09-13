#8-12-2022
rm(list=ls())#reset working directory

#Set working directory.................................................................
setwd("~/Dropbox/1-Dissertation/1-Research/HolyFire/1-HolyFire-Pub1/1c-Qpcr/Year1")

#Load packages ........................................................................
library(MuMIn)#for AIC selection
library(lme4)#for neg binom

# Load data ............................................................................................
MetaBioFun<-read.csv("Analysis/Metadata/MetaFun-Biomass.csv", header = TRUE)



#----
#----
#******************************************************************************************************************----
#  ............ Quality control  ................................................................................
#******************************************************************************************************************----
MetaBioFun$Plot<- as.factor(MetaBioFun$Plot)
MetaBioFun$Treatment<-as.factor(MetaBioFun$Treatment)
MetaBioFun$Treatment <- try(relevel(MetaBioFun$Treatment , "Unburned"))#relevel so that unburned is the base level for comparison



#----
#----
#******************************************************************************************************************----
#  ............ Significance testing--Funkward model selection  ................................................................................
#******************************************************************************************************************----
#Per model selection script, we will use the folling model for testing significance
#Mod<- glmer.nb(CopyNumberPerGramPerGram ~ (1|Plot)+(1|TSFdays), nAGQ=0,data = MetaBioFun);summary(Mod)#

#Rescale variables..................................................................................----
#Can use the center= True and scale=True but that is the normal so do not need to specify
# .....Example scale(MetaBioFun$TotPrecip, center = TRUE, scale = TRUE)
#......Scale, center = gives you a central point along which data is centered around,
#......Scale, scale=True means tart it will subtract the mean and divide by the std deviation
#------

summary(MetaBioFun)#Check to see what to rescale

#Rescale and center the numerical data
MetaBioFun$TSF<-scale(MetaBioFun$TSFdays);summary(MetaBioFun$TSF)
MetaBioFun$Precip<-scale(MetaBioFun$TotPrecip);summary(MetaBioFun$TSF)

#nAGQ=0 #only use when model will not converge

attach(MetaBioFun)

#Run backward model selection.............................................................................
#Trt*ash not run since it makes no sense, initial ash depth in unburned = 0

Full<-glmer.nb(CopyNumberPerGram ~ Treatment*TSF + Treatment*Precip + TSF*Precip + 
                 TSF*InitialAshDepth + Treatment + TSF + Precip + InitialAshDepth +
                 (1|Plot) + (1|TSFdays), nAGQ=0,  data = MetaBioFun);summary(Full)

Null<-glmer.nb(CopyNumberPerGram ~ (1|Plot)+ (1|TSFdays), nAGQ=0, data=MetaBioFun);summary(Null)
anova(Full, Null)#Full model better



#backward selection...start by removing interaction w the lowest p-value................................

#Remove TSF*InitialAshDepth +  (highest P)
Mod1<-glmer.nb(CopyNumberPerGram ~ Treatment*TSF + Treatment*Precip + TSF*Precip + 
                  Treatment + TSF + Precip + InitialAshDepth +
                 (1|Plot) + (1|TSFdays), nAGQ=0,  data = MetaBioFun);summary(Mod1)
anova(Mod1,Full)#Mod 1 better, not signif, remove


#Remove  InitialAshDepth +  (highest P)
Mod2<-glmer.nb(CopyNumberPerGram ~ Treatment*TSF + Treatment*Precip + TSF*Precip + 
                 Treatment + TSF + Precip + 
                 (1|Plot) + (1|TSFdays), nAGQ=0,  data = MetaBioFun);summary(Mod2)
anova(Mod1,Mod2)#Mod 1, Significant, do not remove...............................................


#Remove Treatment*Precip +  (highest P)
Mod3<-glmer.nb(CopyNumberPerGram ~ Treatment*TSF +  TSF*Precip + 
                 Treatment + TSF + Precip + InitialAshDepth +
                 (1|Plot) + (1|TSFdays), nAGQ=0,  data = MetaBioFun);summary(Mod3)
anova(Mod1,Mod3)#Mod 1 better, signif, do not remove...........................................

#Remove Precip +  (highest P)
Mod4<-glmer.nb(CopyNumberPerGram ~ Treatment*TSF +  TSF*Precip + 
                 Treatment + TSF +  InitialAshDepth + Treatment*Precip +
                 (1|Plot) + (1|TSFdays), nAGQ=0,  data = MetaBioFun);summary(Mod4)
anova(Mod1,Mod4)#Same AIC, not signif, remove


#Remove  TSF +  (highest P)
Mod5<-glmer.nb(CopyNumberPerGram ~ Treatment*TSF +  TSF*Precip + 
                 Treatment +  InitialAshDepth + Treatment*Precip +
                 (1|Plot) + (1|TSFdays), nAGQ=0,  data = MetaBioFun);summary(Mod5)
anova(Mod4,Mod5)#Same AIC, Not significant, remove


#Remove  Treatment + (highest P)
Mod6<-glmer.nb(CopyNumberPerGram ~ Treatment*TSF +  TSF*Precip + 
                  InitialAshDepth + Treatment*Precip +
                 (1|Plot) + (1|TSFdays), nAGQ=0,  data = MetaBioFun);summary(Mod6)
anova(Mod5,Mod6)#Same AIC, not signif, remove


#Remove TSF*Precip +  (highest P)
Mod7<-glmer.nb(CopyNumberPerGram ~ Treatment*TSF + InitialAshDepth +  TSF*InitialAshDepth +  
          Treatment*Precip + (1|Plot) + (1|TSFdays), nAGQ=0,  data = MetaBioFun);summary(Mod7)
anova(Mod6,Mod7)#Mod 7 better, not signif, remove

#Remove Treatment*TSF +  (highest P)
Mod8<-glmer.nb(CopyNumberPerGram ~ InitialAshDepth + Treatment*Precip + (1|Plot) + (1|TSFdays), 
               nAGQ=0,  data = MetaBioFun);summary(Mod8)
anova(Mod7,Mod8)#Significant, do not remove............................................................



#Compare all models
anova(Mod1, Mod2, Mod3, Mod4, Mod5, Mod6, Mod7, Mod8) 
AICc(Mod1, Mod2, Mod3, Mod4, Mod5, Mod6, Mod7, Mod8)
#Mod 7 butter,




#----
#----
#******************************************************************************************************************----
#  ............ Final Model -Results  ................................................................................
#******************************************************************************************************************----
FunMod<-glmer.nb(CopyNumberPerGram ~ Treatment*TSF + InitialAshDepth + 
                   Treatment*Precip + (1|Plot) + (1|TSFdays), nAGQ=0,  data = MetaBioFun);summary(FunMod)
R2Full<-r.squaredGLMM(FunMod);R2Full#Delta; R2m = 0.4728572 0.5570665

#Export model results.............................................................................................
capture.output(summary(FunMod),file="Analysis/Tables/Significance/Fun-NegBinom-FullModel.csv")
write.csv(R2Full,"Analysis/Tables/Significance/Fun-NegBinom-R2FullModel.csv")


#Model Diagnostic & Validation plots ...........................................
par(mfrow=c(1,1));plot(FunMod)
Fitted <- fitted(FunMod) #command fitted gives already e^model
Res <- resid(FunMod, type = "pearson")
plot(x = Fitted,  y = Res,xlab = "Fitted values",ylab = "Pearson residuals")
abline(h = 0, lty = 2, col="red")


