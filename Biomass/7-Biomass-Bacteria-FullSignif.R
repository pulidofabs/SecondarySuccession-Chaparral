#8-12-2022
rm(list=ls())#reset working directory

#Set working directory.................................................................
setwd("~/Dropbox/1-Dissertation/1-Research/HolyFire/1-HolyFire-Pub1/1c-Qpcr/Year1")

#Load packages ........................................................................
library(MuMIn)#for AIC selection
library(lme4)#for neg binom

# Load data ............................................................................................
MetaBioBac<-read.csv("Analysis/Metadata/MetaBac-Biomass.csv", header = TRUE)



#----
#----
#******************************************************************************************************************----
#  ............ Quality control  ................................................................................
#******************************************************************************************************************----
MetaBioBac$Plot<- as.factor(MetaBioBac$Plot)
MetaBioBac$Treatment<-as.factor(MetaBioBac$Treatment)
MetaBioBac$Treatment <- try(relevel(MetaBioBac$Treatment , "Unburned"))#relevel so that unburned is the base level for comparison



#----
#----
#******************************************************************************************************************----
#  ............ Significance testing--Backward model selection  ................................................................................
#******************************************************************************************************************----
#Per model selection script, we will use the folling model for testing significance
#Mod<- glmer.nb(CopyNumberPerGramPerGram ~ (1|Plot)+(1|TSFdays), nAGQ=0,data = MetaBioBac);summary(Mod)#

#Rescale variables..................................................................................----
#Can use the center= True and scale=True but that is the normal so do not need to specify
# .....Example scale(MetaBioBac$TotPrecip, center = TRUE, scale = TRUE)
#......Scale, center = gives you a central point along which data is centered around,
#......Scale, scale=True means tart it will subtract the mean and divide by the std deviation
#------

summary(MetaBioBac)#Check to see what to rescale

#Rescale and center the numerical data
MetaBioBac$TSF<-scale(MetaBioBac$TSFdays);summary(MetaBioBac$TSF)
MetaBioBac$Precip<-scale(MetaBioBac$TotPrecip);summary(MetaBioBac$TSF)

#nAGQ=0 #only use when model will not converge

attach(MetaBioBac)

#Run backward model selection.............................................................................
#Did not use trt and initial ash depht, since the interaction does not make sense, unburn ash depth = 0
Full<-glmer.nb(CopyNumberPerGram ~ Treatment*TSF + Treatment*Precip + TSF*Precip + 
                 TSF*InitialAshDepth + Treatment + TSF + Precip + InitialAshDepth +
                 (1|Plot) + (1|TSFdays), nAGQ=0, data = MetaBioBac);summary(Full)

Null<-glmer.nb(CopyNumberPerGram ~ (1|Plot)+ (1|TSFdays), data=MetaBioBac);summary(Null)
anova(Full, Null)#Full model better

#backward selection...start by removing interaction w the lowest p-value.....................


#Remove TSF, highest P
Mod1<-glmer.nb(CopyNumberPerGram ~ Treatment*TSF + Treatment*Precip + TSF*Precip + 
                 TSF*InitialAshDepth + Treatment +  Precip + InitialAshDepth +
                 (1|Plot) + (1|TSFdays), nAGQ=0, data = MetaBioBac);summary(Mod1)
anova(Mod1,Full)#Same AIC, mod 1, not signif remove


#Remove Treatment*TSF +  (highest P)
Mod2<-glmer.nb(CopyNumberPerGram ~  Treatment*Precip + TSF*Precip + 
                 TSF*InitialAshDepth + Treatment +  Precip + InitialAshDepth +
                 (1|Plot) + (1|TSFdays), nAGQ=0, data = MetaBioBac);summary(Mod2)
anova(Mod1,Mod2)#Mod 1 better, not signif, remove


#Remove TSF*Precip +    (highest P)
Mod3<-glmer.nb(CopyNumberPerGram ~  Treatment*Precip + 
            TSF*InitialAshDepth + Treatment +  Precip + InitialAshDepth +
            (1|Plot) + (1|TSFdays), nAGQ=0, data = MetaBioBac);summary(Mod3)
anova(Mod2,Mod3)#Mod 3 better, not signif, remove


#Remove Precip + (highest P)
Mod4<-glmer.nb(CopyNumberPerGram ~  Treatment*Precip + 
                 TSF*InitialAshDepth + Treatment +   InitialAshDepth +
                 (1|Plot) + (1|TSFdays), nAGQ=0, data = MetaBioBac);summary(Mod4)
anova(Mod3,Mod4)#Same AIC, not signif, remove



#Remove Treatment*Precip + (highest P)
Mod5<-glmer.nb(CopyNumberPerGram ~   
                 TSF*InitialAshDepth + Treatment +   InitialAshDepth +
                 (1|Plot) + (1|TSFdays), nAGQ=0, data = MetaBioBac);summary(Mod5)
anova(Mod4,Mod5)#Mod 4 better, significant, do not remove.................................


#Remove  InitialAshDepth + (highest P)
Mod6<-glmer.nb(CopyNumberPerGram ~  TSF*InitialAshDepth + Treatment +  Treatment*Precip +
          (1|Plot) + (1|TSFdays), nAGQ=0, data = MetaBioBac);summary(Mod6)
anova(Mod4,Mod6)#Same AIC, not significant, remove


#Remove  Treatment +  (highest P)
Mod7<-glmer.nb(CopyNumberPerGram ~  TSF*InitialAshDepth +  Treatment*Precip +
                 (1|Plot) + (1|TSFdays), nAGQ=0, data = MetaBioBac);summary(Mod7)
anova(Mod6,Mod7)#Same AIC, not significant, remove


#Remove   TSF*InitialAshDepth +  (highest P)
Mod8<-glmer.nb(CopyNumberPerGram ~  Treatment*Precip +
                 (1|Plot) + (1|TSFdays), nAGQ=0, data = MetaBioBac);summary(Mod8)
anova(Mod7,Mod8)#Significant, do not remove............................................


#Compare all models
anova(Mod1, Mod2, Mod3, Mod4, Mod5, Mod6, Mod7, Mod8) 
AICc(Mod1, Mod2, Mod3, Mod4, Mod5, Mod6, Mod7, Mod8)
#3 models same AIC, but model 7 is more parsimpnious



#----
#----
#******************************************************************************************************************----
#  ............ Final Model -Results  ................................................................................
#******************************************************************************************************************----
FunMod<-glmer.nb(CopyNumberPerGram ~  TSF*InitialAshDepth +  Treatment*Precip +
               (1|Plot) + (1|TSFdays), nAGQ=0, data = MetaBioBac);summary(FunMod)
R2Full<-r.squaredGLMM(FunMod);R2Full#Delta; R2m = 0.3026104; 0.4639148

#Export model results.............................................................................................
capture.output(summary(FunMod),file="Analysis/Tables/Significance/Bac-NegBinom-FullModel.csv")
write.csv(R2Full,"Analysis/Tables/Significance/Bac-NegBinom-R2FullModel.csv")


#Model Diagnostic & Validation plots ...........................................
par(mfrow=c(1,1));plot(FunMod)
Fitted <- fitted(FunMod) #command fitted gives already e^model
Res <- resid(FunMod, type = "pearson")
plot(x = Fitted,  y = Res,xlab = "Fitted values",ylab = "Pearson residuals")
abline(h = 0, lty = 2, col="red")


