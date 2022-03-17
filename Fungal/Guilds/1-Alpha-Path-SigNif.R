#Last edits: 10-6-2021

#reset working directory.....
rm(list=ls())


#Set working directory....................................................................
setwd("C:/Users/fabipc/Dropbox/3-HolyFire/1-HolyFire-Pub1/2-Fungi/")


#Load additional packages.................................................................
library(dplyr)
library(lme4)#glmer.nb
library(lmerTest)#to get P value
library(emmeans)#tukey comparisons
library(rcompanion)#pseudo R2
library(MuMIn)#
library(MASS)


#Load packages.............................................................................
MetaRarePath<- read.csv("Metadata/Path/MetaRareSppPath-Soil.csv", na.strings = "NA", header = TRUE)



#.----
#.----
#**************************************************************************************************----
# .................QUALITY CONTROL.....................................................................
#**************************************************************************************************----

# * Look at data structure ................................................................
MetaRarePath$Plot<- as.factor(MetaRarePath$Plot)
MetaRarePath$TSFdays<-as.numeric(MetaRarePath$TSFdays)
MetaRarePath$Treatment<-as.factor(MetaRarePath$Treatment)
MetaRarePath$Treatment <- try(relevel(MetaRarePath$Treatment,"Unburned")) #change base level
MetaRarePath$AshSeverity<-ordered(MetaRarePath$AshSeverity,c("Unburned","Low","Moderate","High"))
MetaRarePath$Date<-ordered(MetaRarePath$Date, c("9/30/18","10/8/18","10/17/18","11/19/18",
                     "12/17/18","1/22/19","3/19/19","6/26/19","9/24/19"))


# *  * * TEST FOR NORMALITY ............................................................
par(mfrow=c(1,3))
hist(MetaRarePath$S.obs)
shapiro.test(MetaRarePath$S.obs)#< 2.2e-16------Non-normal



#.----
#.----
#**************************************************************************************************----
# .................MODEL SELECTION ....................................................................
#***************************************************************************************************----
attach(MetaRarePath)

#* * * TREATMENT EFFECTS .......................................................
#Checking mean and variance of the data to see if its NegBinom..................
df <- MetaRarePath %>%
  group_by(Treatment) %>%
  summarise(means = mean(S.obs),
            variance = var(S.obs));df #suggest negative binomial


#Start with a poisson model ....................................................
Pois1 <- glm(S.obs ~ Treatment, data = MetaRarePath, family = poisson(link = log))
ods <- Pois1$deviance / Pois1$df.residual;ods#0verdispersion 0.59
par(mfrow=c(2,2))
plot(Pois1)#look at plots, cooks, levenes etc

#Zero inflaited model ..........................................................
sum(S.obs == 0) #0


#Test different identity link to see if disperssion decreases ..................
Pois2 <- glm(S.obs ~ Treatment,family = poisson(link = sqrt))
ods2<- Pois2$deviance / Pois2$df.residual;ods2 #0.59

Pois3 <- glm(S.obs ~ Treatment,family = poisson(link = identity))
ods3<- Pois3$deviance / Pois3$df.residual;ods3##0.59

Gam <- glm(S.obs ~ Treatment,family = "Gamma")
odsGam<-Gam$deviance/Gam$df.residual;odsGam#Does not work

#Negative binomial... ..........................................................
#* mean and variance suggest this is the model to use
nb1 <- glm.nb(S.obs ~ Treatment,data = MetaRarePath)
ods_nb <- (nb1$deviance / nb1$df.residual);ods_nb#0.77 mild= acceptable
plot(nb1)

#Compare models.................................................................
AIC(Pois1, Pois2, Pois3, nb1)#neg binomial slightly better 

#Gam and negative binomial were same, but for consistency
#will use Negative Binomial


#----
#----
#*****************************************************************************************************----
#........... TEST THE LEVEL OF NESTEDNESS .............................................................
#*****************************************************************************************************----
#Load metadata to make it easier to call out variables..............
attach(MetaRarePath)
summary(MetaRarePath)

#Rescale variables..................................................
#Rescale and center the numerical data
MetaRarePath$TSF<-scale(MetaRarePath$TSFdays);summary(MetaRarePath$TSF)
MetaRarePath$Precip<-scale(MetaRarePath$TotPrecip)

#Explore Temporal correlations to ensure that we are capturing as much of the variance 
ASVneg1 <- glmer.nb(S.obs ~ (1|Plot) + (1|Subplot)+(1|TSFdays),data = MetaRarePath);summary(ASVneg1)
ASVneg2 <- glmer.nb(S.obs ~ (1|Plot)+(1|TSFdays),data = MetaRarePath);summary(ASVneg2)
ASVneg3 <- glmer.nb(S.obs ~ (1|Subplot)+(1|TSFdays),data = MetaRarePath);summary(ASVneg3)


#MuMIN for model selection................................
AICc(ASVneg1,ASVneg2,ASVneg3)#ASVneg1 better 292



#.----
#***************************************************************************************************----
# .................SIGNIFICANCE  .....................................................................
#***************************************************************************************************----

# * * FULL MODEL ..............................................................................
# * * * * Scale variables to fit the model............
Full<-glmer.nb(S.obs~Treatment*TSF + Treatment*Precip + Treatment*InitialAshDepth + TSF*Precip + 
               TSF*InitialAshDepth + Treatment + TSF + Precip + InitialAshDepth + (1|Plot) +
              (1|Subplot) + (1|TSFdays),data = MetaRarePath);summary(Full)

Null<-glmer.nb(S.obs ~ (1|Plot) + (1|Subplot) + (1|TSFdays), data=MetaRarePath);summary(Null)
anova(Full,Null)#Full model better



#BACKWARDS MODEL SELECTION........................................................................
Mod1<-glmer.nb(S.obs~Treatment*TSF + Treatment*Precip + Treatment*InitialAshDepth + 
         TSF*Precip + TSF*InitialAshDepth + Treatment + TSF + Precip + InitialAshDepth +
         (1|Plot) + (1|Subplot) + (1|TSFdays),data = MetaRarePath);summary(Mod1)


#Remove Treatment*TSF + 
Mod2<-glmer.nb(S.obs~Treatment*Precip + Treatment*InitialAshDepth + TSF*Precip + 
          TSF*InitialAshDepth + Treatment + TSF + Precip + InitialAshDepth +
          (1|Plot) + (1|Subplot) + (1|TSFdays),data = MetaRarePath);summary(Mod2)
anova(Mod1,Mod2)#Same AIC, not signif, remove


#Remove Treatment*Precip + 
Mod3<-glmer.nb(S.obs~Treatment*InitialAshDepth + TSF*Precip + TSF*InitialAshDepth +
            Treatment + TSF + Precip + InitialAshDepth + (1|Plot) + (1|Subplot) +
            (1|TSFdays),data = MetaRarePath);summary(Mod3)
anova(Mod2,Mod3)#Mod 2 better,  signif, do not remove................................


#Remove Treatment*InitialAshDepth +
Mod4<-glmer.nb(S.obs~ TSF*Precip + TSF*InitialAshDepth + Treatment + TSF + 
            Precip + InitialAshDepth + Treatment*Precip + (1|Plot) +
            (1|Subplot) + (1|TSFdays),data = MetaRarePath);summary(Mod4)
anova(Mod2,Mod4)#Same AIC, not signif, do not remove


#Remove TSF*Precip + 
Mod5<-glmer.nb(S.obs ~ TSF*InitialAshDepth + Treatment + TSF + Precip + 
            InitialAshDepth + Treatment*Precip + (1|Plot) + (1|Subplot) +
            (1|TSFdays),data = MetaRarePath);summary(Mod5)
anova(Mod4,Mod5)#Mod 5 better, not signif, do not remove


#Remove TSF*InitialAshDepth + 
Mod6<-glmer.nb(S.obs ~ Treatment + TSF + Precip +  InitialAshDepth + 
          Treatment*Precip + (1|Plot) + (1|Subplot) +(1|TSFdays),
                 data = MetaRarePath);summary(Mod6)
anova(Mod5,Mod6)#Mod 6 better, not signif, do not remove


#Remove InitialAshDepth + 
Mod7<-glmer.nb(S.obs ~ Treatment + TSF + Precip +  Treatment*Precip + (1|Plot) + 
          (1|Subplot) +(1|TSFdays), data = MetaRarePath);summary(Mod7)
anova(Mod6,Mod7)#Mod 6 better, not signif, do not remove


#Remove Precip +  
Mod8<-glmer.nb(S.obs ~ Treatment + TSF + Treatment*Precip + (1|Plot) + 
            (1|Subplot) +(1|TSFdays), data = MetaRarePath);summary(Mod8)
anova(Mod6,Mod8)#Mod 8 better, not signif, do not remove


#Remove TSF + 
Mod9<-glmer.nb(S.obs ~ Treatment + Treatment*Precip + (1|Plot) + 
          (1|Subplot) +(1|TSFdays), data = MetaRarePath);summary(Mod9)
anova(Mod8,Mod9)#Mod 9 better, not signif, do not remove

#Remove Treatment + 
Mod10<-glmer.nb(S.obs ~ Treatment*Precip + (1|Plot) + (1|Subplot) +(1|TSFdays), 
                 data = MetaRarePath);summary(Mod10)
anova(Mod9,Mod10)#Same AIC, not signif, do not remove


#Compare all models...................................................
anova(Mod1, Mod2, Mod3, Mod4, Mod5, Mod6, Mod7, Mod8, Mod9,Mod10) 
AICc(Mod1, Mod2, Mod3, Mod4, Mod5, Mod6, Mod7, Mod8, Mod9,Mod10)
#Mod 10 better


#----
#----
#*********************************************************************************************----
#Final Model.....................................................................................
#*********************************************************************************************----

FullNB<-glmer.nb(S.obs ~ Treatment*Precip + (1|Plot) + (1|Subplot) +(1|TSFdays), 
                 data = MetaRarePath);summary(FullNB)
R2Full<-r.squaredGLMM(FullNB);R2Full#Delta:0.1134819 0.15549934



#Model Diagnostic & Validation plots ...........................................
par(mfrow=c(1,1));plot(FullNB)
Fitted <- fitted(FullNB) #command fitted gives already e^model
Res <- resid(FullNB, type = "pearson")

plot(x = Fitted,  y = Res,xlab = "Fitted values",ylab = "Pearson residuals")
abline(h = 0, lty = 2, col="red")


#Export model results.............................................................................
capture.output(summary(FullNB),file="1-Analysis/Diversity/Alpha/Guilds/Tables/Path-NegBinom-FullModel.csv")
write.csv(R2Full,"1-Analysis/Diversity/Alpha/Guilds/Tables/Path-NegBinom-R2FullModel.csv")







#.----
#.----
#*************************************************************************************************************----
# .................SIGNIFICANCE PER VARIABLE, INDEPENDENTLY ..................................................
#*************************************************************************************************************----

#Test treatment * tsf significance..................................................................................
TrTSFNB<-glmer.nb(S.obs ~ Treatment * TSF2 + (1|Plot) + (1|Subplot) + (1|TSFdays),data = MetaRarePath);summary(TrTSFNB)
R2TrTSF<-r.squaredGLMM(TrTSFNB);R2TrTSF#Delta:0.2323200, 0.5385855
plot(TrTSFNB)

#Get pairwise comparisons w differences indicated as letters.......
UTrtTSF<- emmeans(TrTSFNB, ~ Treatment*TSF2)
TrtTSFPost<-pairs(UTrtTSF);TrtTSFPost

#Export model results..........................................................................................
capture.output(summary(TrTSFNB),file="1-Analysis/Diversity/Alpha/Tables/Significance/Trt-TSF-NegBinom.csv")
write.csv(R2TrTSF,"1-Analysis/Diversity/Alpha/Tables/Significance/Trt-TSF-R2-NegBinom.csv")
write.csv(TrtTSFPost,"1-Analysis/Diversity/Alpha/Tables/Significance/PostHoc-Trt-TSF-NegBinom-.csv")




#Test treatment * Precip significance..........................................................................................
PrecipNB<-glmer.nb(S.obs ~ Treatment * Precip + TSF * Precip +(1|Plot)+(1|Subplot)+(1|TSFdays),data = MetaRarePath);summary(PrecipNB)
R2Precip<-r.squaredGLMM(PrecipNB);R2Precip#Delta:0.2856056 0.5611148
plot(TrPrecipNB)

#Get pairwise comparisons w differences indicated as letters.......
UPrecip<- emmeans(PrecipNB, ~ Treatment*Precip);UPrecip
pairs(UPrecip)

#Export model results..........................................................................................
capture.output(summary(PrecipNB),file="1-Analysis/Diversity/Alpha/Tables/Significance/Precip-Trt-TSF-NegBinom.csv")
write.csv(R2Precip,"1-Analysis/Diversity/Alpha/Tables/Significance/Precip-R2-NegBinom.csv")



#Test treatment * AshDepth significance...................................................................................................
AshDepthNB<-glmer.nb(S.obs ~ Treatment * AshDepth + TSF * AshDepth +(1|Plot)+(1|Subplot)+(1|TSFdays),data = MetaRarePath);summary(AshDepthNB)
R2AshDepth<-r.squaredGLMM(AshDepthNB);R2AshDepth#Delta:0.2887868 0.5510319
plot(TrAshDepthNB)

#Export model results..........................................................................................
capture.output(summary(AshDepthNB),file="1-Analysis/Diversity/Alpha/Tables/Significance/AshDepth-Trt-TSF-NegBinom.csv")
write.csv(R2AshDepth,"1-Analysis/Diversity/Alpha/Tables/Significance/AshDepth-R2-NegBinom.csv")


























#-----
#Model to compare S.obs at each timepoint, between treatments ....................................
#Model requires TSF as a factor for comparison
TrtTSF<-glmer.nb(S.obs ~ Treatment + TSF2 + (1|Plot) +  (1|Subplot) + (1|TSFdays),
                 data = MetaRarePath);summary(TrtTSF)
R2TrtTSF<-r.squaredGLMM(TrtTSF);R2TrtTSF#


#Get pairwise comparisons w differences indicated as letters.......
m_means2 <- emmeans(TrtTSF, ~ Treatment*TSF2)
pairs(m_means2, Letters = letters)
pwpp(m_means2)

library(ggeffects)
pr1 <- ggpredict(FullNB, c("AshDepth"))
plot(pr1)

#Get the confidence interval levels
(est <- cbind(Estimate = coef(FullNB), confint(FullNB)))
exp(est)#Incidence ratio



