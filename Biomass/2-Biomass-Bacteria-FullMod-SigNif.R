#10-6-2021

#Reset R's Brain
rm(list=ls())

#set working directory...........................................................
setwd("C:/Users/fabipc/Dropbox/3-HolyFire/1-HolyFire-Pub1/4-Qpcr")


#load libraries..................................................................
library("ggplot2")
library("dplyr")
library(lme4)
library(lmerTest)
library(MASS) #for negative binomial
library(dplyr)

#load in data....................................................................................................
MetaBac<-read.csv("Analysis/Tables/Bacteria-Biomass.csv",na.strings = "NA", header = TRUE)


#----
#----
#****************************************************************************************************************-----
#---------------------- QUALITY CONTROL ------------------------------------------------------------------------------
#****************************************************************************************************************-----
#Bacteria...............................................................
MetaBac$TSFdays <- as.numeric(MetaBac$TSFdays)
MetaBac$TotPrecip <- as.numeric(MetaBac$TotPrecip)
MetaBac$Treatment <- as.factor(MetaBac$Treatment)
MetaBac$Treatment <- try(relevel(MetaBac$Treatment , "Unburned"))
str(MetaBac)


#----
#****************************************************************************************************************-----
#---------------------- FUNGI-MODEL SELECTION ------------------------------------------------------------------------------
#****************************************************************************************************************-----

df <- MetaBac %>%
  group_by(Treatment) %>%
  summarise(means = mean(CopyNumber),
            variance = var(CopyNumber));df #NegNB 
#will compare different models just in case

attach(MetaBac)
#Start with a poisson model ....................................................
Pois1 <- glm(CopyNumber ~ Treatment,data = MetaBac,family = poisson(link = log))
ods <- Pois1$deviance / Pois1$df.residual;ods#0verdispersion 198079.9
par(mfrow=c(2,2))
plot(Pois1)#look at plots, cooks, levenes etc

.
#Zero inflaited model ..........................................................
sum(CopyNumber == 0) #2


#Test different identity link to see if disperssion decreases ..................
Pois2 <- glm(CopyNumber ~ Treatment,family = poisson(link = sqrt))
ods2<- Pois2$deviance / Pois2$df.residual;ods2 #198079.9

Pois3 <- glm(CopyNumber ~ Treatment,family = poisson(link = identity))
ods3<- Pois3$deviance / Pois3$df.residual;ods3#198079.9

Gam <- glm(CopyNumber ~ Treatment,family = "Gamma")#0.93
odsGam<-Gam$deviance/Gam$df.residual;odsGam#

#Negative binomial... ..........................................................
#* mean and variance suggest this is the model to use
nb1 <- glm.nb(CopyNumber ~Treatment,data = MetaBac)
ods_nb <- (nb1$deviance / nb1$df.residual);ods_nb#1.13, mild= acceptable
plot(nb1)

#Compare models.................................................................
AICc(Pois1, Pois2, Pois3,nb1,Gam)#neg binomial better 




#----
#--------
#******************************************************************************************************-----
#------------------ TEST OF SIGNIFICANCE-------------------------------------------------------
#*****************************************************************************************************-----
#Will perform a negative binomial model w glmer to control for random effect
#r2m=marginal (proportion of variance explained by the fixed factor(s) alone)
#r2c=conditional(which describes the proportion of variance explained by both 
#the fixed and random factors)

attach(MetaBac) 
summary(MetaBac)
str(MetaBac)

#Rescale variables..............................................................................
#Can use the center= True and scale=True but that is the normal so do not need to specify
# .....Example scale(MetaBac$TotPrecip, center = TRUE, scale = TRUE)
#......Scale, center = gives you a central point along which data is centered around,
#......Scale, scale=True means tart it will subtract the mean and divide by the std deviation

#Rescale and center the numerical data
#MetaBac$TSF2<-as.factor(MetaBac$TSFdays)
MetaBac$TSF<-scale(MetaBac$TSFdays);summary(MetaBac$TSF)
MetaBac$Precip<-scale(MetaBac$TotPrecip)



#Explore Temporal correlations to ensure that we are capturing as much of the variance 
ASVneg1 <- glmer.nb(CopyNumber ~ (1|Plot)+(1|Subplot)+(1|TSFdays),data = MetaBac);summary(ASVneg1)
ASVneg2 <- glmer.nb(CopyNumber ~ (1|Plot)+(1|TSFdays),data = MetaBac);summary(ASVneg2)
ASVneg3 <- glmer.nb(CopyNumber ~ (1|Subplot)+(1|TSFdays),data = MetaBac);summary(ASVneg3)#does not converge

#Used only to trouble shoot to check optomizer
#Mod1@optinfo[c("optimizer","control")]


#MuMIN for model selection.............................
library(MuMIn)
AICc(ASVneg1,ASVneg2)#ASVneg2 better, lower AIC



#Run Full and Null model.............................................................................
Full<-glmer.nb(CopyNumber~Treatment*TSF + Treatment*InitialAshDepth + Treatment*Precip +
                 TSF*InitialAshDepth + TSF*Precip + Treatment + TSF + Precip + InitialAshDepth +
                 (1|Plot) + (1|TSFdays),data = MetaBac);summary(Full)

Null<-glmer.nb(CopyNumber ~ (1|Plot)+ (1|TSFdays), data=MetaBac);summary(Null)
anova(Full, Null)#Full model better

attach(MetaBac)
str(MetaBac)
#Begin backward model selection.............................................................
Mod1<-glmer.nb(CopyNumber ~ Treatment*TSF + Treatment*InitialAshDepth + Treatment*Precip +
                 TSF*InitialAshDepth + TSF*Precip + Treatment + TSF + InitialAshDepth + 
                 Precip + (1|Plot) + (1|TSFdays), MetaBac); summary(Mod1)


#Remove Treatment*TSF + 
Mod2<-glmer.nb(CopyNumber ~ Treatment*InitialAshDepth + Treatment*Precip +
            TSF*InitialAshDepth + TSF*Precip + Treatment + TSF + InitialAshDepth + 
            Precip + (1|Plot) + (1|TSFdays), MetaBac); summary(Mod2)
anova(Mod1,Mod2)#Mod2 better, not significant, remove


#Remove Treatment*InitialAshDepth + 
Mod3<-glmer.nb(CopyNumber ~ Treatment*Precip + TSF*InitialAshDepth + TSF*Precip +
              Treatment + TSF + InitialAshDepth + Precip + (1|Plot) + (1|TSFdays), 
              MetaBac); summary(Mod3)
anova(Mod2,Mod3)#Mod2 better, not significant, remove


#Remove Treatment*Precip + 
Mod4<-glmer.nb(CopyNumber ~ TSF*InitialAshDepth + TSF*Precip + Treatment + TSF + 
          InitialAshDepth + Precip + (1|Plot) + (1|TSFdays), MetaBac); summary(Mod4)
anova(Mod2,Mod4)#Mod2 better, significant, do not remove ..................................


#Remove TSF*InitialAshDepth 
Mod5<-glmer.nb(CopyNumber ~ TSF*Precip + Treatment + TSF + 
                 InitialAshDepth + Precip + Treatment*Precip + (1|Plot) +
                 (1|TSFdays), MetaBac); summary(Mod5)
anova(Mod2,Mod5)#Mod 2 better, significant, do not remove. ...............................


#Remove TSF*Precip + 
Mod6<-glmer.nb(CopyNumber ~ Treatment + TSF + InitialAshDepth + Precip + 
                 Treatment*Precip +TSF*InitialAshDepth + (1|Plot) +
                 (1|TSFdays), MetaBac); summary(Mod6)
anova(Mod2,Mod6)#Mod 2 better, not significant, remove


#Remove Treatment + 
Mod7<-glmer.nb(CopyNumber ~ TSF + InitialAshDepth + Precip + 
                 Treatment*Precip +TSF*InitialAshDepth + (1|Plot) +
                 (1|TSFdays), MetaBac); summary(Mod7)
anova(Mod2,Mod7)#Mod 2 better, not significant, remove


#Remove Precip + 
Mod8<-glmer.nb(CopyNumber ~  InitialAshDepth + TSF + Treatment*Precip +
          TSF*InitialAshDepth + (1|Plot) + (1|TSFdays), MetaBac); summary(Mod8)
anova(Mod2,Mod8)#Mod 2 better, not significant, remove


#Remove  TSF +
Mod9<-glmer.nb(CopyNumber ~ InitialAshDepth +Treatment*Precip + TSF*InitialAshDepth +
           (1|Plot) + (1|TSFdays), MetaBac); summary(Mod9)
anova(Mod2,Mod9)#Mod 2 better, not significant, remove


#Remove InitialAshDepth +
Mod10<-glmer.nb(CopyNumber ~ Treatment*Precip + TSF*InitialAshDepth +
                 (1|Plot) + (1|TSFdays), MetaBac); summary(Mod10)
anova(Mod2,Mod10)#Mod 2 better, not significant, remove



#Compare all models.................................................
anova(Mod1, Mod2, Mod3, Mod4, Mod5, Mod6, Mod7, Mod8, Mod9,Mod10)
AICc(Mod1, Mod2, Mod3, Mod4, Mod5, Mod6, Mod7, Mod8, Mod9,Mod10)
#Mod 2 Better, lower AIC


#----

#***************************************************************************************----
#Final Model............................................................................................................
#***************************************************************************************----

FinalNB<-glmer.nb(CopyNumber ~ Treatment*InitialAshDepth + Treatment*Precip +
            TSF*InitialAshDepth + TSF*Precip + Treatment + TSF + InitialAshDepth + 
            Precip + (1|Plot) + (1|TSFdays), MetaBac);summary(FinalNB)
R2Full<-r.squaredGLMM(FinalNB);R2Full#Delta;0.3030233 0.4646893 



#Model Diagnostic & Validation plots ...........................................
par(mfrow=c(1,1));plot(FullNB)
Fitted <- fitted(FullNB) #command fitted gives already e^model
Res <- resid(FullNB, type = "pearson")

plot(x = Fitted,  y = Res,xlab = "Fitted values",ylab = "Pearson residuals")
abline(h = 0, lty = 2, col="red")


#Export model results.............................................................................
capture.output(summary(FinalNB),file="Analysis/Tables/NegBinom-FullModel-Bacteria.csv")
write.csv(R2Full,"Analysis/Tables/NegBinom-R2FullModel-Bacteria.csv")




