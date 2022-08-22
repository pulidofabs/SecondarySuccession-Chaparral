#8-12-2022
rm(list=ls())#reset working directory

#Set working directory.................................................................
setwd("~/Dropbox/1-Dissertation/1-Research/HolyFire/1-HolyFire-Pub1/1c-Qpcr/Year1")

#Load packages ........................................................................
library(MuMIn)#for AIC selection
library(lme4)#for neg binom


# Load data ............................................................................................
MetaBioBac<-read.csv("Analysis/Metadata/MetaBac-Biomass.csv", header = TRUE)
MetaBioFun<-read.csv("Analysis/Metadata/MetaBioFun-Biomass.csv", header = TRUE)

#----
#----
#******************************************************************************************************************----
#  ............ Model Selection Fungi  ................................................................................
#******************************************************************************************************************----
attach(MetaBioFun)

VarFun <- MetaBioFun %>%
  group_by(Treatment) %>%
  summarise(means = mean(CopyNumberPerGram),
            variance = var(CopyNumberPerGram));VarFun #Variance > mean, will compare diff models just in case


#Overdisperssion needs to be ~1
#Start with a poisson model ....................................................................................
Pois1 <- glm(CopyNumberPerGram ~ Treatment,data = MetaBioFun,family = poisson(link = log))
ods <- Pois1$deviance / Pois1$df.residual;ods#0verdispersion 18034276
par(mfrow=c(2,2))
plot(Pois1)#look at plots, cooks, levenes etc

#Zero inflaited model ..................................................................................
sum(CopyNumberPerGram == 0) #2

#Test different identity link to see if disperssion decreases ..........................................
Pois2 <- glm(CopyNumberPerGram ~ Treatment,family = poisson(link = sqrt))
ods2<- Pois2$deviance / Pois2$df.residual;ods2 #18034276

Pois3 <- glm(CopyNumberPerGram ~ Treatment,family = poisson(link = identity))
ods3<- Pois3$deviance / Pois3$df.residual;ods3#18034276

Gam <- glm(CopyNumberPerGram ~ Treatment,family = "Gamma")
odsGam<-Gam$deviance/Gam$df.residual;odsGam#non-positive values not allowed ..does not work

#Negative binomial... ..................................................................................
#* mean and variance suggest this is the model to use
nb1 <- glm.nb(CopyNumberPerGram ~Treatment,data = MetaBioFun)
ods_nb <- (nb1$deviance / nb1$df.residual);ods_nb#1.19, mild= acceptable
plot(nb1)

#Compare models...........................................................................
AICc(Pois1, Pois2, Pois3,nb1)#nb1=1.086222e+04 best model

#We will proceed with a negative binomial model but need to test the nestedness level
#to make sure we have the correct model before doin any statistics 

#----
#----
#**************NESTEDNESS LEVEL**************************************************************************************************----
#Test null model with different nestedness level

#Explore Temporal correlations to ensure that we are capturing as much of the variance ....................................
ASVneg1 <- glmer.nb(CopyNumberPerGram ~ (1|Plot)+(1|Subplot)+(1|TSFdays), nAGQ=0,data = MetaBioFun);summary(ASVneg1)
ASVneg2 <- glmer.nb(CopyNumberPerGram ~ (1|Plot)+(1|Subplot),nAGQ=0,data = MetaBioFun);summary(ASVneg2)#
ASVneg3 <- glmer.nb(CopyNumberPerGram ~ (1|Plot)+(1|TSFdays), nAGQ=0,data = MetaBioFun);summary(ASVneg3)#
ASVneg4 <- glmer.nb(CopyNumberPerGram ~ (1|Subplot)+(1|TSFdays), nAGQ=0, data = MetaBioFun);summary(ASVneg4)#did not converge
ASVneg5 <- glmer.nb(CopyNumberPerGram ~ (1|Plot), nAGQ=0, data = MetaBioFun);summary(ASVneg5)
ASVneg6 <- glmer.nb(CopyNumberPerGram ~ (1|Subplot),  nAGQ=0,data = MetaBioFun);summary(ASVneg6)#

#MuMIN for model selection....................................................................
AICc(ASVneg1,ASVneg2,ASVneg3,ASVneg4,ASVneg5,ASVneg6)#ASVneg3 better, lower AIC
#Mod 3 better, will continue with Plot and TSF
#ASVneg1  5 10861.68
#ASVneg3  4 10859.78

#Model we will be using for fungal biomass will be 
#Mod<- glmer.nb(CopyNumberPerGram ~ (1|Plot)+(1|TSFdays), nAGQ=0,data = MetaBioFun);summary(Mod)#


#----
#----
#******************************************************************************************************************----
#  ............ Model Selection Bacteria  ................................................................................
#******************************************************************************************************************----
attach(MetaBioBac)

VarBac <- MetaBioBac %>%
  group_by(Treatment) %>%
  summarise(means = mean(CopyNumberPerGram),
            variance = var(CopyNumberPerGram));VarBac #Variance > mean, will compare diff models just in case


#Overdisperssion needs to be ~1
#Start with a poisson model ....................................................................................
Pois1 <- glm(CopyNumberPerGram ~ Treatment,data = MetaBioBac,family = poisson(link = log))
ods <- Pois1$deviance / Pois1$df.residual;ods#0verdispersion 79231978
par(mfrow=c(2,2))
plot(Pois1)#look at plots, cooks, levenes etc

#Zero inflaited model ..................................................................................
sum(CopyNumberPerGram == 0) #0

#Test different identity link to see if disperssion decreases ..........................................
Pois2 <- glm(CopyNumberPerGram ~ Treatment,family = poisson(link = sqrt))
ods2<- Pois2$deviance / Pois2$df.residual;ods2 #79231978

Pois3 <- glm(CopyNumberPerGram ~ Treatment,family = poisson(link = identity))
ods3<- Pois3$deviance / Pois3$df.residual;ods3#79231978

Gam <- glm(CopyNumberPerGram ~ Treatment,family = "Gamma")
odsGam<-Gam$deviance/Gam$df.residual;odsGam#0.93
plot(Gam)

#Negative binomial... ..................................................................................
#* mean and variance suggest this is the model to use
nb1 <- glm.nb(CopyNumberPerGram ~Treatment,data = MetaBioBac)
ods_nb <- (nb1$deviance / nb1$df.residual);ods_nb#1.13, mild= acceptable
plot(nb1)

#Compare models...........................................................................
AICc(Pois1, Pois2, Pois3, Gam, nb1)#nb1 and Gam are the same, but for 
                                  #consistency will proceed with nb1

#We will proceed with a negative binomial model but need to test the nestedness level
#to make sure we have the correct model before doin any statistics 

#----
#----
#**************NESTEDNESS LEVEL**************************************************************************************************----
#Test null model with different nestedness level

#Explore Temporal correlations to ensure that we are capturing as much of the variance ....................................
ASVneg1 <- glmer.nb(CopyNumberPerGram ~ (1|Plot)+(1|Subplot)+(1|TSFdays), nAGQ=0,data = MetaBioBac);summary(ASVneg1)
ASVneg2 <- glmer.nb(CopyNumberPerGram ~ (1|Plot)+(1|Subplot),nAGQ=0,data = MetaBioBac);summary(ASVneg2)#
ASVneg3 <- glmer.nb(CopyNumberPerGram ~ (1|Plot)+(1|TSFdays), nAGQ=0,data = MetaBioBac);summary(ASVneg3)#
ASVneg4 <- glmer.nb(CopyNumberPerGram ~ (1|Subplot)+(1|TSFdays), nAGQ=0, data = MetaBioBac);summary(ASVneg4)#did not converge
ASVneg5 <- glmer.nb(CopyNumberPerGram ~ (1|Plot), nAGQ=0, data = MetaBioBac);summary(ASVneg5)
ASVneg6 <- glmer.nb(CopyNumberPerGram ~ (1|Subplot),  nAGQ=0,data = MetaBioBac);summary(ASVneg6)#

#MuMIN for model selection....................................................................
AICc(ASVneg1,ASVneg2,ASVneg3,ASVneg4,ASVneg5,ASVneg6)#ASVneg3 better, lower AIC
#Mod 3 better, will continue with Plot and TSF
#ASVneg1  5 12098.61
#ASVneg3  4 12097.35

#Model we will be using for Bacgal biomass will be 
#Mod<- glmer.nb(CopyNumberPerGram ~ (1|Plot)+(1|TSFdays), nAGQ=0,data = MetaBioBac);summary(Mod)#


