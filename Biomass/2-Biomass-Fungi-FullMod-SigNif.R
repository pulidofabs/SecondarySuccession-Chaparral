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
MetaFun<-read.csv("Analysis/Tables/Fungi-Biomass.csv",na.strings = "NA", header = TRUE)


#----
#----
#****************************************************************************************************************-----
#---------------------- QUALITY CONTROL ------------------------------------------------------------------------------
#****************************************************************************************************************-----
MetaFun$TSFdays <- as.numeric(MetaFun$TSFdays)
MetaFun$TotPrecip <- as.numeric(MetaFun$TotPrecip)
MetaFun$Treatment <- as.factor(MetaFun$Treatment)
MetaFun$Treatment <- try(relevel(MetaFun$Treatment , "Unburned"))
str(MetaFun)


#----
#****************************************************************************************************************-----
#---------------------- FUNGI-MODEL SELECTION ------------------------------------------------------------------------------
#****************************************************************************************************************-----

df <- MetaFun %>%
  group_by(Treatment) %>%
  summarise(means = mean(CopyNumber),
            variance = var(CopyNumber));df #NegNB 
#will compare different models just in case

attach(MetaFun)
#Start with a poisson model ....................................................
Pois1 <- glm(CopyNumber ~ Treatment,data = MetaFun,family = poisson(link = log))
ods <- Pois1$deviance / Pois1$df.residual;ods#0verdispersion 45085.69
par(mfrow=c(2,2))
plot(Pois1)#look at plots, cooks, levenes etc

.
#Zero inflaited model ..........................................................
sum(CopyNumber == 0) #2


#Test different identity link to see if disperssion decreases ..................
Pois2 <- glm(CopyNumber ~ Treatment,family = poisson(link = sqrt))
ods2<- Pois2$deviance / Pois2$df.residual;ods2 #45085.69

Pois3 <- glm(CopyNumber ~ Treatment,family = poisson(link = identity))
ods3<- Pois3$deviance / Pois3$df.residual;ods3#45085.69

Gam <- glm(CopyNumber ~ Treatment,family = "Gamma")
odsGam<-Gam$deviance/Gam$df.residual;odsGam#

#Negative binomial... ..........................................................
#* mean and variance suggest this is the model to use
nb1 <- glm.nb(CopyNumber ~Treatment,data = MetaFun)
ods_nb <- (nb1$deviance / nb1$df.residual);ods_nb#1.18, mild= acceptable
plot(nb1)

#Compare models.................................................................
AICc(Pois1, Pois2, Pois3,nb1)#neg binomial better 




#----
#--------
#******************************************************************************************************-----
#------------------ TEST OF SIGNIFICANCE-------------------------------------------------------
#*****************************************************************************************************-----
#Will perform a negative binomial model w glmer to control for random effect
#r2m=marginal (proportion of variance explained by the fixed factor(s) alone)
#r2c=conditional(which describes the proportion of variance explained by both 
#the fixed and random factors)

attach(MetaFun) 
summary(MetaFun)
str(MetaFun)

#Rescale variables..............................................................................
#Can use the center= True and scale=True but that is the normal so do not need to specify
# .....Example scale(MetaFun$TotPrecip, center = TRUE, scale = TRUE)
#......Scale, center = gives you a central point along which data is centered around,
#......Scale, scale=True means tart it will subtract the mean and divide by the std deviation

#Rescale and center the numerical data
#MetaFun$TSF2<-as.factor(MetaFun$TSFdays)
MetaFun$TSF<-scale(MetaFun$TSFdays);summary(MetaFun$TSF)
MetaFun$Precip<-scale(MetaFun$TotPrecip)



#Explore Temporal correlations to ensure that we are capturing as much of the variance 
ASVneg1 <- glmer.nb(CopyNumber ~ (1|Plot)+(1|Subplot)+(1|TSFdays),data = MetaFun);summary(ASVneg1)
ASVneg2 <- glmer.nb(CopyNumber ~ (1|Plot)+(1|TSFdays),data = MetaFun);summary(ASVneg2)
ASVneg3 <- glmer.nb(CopyNumber ~ (1|Subplot)+(1|TSFdays),data = MetaFun);summary(ASVneg3)#does not converge


#MuMIN for model selection.............................
library(MuMIn)
AICc(ASVneg1,ASVneg2)#ASVneg2 better, lower AIC



#Run Full and Null model.............................................................................
Full<-glmer.nb(CopyNumber~Treatment*TSF + Treatment*InitialAshDepth + Treatment*Precip +
            TSF*InitialAshDepth + TSF*Precip + Treatment + TSF + Precip + InitialAshDepth +
            (1|Plot) + (1|TSFdays),data = MetaFun);summary(Full)

Null<-glmer.nb(CopyNumber ~ (1|Plot)+ (1|TSFdays), data=MetaFun);summary(Null)
anova(Full, Null)#Full model better



#Remove InitialAshDepth +
Mod2<-glmer.nb(CopyNumber~Treatment*TSF + Treatment*InitialAshDepth + Treatment*Precip +
         TSF*InitialAshDepth + TSF*Precip + Treatment + TSF + Precip + 
        (1|Plot) + (1|TSFdays),data = MetaFun);summary(Mod2)
anova(Full, Mod2)#Full better, signif, do not remove


#Remove TSF*Precip + 
Mod3<-glmer.nb(CopyNumber~InitialAshDepth +Treatment*TSF + Treatment*InitialAshDepth +
              Treatment*Precip + TSF*InitialAshDepth + Treatment +  TSF +
              Precip + (1|Plot) + (1|TSFdays),data = MetaFun);summary(Mod3)
anova(Full, Mod3)#Mod 3 better,not signif, remove



#Remove  Precip + 
Mod4<-glmer.nb(CopyNumber~InitialAshDepth +Treatment*TSF + Treatment*InitialAshDepth +
                 Treatment*Precip + TSF*InitialAshDepth + Treatment + TSF +
                 (1|Plot) + (1|TSFdays),data = MetaFun);summary(Mod4)
anova(Mod3, Mod4)#Same AIC, not signif, remove


#Remove Treatment*InitialAshDepth +
Mod5<-glmer.nb(CopyNumber~InitialAshDepth +Treatment*TSF + 
                 Treatment*Precip + TSF*InitialAshDepth + Treatment + TSF +
                 (1|Plot) + (1|TSFdays),data = MetaFun);summary(Mod5)
anova(Mod4, Mod5)#Same AIC, not signif, remove


#Remove Treatment*TSF + 
Mod6<-glmer.nb(CopyNumber~InitialAshDepth +Treatment*Precip +
                  TSF*InitialAshDepth + Treatment + TSF +
                 (1|Plot) + (1|TSFdays),data = MetaFun);summary(Mod6)
anova(Mod5, Mod6)#Mod 5 better not significant, remove


#Remove Treatment*Precip +
Mod7<-glmer.nb(CopyNumber~InitialAshDepth +
                 TSF*InitialAshDepth + Treatment + TSF +
                 (1|Plot) + (1|TSFdays),data = MetaFun);summary(Mod7)
anova(Mod6, Mod7)#Mod 7 better, not significant, remove


#Remove TSF +
Mod8<-glmer.nb(CopyNumber~InitialAshDepth +
                 TSF*InitialAshDepth + Treatment + 
                 (1|Plot) + (1|TSFdays),data = MetaFun);summary(Mod8)
anova(Mod7, Mod8)#Same AIC, not signif, remove


#Remove Treatment +
Mod9<-glmer.nb(CopyNumber~InitialAshDepth +TSF*InitialAshDepth +  
                 (1|Plot) + (1|TSFdays),data = MetaFun);summary(Mod9)
anova(Mod8, Mod9)#Mod 8 bette, Significant, do not remove...............................

#Remove InitialAshDepth +
Mod10<-glmer.nb(CopyNumber~TSF*InitialAshDepth + Treatment + (1|Plot) +
                (1|TSFdays),data = MetaFun);summary(Mod10)
anova(Mod8, Mod10)#Same AIC, not signif, remove

#Remove TSF*InitialAshDepth + 
Mod11<-glmer.nb(CopyNumber~Treatment + (1|Plot) +
                  (1|TSFdays),data = MetaFun);summary(Mod11)
anova(Mod10, Mod11)#Mod 10 bette, Significant, do not remove...............................



#Compare all models..........................................................................
anova( Mod2, Mod3, Mod4, Mod5, Mod6, Mod7, Mod8, Mod9, Mod10, Mod11)
AICc(Mod1, Mod2, Mod3, Mod4, Mod5, Mod6, Mod7, Mod8)
#Mod 10 suggested to be the best model, lower AIC values


#Final Model............................................................................................................
FinalNB<-glmer.nb(CopyNumber~TSF*InitialAshDepth + Treatment + (1|Plot) +
                    (1|TSFdays),data = MetaFun);summary(FinalNB)
R2Full<-r.squaredGLMM(FinalNB);R2Full#Delta;0.4815221 0.5728228



#Model Diagnostic & Validation plots ...........................................
par(mfrow=c(1,1));plot(FinalNB)
Fitted <- fitted(FinalNB) #command fitted gives already e^model
Res <- resid(FinalNB, type = "pearson")

plot(x = Fitted,  y = Res,xlab = "Fitted values",ylab = "Pearson residuals")
abline(h = 0, lty = 2, col="red")


#Export model results.............................................................................
capture.output(summary(FinalNB),file="Analysis/Tables/NegBinom-FullModel-Fungi.csv")
write.csv(R2Full,"Analysis/Tables/NegBinom-R2FullModel-Fungi.csv")




