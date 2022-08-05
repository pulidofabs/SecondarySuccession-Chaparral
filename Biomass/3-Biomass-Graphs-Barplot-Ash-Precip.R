setwd("~/Dropbox/1-Dissertation/1-Research/HolyFire/1-HolyFire-Pub1")

library(ggpubr) 

#Fungal Biomass
BioFun<-read.csv("Metadata/MetaFun-Biomass.csv")
BioFun$Treatment<-as.factor(BioFun$Treatment)
BioFun$Treatment <- try(relevel(BioFun$Treatment , "Unburned"))


#Bacteria Biomass
BioBac<-read.csv("Metadata/MetaBac-Biomass.csv")
BioBac$Treatment<-as.factor(BioBac$Treatment)
BioBac$Treatment <- try(relevel(BioBac$Treatment , "Unburned"))



#----
#------
#****************************************************************************************************************----
#*----------------- GRAPHS BIOMASS ---------------------------------------------------------------------
#*#**************************************************************************************************************----

#Fungal biomass ...............................................................
BioFunTrt<-ggboxplot(BioFun, x="Treatment", y="CopyNumber") +   
  geom_boxplot(aes(fill=Treatment))+ 
  theme_bw() +  ylab("") + 
  theme(panel.background = element_blank(),
        panel.grid= element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.position = "",
        text = element_text(size=25),
        axis.title.x=element_blank(),
        axis.text.y = element_text(size=20, angle=90, 
              hjust=.5,vjust=.5,colour = "black"),
        axis.text.x = element_text(size=25, colour = "black"))+ 
  scale_fill_manual(values=c("#45877f","#a2673f"));BioFunTrt



#Bacteria biomass ...............................................................
BioBacTrt<-ggboxplot(BioBac, x="Treatment", y="CopyNumber") +   
  geom_boxplot(aes(fill=Treatment))+ 
  theme_bw() +  ylab("Gene Copy Number") + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.position = "",
        text = element_text(size=25),
        axis.title.x=element_blank(),
        axis.text.y = element_text(size=20, angle=90, 
              hjust=.7,vjust=.7,colour = "black"),
        axis.text.x = element_text(size=25, colour = "black"))+ 
  scale_y_continuous(labels = scales::scientific)+
  scale_fill_manual(values=c("#45877f","#a2673f"));BioBacTrt


pdf("1a-PanelsBac-Fun/Graphs/Alpha/Biomass-BacFun.pdf", height=10, width=6,onefile=FALSE)
ggarrange(BioBacTrt,BioFunTrt,ncol=2, nrow=1, common.legend = TRUE, 
          legend="bottom", labels = c("a", "b"), hjust = c(-.8,-1),
          align="hv",font.label = list(size=16,face="bold"))
dev.off()


#----
#----
#**************************************************************************************************----
#------------ ASH DEPTH
#**************************************************************************************************----
#Creat color scheme
colorTP<-c("#de0025","#b05800","#ff6a00","#ffd000","#fac2e3","black","#0a2bff","#00b5f7","#007d2a")

#Maintain same color for all graphs .....................................
colorPlot<-c("black","#5e4839","#947f70","#c49f68","#bfb3aa","#ebddc7","#235978","#45877f","#aed4ca")

attach(BioBac)
AshBioBac<-ggplot(BioBac, aes(x=InitialAshDepth, y=CopyNumber))  +
  geom_smooth(method = MASS::glm.nb, se = TRUE)+ 
  scale_shape_manual(values=c(15,24,19,1,18,25,8,17,13))+
  geom_point(size=2,aes(shape = as.factor(TSFdays), colour = Treatment))+
  stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label..,
                                           sep = "~~~~")), formula = y~x, label.x = 3, label.y =1550000, size=6)+
  scale_color_manual(values=  c("#45877f","#a2673f")) +
  scale_fill_manual(values=colorTP) +
  theme_bw()+ 
  theme(panel.grid = element_blank(),
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
  labs(x = "Ash Depth (cm)",y="Gene Copy Number");AshBioBac
detach(BioBac)



attach(BioFun)
AshBioFun<-ggplot(BioFun, aes(x=InitialAshDepth, y=CopyNumber))  +
  geom_smooth(method = MASS::glm.nb, se = TRUE)+ 
  scale_shape_manual(values=c(15,24,19,1,18,25,8,17,13))+
  geom_point(size=2,aes(shape = as.factor(TSFdays), colour = Treatment))+
  stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label..,
      sep = "~~~~")), formula = y~x, label.x = 3, label.y =600000, size=6)+
  scale_color_manual(values=  c("#45877f","#a2673f"))+
  scale_fill_manual(values=colorTP) +
  theme_bw()+ 
  theme(panel.grid = element_blank(),
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
detach(BioFun)




#Create panel...#A and BLabels are not displayin properties.............................
AshPlots<-ggarrange(AshBioBac,AshBioFun,ncol=2, nrow=1, common.legend = TRUE,
             legend="bottom", align="hv", labels = c("A Bacteria","B Fungi"), 
             font.label = list(size=24,face="bold"));AshPlots

pdf("Analysis/Graphs/AhDepthBiomass-Bac-Fun.pdf", height=6, width=10)
AshPlots
dev.off()


#----
#----
#***********************************************************************************************************----
#------------------ TPRECIPITATION GRAPHS  ----------------------------------------------------------------------
#***********************************************************************************************************----
#Put two separate graphs together but apart....................................................
library(scales)#changes legend from sci to comma
library(patchwork)#to put graphs together

#Quality control
BioBac$Plot<- as.factor(BioBac$Plot)
BioBac$TSFdays<- as.factor(BioBac$TSFdays)
BioBac$Treatment<-as.factor(BioBac$Treatment)
BioBac$InitialAshDepth<-as.factor(BioBac$InitialAshDepth)
BioBac$Treatment <- try(relevel(BioBac$Treatment , "Unburned"))#relevel so that unburned is the base level for comparison

BioFun$Plot<- as.factor(BioFun$Plot)
BioFun$TSFdays<- as.factor(BioFun$TSFdays)
BioFun$Treatment<-as.factor(BioFun$Treatment)
BioFun$InitialAshDepth<-as.factor(BioFun$InitialAshDepth)
BioFun$Treatment <- try(relevel(BioFun$Treatment , "Unburned"))#relevel so that unburned is the base level for comparison


#FUNGI ........................................................................................
P1<-ggplot(BioFun, aes(x=TSFdays,y=CopyNumber, group=Treatment, col=Treatment)) +
  stat_summary(fun = "mean", size = 2, geom = "point")+
  stat_summary(fun=mean,geom="line", size=1.5, aes(linetype = Treatment),show.legend=TRUE)+
  stat_summary(fun.data = mean_se,geom = "errorbar",size=.85, 
               alpha=0.7,position = position_dodge(0.01))+
  scale_linetype_manual(values=c("solid","twodash"))+
  scale_color_manual(values=c("#45877f","#a2673f"))+
  theme_bw()+  
theme(panel.grid  = element_blank(),
     axis.title.y = element_text(size = 19, colour = "black"),
     axis.text.y = element_text(size = 18, angle=90, hjust=.5,colour = "black"), 
     #axis.title.x =  element_text(size =20,  colour = "black"), 
     #axis.text.x = element_text(size = 19, colour = "black"), 
     axis.ticks.x=element_blank(),axis.text.x = element_blank(),
     axis.title.x = element_blank(),
     legend.title = element_text(size=20),
     legend.text=element_text(size = 19),
     legend.position = "botttom")+ ylab("")+
  scale_y_continuous(labels = comma);P1

#Just quick and easy way to look but does not account for nestedness..look at 
#signif scrip
P2 <-P1 + stat_compare_means(method = "kruskal.test", 
           label="p.signif", label.y=300000,show.legend = FALSE );P2





#BACTERIA..........................................................................................
P1B<-ggplot(BioBac, aes(x=TSFdays,y=CopyNumber, group=Treatment, col=Treatment)) +
  stat_summary(fun = "mean", size = 2, geom = "point")+
  stat_summary(fun=mean,geom="line", size=1.5, aes(linetype = Treatment),show.legend=TRUE)+
  stat_summary(fun.data = mean_se,geom = "errorbar",size=.85, 
               alpha=0.7,position = position_dodge(0.01))+
  scale_linetype_manual(values=c("solid","twodash"))+
  scale_color_manual(values=c("#45877f","#a2673f"))+
  theme_bw()+  
  theme(panel.grid  = element_blank(),
    axis.title.y = element_text(size = 19, colour = "black"),
    axis.text.y = element_text(size = 18, angle=90, hjust=.5,colour = "black"), 
    #axis.title.x =  element_text(size =20,  colour = "black"), 
    #axis.text.x = element_text(size = 19, colour = "black"), 
    axis.ticks.x=element_blank(),axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_text(size=20),
    legend.text=element_text(size = 19),
    legend.position = "botttom")+
  ylab("Gene Copy Number");P1B


P2B <-P1B + stat_compare_means(method = "kruskal.test", 
          label="p.signif", label.y=1000000,show.legend = FALSE );P2B


# RAIN GRAPH.........................................................................
#for adding below graphs..............................................
Rain<-ggplot(BioFun, aes(x=TSFdays,y=CopyNumber,group=Treatment, col=Treatment))+
  geom_line(aes(y=TotPrecip),col="#94C5E0",  size=1) +
  theme_bw()+ 
  theme(panel.grid  = element_blank(),
    axis.text.y = element_text(size = 15,angle=90, hjust=1,color="black"),  
    axis.title.x =  element_text(size =20,  colour = "black"), 
    axis.text.x = element_text(size = 22, colour = "black"), 
    axis.title = element_text(size = 18),
    legend.title = element_text(size=20),
    legend.text=element_text(size = 19),
    legend.position = "none")+
  scale_fill_manual(values ="red")+
  ylab("")+  xlab ("Time Since Fire (days)");Rain


#Add the precip graphs..............................................
Bac<-BacRain2/Rain + plot_layout(nrow = 2, heights = c(4,1),
            guides = 'collect') & theme(legend.position = 'bottom');Bac


Fun<-FunRain2/Rain + plot_layout(nrow = 2, heights = c(4,1),
           guides = 'collect') & theme(legend.position = 'bottom');Fun


#Put the graphs together...........................................
All<-P2B + Fun + plot_layout(nrow = 2, ncol=1,heights = c(1.6,2),
           guides = 'collect') & theme(legend.position = 'bottom');All

#Export graph..............................................................
pdf("Analysis/Graphs/Rain-Bio-Biomass.pdf", height = 10.5, width = 6.6)
All
dev.off()







#----
#----  
#*********************************************************************************************************----
#------------- SIGNIFICANCE -------------------------------------------------------------------------
#*********************************************************************************************************----
#Start with model validation

#Load required libraries
library(MASS)
library(lme4)
library(rcompanion)
library(MuMIn)#pseudo R2
library(multcomp)#post-hoc tukey
library(pscl)#zero inflaited

# BACTERIA ..............................................................................

# 1. Poisson
Pois1 <- glm(CopyNumber ~ ,data = BioBac, family = poisson(link = log))

#check overdispersion (should be 1 anything > = overdips)
ods <- Pois1$deviance / Pois1$df.residual;ods #198079.9

#look at plots, cooks, levenes etc
par(mfrow=c(2,2))
plot(Pois1)


#Test different identity link to see if disperssion decreases ...........
Pois2 <- glm(CopyNumber~,family = poisson(link = sqrt), data = BioBac)
ods2<- Pois2$deviance / Pois2$df.residual; ods2 #198079.9

Pois3 <- glm(CopyNumber~,family = poisson(link = identity), data = BioBac)
ods3<- Pois3$deviance / Pois3$df.residual;ods3 #198079.9

Gam <- glm(CopyNumber~,family = "Gamma", data = BioBac)#non positive values cannot use
odsGam<-Gam$deviance/Gam$df.residual;odsGam #0.93

#Could it be a Zero inflated model, check how many zero's in response
attach(MetaBioBac)
sum(CopyNumber == 0)#0 not a zero inflaited model
#100*sum(CopyNumber == 0)/nrow(MetaSoil)#--shows that 70% of your data are zeroes


#** none of above models work, will test negative binomial as 
#* mean and variance sugges this is the model to use
nb1 <- glm.nb(CopyNumber~, data = BioBac)
ods_nb <- (nb1$deviance / nb1$df.residual);ods_nb#dispersion 1.13
plot(nb1)#look at residuals, leverage etc

#Compare models
AIC(nb1, Pois1, Pois2,Pois3, Gam)#Negative binomial better for bacteria 



# FUNGI..............................................................................
# 1. Poisson
Pois1 <- glm(CopyNumber ~ ,data = BioFun, family = poisson(link = log))


#check overdispersion (should be 1 anything > = overdips)
ods <- Pois1$deviance / Pois1$df.residual;ods #45085.69

#look at plots, cooks, levenes etc
par(mfrow=c(2,2))
plot(Pois1)


#Test different identity link to see if disperssion decreases ...........
Pois2 <- glm(CopyNumber~,family = poisson(link = sqrt), data = BioFun)
ods2<- Pois2$deviance / Pois2$df.residual; ods2 #45085.69

Pois3 <- glm(CopyNumber~,family = poisson(link = identity), data = BioFun)
ods3<- Pois3$deviance / Pois3$df.residual;ods3 #45085.69

#Gam <- glm(CopyNumber~Treatment,family = "Gamma",MetaBioFun)# Non positive num cannot use
#odsGam<-Gam$deviance/Gam$df.residual;odsGam #0.93

#Could it be a Zero inflated model, check how many zero's in response
attach(MetaBioFun)
sum(CopyNumber == 0)#2 not a zero inflaited model
100*sum(CopyNumber == 0)/nrow(MetaBioFun)#--shows that 70% of your data are zeroes


#** none of above models work, will test negative binomial as 
#* mean and variance sugges this is the model to use
nb1 <- glm.nb(CopyNumber~,data = BioFun)
ods_nb <- (nb1$deviance / nb1$df.residual);ods_nb#dispersion 1.18
plot(nb1)#look at residuals, leverage etc

#Compare models
AIC(nb1, Pois1, Pois2,Pois3)#Negative binomial better for Fungi



#******TEST NESTEDNESS LEVEL***************************************************************************-----
#Rescale variables..............................................................................
#Can use the center= True and scale=True but that is the normal so do not need to specify
# .....Example scale(MetaBac$TotPrecip, center = TRUE, scale = TRUE)
#......Scale, center = gives you a central point along which data is centered around,
#......Scale, scale=True means tart it will subtract the mean and divide by the std deviation

#Rescale and center the numerical data
#BioBac$TSF2<-as.factor(BioBac$TSFdays)
Bac$TSF<-scale(BioBac$TSFdays);summary(BioBac$TSF)
BioBac$Precip<-scale(BioBac$TotPrecip)



#Explore Temporal correlations to ensure that we are capturing as much of the variance 
ASVneg1 <- glmer.nb(CopyNumber ~ (1|Plot)+(1|Subplot)+(1|TSFdays),data = BioBac);summary(ASVneg1)
ASVneg2 <- glmer.nb(CopyNumber ~ (1|Plot)+(1|Subplot),data = BioBac);summary(ASVneg2)
ASVneg3 <- glmer.nb(CopyNumber ~ (1|Plot)+(1|TSFdays),data = BioBac);summary(ASVneg3)
ASVneg4 <- glmer.nb(CopyNumber ~ (1|Subplot)+(1|TSFdays),data = BioBac);summary(ASVneg4)#does not converge

#Used only to trouble shoot to check optomizer
#Mod1@optinfo[c("optimizer","control")]


#MuMIN for model selection....................................................
library(MuMIn)
AICc(ASVneg1,ASVneg2,ASVneg3)#ASVneg3 better, lower AIC (plot & TSF)
#ASV neg 1 very similar, so for consisteny with other analsys
#I will go with neg1

#Fungal nestedness level.....................
#Explore Temporal correlations to ensure that we are capturing as much of the variance 
ASVneg1F <- glmer.nb(CopyNumber ~ (1|Plot)+(1|Subplot)+(1|TSFdays),data = BioFun);summary(ASVneg1F)
ASVneg2F <- glmer.nb(CopyNumber ~ (1|Plot)+(1|Subplot),data = BioFun);summary(ASVneg2F)
ASVneg3F <- glmer.nb(CopyNumber ~ (1|Plot)+(1|TSFdays),data = BioFun);summary(ASVneg3F)
ASVneg4F <- glmer.nb(CopyNumber ~ (1|Subplot)+(1|TSFdays),data = BioFun);summary(ASVneg4F)#does not converge

#Used only to trouble shoot to check optomizer
#Mod1@optinfo[c("optimizer","control")]


#MuMIN for model selection....................................................
AICc(ASVneg1F,ASVneg2F,ASVneg3F,ASVneg4F)

#----
#----
#*************************************************************************************************************************-----
#........ACTUAL TEST FOR SIGNIFICANCE ..............................................................................................
#*************************************************************************************************************************-----

#Time since fire-* ash depth- biomass
BacBioAsh<-glmer.nb(CopyNumber ~ InitialAshDepth*TSFdays*Treatment+ (1|Plot)+(1|Subplot)+(1|TSFdays), BioBac)
summary(BacBioAsh)#
R2BacBioAsh<-r.squaredGLMM(BacBioAsh);R2BacBioAsh#delta 0.5079132 0.5889870

#Fungi
BacBioAshF<-glmer.nb(CopyNumber ~ InitialAshDepth*TSFdays*Treatment+ (1|Plot)+ (1|TSFdays), BioFun)
summary(BacBioAshF)#
R2BacBioAshF<-r.squaredGLMM(BacBioAshF);R2BacBioAshF#0.6513008 0.6968298




#EXPORT RESULTS.......................................................................................
#Bacteria...........................................................................................................
dir.create("Analysis/Tables/Significance")
capture.output(BacBioTrt, file="Analysis/Tables/Significance/GLMNB-BacBiomass-TrtAshTSF.csv")
write.csv(R2BacBioTrt, "Analysis/Tables/Significance/R2-GLMNB-BacBiomass-TrtAshTSF.csv")

#Fungi...........................................................................................................
capture.output(BacBioTrt, file="Analysis/Tables/Significance/GLMNB-FunBiomass-TrtAshTSF.csv")
write.csv(R2BacBioTrt, "Analysis/Tables/Significance/R2-GLMNB-FunBiomass-TrtAshTSF.csv")

