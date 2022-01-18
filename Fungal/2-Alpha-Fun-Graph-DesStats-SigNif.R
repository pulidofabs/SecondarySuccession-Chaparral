#Last edits:10-6-2021

#reset working directory.....
rm(list=ls())


#Set working directory....................................................................
setwd("C:/Users/fabipc/Dropbox/3-HolyFire/1-HolyFire-Pub1/2-Fungi")


#Load additional packages.................................................................
library(dplyr)
library(lme4)#glmer.nb
library(lmerTest)#to get P value
library(emmeans)#tukey comparisons
library(rcompanion)#pseudo R2
library(MuMIn)#Aic
library(MASS)


#Load packages.............................................................................
MetaSoil<- read.csv("Metadata/Fungi/MetaSoilSpp.csv", na.strings="NA",header = TRUE)



#.----
#.----
#**************************************************************************************************----
# .................QUALITY CONTROL.....................................................................
#**************************************************************************************************----
# * Look at data structure ................................................................
MetaSoil$Plot<- as.factor(MetaSoil$Plot)
MetaSoil$TSFdays<-as.numeric(MetaSoil$TSFdays)
MetaSoil$TSF2<-as.factor(MetaSoil$TSFdays)
MetaSiuk$InitialAshDepth<-as.numeric(MetaSoil$InitialAshDepth)

MetaSoil$Treatment<-as.factor(MetaSoil$Treatment)
MetaSoil$Treatment <- try(relevel(MetaSoil$Treatment,"Unburned")) #change base level

MetaSoil$AshSeverity<-ordered(MetaSoil$AshSeverity,c("Unburned","Low","Moderate","High"))
MetaSoil$Date<-ordered(MetaSoil$Date, c("9/30/18","10/8/18","10/17/18","11/19/18",
                     "12/17/18","1/22/19","3/19/19","6/26/19","9/24/19"))


# *  * * TEST FOR NORMALITY ............................................................
par(mfrow=c(1,3))
hist(MetaSoil$S.obs)
hist(Burned$S.obs)
hist(Unburned$S.obs)


#shapiro test .....................................................
shapiro.test(MetaSoil$S.obs)#2.161e-11------Non-normal
shapiro.test(Burned$S.obs)#2.383e-06--------Non-normal
shapiro.test(Unburned$S.obs)#0.3859---------Normal


#tranform to see if it makes a difference..........................
#None worked--log, log10, sqrt, cubed
shapiro.test(log10(MetaSoil$S.obs))#--------Not normal
shapiro.test(sqrt(Burned$S.obs))#-----------Normalized

#TRANSFORMATIONS...................................................
OtuB<-sqrt(Burned$S.obs)

#Unburned richness is normally distributed
#Burned samples transformed using Sqrt transformation
#Total could not transformed



#-----
#----
#-------
#**********************************************************************************************************----
#DESRIPTIVE STATISTICS-----------------------------------------------------------------------------------------
#**********************************************************************************************************----
library(dplyr)
library(srvyr)#prevents error for dplyr

#Obs
TrtStatsObs<-MetaSoil %>%
  filter(!is.na(S.obs)) %>%
  group_by(Treatment) %>%
  summarise(N = n(),
            Avg_ASV = mean(S.obs, na.rm = TRUE), 
            Variance = var(S.obs, na.rm = TRUE), 
            sd = sd(S.obs), se =sd / sqrt(N)) %>%
  filter(N > 1) %>% as.data.frame();TrtStatsObs

TSFTrTstatsObs<-MetaSoil %>%
  filter(!is.na(S.obs)) %>%
  group_by(Treatment, TSFdays) %>%
  summarise(N = n(),Avg_ASV = mean(S.obs, na.rm = TRUE), 
            Variance = var(S.obs, na.rm = TRUE), sd = sd(S.obs), 
            se =sd / sqrt(N)) %>%
  filter(N > 1) %>% as.data.frame();TSFTrTstatsObs



#Simpson
TrtStatSimp<-MetaSoil %>%
  filter(!is.na(simpson)) %>%
  group_by(Treatment) %>%
  summarise(N = n(),
      Avg_ASV = mean(simpson, na.rm = TRUE), 
      Variance = var(simpson, na.rm = TRUE), 
      sd = sd(simpson), se =sd / sqrt(N)) %>%
  filter(N > 1) %>% as.data.frame();TrtStatSimp

TSFTrTstatSimp<-MetaSoil %>%
  filter(!is.na(simpson)) %>%
  group_by(Treatment, TSFdays) %>%
  summarise(N = n(),Avg_ASV = mean(simpson, na.rm = TRUE), 
      Variance = var(simpson, na.rm = TRUE), sd = sd(simpson), 
      se =sd / sqrt(N)) %>%
  filter(N > 1) %>% as.data.frame();TSFTrTstatSimp



#SimpsonEvenness
TrtStatSimpE<-MetaSoil %>%
  filter(!is.na(simpEven)) %>%
  group_by(Treatment) %>%
  summarise(N = n(),
            Avg_ASV = mean(simpEven, na.rm = TRUE), 
            Variance = var(simpEven, na.rm = TRUE), 
            sd = sd(simpEven), se =sd / sqrt(N)) %>%
  filter(N > 1) %>% as.data.frame();TrtStatSimpE

TSFTrTstatSimpE<-MetaSoil %>%
  filter(!is.na(simpEven)) %>%
  group_by(Treatment, TSFdays) %>%
  summarise(N = n(),Avg_ASV = mean(simpEven, na.rm = TRUE), 
            Variance = var(simpEven, na.rm = TRUE),
            sd = sd(simpEven), se =sd / sqrt(N)) %>%
  filter(N > 1) %>% as.data.frame();TSFTrTstatSimpE



#S.chao1
TrtStatChao1<-MetaSoil %>%
  filter(!is.na(S.chao1)) %>%
  group_by(Treatment) %>%
  summarise(N = n(),
            Avg_ASV = mean(S.chao1, na.rm = TRUE), 
            Variance = var(S.chao1, na.rm = TRUE), 
            sd = sd(S.chao1), se =sd / sqrt(N)) %>%
  filter(N > 1) %>% as.data.frame();TrtStatChao1

TSFTrTstatChao1<-MetaSoil %>%
  filter(!is.na(S.chao1)) %>%
  group_by(Treatment, TSFdays) %>%
  summarise(N = n(),Avg_ASV = mean(S.chao1, na.rm = TRUE), 
            Variance = var(S.chao1, na.rm = TRUE),
            sd = sd(S.chao1), se =sd / sqrt(N)) %>%
  filter(N > 1) %>% as.data.frame();TSFTrTstatChao1


#ACE
TrtStatsACE<-MetaSoil %>%
  filter(!is.na(S.ACE)) %>%
  group_by(Treatment) %>%
  summarise(N = n(),
            Avg_ASV = mean(S.ACE, na.rm = TRUE), 
            Variance = var(S.ACE, na.rm = TRUE), 
            sd = sd(S.ACE), se =sd / sqrt(N)) %>%
  filter(N > 1) %>% as.data.frame();TrtStatsACE

TSFTrTstatsACE<-MetaSoil %>%
  filter(!is.na(S.ACE)) %>%
  group_by(Treatment, TSFdays) %>%
  summarise(N = n(),Avg_ASV = mean(S.ACE, na.rm = TRUE), 
            Variance = var(S.ACE, na.rm = TRUE), sd = sd(S.ACE), 
            se =sd / sqrt(N)) %>%
  filter(N > 1) %>% as.data.frame();TSFTrTstatsACE


#Shan
TrtStatsShan<-MetaSoil %>%
  filter(!is.na(ShannonH)) %>%
  group_by(Treatment) %>%
  summarise(N = n(),Avg_ASV = mean(ShannonH, na.rm = TRUE), 
            Variance = var(ShannonH, na.rm = TRUE), 
            sd = sd(ShannonH), se =sd / sqrt(N)) %>%
  filter(N > 1) %>% as.data.frame();TrtStatsShan

TSFTrTstatsShan<-MetaSoil %>%
  filter(!is.na(ShannonH)) %>%
  group_by(Treatment, TSFdays) %>%
  summarise(N = n(),Avg_ASV = mean(ShannonH, na.rm = TRUE), 
            Variance = var(ShannonH, na.rm = TRUE), sd = sd(ShannonH), 
            se =sd / sqrt(N)) %>%
  filter(N > 1) %>% as.data.frame();TSFTrTstatsShan



#----
#--------
#**************************************************************************************----
#CALCULATE PERCENT CHANGE BETWEEN TREATMENTS AND TIMEPOINTS--------------------------------
#**************************************************************************************----
# * Percent change all samples .................................................
TrtSobs <-  MetaSoil%>% 
  group_by(Treatment) %>% 
  summarise(mean = mean(S.obs,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();TrtSobs


TrtSimp<-  MetaSoil%>% 
  group_by(Treatment) %>% 
  summarise(mean = mean(simpson,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();TrtSimp


TrtSimpEven<-  MetaSoil%>% 
  group_by(Treatment) %>% 
  summarise(mean = mean(simpEven,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();TrtSimpEven


TrtChao1<-  MetaSoil%>% 
  group_by(Treatment) %>% 
  summarise(mean = mean(S.chao1,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();TrtChao1

TrtACE<-  MetaSoil%>% 
  group_by(Treatment) %>% 
  summarise(mean = mean(S.ACE,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();TrtACE

TrtShan<-  MetaSoil%>% 
  group_by(Treatment) %>% 
  summarise(mean = mean(ShannonH,na.rm = T)) %>% 
  mutate(percent = (mean - first(mean))/first(mean)*100)%>%
  as.data.frame();TrtShan






#-----
#--------
#**************************************************************************************************----
#..................PLOTS SPECIES RICHNESS (ASV's)..............................................
#*************************************************************************************************----

# TREATMENT....................................................................................
Obs<-ggboxplot(MetaSoil, x="Treatment", y="S.obs",add = "jitter",color = "Treatment") +   
  scale_colour_manual(values=c("#45877f","#a2673f"))+
  stat_compare_means(label = "p.format", size=4.5, label.y = 590,label.x = 2.1)+
  annotate("text", x = 2.47, y=600, label = c("; -68%"), size=4.5)+
  theme_bw()+ ylab("Species Richness") + 
  theme(panel.grid = element_blank(),
        plot.margin = margin(0.3, 1.5, 0.3, 0.3, "mm"),
        axis.text.y = element_text(size=18, angle=90, vjust=.5,
                      hjust=.5, colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text.x = element_text(size=18, colour = "black"),
        legend.position = "bottom",text = element_text(size=20),
        axis.title.x=element_blank())+ ylim(0,600);Obs# Use only p.format as label. Remove method name.


# SEVERITY ...........................................................................................
SppSev<-ggboxplot(MetaSoil, x = "AshSeverity", y = "S.obs",color = "AshSeverity",add = "jitter")+
  scale_colour_manual(values=c("#027300","#93b396","#d9cd2c","#8B1919"))+
  stat_compare_means(label = "p.format", size=4.5, label.y = 594,label.x = 4)+
  #annotate("text", x = 4.47, y=600, label = c("; -68%"), size=4.5)+
  theme(plot.margin = margin(0.3, 1.5, 0.3, 0.3, "mm"),
        panel.border = element_rect(color="black", fill = NA),
        axis.text.y = element_text(size=16,angle=90, hjust=.5),
        text = element_text(size=16),
        legend.position = "")+ ylim(0,600)+
  ylab("Species Richness");SppSev

#SITES ...............................................................................................
SppSite<-ggboxplot(MetaSoil, x = "Plot", y = "S.obs",color = "Treatment",add = "jitter")+
  scale_colour_manual(values=c("#45877f","#a2673f"))+
  stat_compare_means(label = "p.format", size=4.5, label.y = 594,label.x = 8.5)+
  theme(plot.margin = margin(0.3, 1.5, 0.3, 0.3, "mm"),
        panel.border = element_rect(color="black", fill = NA),
        axis.text.y = element_text(size=16,angle=90, hjust=.5),
        text = element_text(size=16),
        legend.position = "")+ ylim(0,600)+
  ylab("Species Richness");SppSite


# TSF per TREATMENT ...................................................................................
mycolors<-c("#65000b","#a50F15","#EF3B2C","#FC9272","#FCBBA1","#fac2e3","#001499","#0a2bff","#9acfe4")
mycol2<-c("#40271f","#5e4839","#947f70","#c49f68","#bfb3aa","#ebddc7",
          "#235978","#45877f","#aed4ca")

ObsTSF<-ggplot(MetaSoil, aes(x=TSF2, y=S.obs, group=Treatment, col=Treatment))+
  stat_summary(fun=mean,geom="line", size=1, aes(linetype = Treatment),show.legend=TRUE)+
  stat_summary(fun.data = mean_se,geom = "errorbar",size=.85,alpha=0.7,position = position_dodge(0.01))+
  stat_compare_means(method = "kruskal.test", label.y = 410, label="p.signif", show.legend = FALSE, size=7)+
  theme_bw()+ theme(plot.margin = margin(0.3, 1.5, 0.3, 0.3, "mm"),
                    panel.grid = element_blank(),legend.key=element_blank(),
                    axis.text.y = element_text(size =18, angle=90,colour = "black", hjust=0.7), 
                    axis.text = element_text(size = 18, colour = "black"), 
                    axis.title = element_text(size = 18),
                    legend.title = element_text(size=18),
                    legend.text=element_text(size = 18),
                    legend.position = "bottom")+ 
  scale_linetype_manual(values=c("solid","twodash"))+
  scale_color_manual(values=c("#45877f","#a2673f"))+
  labs(x = "Time Since Fire (days)",y="Species Richness");ObsTSF


# TSF per PLOT ...................................................................................
ObsTSF2<-ggplot(MetaSoil, aes(x=TSF2, y=S.obs, group=Plot, col=Plot, shape=Treatment))+
  stat_summary(fun=mean,geom="line", size=1, aes(linetype = Treatment),show.legend=TRUE)+
  stat_summary(fun.data = mean_se,geom = "errorbar",size=.85, 
               alpha=0.7,position = position_dodge(0.01))+
  scale_linetype_manual(values=c("solid","twodash"))+
  scale_color_manual(values=mycol2) + 
  theme_bw()+ theme(plot.margin = margin(0.3, 1.5, 0.3, 0.3, "mm"),
                    panel.grid = element_blank(),
                    legend.key=element_blank(),
                    axis.text.y = element_text(size =18, angle=90,colour = "black", hjust=0.7), 
                    axis.text = element_text(size = 18, colour = "black"), 
                    axis.title = element_text(size = 18),
                    legend.title = element_text(size=18),
                    legend.text=element_text(size = 18),
                    legend.position = "right")+ 
  labs(x = "Time Since Fire (days)",y="Species Richness") + 
  stat_compare_means(method = "kruskal.test", label.y = 510,
                     label="p.signif", show.legend = FALSE, size=7);ObsTSF2



# TSF per SEVERITY ...................................................................................
SppTSF.Sev<-ggplot(MetaSoil, aes(x=TSF2, y=S.obs, group=AshSeverity, col=AshSeverity))+
  stat_summary(fun=mean,geom="line", size=1)+
  stat_summary(fun.data = mean_se,geom = "errorbar", size=.85, alpha=0.7,
               position = position_dodge(0.01))+
  scale_color_manual(values=c("#027300","#93b396","#d9cd2c","#8B1919")) + 
  theme_bw()+ theme(plot.margin = margin(0.3, 1.5, 0.3, 0.3, "mm"),
                    panel.grid = element_blank(),
                    legend.key=element_blank(),
                    axis.text.y = element_text(size =18, angle=90,colour = "black", hjust=0.7), 
                    axis.text = element_text(size = 18, colour = "black"), 
                    axis.title = element_text(size = 18),
                    legend.title = element_text(size=18),
                    legend.text=element_text(size = 18),
                    legend.position = "bottom")+ 
  labs(x = "Time Since Fire (days)",y="Species Richness");SppTSF.Sev



#-----
#-------
#****************************************************************************************************----
# .................ALPHA DIVERSITY ALL METRICS ..........................................................................
#****************************************************************************************************----
#margin(t, r, l, b)

SimpETrt<-ggboxplot(MetaSoil, x="Treatment", y="simpEven",add = "jitter",color = "Treatment") +   
  scale_colour_manual(values=c("#45877f","#a2673f"))+
  stat_compare_means(label = "p.format", size=4.5, label.y = 88,label.x = 2)+
  annotate("text", x = 2.47, y=89, label = c("; -78%"), size=4.5)+
  theme_bw()+ ylab("Simpson Evenness") + 
  theme(panel.grid = element_blank(),
        plot.margin = margin(0.3, 1.5, 0.3, 0.1, "mm"),
        axis.text.y = element_text(size=18, angle=90, vjust=.5,
                                   hjust=.5, colour = "black"),
        axis.text.x = element_text(size=18, colour = c("#45877f","#a2673f")),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",text = element_text(size=20),
        axis.title.x=element_blank()+ ylim(0,90));SimpETrt


SimpTrt<-ggboxplot(MetaSoil, x="Treatment", y="simpson",
                   add = "jitter",color = "Treatment") +
  scale_colour_manual(values=c("#45877f","#a2673f"))+
  stat_compare_means(label = "p.format", size=4.5, label.y = .99,label.x = 2)+
  annotate("text", x = 2.47, y=1, label = c("; -24%"), size=4.5)+
  theme_bw()+ ylab("Simpson") + 
  theme(panel.grid = element_blank(),
        plot.margin = margin(0.3, 1.5, 0.3, 0.3, "mm"),
        axis.text.y = element_text(size=18, angle=90, vjust=.5, hjust=.5, colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text.x = element_text(size=18, colour = c("#45877f","#a2673f")),
        legend.position = "bottom", text = element_text(size=20),
        axis.title.x=element_blank());SimpTrt


ChaoTrt<-ggboxplot(MetaSoil, x="Treatment", y="S.chao1",
                   add = "jitter",color = "Treatment") +
  scale_colour_manual(values=c("#45877f","#a2673f"))+
  stat_compare_means(label = "p.format", size=4.5, label.y = 594,label.x = 2)+
  annotate("text", x = 2.47, y=600, label = c("; -67%"), size=4.5)+
  theme_bw()+ ylab("Chao") + 
  theme(panel.grid = element_blank(),
        plot.margin = margin(0.3, 1.5, 0.3, 0.1, "mm"),
        axis.text.y = element_text(size=18, angle=90, vjust=.5,
                                   hjust=.5, colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text.x = element_text(size=18, colour = c("#45877f","#a2673f")),
        legend.position = "bottom", text = element_text(size=20),
        axis.title.x=element_blank());ChaoTrt


AceTrt<-ggboxplot(MetaSoil, x="Treatment", y="S.ACE",
                  add = "jitter",color = "Treatment") +
  scale_colour_manual(values=c("#45877f","#a2673f"))+
  stat_compare_means(label = "p.format", size=4.5, label.y = 594,label.x = 2)+
  annotate("text", x = 2.47, y=600, label = c("; -66%"), size=4.5)+
  theme_bw()+ ylab("Ace") + 
  theme(panel.grid = element_blank(),
        plot.margin = margin(0.3, 1.5, 0.3, 0.1, "mm"),
        axis.text.y = element_text(size=18, angle=90, vjust=.5,
                                   hjust=.5, colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text.x = element_text(size=18, colour = c("#45877f","#a2673f")),
        legend.position = "bottom", text = element_text(size=20),
        axis.title.x=element_blank());AceTrt


ShannonTrt<-ggboxplot(MetaSoil, x="Treatment", y="ShannonH",
                      add = "jitter",color = "Treatment") +
  scale_colour_manual(values=c("#45877f","#a2673f"))+
  stat_compare_means(label = "p.format", size=4.5, label.y = 5.95,label.x = 2)+
  annotate("text", x = 2.47, y=6, label = c("; -47%"), size=4.5)+
  theme_bw()+ ylab("Shannon") + 
  theme(panel.grid = element_blank(),
        plot.margin = margin(0.3, 1.5, 0.3, 0.3, "mm"),
        axis.text.y = element_text(size=18, angle=90, vjust=.5,
            hjust=.5, colour = "black"),
        axis.text.x = element_text(size=18, colour = c("#45877f","#a2673f")),
        legend.position = "bottom", 
        text = element_text(size=20));ShannonTrt


Alpha<-ggarrange(Obs,SimpTrt,ChaoTrt, AceTrt,ShannonTrt,SimpETrt, 
                 ncol=2, nrow=3,common.legend = TRUE, align="hv",
                 labels = c("A", "B","C","D","E","F"),legend="bottom",
                 hjust = c(-.5),font.label=list(size=18,face="bold"));Alpha
#-----
#-----
#****************************************************************************************************-----
#..................EXPORT PLOTS ...............................................................
#**************************************************************************************************-----
dir.create("1-Analysis/Diversity/Alpha/Fungi/Graphs")

pdf("1-Analysis/Diversity/Alpha/Fungi/Graphs/Obs-Trt.pdf", height=4.5, width=6)
Obs
dev.off()

pdf("1-Analysis/Diversity/Alpha/Fungi/Graphs/Obs-Severity.pdf", height=4.5, width=6)
SppSev
dev.off()

pdf("1-Analysis/Diversity/Alpha/Fungi/Graphs/Obs-Site.pdf", height=4.5, width=6)
SppSite
dev.off()

pdf("1-Analysis/Diversity/Alpha/Fungi/Graphs/Obs-TSF.pdf", height=5, width=6.5)
ObsTSF
dev.off()

pdf("1-Analysis/Diversity/Alpha/Fungi/Graphs/Obs-TSF-Per-Site.pdf", height=5, width=8)
ObsTSF2
dev.off()

pdf("1-Analysis/Diversity/Alpha/Fungi/Graphs/Obs-TSF-Per-Site.pdf", height=5, width=8)
ObsTSF2
dev.off()

pdf("1-Analysis/Diversity/Alpha/Fungi/Graphs/Alpha-Trt-All.pdf", height=13, width=12)
Alpha
dev.off()



#.----
#.----
#**************************************************************************************************----
# .................SIGNIFICANCE TESTING-MODEL SELECTION ....................................................................
#***************************************************************************************************----
attach(MetaSoil)

#* * * TREATMENT EFFECTS .......................................................
#Checking mean and variance of the data to see if its NegBinom..................
df <- MetaSoil %>%
  group_by(Treatment) %>%
  summarise(means = mean(S.obs),
            variance = var(S.obs));df #suggest negative binomial


#Start with a poisson model ....................................................
Pois1 <- glm(S.obs ~ Treatment, data = MetaSoil, family = poisson(link = log))
ods <- Pois1$deviance / Pois1$df.residual;ods#0verdispersion 29.69
par(mfrow=c(2,2))
plot(Pois1)#look at plots, cooks, levenes etc

#Zero inflaited model ..........................................................
sum(S.obs == 0) #0


#Test different identity link to see if disperssion decreases ..................
Pois2 <- glm(S.obs ~ Treatment,family = poisson(link = sqrt))
ods2<- Pois2$deviance / Pois2$df.residual;ods2 #29.69

Pois3 <- glm(S.obs ~ Treatment,family = poisson(link = identity))
ods3<- Pois3$deviance / Pois3$df.residual;ods3#29.69

Gam <- glm(S.obs ~ Treatment,family = "Gamma")
odsGam<-Gam$deviance/Gam$df.residual;odsGam#0.28

#Negative binomial... ..........................................................
#* mean and variance suggest this is the model to use
library(MASS)
nb1 <- glm.nb(S.obs ~Treatment,data = MetaSoil)
ods_nb <- (nb1$deviance / nb1$df.residual);ods_nb#1.05, mild= acceptable
plot(nb1)

#Compare models.................................................................
AIC(Pois1, Pois2, Pois3, Gam, nb1)#neg binomial better 

#Not needed for burned and unburnend data as data was transformed 
#to meet normality assumptions





#.----
#.-----
#***************************************************************************************************----
# .................TEST FOR SIGNIFICANCE  .....................................................................
#***************************************************************************************************----
#Load metadata to make it easier to call out variables..............
attach(MetaSoil)
summary(MetaSoil)

#Rescale variables..................................................
#Rescale and center the numerical data
MetaSoil$TSF<-scale(MetaSoil$TSFdays);summary(MetaSoil$TSF)
MetaSoil$Precip<-scale(MetaSoil$TotPrecip)



#Explore Temporal correlations to ensure that we are capturing as much of the variance............... 
ASVneg1 <- glmer.nb(S.obs ~ (1|Plot)+(1|Subplot)+(1|TSFdays), data = MetaSoil);summary(ASVneg1)
ASVneg2 <- glmer.nb(S.obs ~ (1|Plot)+(1|TSFdays), data = MetaSoil);summary(ASVneg2)
ASVneg3 <- glmer.nb(S.obs ~ (1|Subplot)+(1|TSFdays), data = MetaSoil);summary(ASVneg3)


#MuMIN for model selection..................................................
AICc(ASVneg1,ASVneg2,ASVneg3)#ASVneg1 better 3379.203



# * * FULL MODEL ..............................................................................
Full<-glmer.nb(S.obs~Treatment*TSF + Treatment*Precip + Treatment*InitialAshDepth + 
             TSF*Precip + TSF*InitialAshDepth + Treatment + TSF + Precip + InitialAshDepth +
             (1|Plot) + (1|Subplot) + (1|TSFdays),data = MetaSoil);summary(Full)

Null<-glmer.nb(S.obs ~ (1|Plot) + (1|Subplot) + (1|TSFdays), data=MetaSoil);summary(Null)
anova(Full,Null)#Full model better


#Start backwards model selection..............................................................

#Remove Treatment*TSF + 
Mod2<-glmer.nb(S.obs~Treatment*Precip + Treatment*InitialAshDepth + TSF*Precip + 
              TSF*InitialAshDepth + Treatment + TSF + Precip + InitialAshDepth +
              (1|Plot) + (1|Subplot) + (1|TSFdays),data = MetaSoil);summary(Mod2)
anova(Full, Mod2)#Mod2 better, not significant. remove


#Remove Treatment*Precip + 
Mod3<-glmer.nb(S.obs~Treatment*InitialAshDepth + TSF*Precip + TSF*InitialAshDepth +
            Treatment + TSF + Precip + InitialAshDepth + (1|Plot) + (1|Subplot) +
            (1|TSFdays),data = MetaSoil);summary(Mod3)
anova(Mod2,Mod3)#Mod 3 better, not signif, remove


#Remove Treatment*InitialAshDepth +
Mod4<-glmer.nb(S.obs ~ TSF*Precip + TSF*InitialAshDepth + Treatment + TSF + Precip + 
            InitialAshDepth + (1|Plot) + (1|Subplot) + (1|TSFdays),
            data = MetaSoil);summary(Mod4)
anova(Mod3,Mod4)#Same AIC, not significant, remove


#Remove TSF*Precip + 
Mod5<-glmer.nb(S.obs ~ TSF*InitialAshDepth + Treatment + TSF + Precip + 
                 InitialAshDepth + (1|Plot) + (1|Subplot) + (1|TSFdays),
               data = MetaSoil);summary(Mod5)
anova(Mod4,Mod5)#Mod 4 better, significant, do not remove......................................


#Remove TSF*InitialAshDepth +
Mod6<-glmer.nb(S.obs ~ Treatment + TSF + Precip +  InitialAshDepth + TSF*Precip + 
            (1|Plot) + (1|Subplot) + (1|TSFdays),data = MetaSoil);summary(Mod6)
anova(Mod4,Mod6)#Mod 4 better, signif, do not remove...........................................


#Remove Treatment + 
Mod7<-glmer.nb(S.obs ~ TSF + Precip +  InitialAshDepth + TSF*Precip + TSF*InitialAshDepth +
               (1|Plot) + (1|Subplot) + (1|TSFdays),data = MetaSoil);summary(Mod7)
anova(Mod4,Mod7)#Mod 4 better, signif, do not remove...........................................

#Remove TSF + 
Mod8<-glmer.nb(S.obs ~ Precip + InitialAshDepth + TSF*Precip + TSF*InitialAshDepth +
        Treatment + (1|Plot) + (1|Subplot) + (1|TSFdays),data = MetaSoil);summary(Mod8)
anova(Mod4,Mod8)#Same AIC, not signif, remove


#Precip + 
Mod9<-glmer.nb(S.obs ~ InitialAshDepth + TSF*Precip + TSF*InitialAshDepth +Treatment + 
         (1|Plot) + (1|Subplot) + (1|TSFdays),data = MetaSoil);summary(Mod9)
anova(Mod8,Mod9)#Same AIC, not signif, remove


#Remove InitialAshDepth + 
Mod10<-glmer.nb(S.obs ~ TSF*Precip + TSF*InitialAshDepth +Treatment + (1|Plot) + 
               (1|Subplot) + (1|TSFdays),data = MetaSoil);summary(Mod10)
anova(Mod9,Mod10)#Same AIC, not signif, remove



#Compare all models....................................................................
anova(Full, Mod2, Mod3, Mod4, Mod5, Mod6, Mod7, Mod8, Mod9,Mod10)
AICc(Full, Mod2, Mod3, Mod4, Mod5, Mod6, Mod7, Mod8, Mod9,Mod10)

#Mod 4 suggested to be the best model, lower AIC values
summary(Mod4)


#Final Model...........................................................................
FullNB<-glmer.nb(S.obs ~ TSF*Precip + TSF*InitialAshDepth +Treatment + (1|Plot) + 
              (1|Subplot) + (1|TSFdays),data = MetaSoil);summary(FullNB)
R2Full<-r.squaredGLMM(FullNB);R2Full#Delta;0.6337314,0.7391477 




#Model Diagnostic & Validation plots ...........................................
par(mfrow=c(1,1));plot(FullNB)
Fitted <- fitted(FullNB) #command fitted gives already e^model
Res <- resid(FullNB, type = "pearson")

plot(x = Fitted,  y = Res,xlab = "Fitted values",ylab = "Pearson residuals")
abline(h = 0, lty = 2, col="red")


#Export model results.............................................................................
capture.output(summary(FullNB),file="1-Analysis/Diversity/Alpha/Fungi/Tables/Significance/NegBinom-FullModel.csv")
write.csv(R2Full,"1-Analysis/Diversity/Alpha/Fungi/Tables/Significance/NegBinom-R2FullModel.csv")




#.----
#.----
#*************************************************************************************************************----
# .................SIGNIFICANCE PER VARIABLE, INDEPENDENTLY ..................................................
#*************************************************************************************************************----

#Test treatment * tsf significance..................................................................................
TrTSFNB<-glmer.nb(S.obs ~ Treatment * TSF2 + (1|Plot) + (1|Subplot) + (1|TSFdays),data = MetaSoil);summary(TrTSFNB)
R2TrTSF<-r.squaredGLMM(TrTSFNB);R2TrTSF#Delta:0.6466322 0.7728300
plot(TrTSFNB)

#Get pairwise comparisons w differences indicated as letters.......
UTrtTSF<- emmeans(TrTSFNB, ~ Treatment*TSF2)
TrtTSFPost<-pairs(UTrtTSF);TrtTSFPost


#Export model results..........................................................................................
capture.output(summary(TrTSFNB),file="1-Analysis/Diversity/Alpha/Fungi/Tables/Significance/Trt-TSF-NegBinom.csv")
write.csv(R2TrTSF,"1-Analysis/Diversity/Alpha/Fungi/Tables/Significance/Trt-TSF-R2-NegBinom.csv")
write.csv(TrtTSFPost,"1-Analysis/Diversity/Alpha/Fungi/Tables/Significance/PostHoc-Trt-TSF-NegBinom-.csv")



#Test treatment * Precip significance..........................................................................................
PrecipNB<-glmer.nb(S.obs ~ Treatment * Precip + TSF * Precip +(1|Plot)+(1|Subplot)+(1|TSFdays),data = MetaSoil);summary(PrecipNB)
R2Precip<-r.squaredGLMM(PrecipNB);R2Precip#Delta:0.2856056 0.5611148
plot(TrPrecipNB)

#Get pairwise comparisons w differences indicated as letters.......
UPrecip<- emmeans(PrecipNB, ~ Treatment*Precip);UPrecip
pairs(UPrecip)

#Export model results..........................................................................................
capture.output(summary(PrecipNB),file="1-Analysis/Diversity/Alpha/Fungi/Tables/Significance/Precip-Trt-TSF-NegBinom.csv")
write.csv(R2Precip,"1-Analysis/Diversity/Alpha/Fungi/Tables/Significance/Precip-R2-NegBinom.csv")



#Test treatment * AshDepth significance...................................................................................................
AshDepthNB<-glmer.nb(S.obs ~ Treatment * AshDepth + TSF * AshDepth +(1|Plot)+(1|Subplot)+(1|TSFdays),data = MetaSoil);summary(AshDepthNB)
R2AshDepth<-r.squaredGLMM(AshDepthNB);R2AshDepth#Delta:0.2887868 0.5510319
plot(TrAshDepthNB)

#Export model results..........................................................................................
capture.output(summary(AshDepthNB),file="1-Analysis/Diversity/Alpha/Fungi/Tables/Significance/AshDepth-Trt-TSF-NegBinom.csv")
write.csv(R2AshDepth,"1-Analysis/Diversity/Alpha/Fungi/Tables/Significance/AshDepth-R2-NegBinom.csv")




