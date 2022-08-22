#8-12-2022
rm(list=ls())#reset working directory

#Set working directory.................................................................
setwd("~/Dropbox/1-Dissertation/1-Research/HolyFire/1-HolyFire-Pub1/1c-Qpcr/Year1")

#Load packages 
library(ggpubr) 
library(ggfortify) #Manipulate and graph plots
library(ggplot2)


# Load data ............................................................................................
MetaBioBac<-read.csv("Analysis/Metadata/MetaBac-Biomass.csv", header = TRUE)
MetaBioFun<-read.csv("Analysis/Metadata/MetaFun-Biomass.csv", header = TRUE)
RainData<-read.csv("Analysis/Metadata/MetaSoilSpp.csv", header=TRUE)

#.----
#.----
#**********************************************************************************************************----
# .................. Quality Control ...................................................
#**********************************************************************************************************----
MetaBioBac$Plot<- as.factor(MetaBioBac$Plot)
MetaBioBac$TSFdays<- as.factor(MetaBioBac$TSFdays)
MetaBioBac$Treatment<-as.factor(MetaBioBac$Treatment)
MetaBioBac$InitialAshDepth<-as.factor(MetaBioBac$InitialAshDepth)
MetaBioBac$Treatment <- try(relevel(MetaBioBac$Treatment , "Unburned"))#relevel so that unburned is the base level for comparison

MetaBioFun$Plot<- as.factor(MetaBioFun$Plot)
MetaBioFun$TSFdays<- as.factor(MetaBioFun$TSFdays)
MetaBioFun$Treatment<-as.factor(MetaBioFun$Treatment)
MetaBioFun$InitialAshDepth<-as.factor(MetaBioFun$InitialAshDepth)
MetaBioFun$Treatment <- try(relevel(MetaBioFun$Treatment , "Unburned"))#relevel so that unburned is the base level for comparison


#----
#----
#*************************************************************************************************************************************----
#-------------BIOMASS GRAPHS --------------------------------------------------------------------------------
#*************************************************************************************************************************************----
library(scales)#changes legend from sci to comma
library(patchwork)#to put graphs together

#BACTERIA..........................................................................................
Bac<-ggplot(MetaBioBac, aes(x=TSFdays,y=CopyNumberPerGram, group=Treatment, col=Treatment)) +
  stat_summary(fun = "mean", size = 1.6, geom = "point")+
  stat_summary(fun=mean,geom="line", size=.75, aes(linetype = Treatment),show.legend=TRUE)+
  stat_summary(fun.data = mean_se,geom = "errorbar",size=.85, 
               alpha=0.7,position = position_dodge(0.01))+
  scale_linetype_manual(values=c("solid","twodash"))+
  scale_color_manual(values=c("#45877f","#a2673f"))+
  labs(title = "Biomass", subtitle = "Bacteria")+
  theme_bw()+  
  theme(panel.grid  = element_blank(),
        axis.title.y = element_text(size = 17, colour = "black"),
        axis.text.y = element_text(size = 16, angle=90, hjust=.5,colour = "black"), 
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_text(size=18),
        legend.text=element_text(size = 18),
        legend.position = "none")+
  ylab("Gene Copy Number (g-1 soil)")+
  coord_cartesian(ylim = c(0, NA))+
  scale_y_continuous(labels = label_scientific(digits = 1));Bac

Bac2 <-Bac + stat_compare_means(method = "kruskal.test", 
       label="p.signif", label.y=4e+08,show.legend = FALSE );Bac2


#FUNGI ........................................................................................
Fun1<-ggplot(MetaBioFun, aes(x=TSFdays,y=CopyNumberPerGram, group=Treatment, col=Treatment)) +
  stat_summary(fun = "mean", size = 1.6, geom = "point")+
  stat_summary(fun=mean,geom="line", size=.75, aes(linetype = Treatment),show.legend=TRUE)+
  stat_summary(fun.data = mean_se,geom = "errorbar",size=.85, 
               alpha=0.7,position = position_dodge(0.01))+
  scale_linetype_manual(values=c("solid","twodash"))+
  scale_color_manual(values=c("#45877f","#a2673f"))+
  labs(subtitle = "Fungi")+
  theme_bw()+  
  theme(panel.grid  = element_blank(),
        axis.text.y = element_text(size = 16, angle=90, hjust=.5,colour = "black"), 
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),
        axis.title = element_blank(),
        legend.title = element_text(size=18),
        legend.text=element_text(size = 18),
        legend.position = "none")+
  coord_cartesian(ylim = c(0, NA))+
  scale_y_continuous(labels = label_scientific(digits = 1));Fun1

Fun2 <-Fun1 + stat_compare_means(method = "kruskal.test", 
        label="p.signif",label.y=11e07,show.legend = FALSE );Fun2


# RAIN GRAPH...........................................................................
Rain<-ggplot(RainData, aes(x=as.factor(TSFdays), y=TotPrecip, group=Treatment, col=Treatment))+
  stat_summary(fun = "mean", size = 1.6, geom = "point", col="#94C5E0")+
  geom_line(aes(y=TotPrecip), col="#94C5E0",  size=1, show.legend = T) +
  theme_bw()+ 
  theme(panel.grid  = element_blank(),
        axis.title.x =  element_text(size =17,  colour = "black"), 
        axis.text.y = element_text(size = 12,angle=90, hjust=1,color="black"),  
        axis.text.x = element_text(size = 16, colour = "black"), 
        axis.title = element_text(size = 18),
        legend.title = element_text(size=18),
        legend.text=element_text(size = 18),
        legend.position = "bottom")+
  ylab("")+  xlab ("Time Since Fire (days)");Rain

Bac<-Bac2/Rain + plot_layout(nrow = 2, heights = c(4,1),
       guides = 'collect') & theme(legend.position = 'none');Bac


Fun<-Fun2/Rain + plot_layout(nrow = 2, heights = c(4,1),
        guides = 'collect') & theme(legend.position = 'none');Fun


Plots<-Bac2 + Fun + plot_layout(ncol = 1, nrow=3,guides = 'collect',heights = c(1.6,2)) &
  theme(legend.position = 'bottom');Plots

pdf("Analysis/Graphs/Trt/Rain-BiomassBacFun.pdf", height = 12, width = 8)
Plots
dev.off()



