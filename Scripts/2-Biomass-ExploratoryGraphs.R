#Reset R's Brain
rm(list=ls())

#set working directory
setwd("~/Dropbox/1-Dissertation/1-Research/HolyFire/1-HolyFire-Pub1/1c-Qpcr/Year1")

#load libraries
library("tidyverse")
library("scales")
library("gridExtra")
library("ggpubr")
library("ggplot2")

#Export tables for graphing purposes in future, clean graphs.............................
Metadata<-read.csv("Analysis/Metadata/MetadataFilter.csv",row.names = 1)
MetaFun<-read.csv("Analysis/Metadata/MetaFun-Biomass.csv", row.names = 1)
MetaBac<-read.csv("Analysis/Metadata/MetaBac-Biomass.csv", row.names = 1)
BacBurn<-read.csv("Analysis/Metadata/MetaBac-Biomass-Burn.csv", row.names = 1)
BacUn<-read.csv("Analysis/Metadata/MetaBac-Biomass-Unburn.csv", row.names = 1)
FunBurn<-read.csv("Analysis/Metadata/MetaFun-Biomass-Burn.csv", row.names = 1)
FunUn<-read.csv("Analysis/Metadata/MetaFun-Biomass-Unburn.csv", row.names = 1)

#----
#----
#***********************************************************************************************************----
#-------------------------- Quality control ---------------------------------------------------------------------
#***********************************************************************************************************----
Metadata$TSFdays <- as.factor(Metadata$TSFdays)
Metadata$Treatment<-as.factor(Metadata$Treatment)
Metadata$Treatment <- try(relevel(Metadata$Treatment , "Unburned"))
Metadata$AshSeverity<-ordered(Metadata$AshSeverity, c("Unburned","Low","Moderate", "High"))



#----
#-----
#*******************************************************************************************************************************----
#................... Create graphs................................................................................................
#*******************************************************************************************************************************----
mycol<-c("#40271f","#5e4839","#947f70","#c49f68","#bfb3aa","#ebddc7","#235978","#45877f","#aed4ca")

#TSF and treatment....................................................................................................
TSFtrt<-ggplot(Metadata, aes(x=TSFdays, y=CopyNumberPerGram, group=Treatment, col=Treatment)) +
  stat_summary(fun = "mean", size = 1.6, geom = "point")+
  stat_summary(fun=mean,geom="line", size=.75, aes(linetype = Treatment),show.legend=TRUE)+
  stat_summary(fun.data = mean_se,geom = "errorbar",size=.85, alpha=0.7,position = position_dodge(0.01))+
  scale_linetype_manual(values=c("solid","twodash"))+
  scale_color_manual(values=c("#45877f","#a2673f"))+ 
  theme_bw()+ 
  theme(panel.grid = element_blank(),
        legend.key=element_blank(),
        axis.text.y = element_text(size =14,angle=90,hjust=0.5,colour = "black"), 
        axis.text.x = element_text(size = 14,colour = "black"), 
        axis.title.y = element_text(size = 16, colour = "black"), 
        axis.title.x = element_text(size = 16, colour = "black"),
        legend.title = element_text(size=15),
        legend.text=element_text(size = 14),legend.position = "bottom")+ 
  labs(x = "Time Since Fire (days)",y="Gene Copy Number (g-1 soil) ")+
  scale_y_continuous(labels = label_scientific(digits = 1))+
  coord_cartesian(ylim = c(0, NA));TSFtrt

#Add  kruskal.wallis signfiicance at each time point
TSFtrtStat<- TSFtrt + stat_compare_means(method = "kruskal.test", 
               label="p.signif", label.y = -1000000,
                show.legend = FALSE );TSFtrtStat 

(TSFtrtStat2<- TSFtrtStat + facet_wrap( ~ Kingdom, scale = "free")+
    theme(strip.text = element_text(size=20)));TSFtrtStat2

pdf("Analysis/Graphs/Trt/TSF-Trt-BacFun-Facet.pdf", height=7, width=12)
TSFtrtStat2
dev.off()

#Look at variance per plot, just to look...............................................................
Sub<-ggboxplot(Metadata, x="Subplot", y="CopyNumberPerGram") +   
  geom_boxplot(aes(fill=Treatment))+ 
  theme_bw() +  ylab("Species Richness") + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    text = element_text(size=16),
    axis.title.x=element_blank(),
    axis.text.y = element_text(size=12,angle=90, hjust=.5,colour = "black"),
    axis.text.x = element_text(size=12,angle=90, hjust=.5, colour = "black"))+ 
  scale_fill_manual(values=c("#45877f","#a2673f"));Sub

(Sub2<- Sub  + facet_wrap( ~ Timepoint, scale = "free")+
    theme(strip.text = element_text(size=20)));Sub2

pdf("Analysis/Graphs/Trt/Trt-Subplot-CopyNum.pdf", height=12, width=16)
Sub2
dev.off()



#Soil burn severity......................................................................................................
TSFSevAll<-ggplot(Metadata, aes(x=TSFdays, y=CopyNumberPerGram,group=AshSeverity, col=AshSeverity, label="Sample")) +
  stat_summary(fun = "mean", size = 1.6, geom = "point")+
  stat_summary(fun=mean,geom="line", size=.75, aes(linetype = Treatment),show.legend=TRUE)+
  stat_summary(fun.data = mean_se,geom = "errorbar",size=.85, alpha=0.7,position = position_dodge(0.01))+
  scale_linetype_manual(values=c("solid","twodash"))+
  scale_color_manual(values=c("#027300","#93b396","#d9cd2c","#8B1919"))+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),  
        legend.key=element_blank(),
        axis.text.y = element_text(size =14, angle=90, hjust=0.5, colour = "black"), 
        axis.text.x = element_text(size = 14, colour = "black"), 
        axis.title.y = element_text(size = 18, colour = "black"), 
        axis.title.x = element_text(size = 18, colour = "black"),
        legend.title = element_text(size=18),
        legend.text=element_text(size = 14),
        legend.position = "bottom")+ 
  labs(x = "Time Since Fire (days)",y="Gene Copy Number") +
  scale_y_continuous(labels = label_scientific(digits = 1))+
  coord_cartesian(ylim = c(0, NA));TSFSevAll

(TSFSevAll2<- TSFSevAll  + facet_wrap( ~ Kingdom, scale = "free")+
    theme(strip.text = element_text(size=18)));TSFSevAll2

pdf("Analysis/Graphs/Trt/TSF-Severity-CopyNum.pdf", height=5, width=8)
TSFSevAll2
dev.off()

#Severity by boxplot .........................................................................
SevBox<-ggboxplot(Metadata, x="AshSeverity", y="CopyNumberPerGram") +   
  geom_boxplot(aes(fill=AshSeverity))+ 
  theme_bw() +  
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
    panel.background = element_blank(),panel.grid = element_blank(),
    legend.position = "bottom",
    text = element_text(size=14),
    axis.title.x=element_blank(),
    axis.text.y = element_text(size=18, angle=90, hjust=.5, colour = "black"),
    axis.text.x = element_text(size=18, colour = "black"))+ 
  scale_fill_manual(values=c("#027300","#93b396","#d9cd2c","#8B1919"))+
  ylab("Gene Copy Number")+
  scale_y_continuous(labels = label_scientific(digits = 1))+
  coord_cartesian(ylim = c(0, NA));SevBox
                  
(SevBox2<-SevBox + facet_wrap( ~ Kingdom, scale = "free")+
    theme(strip.text = element_text(size=18)));SevBox2


pdf("Analysis/Graphs/Trt/Severity-CopyNum-BoxPlot.pdf", height=5, width=9)
SevBox2
dev.off()


TSFPt <- ggplot(Metadata, aes(x=TSFdays, y=CopyNumberPerGram, group=Treatment, col=Treatment, label="Sample")) +
  stat_summary(fun=mean,geom="line", size=1, aes(linetype = Treatment),show.legend=TRUE)+
  stat_summary(fun.data = mean_se,geom = "errorbar",size=.85, 
               alpha=0.7,position = position_dodge(0.01))+
  geom_point(size=1, alpha=0.5)+ #Adds in the data points
  scale_linetype_manual(values=c("solid","twodash"))+
  scale_color_manual(values=c("#45877f","#a2673f"))+
  theme_bw()+ 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1.2),  
    panel.grid= element_blank(),
        legend.key=element_blank(),
        axis.text.y = element_text(size =14, angle=90, hjust=0.5, colour = "black"), 
        axis.text.x = element_text(size = 14, colour = "black"), 
        axis.title.y = element_text(size = 18, colour = "black"), 
        axis.title.x = element_text(size = 18, colour = "black"),
        legend.title = element_text(size=18),
        legend.text=element_text(size = 14),
        legend.position = "bottom")+ 
  labs(x = "Time Since Fire (days)", y="Gene Copy Number (g-1 soil)");TSFPt 

#add in the significance
TSFPt2 <- TSFPt+ 
  stat_compare_means(method = "kruskal.test", 
            label="p.signif");TSFPt2

# This shows the variation between the burned and unburnded plots, free scales
TSFPt3 <- TSFPt2 +facet_wrap( ~ Kingdom, scale = "free")+
  theme(strip.text = element_text(size=20));TSFPt3

pdf("Analysis/Graphs/Trt/Trt-TSF-CopyNum-Points.pdf", height=6, width=11)
TSFPt3
dev.off()

#----  
#----
#***************************************************************************************************************----
#MetaFun=-------------------------------------------------------------------------------------
#***************************************************************************************************************----
#Fungi per timesfire per plot to look at variance
TSFtrtFun<-ggplot(MetaFun, aes(x=TSFdays, y=CopyNumberPerGram, group=Plot, col=Plot, label="Sample")) +
  stat_summary(fun = "mean", size = 1.6, geom = "point")+
  stat_summary(fun=mean,geom="line", size=1, aes(linetype = Treatment),show.legend=TRUE)+
  stat_summary(fun.data = mean_se,geom = "errorbar",size=.85,alpha=0.7,position = position_dodge(0.01))+
  scale_linetype_manual(values=c("solid","twodash"))+
  scale_color_manual(values=c("#40271f","#5e4839","#947f70",
     "#c49f68","#bfb3aa","#ebddc7","#235978","#45877f","#aed4ca"))+
  labs(subtitle = "Fungi")+
  theme_bw()+ 
  theme(panel.grid= element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),  
        legend.key=element_blank(),
        axis.text.y = element_text(size =14, angle=90, hjust=0.5, colour = "black"), 
        axis.text.x = element_text(size = 14, colour = "black"), 
        axis.title.y = element_text(size = 18, colour = "black"), 
        axis.title.x = element_text(size = 18, colour = "black"),
        legend.title = element_text(size=18),
        legend.text=element_text(size = 14),
        legend.position = "right")+ labs(x = "",y="")+
  scale_y_continuous(labels = label_scientific(digits = 1))+
  coord_cartesian(ylim = c(0, NA));TSFtrtFun

#MetaBac
TSFtrtBac<-ggplot(MetaBac, aes(x=TSFdays, y=CopyNumberPerGram, group=Plot, col=Plot, label="Sample")) +
  stat_summary(fun = "mean", size = 1.6, geom = "point")+
  stat_summary(fun=mean,geom="line", size=1, aes(linetype = Treatment),show.legend=TRUE)+
  stat_summary(fun.data = mean_se,geom = "errorbar",size=.85,alpha=0.7,position = position_dodge(0.01))+
  scale_linetype_manual(values=c("solid","twodash"))+
  scale_color_manual(values=c("#40271f","#5e4839","#947f70",
     "#c49f68","#bfb3aa","#ebddc7","#235978","#45877f","#aed4ca"))+
  labs(subtitle = "Bacteria")+
  theme_bw()+ 
  theme(panel.grid= element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),  
        legend.key=element_blank(),
        axis.text.y = element_text(size =14, angle=90, hjust=0.5, colour = "black"), 
        axis.text.x = element_text(size = 14, colour = "black"), 
        axis.title.y = element_text(size = 18, colour = "black"), 
        axis.title.x = element_text(size = 18, colour = "black"),
        legend.title = element_text(size=18),
        legend.text=element_text(size = 14),
        legend.position = "right")+ labs(x = "",y="Gene Copy Number (g-1 soil)")+
  scale_y_continuous(labels = label_scientific(digits = 1))+
  coord_cartesian(ylim = c(0, NA));TSFtrtBac


TSFtrt<-ggarrange(TSFtrtBac,TSFtrtFun, ncol=2, nrow=1, common.legend = TRUE, 
        legend="right", hjust = c(-.8,-1),
         align="hv",font.label = list(size=16,face="bold"));TSFtrt

pdf("Analysis/Graphs/Trt/Biomass-BacFun-Panels-TSFPlot-trt.pdf", height=5, width=10)
TSFtrt
dev.off()
