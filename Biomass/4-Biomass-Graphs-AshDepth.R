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

#set as factor first to set base level..............................................
MetaBioFun$Treatment<-as.factor(MetaBioFun$Treatment)
MetaBioBac$Treatment<-as.factor(MetaBioBac$Treatment)
MetaBioFun$Treatment <- try(relevel(MetaBioFun$Treatment , "Unburned"))
MetaBioBac$Treatment <- try(relevel(MetaBioBac$Treatment , "Unburned"))


#----
#----
#*****************************************************************************************************-----
# ----------------------------- CREATE PLOTS --------------------------------------------------------------
#*****************************************************************************************************-----
# * * CREATE Color scheme ...............................................
#Scheme for timepoitns ..................................................
colorTP<-c("#de0025","#b05800","#ff6a00","#ffd000","#fac2e3","black","#0a2bff","#00b5f7","#007d2a")

#Maintain same color for all graphs .....................................
colorPlot<-c("black","#5e4839","#947f70","#c49f68","#bfb3aa","#ebddc7","#235978","#45877f","#aed4ca")


#Bacteria ............................................................................
BioBac<-ggplot(MetaBioBac, aes(x=InitialAshDepth, y=CopyNumberPerGram))  +
  geom_smooth(method = MASS::glm.nb, se = TRUE)+ 
  scale_shape_manual(values=c(15,24,19,1,18,25,8,17,13))+
  geom_point(size=2,aes(shape = as.factor(TSFdays), colour = Treatment))+
  stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label..,
           sep = "~~~~")), formula = y~x, label.y =6e08,label.x = 3,  size=3)+
  scale_color_manual(values=  c("#45877f","#a2673f")) +
  scale_fill_manual(values=colorTP) +
  labs(title = "A Bacteria", size=14)+
  theme_bw()+ 
  theme( panel.grid= element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),  
        legend.key=element_blank(),
        axis.title.y = element_text(size = 18, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size =16, angle=90,hjust=0.5,colour = "black"), 
        axis.text.x = element_text(size = 16,colour = "black"), 
        legend.title = element_text(size=18),
        legend.text=element_text(size = 16),
        legend.position = "bottom")+ 
  labs(colour="Treatment", shape="Time Since Fire (days)")+
  labs(x = "Ash Depth (cm)",y="Gene Copy Number (g-1 soil)")+
  coord_cartesian(ylim = c(0, NA))+
  scale_y_continuous(labels = label_scientific(digits = 1));BioBac

#Fungi............................................................................
BioFun<-ggplot(MetaBioFun, aes(x=InitialAshDepth, y=CopyNumberPerGram))  +
  geom_smooth(method = MASS::glm.nb, se = TRUE)+ 
  scale_shape_manual(values=c(15,24,19,1,18,25,8,17,13))+
  geom_point(size=2,aes(shape = as.factor(TSFdays), colour = Treatment))+
  stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label..,
       sep = "~~~~")), formula = y~x, label.x = 3, size=3)+
  scale_color_manual(values=  c("#45877f","#a2673f"))+
  scale_fill_manual(values=colorTP) +
  labs(title = "B Fungi", size=14)+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),  
        legend.key=element_blank(),
        axis.text.y = element_text(size =16, angle=90,hjust=0.5,colour = "black"), 
        axis.text.x = element_text(size = 16,colour = "black"), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_text(size=18),
        legend.text=element_text(size = 16),
        legend.position = "bottom")+ 
  labs(colour="Treatment", shape="Time Since Fire (days)")+
  #labs(x = "Ash Depth (cm)",y="")+
  scale_y_continuous(labels = label_scientific(digits = 1));BioFun


#Create panel.........................................................................................
AshPlots<-ggarrange(BioBac,BioFun, ncol=1, nrow=2, 
              common.legend = TRUE, legend="bottom", align="hv", 
              font.label = list(size=16,face="bold"));AshPlots

pdf("Analysis/Graphs/Trt/AhDepth-BacFun.pdf", height=9, width=6)
AshPlots
dev.off()

