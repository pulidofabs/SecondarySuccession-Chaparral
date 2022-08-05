#10-7-2021

#Reset R's Brain
rm(list=ls())

#Set working directory.............................................................
setwd("C:/Users/fabipc/Dropbox/3-HolyFire/1-HolyFire-Pub1/1-Bacteria/")

#Load librarires...................................................................
#Ape (import tree-qiime), codyn (turnover), tidy (load taxtable)
library(phyloseq,qiime2R,ape,tidyverse,codyn)
library(phyloseq)#RelAbund
library(qiime2R)#use qiime files
library(ape)#
library(tidyverse)
library(codyn)

#Load data .........................................................................................
Metadata<-read_tsv("1-Analysis/Metadata/MetaRareSoilPhylo.tsv")#Rare metadata
Table<-read_qza("Qiime/core-metrics-results/rarefied_table.qza", tmp = "C:/tmp")#Rarefied table
Tree<-read_qza("Qiime/Phylo/Bacteria-rooted-tree-NC.qza",tmp = "C:/tmp" )
Taxonomy<-read_qza("Qiime/Taxonomy/Bacteria-taxonomy-paired-Clean-2020-132.qza",tmp = "C:/tmp")
Taxtable<-Taxonomy$data %>% as_tibble() %>% separate(Taxon, sep = ";",
              c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))


#Create Phyloseq object ............................................................................
physeq<-phyloseq(
  otu_table(Table$data, taxa_are_rows = TRUE), 
  phy_tree(Tree$data), 
  tax_table(as.data.frame(Taxtable) %>% column_to_rownames("Feature.ID") %>% as.matrix()), 
  sample_data(Metadata%>% as.data.frame() %>% column_to_rownames("SampleID")))


#********************************************************************************************************************----
#------------------------------------     QUALITY CONTROL        ---------------------------------------------
#********************************************************************************************************************----
rank_names(physeq)# Look at rank names


#Quality control: Remove the g__ from each rank number..............................................................
colnames(tax_table(physeq))= c("Kingdom","Phylum","Class","Order","Family","Genus","Species", "Confidence")
tax_table(physeq)[, "Kingdom"] <- gsub("D_0__", "", tax_table(physeq)[, "Kingdom"])
tax_table(physeq)[, "Phylum"] <- gsub("D_1__", "", tax_table(physeq)[, "Phylum"])
tax_table(physeq)[, "Class"] <- gsub("D_2__", "", tax_table(physeq)[, "Class"])
tax_table(physeq)[, "Order"] <- gsub("D_3__", "", tax_table(physeq)[, "Order"])
tax_table(physeq)[, "Family"] <- gsub("D_4__", "", tax_table(physeq)[, "Family"])
tax_table(physeq)[, "Genus"] <- gsub("D_5__", "", tax_table(physeq)[, "Genus"])
tax_table(physeq)[, "Species"] <- gsub("D_6__", "", tax_table(physeq)[, "Species"])


#First convert TSF to factor .....................................................
#sample_data(physeq)$TSFdays<-factor(sample_data(physeq)$TSFdays)

#SUBSET DATA
Burned<-subset_samples(physeq,Treatment=="Burned");Burned#29819x229
sample_names(Burned);sample_data(Burned)$TSFdays

Unburned<-subset_samples(physeq,Treatment=="Unburned");Unburned#29819x110
sample_names(Unburned);sample_data(Unburned)$TSFdays



#---
#----
#*************************************************************************************************----
#.................... CALCULATE ABUNDANCE ............................................................
#*************************************************************************************************----
#Will not agglomorate since I want all the data so that I can run turnover rates based
#on whatever variable of interest

#Total relative abundance ............................................................
Genus1 <- tax_glom(physeq, taxrank = 'Genus') # agglomerate taxa
RelGen <- transform_sample_counts(Genus1, function(x) x/sum(x)) #get abundance in %
RelGen1 <- psmelt(RelGen) # create dataframe from phyloseq object
RelGen1$Genus <- as.character(RelGen1$Genus) #convert to character

#3 percent calculation.............................................
Genus1a <- tax_glom(physeq, taxrank = 'Genus')
RelGena <- transform_sample_counts(Genus1a, function(x) x/sum(x))
RelGen1a <- psmelt(RelGena)
RelGen1a$Genus <- as.character(RelGen1a$Genus)
RelGen1a$Genus[RelGen1a$Abundance < 0.03] <- "< 3% abund."

#Burned...........................................................
GenB1 <- tax_glom(Burned, taxrank = 'Genus')
RelGenB <- transform_sample_counts(GenB1, function(x) x/sum(x))
RelGenB1 <- psmelt(RelGenB)
RelGenB1$Genus <- as.character(RelGenB1$Genus) 

#Unburned.........................................................
GenU1 <- tax_glom(Unburned, taxrank = 'Genus')
RelGenU <- transform_sample_counts(GenU1, function(x) x/sum(x))
RelGenU1 <- psmelt(RelGenU)
RelGenU1$Genus <- as.character(RelGenU1$Genus) 

#Burned...........................................................
GenB1a <- tax_glom(Burned, taxrank = 'Genus')
RelGenBa <- transform_sample_counts(GenB1a, function(x) x/sum(x))
RelGenB1a <- psmelt(RelGenBa)
RelGenB1a$Genus <- as.character(RelGenB1a$Genus) 
RelGenB1a$Genus[RelGenB1a$Abundance < 0.03] <- "< 3% abund."

#Unburned.........................................................
GenU1a <- tax_glom(Unburned, taxrank = 'Genus')
RelGenUa <- transform_sample_counts(GenU1a, function(x) x/sum(x))
RelGenU1a <- psmelt(RelGenUa)
RelGenU1a$Genus <- as.character(RelGenU1a$Genus) 
RelGenU1a$Genus[RelGenU1a$Abundance < 0.03] <- "< 3% abund."


#Crete directory for storage ........................................................................
dir.create(file.path("1-Analysis/RelativeAbundance/TurnOverRates/Data"), recursive = TRUE)

write.csv(RelGen1,"1-Analysis/RelativeAbundance/TurnOverRates/Data/RelativeAbund-TurnoverSoil.csv")
write.csv(RelGen1a,"1-Analysis/RelativeAbundance/TurnOverRates/Data/RelativeAbund-TurnoverSoil3per.csv")
write.csv(RelGenB1,"1-Analysis/RelativeAbundance/TurnOverRates/Data/RelativeAbund-TurnoverSoil-Burned.csv")
write.csv(RelGenU1,"1-Analysis/RelativeAbundance/TurnOverRates/Data/RelativeAbund-TurnoverSoil-Unburned.csv")
write.csv(RelGenB1a,"1-Analysis/RelativeAbundance/TurnOverRates/Data/RelativeAbund-TurnoverSoil-Burned3Per.csv")
write.csv(RelGenU1a,"1-Analysis/RelativeAbundance/TurnOverRates/Data/RelativeAbund-TurnoverSoil-Unburned3Per.csv")



#---
#----
#*******************************************************************************************************----
#............AGGLOMORATE DATA BY TAXA OF INTEREST ....................................................
#*******************************************************************************************************----
#Calculate rel abun per variables of interest to calculate turnover

#Calculate abundance per treatment, tsf and genus ................
All <- RelGen1 %>%
  filter(!is.na(Genus)) %>%
  group_by(Treatment,TSFdays,Genus) %>%
  summarise(n_obs = n(),
      TotalAbund= sum(Abundance)) %>%
  as.data.frame();All

Filtered<- RelGen1a %>%
  filter(!is.na(Genus)) %>%
  group_by(Treatment,TSFdays,Genus) %>%
  summarise(n_obs = n(),
    TotalAbund= sum(Abundance)) %>%
  as.data.frame();Filtered


Burned<- RelGenB1%>%
  filter(!is.na(Genus)) %>%
  group_by(Plot,TSFdays,Genus) %>%
  summarise(n_obs = n(),
    TotalAbund= sum(Abundance)) %>%
  as.data.frame();Burned

Unburned<- RelGenU1%>%
  filter(!is.na(Genus)) %>%
  group_by(Treatment,TSFdays,Genus) %>%
  summarise(n_obs = n(),
    TotalAbund= sum(Abundance)) %>%
  as.data.frame();Unburned

Burned1<- RelGenB1a%>%
  filter(!is.na(Genus)) %>%
  group_by(Treatment,TSFdays,Genus) %>%
  summarise(n_obs = n(),
            TotalAbund= sum(Abundance)) %>%
  as.data.frame();Burned1

Unburned1<- RelGenU1a%>%
  filter(!is.na(Genus)) %>%
  group_by(Treatment,TSFdays,Genus) %>%
  summarise(n_obs = n(),
            TotalAbund= sum(Abundance)) %>%
  as.data.frame();Unburned1



#Crete directory for storage ........................................................................
write.csv(All, "1-Analysis/RelativeAbundance/TurnOverRates/Data/TurnoverSoilAll.csv")
write.csv(Filtered, "1-Analysis/RelativeAbundance/TurnOverRates/Data/TurnoverSoil3per.csv")
write.csv(Burned, "1-Analysis/RelativeAbundance/TurnOverRates/Data/TurnoverSoil-Burned.csv")
write.csv(Burned1, "1-Analysis/RelativeAbundance/TurnOverRates/Data/TurnoverSoil-Burned3Per.csv")
write.csv(Unburned, "1-Analysis/RelativeAbundance/TurnOverRates/Data/TurnoverSoil-Unburned.csv")
write.csv(Unburned1, "1-Analysis/RelativeAbundance/TurnOverRates/Data/TurnoverSoil-Unburned3Per.csv")





#---
#----
#************************************************************************************************************----
#.................... DETERMINE TURNOVER RATES....................................................................
#************************************************************************************************************----
#NOTES..........-----------
#turnover = total turnover as well as the proportion of species that either 
            #appear or disappear between time points.

#rank_shift = quantifies relative changes in sp. rank abundances by taking 
            #the sum difference of species ranks in consecutive time points. 
            #This metric goes hand-in-hand with "rank clocks," visualization 
            #tool of shifts in species ranks.

#rate_change analyzes differences in species composition between samples at 
            #increasing time lags. It reflects the rate of directional change 
            #in community composition. Associated function rate_change_interval 
            #returns the full set of community distance values and associated 
            #time lag intervals to use in visualization.


#appearance = will return the proportion of species that appeared relative to 
            #the total number of species observed in both time points. 

#"disappearance" will return the proportion of species that disappeared relative
            #to the total number of species observed in both time points.

#Close notes----

#----------------------------------------------------------------------------------------------

#Load data for calculations.........................................................................
library(codyn)

Abund<-read.csv("1-Analysis/RelativeAbundance/TurnOverRates/Data/TurnoverSoil-TRT-TSF.csv")
AbundB<-read.csv("1-Analysis/RelativeAbundance/TurnOverRates/Data/RelativeAbund-TurnoverSoil-Burned3Per.csv")

head(AbundB)

#Turnover rate in burned soils.......................................................
#Calculate total turnover
TurnO<- turnover(df = Abund, time.var = "TSFdays", species.var = "Genus", 
            abundance.var = "TotalAbund",replicate.var = "Treatment");TurnO


#Appearance:returns proportion of sp. that appeared relative to the total # of sp. 
#observed in both time points.
TurnOA<- turnover(df = Abund, time.var = "TSFdays", species.var = "Genus", 
                  abundance.var = "TotalAbund",replicate.var = "Treatment",
                  metric = "appearance");TurnOA

#Disappearance:returns proportion of sp. that appeared relative to the total # of sp. 
#observed in both time points.
TurnOD<- turnover(df = Abund, time.var = "TSFdays", species.var = "Genus", 
                  abundance.var = "TotalAbund",replicate.var = "Treatment",
                  metric = "disappearance");TurnOD


ROC<-rate_change(df = Abund, time.var = "TSFdays", species.var = "Genus", 
                  abundance.var = "TotalAbund",replicate.var = "Treatment"); ROC


Stability<-community_stability(df = Abund, time.var = "TSFdays",
                abundance.var = "TotalAbund",replicate.var = "Treatment");Stability

Sync<-synchrony(df = Abund, time.var = "TSFdays", abundance.var = "TotalAbund",
          replicate.var = "Treatment", species.var = "Genus",metric = "Loreau");Sync


#Unburned plots which have higher richness and abundance also displayed higher stability and 
#and synchrony


#Export results ....................................................................................
write.csv(TurnO,"1-Analysis/RelativeAbundance/TurnOverRates/Results/TurnoverRate.csv")
write.csv(TurnOA,"1-Analysis/RelativeAbundance/TurnOverRates/Results/Turnover-Appearance.csv")
write.csv(TurnOD,"1-Analysis/RelativeAbundance/TurnOverRates/Results/Turnover-Disappearance.csv")
write.csv(Stability,"1-Analysis/RelativeAbundance/TurnOverRates/Results/Turnover-Stability.csv")
write.csv(Sync,"1-Analysis/RelativeAbundance/TurnOverRates/Results/Turnover-Synchrony.csv")




#RANKSHIFT...............................................................................----
#relative changes in sp rank abundances, which indicate the degree of species reordering
#between two time points. 

#Rankshift <- rank_shift(df = Abund, time.var = "TSFdays", species.var = "Genus", 
           #abundance.var = "TotalAbund",replicate.var = "Treatment");Rankshift
#Select the final time point from the returned time.var_pair
#Rankshift$Year <- as.numeric(substr(Rankshift$year_pair, 1,16))
#ggplot(Rankshift, aes(year, MRS, color=replicate)) + geom_line(size= 2) + theme_bw() 

#----


#----
#-----
#***********************************************************************************************************----
#RATE OF CHAGE..............................................................................................
#***********************************************************************************************************----
#NOTES.........----
#Rate and pattern of variability within a community, which indicates whether species reordering 
#over time is resulting in directional change. This function calculates differences in species 
#composition between samples at increasing time intervals. Differences in species composition are 
#characterized by Euclidean distances, which are calculated on pair-wise communities across the 
#entire time series. For example, a data set with 6 time intervals will have distance values for 
#five one-interval time lags (e.g., time 1 vs time 2, time 2 vs time 3 .), 4 two-interval time 
#lags (e.g., time 1 vs time 3, time 2 vs time 4 .) and so forth. These distance values are 
#regressed against the time lag interval. The slope of the regression line is an indication of 
#the rate and direction of compositional change in the community (Collins, Micheli, and Hartt 2000)

#Calculates the slope of the differences in species composition within a community over 
#increasing time intervals, which provides a measures of the rate of directional change
#in community composition. Slope of regression = indication of the rate and direction of 
#compositional change in the community (Collins, Micheli, and Hartt 2000).


#END NOTES ----


#Community appears to have dirctional change, positive valie 0.4853 vs 0.08 unburned
RateRes <- rate_change(df = Abund, time.var = "TSFdays", species.var = "Genus", 
               abundance.var = "TotalAbund",replicate.var = "Treatment");RateRes


CommRes <- rate_change_interval(df = Abund, time.var = "TSFdays", species.var = "Genus", 
               abundance.var = "TotalAbund",replicate.var = "Treatment");CommRes

#Export tables..................................................................,,,,,,,,,,,,,,,,,,,,
write.csv(CommRes,"1-Analysis/RelativeAbundance/TurnOverRates/Results/RateChangeInterval.csv")
write.csv(RateRes,"1-Analysis/RelativeAbundance/TurnOverRates/Results/RateChangeInterval-Slope.csv")



RateBac<-ggplot(CommRes, aes(interval, distance, color = Treatment)) +
  geom_point(size=3) +
  stat_smooth(method = "lm", se = TRUE, size = 2, alpha=.3)+
  stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label..,
      sep = "~~~~")), formula = y~x, size=4.5, label.x = 5)+
  scale_color_manual(values=  c("#a2673f","#45877f"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size =16, angle=90,hjust=0.5,colour = "black"), 
        axis.text = element_text(size =16, colour = "black"),
        axis.title = element_text(size = 20, colour = "black"),
        legend.title = element_text(size=22),
        legend.text=element_text(size = 20),
  legend.position = "bottom")+
  labs(x = "Time Interval",y="Distance");RateBac
      
RateBac2<-RateBac+ facet_wrap(~Treatment) +
  theme(strip.text.x = element_text(size = 17),
        strip.background =element_rect(fill="white"));RateBac2




#Export graphs and table to calculate significance......................................................
pdf("1-Analysis/RelativeAbundance/TurnOverRates/Graphs/RateChangeInterval-Bac.pdf", height=8, width=8)
RateBac
dev.off()


pdf("1-Analysis/RelativeAbundance/TurnOverRates/Graphs/RateChangeInterval-Bac-Panel.pdf", height=8, width=8)
RateBac2
dev.off()



