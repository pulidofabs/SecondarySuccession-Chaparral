Note: All files here are for the bacterial analysis, but the graph for bacterial biomass and richness and their interaction with ash depth and time since fire is located in the fungal folder, file named "#9-RichnessBiomass-AshDepth-Bac-Fun", this file is to create the graph for both bacteria and fungi for both biomass and species richness. 

FILES: BASED ON NUMBER #-Name
1. Quality control: takes a raw AV table exported from qiime2, rarefies the dataset, removes controls, and calculates alpha diversity metrics. All files from here are saved for downstream analysis 
2. Alpha-fun-graphs-DesStats-SigNif: file used to create the alpha diversity graphs and calculate the descriptive statistics using dplyr (percent changes, mean, variance, etc.). Significance analysis is based on a full/global model using a negative binomial model, including the nested levels. 
3. Beta-Graphs-TP-Signif: File used to create beta diversity graphs and determine significance using PERMANOVA. TP refers to facet wrapping the NMDS plot by time point (time since fire) to see the diversity at each sampling time point. 
4. TO-Rates-Stability-Sync-Fun: R script is used to calculate the turnover rates, community stability, and synchrony for the bacterial datasets as a measure of succession
5. Relative abundance bar plots, TSF: Uses phyloseq and the already rarefied exported qiime2 outputs (.qza files) to calculate the relative abundance for fungi per treatment (burned vs. unburned) and per time since fire for the burned communities and unburned communities independently) to visualize successional patterns. 
6. Rel-Abund-PerChang-Trt-TSF: Uses the files created above to calculate each taxon of interest's percent change and total sequence at each specific time point or treatment (burned vs. unburned) to quantify succession. 
7. RelAbund-LineGraphs-Top5: uses phyloseq to calculate the top5 taxa in the community and create a line graph over time of their respective changes (At the genus level). In this case, relative abundance is in total relative abundance as we do not select the taxa above 3 percent relative abundance. 
8. Rarefaction curves: Create rarefaction curves using phyloseq and the unrarefied qiime2 files; in other words, we are using the raw files from qiime2, rarefying using phyloseq, and creating rarefaction curves per treatment, time point, sites (graphs).

**See note above

