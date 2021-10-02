## Packages
#install.packages(c("ape","dplyr","ggplot2","gplots","lme4","phangorn","phyloseq","plotly","tidyr","vegan","VennDiagram")) 
library(ape)
library(dplyr)
library(ggplot2)
library(gplots)
library(lme4)
library(phangorn)
library(phyloseq)
library(plotly)
library(tidyr)
library(vegan)
library(VennDiagram)

rm(list = ls())

setwd("~/Box/Mendota_16S/dataEdited/16S_filtering") 
##### Loading in data #####
#OTU table
OTU = read.table("5M_cleaner.count", header=TRUE, sep="\t", stringsAsFactors = FALSE)

#Metadata. Since we made this in Excel, not mothur, we can use the "row.names" modifier to automatically name the rows by the values in the first column (sample names)
meta = read.csv("16s_metadata.csv", sep=",",header=TRUE, stringsAsFactors = FALSE)

colnames(OTU) <-  gsub("X5M", "5M", colnames(OTU))

#####Beta diversity #####

# Remove the singletons
OTU$total ==  1
OTU[(OTU$total==1), ]
OTU %>% filter(total == 1)
OTUs.filtered <- OTU %>% filter(total > 1)

# Name the OTU table's rows as the Representative sequences
row.names(OTUs.filtered) = OTUs.filtered$Representative_Sequence

# Remove the entire Representative Sequence column (the rep seqs are now the row names)
OTUs.filtered <- OTUs.filtered %>% select(-Representative_Sequence)
OTUs.filtered <- OTUs.filtered %>% select(meta[,"sampleID"])

# Transpose the OTU table so that the samples are the rows instead of the columns
OTUs.transposed <- t(OTUs.filtered)

#Normalize the OTU table
rel.OTUs.filtered <- OTUs.transposed /rowSums(OTUs.transposed)

###### Bray-Curtis metaMDS #####
BC.nmds = metaMDS(rel.OTUs.filtered, distance= "bray", k=2, trymax=1000)

# Make Bray-Curtis graph
par(mfrow = c(1, 1))

#### Robert's ploting code ####
#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(BC.nmds),stringsAsFactors = FALSE)

#add columns to data frame these will serve as groups
data.scores$redoxState = meta[,"redoxState"]
data.scores$Year <- as.character(meta$year)

brayPlot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2, col = redoxState, shape = Year)) + 
  #place groups in aes!
  geom_point(size = 1, aes(colour = redoxState))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  stat_ellipse(aes(lty = Year)) + #stat_ellipse is adding the elipses
  labs(x = "NMDS1", colour = "redoxState", y = "NMDS2")  + 
  #colors here, make sure correct amount
  scale_colour_manual(values = c("#8A2BE2","#009E73","#E69F00")) 

brayPlot

