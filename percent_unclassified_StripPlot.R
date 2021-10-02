library(ape)
library(dplyr)
library(ggplot2)
library(lme4)
library(phangorn)
library(phyloseq)
library(plotly)
library(tidyr)
library(vegan)
library(VennDiagram)
library(gridExtra)
library(RColorBrewer)

rm(list = ls())
setwd("/Users/robertmarick/Box/Mendota_16S/dataEdited/16S_filtering") 

##SEE LINE 58 FOR VARIBLES TO CHANGE ANALYSIS

#OTU table
OTU = read.table("5M_cleaner.count", header=TRUE, sep="\t", stringsAsFactors = FALSE)
#Cleaned tax data
#taxaData = read.csv("cleanedTax.csv", header=TRUE, stringsAsFactors = FALSE, row.names=1, sep=",")
taxaData <- read.table(file = "5M_cleaner.taxonomy",
                       sep = ",", 
                       stringsAsFactors = FALSE)

# Name the OTU and taxonomy columns
colnames(taxaData) <- c("SeqID", "kingdom", "phylum", "class", 
                        "order", "lineage", "clade", "tribe")

#metadata
meta = read.csv("16s_metadata.csv", sep=",",header=TRUE, stringsAsFactors = FALSE)

#To start need to create filter all the OTUs 
# Remove the singletons
OTU$total ==  1
OTU[(OTU$total==1), ]
OTU %>% filter(total == 1)
OTUs.filtered <- OTU %>% filter(total > 1)

colnames(OTU) <-  gsub("X5M", "5M", colnames(OTU))

# Name the OTU table's rows as the Representative sequences
row.names(OTUs.filtered) = OTUs.filtered$Representative_Sequence
# Remove the entire Representative Sequence column (the rep seqs are now the row names)
OTUs.filtered <- OTUs.filtered %>% select(-Representative_Sequence)
OTUs.filtered <- OTUs.filtered %>% select(meta[,"sampleID"])

# Remove the entire total column
OTUs.filtered <- OTUs.filtered %>% select(-total)

#Normalize the OTU table
OTUs.transposed <- t(OTUs.filtered)
rel.OTUs.filtered <- OTUs.transposed /rowSums(OTUs.transposed) 
rel.OTUs.filtered <- t(rel.OTUs.filtered)
rel.OTUs.filtered <- as.data.frame(rel.OTUs.filtered) 
SeqIDs <- row.names(rel.OTUs.filtered)
row.names(rel.OTUs.filtered) <- NULL
rel.OTUs.filtered$SeqID <- SeqIDs

figureMaker <- function(percentWanted = 0.00,taxaWanted = "phylum", leftText = FALSE){
  taxData <- taxaData
  taxData$taxa <- taxData[, taxaWanted]
  wanted <- colnames(rel.OTUs.filtered)[1:length(colnames(rel.OTUs.filtered))-1]
  rel.OTUs.filtered <- filter(rel.OTUs.filtered, rel.OTUs.filtered$SeqID %in% taxData$SeqID)
  
  allData <- left_join(rel.OTUs.filtered, taxData) 
  
  #filter out unclassified
  allData <- allData[grepl("unclassified", allData$taxa),]
  
  allData <-  gather(allData, key = sampleID, value = Abundance, wanted) %>%
    group_by(sampleID) %>%
    summarize(Abundance = sum(Abundance)) %>%
    ungroup()
  
  #add meta data: day, month, year, depth and rep to allData
  allData$sampleID <- gsub("X5M", "5M", allData$sampleID)
  allData <- left_join(allData, meta)
  
  stripData <- allData[,c("sampleID", "Abundance", "redoxState")] 
  stripData <- filter(stripData, !is.na(stripData$redoxState))
  stripData$percentUnclassified <- stripData$Abundance * 100 
  
  if(leftText == TRUE){
    stripPlot <- ggplot(stripData, aes(x = redoxState, y = percentUnclassified, col = redoxState)) +
      geom_jitter() +
      theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
            axis.text.x = element_text(colour = "black", face = "bold", size = 12),
            legend.text = element_text(size = 12, face ="bold", colour ="black"), 
            legend.position = "none", axis.title.y = element_text(face = "bold", size = 14), 
            axis.title.x.bottom = element_blank(), 
            legend.title = element_text(size = 14, colour = "black", face = "bold"), 
            panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
            legend.key=element_blank(), plot.title = element_text(face = "bold", size = 12, colour = "black", hjust = 0.5)) +
      labs(x = "Redox State",y = "% unclassified", colour = "Redox State") + 
      ggtitle(paste(taxaWanted, sep = "")) +
      #colors here, make sure correct amount
      scale_colour_manual(values = c("blue", "green", "orange")) +
      ylim(0,30)
  }
  else{
    stripPlot <- ggplot(stripData, aes(x = redoxState, y = percentUnclassified, col = redoxState)) +
      geom_jitter() +
      theme(axis.text.y.left = element_blank(), 
            axis.text.x = element_text(colour = "black", face = "bold", size = 12),
            legend.text = element_text(size = 12, face ="bold", colour ="black"), 
            legend.position = "none", axis.title.y = element_text(face = "bold", size = 14), 
            axis.title.x.bottom = element_blank(), 
            axis.title.y.left = element_blank(),
            legend.title = element_text(size = 14, colour = "black", face = "bold"), 
            panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
            legend.key=element_blank(), plot.title = element_text(face = "bold", size = 12, colour = "black", hjust = 0.5)) +
      labs(x = "Redox State", colour = "Redox State") + 
      ggtitle(paste(taxaWanted, sep = "")) +
      #colors here, make sure correct amount
      scale_colour_manual(values = c("blue", "green", "orange")) +
      ylim(0,30)
  }
  
  return(stripPlot)
}

phylaPlot <- figureMaker(leftText = TRUE)
lineagePlot <- figureMaker(taxaWanted = "lineage")
tribePlot <- figureMaker(taxaWanted = "tribe")

grid.arrange(phylaPlot, lineagePlot, tribePlot, nrow = 1, ncol = 3, 
             top = text_grob("Percent Unclassified of Different Taxonomic Levels", face = "bold", size = 14))
