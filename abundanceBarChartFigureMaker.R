library(ape)
library(dplyr)
library(ggplot2)
library(lme4)
library(plotly)
library(tidyr)
library(vegan)
library(VennDiagram)
library(gridExtra)
library(RColorBrewer)
library(ggpubr)

#"ape", "dplyr", "lme4", "phangorn", "phyloseq", "plotly", "tidyr", "vegan", "VennDiagram", "gridExtra", "RColorBrewer", "ggpubr", "indicspecies"
rm(list = ls())
setwd("~/Box/Mendota_16S/dataEdited/16S_filtering") 

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
meta = read.csv("metaData_withNOx.csv", sep=",",header=TRUE, stringsAsFactors = FALSE)

#To start need to create filter all the OTUs 
# Remove the singletons
OTU$total ==  1
OTU[(OTU$total==1), ]
OTU %>% filter(total == 1)
OTUs.filtered <- OTU %>% filter(total > 1)

# Name the OTU table's rows as the Representative sequences
row.names(OTUs.filtered) = OTUs.filtered$Representative_Sequence
# Remove the entire Representative Sequence column (the rep seqs are now the row names)
OTUs.filtered <- OTUs.filtered %>% select(-Representative_Sequence)

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

figureMaker <- function(percentWanted = 1, taxaWanted = "phylum", specificTaxa = FALSE, composition = FALSE, 
                        specificAnalyte = "Proteobacteria", oneTaxaUp = "phylum", abundanceScale = 100) {
  
  #Use this next two lines if you want a specific taxonomic unit ie all the orders in deltaproteobacteria
  if (composition == TRUE) {
    taxData$taxa <- taxData[, taxaWanted] #for column put one level above taxaWanted
    taxData$specific <-  taxData[, oneTaxaUp]
    #taxData <- filter(taxData, taxData$taxa == specificAnalyte)
    #taxData$taxa <- taxData[, taxaWanted]
    wanted <- colnames(rel.OTUs.filtered)[1:length(colnames(rel.OTUs.filtered))-1]
    rel.OTUs.filtered <- filter(rel.OTUs.filtered, rel.OTUs.filtered$SeqID %in% taxData$SeqID)
    allData <- left_join(rel.OTUs.filtered, taxData) %>% 
      gather(key = sampleID, value = Abundance, wanted) %>%
      group_by(taxa, specific, sampleID) %>%
      summarize(Abundance = sum(Abundance)) 
  }
  else {
    taxData$taxa <- taxData[, taxaWanted]
    wanted <- colnames(rel.OTUs.filtered)[1:length(colnames(rel.OTUs.filtered))-1]
    rel.OTUs.filtered <- filter(rel.OTUs.filtered, rel.OTUs.filtered$SeqID %in% taxData$SeqID)
    allData <- left_join(rel.OTUs.filtered, taxData) %>% 
      gather(key = sampleID, value = Abundance, wanted) %>%
      group_by(taxa, sampleID) %>%
      summarize(Abundance = sum(Abundance)) 
  }
  
  #add meta data: day, month, year, depth and rep to allData
  allData$sampleID <- gsub("X5M", "5M", allData$sampleID)
  allData <- right_join(allData, meta)
  
  #get data from of a specfifc parameter and normalize by month
  desired816 <- filter(allData, month == 8 & year == 2016 & day == 9) %>%
    group_by(depth) %>%
    mutate(relAbsDepth = Abundance / sum(Abundance)) %>%
    ungroup()
  if (specificTaxa == TRUE) {desired816 <- filter(desired816, taxa == specificAnalyte)}  #use this line to see just 1 type ie just see Geobacter
  
  desired916 <- filter(allData,month == 9 & year == 2016 & day == 22) %>%
    group_by(depth) %>%
    mutate(relAbsDepth = Abundance / sum(Abundance)) %>%
    #filter(taxa == "Methylococcales") %>% #use this line to see just 1 type ie just see Geobacter
    ungroup()
  if (specificTaxa == TRUE) {desired916 <- filter(desired916, taxa == specificAnalyte)}
  
  desired1016 <- filter(allData,month == 10 & year == 2016 & day == 4) %>%
    group_by(depth) %>%
    mutate(relAbsDepth = Abundance / sum(Abundance)) %>%
    #filter(taxa == "Methylococcales") %>% #use this line to see just 1 type ie just see Geobacter
    ungroup()
  if (specificTaxa == TRUE) {desired1016 <- filter(desired1016, taxa == specificAnalyte)}
  
  desired817 <- filter(allData,month == 8 & year == 2017 & day == 11) %>%
    group_by(depth) %>%
    mutate(relAbsDepth = Abundance / sum(Abundance)) %>%
    #filter(taxa == "Methylococcales") %>% #use this line to see just 1 type ie just see Geobacter
    ungroup()
  if (specificTaxa == TRUE) {desired817 <- filter(desired817, taxa == specificAnalyte)}
  
  desired917 <- filter(allData,month == 9 & year == 2017 & day == 8) %>%
    group_by(depth) %>%
    mutate(relAbsDepth = Abundance / sum(Abundance)) %>%
    #filter(taxa == "Methylococcales") %>% #use this line to see just 1 type ie just see Geobacter
    ungroup()
  if (specificTaxa == TRUE) {desired917 <- filter(desired917, taxa == specificAnalyte)}
  
  desired1017 <- filter(allData,month == 10 & year == 2017 & day == 4) %>%
    group_by(depth) %>%
    mutate(relAbsDepth = Abundance / sum(Abundance)) %>%
    #filter(taxa == "Methylococcales") %>% #use this line to see just 1 type ie just see Geobacter
    ungroup()
  if (specificTaxa == TRUE) {desired1017 <- filter(desired1017, taxa == specificAnalyte)}
  
  desired816$relAbsDepth <- desired816$relAbsDepth * 100
  desired816 <- filter(desired816, relAbsDepth >= percentWanted)
  
  desired916$relAbsDepth <- desired916$relAbsDepth * 100
  desired916 <- filter(desired916, relAbsDepth >= percentWanted)
  
  desired1016$relAbsDepth <- desired1016$relAbsDepth * 100
  desired1016 <- filter(desired1016, relAbsDepth >= percentWanted)
  
  desired817$relAbsDepth <- desired817$relAbsDepth * 100
  desired817 <- filter(desired817, relAbsDepth >= percentWanted)
  
  desired917$relAbsDepth <- desired917$relAbsDepth * 100
  desired917 <- filter(desired917, relAbsDepth >= percentWanted)
  
  desired1017$relAbsDepth <- desired1017$relAbsDepth * 100
  desired1017 <- filter(desired1017, relAbsDepth >= percentWanted)
  
  if (composition == TRUE) {
    desired816 <- filter(desired816, specific == specificAnalyte)
    desired916 <- filter(desired916, specific == specificAnalyte)
    desired1016 <- filter(desired1016, specific == specificAnalyte)
    desired817 <- filter(desired817, specific == specificAnalyte)
    desired917 <- filter(desired917, specific == specificAnalyte)
    desired1017 <- filter(desired1017, specific == specificAnalyte)
  }
  
  #creates a color palette based on # of unique groups
  taxaColor <- filter(allData, Abundance >= percentWanted/100)
  colorCount <- length(unique(taxaColor$taxa)) 
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  my_colors <- getPalette(colorCount)
  set.seed(3)
  my_colors <- sample(my_colors)
  names(my_colors) <-  levels(factor(unique(taxaColor$taxa)))
  cols <- c("Oxy"="green","Sub"="orange","Eux"="blue")
  
  #Abundance bar charts 
  g816 <- ggplot(desired816) +
    geom_bar(aes(x=depth, y=relAbsDepth, fill = taxa), stat= "identity") +
    #geom_text(aes(x = depth, y = (abundanceScale*1.02), label = redoxState), size = 3) +
    scale_fill_manual(name = c("Phyla", "RedoxState"), values = c(my_colors, cols)) +
    labs(y="% Abundance", x = "Depth (m)", fill = taxaWanted) +
    theme(axis.text.y = element_text(colour = "black", face = "bold", size = 10),
          axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
          legend.text = element_text(size = 6, face ="bold", colour ="black"), 
          legend.position = "none", axis.title.y = element_text(face = "bold", size = 8), 
          axis.title.x = element_text(colour = "black", face = "bold", size = 8), 
          axis.line = element_line(colour = "black", size = 0.5),
          legend.title = element_text(size = 8, colour = "black", face = "bold"), 
          panel.background = element_blank(), 
          plot.title = element_text(colour = "black", size = 10, face = "bold", hjust = 0.5)) +
    coord_flip(xlim = c(24,0), ylim = c(0,abundanceScale)) +
    ggtitle(paste(8, "/", 9, "/", 2016, sep = "")) +
    scale_x_continuous(trans = "reverse") +
    scale_y_continuous(sec.axis = sec_axis(~.)) 
  
  g916 <- ggplot(desired916) +
    geom_bar(aes(x=depth, y=relAbsDepth, fill = taxa), stat= "identity") +
    #geom_text(aes(x = depth, y = (abundanceScale*1.02), label = redoxState), size = 3) +
    scale_fill_manual(name = c("Phyla", "RedoxState"), values = c(my_colors, cols)) +
    labs(y="% Abundance", x = "Depth (m)", fill = taxaWanted) +
    theme(axis.text.y = element_text(colour = "black", face = "bold", size = 10),
          axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
          legend.text = element_text(size = 6, face ="bold", colour ="black"), 
          legend.position = "none", axis.title.y = element_text(face = "bold", size = 8), 
          axis.title.x = element_text(colour = "black", face = "bold", size = 8), 
          axis.line = element_line(colour = "black", size = 0.5),
          legend.title = element_text(size = 8, colour = "black", face = "bold"), 
          panel.background = element_blank(), 
          plot.title = element_text(colour = "black", size = 10, face = "bold", hjust = 0.5)) +
    coord_flip(xlim = c(24,0), ylim = c(0,abundanceScale)) +
    ggtitle(paste(9, "/", 22, "/", 2016, sep = "")) +
    scale_x_continuous(trans = "reverse") +
    scale_y_continuous(sec.axis = sec_axis(~.)) 
  
  g1016 <- ggplot(desired1016) +
    geom_bar(aes(x=depth, y=relAbsDepth, fill = taxa), stat= "identity") +
    #geom_text(aes(x = depth, y = (abundanceScale*1.02), label = redoxState), size = 3) +
    scale_fill_manual(name = c("Phyla", "RedoxState"), values = c(my_colors, cols)) +
    labs(y="% Abundance", x = "Depth (m)", fill = taxaWanted) +
    theme(axis.text.y = element_text(colour = "black", face = "bold", size = 10),
          axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
          legend.text = element_text(size = 6, face ="bold", colour ="black"), 
          legend.position = "none", axis.title.y = element_text(face = "bold", size = 8), 
          axis.title.x = element_text(colour = "black", face = "bold", size = 8), 
          axis.line = element_line(colour = "black", size = 0.5),
          legend.title = element_text(size = 8, colour = "black", face = "bold"), 
          panel.background = element_blank(), 
          plot.title = element_text(colour = "black", size = 10, face = "bold", hjust = 0.5)) +
    coord_flip(xlim = c(24,0), ylim = c(0,abundanceScale)) +
    ggtitle(paste(10, "/", 4, "/", 2016, sep = "")) +
    scale_x_continuous(trans = "reverse") +
    scale_y_continuous(sec.axis = sec_axis(~.)) 
  
  g817 <- ggplot(desired817) +
    geom_bar(aes(x=depth, y=relAbsDepth, fill = taxa), stat= "identity") +
    #geom_text(aes(x = depth, y = (abundanceScale*1.02), label = redoxState), size = 3) +
    scale_fill_manual(name = c("Phyla", "RedoxState"), values = c(my_colors, cols)) +
    labs(y="% Abundance", x = "Depth (m)", fill = taxaWanted) +
    theme(axis.text.y = element_text(colour = "black", face = "bold", size = 10),
          axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
          legend.text = element_text(size = 6, face ="bold", colour ="black"), 
          legend.position = "none", axis.title.y = element_text(face = "bold", size = 8), 
          axis.title.x = element_text(colour = "black", face = "bold", size = 8), 
          axis.line = element_line(colour = "black", size = 0.5),
          legend.title = element_text(size = 8, colour = "black", face = "bold"), 
          panel.background = element_blank(), 
          plot.title = element_text(colour = "black", size = 10, face = "bold", hjust = 0.5)) +
    coord_flip(xlim = c(24,0), ylim = c(0,abundanceScale)) +
    ggtitle(paste(8, "/", 11, "/", 2017, sep = "")) +
    scale_x_continuous(trans = "reverse") +
    scale_y_continuous(sec.axis = sec_axis(~.)) 
  
  g917 <- ggplot(desired917) +
    geom_bar(aes(x=depth, y=relAbsDepth, fill = taxa), stat= "identity") +
    #geom_text(aes(x = depth, y = (abundanceScale*1.02), label = redoxState), size = 3) +
    scale_fill_manual(name = c("Phyla", "RedoxState"), values = c(my_colors, cols)) +
    labs(y="% Abundance", x = "Depth (m)", fill = taxaWanted) +
    theme(axis.text.y = element_text(colour = "black", face = "bold", size = 10),
          axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
          legend.text = element_text(size = 6, face ="bold", colour ="black"), 
          legend.position = "none", axis.title.y = element_text(face = "bold", size = 8), 
          axis.title.x = element_text(colour = "black", face = "bold", size = 8), 
          axis.line = element_line(colour = "black", size = 0.5),
          legend.title = element_text(size = 8, colour = "black", face = "bold"), 
          panel.background = element_blank(), 
          plot.title = element_text(colour = "black", size = 10, face = "bold", hjust = 0.5)) +
    coord_flip(xlim = c(24,0), ylim = c(0,abundanceScale)) +
    ggtitle(paste(9, "/", 8, "/", 2017, sep = "")) +
    scale_y_continuous(sec.axis = sec_axis(~.)) 
  
  g1017 <- ggplot(desired1017) +
    geom_bar(aes(x=depth, y=relAbsDepth, fill = taxa), stat= "identity") +
    #geom_text(aes(x = depth, y = (abundanceScale*1.02), label = redoxState), size = 3) +
    scale_fill_manual(name = c("Phyla", "RedoxState"), values = c(my_colors, cols)) +
    labs(y="% Abundance", x = "Depth (m)", fill = taxaWanted) +
    theme(axis.text.y = element_text(colour = "black", face = "bold", size = 10),
          axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
          legend.text = element_text(size = 6, face ="bold", colour ="black"), 
          legend.position = "none", axis.title.y = element_text(face = "bold", size = 8), 
          axis.title.x = element_text(colour = "black", face = "bold", size = 8), 
          axis.line = element_line(colour = "black", size = 0.5),
          legend.title = element_text(size = 8, colour = "black", face = "bold"), 
          panel.background = element_blank(), 
          plot.title = element_text(colour = "black", size = 10, face = "bold", hjust = 0.5)) +
    ggtitle(paste(10, "/", 4, "/", 2017, sep = "")) +
    coord_flip(xlim = c(24,0), ylim = c(0,abundanceScale)) +
    scale_y_continuous(sec.axis = sec_axis(~.)) 
  
  legendMaker <- ggplot(desired816) +
    geom_bar(aes(x=depth, y=relAbsDepth, fill = taxa), stat= "identity") +
    scale_fill_manual(name = c(taxaWanted, "RedoxState"), values = c(my_colors, cols))
  
  figLegend <- get_legend(legendMaker)
  noLegend <- grid.arrange(g816, g916, g1016, g817, g917, g1017, ncol=3, nrow = 2)
  if (specificTaxa == TRUE) {
    return(grid.arrange(g816, g916, g1016, g817, g917, g1017, ncol=3, nrow = 2, top = text_grob(figureTitle, face = "bold", size = 14)))
  }
  else if (composition == TRUE) {
    return(grid.arrange(as_ggplot(figLegend), noLegend, ncol = 2, nrow = 1, widths = c(3,12), top = text_grob(figureTitle, face = "bold", size = 14)))
  }
  else {
    return(grid.arrange(as_ggplot(figLegend), noLegend, ncol = 2, nrow = 1, widths = c(1,6), top = text_grob(figureTitle, face = "bold", size = 14)))
  }
}

figureTitle <- "Composition of Orders of Gammaproteobacteria"
figureTime <- figureMaker(percentWanted = 1, taxaWanted = "phylum", specificTaxa = FALSE, composition = FALSE, 
                          specificAnalyte = "Proteobacteria", oneTaxaUp = "phylum", abundanceScale = 100)
