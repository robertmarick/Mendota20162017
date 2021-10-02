## Packages
#install.packages(c("ape","dplyr","ggplot2","gplots","lme4","phangorn","phyloseq","plotly","tidyr","vegan","VennDiagram"))


#### Clean up ####
rm(list = ls())
setwd("~/Box/Mendota_16S/dataEdited/16S_filtering")
library(ape)
library(dplyr)
library(ggplot2)
library(gplots)
library(lme4)
library(patchwork)
library(phangorn)
library(phyloseq)
library(plotly)
library(tidyr)
library(vegan)
library(VennDiagram)


#### Read in data ####
#OTU table
OTU.counts.data = read.table("5M_cleaner.count", header=TRUE, sep="\t", stringsAsFactors = FALSE)

# Read in taxonomy data
taxaData.needed <- read.table(file = "5M_cleaner.taxonomy",
                              sep = ",",
                              stringsAsFactors = FALSE)
# Name the OTU and taxonomy columns
colnames(taxaData.needed) <- c("Representative_Sequence", "kingdom", "phylum", "class",
                        "order", "lineage", "clade", "tribe")

#Metadata. Since we made this in Excel, not mothur, we can use the "row.names" modifier to automatically name the rows by the values in the first column (sample names)
meta.df = read.csv("16s_metadata.csv",
                   header = TRUE,
                   stringsAsFactors = FALSE) %>%
  select(sampleID, year, redoxState, sulfide_uM, Mn_uM)



#### Variables to set ####
percentWanted <- 0.01
taxaWanted <- "lineage"
color.vector <- c("#009E73", "#E69F00", "#8A2BE2")
names(color.vector) <- c("oxy", "sub", "eux")



#### Function to plot ordination ####
plot.ordination <- function(OTU.counts = OTU.counts.data,
                            taxaData = taxaData.needed,
                            meta = meta.df,
                            group.by.taxa = TRUE,
                            percentWanted = 0.01,
                            taxaWanted = "lineage") {
  
  if (group.by.taxa == TRUE) {
    
    #### Group by taxonomy ####
    taxaData$taxLevelWanted <- taxaData[, taxaWanted]
    taxaData <- taxaData %>%
      select(Representative_Sequence, taxLevelWanted)
    # wanted <- colnames(OTU)[1:length(colnames(OTU))-1]
    # rel.OTUs.filtered <- filter(rel.OTUs.filtered, rel.OTUs.filtered$SeqID %in% taxData$SeqID)
    
    OTU <- left_join(taxaData,
                     OTU.counts) %>%
      select(-Representative_Sequence) %>%
      #filter out unclassified
      filter(!grepl("unclassified", taxLevelWanted)) %>%
      group_by(taxLevelWanted) %>%
      summarise_all(sum) %>%
      ungroup() %>%
      rename(Representative_Sequence = taxLevelWanted)
    
    title.to.use = paste("CCA ordination, aggregated by ", taxaWanted,
                         ",\nfiltered at ", percentWanted*100, "%",
                         sep = "")
    
  } else {
    
    OTU <- OTU.counts
    title.to.use = paste("CCA ordination of all ASVs,\nfiltered at ",
                         percentWanted*100, "%",
                         sep = "")
  }

  
  ##### Loading in OTU data (not needed if grouping by taxonomy) #####
  #OTU table
  #OTU = read.table("5M_cleaner.count", header=TRUE, sep="\t", stringsAsFactors = FALSE)
  
  colnames(OTU) <-  gsub("X5M", "5M", colnames(OTU))
  
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
  
  # Remove the entire total column
  #OTUs.filtered <- OTUs.filtered %>% select(-total)
  
  # Transpose the OTU table so that the samples are the rows instead of the columns
  OTUs.transposed <- t(OTUs.filtered)
  
  #Normalize the OTU table
  rel.OTUs.filtered <- OTUs.transposed / rowSums(OTUs.transposed)
  rel.OTUs.filtered <- as.data.frame(rel.OTUs.filtered)
  
  # Add sample names as column
  sampleIDs <- row.names(rel.OTUs.filtered)
  row.names(rel.OTUs.filtered) <- NULL
  rel.OTUs.filtered$sampleID <- sampleIDs
  
  
  #### Add metadata and filter out oxic samples ####
  allData <- left_join(meta, rel.OTUs.filtered)
  # filter out oxy samples
  allData <- filter(allData, redoxState != "oxy")
  meta <- filter(meta, redoxState != "oxy")
  
  
  #### Find names of columns needed ####
  all.unnamed.taxa.ids <- colnames(allData)[grep("V", colnames(allData))]
  
  
  #### Filter out by cutoff ####
  unnamed.taxa.to.include <- gather(allData,
                                    key = unnamed_taxa,
                                    value = Abundance,
                                    all_of(all.unnamed.taxa.ids)) %>%
    filter(Abundance > percentWanted) %>%
    select(unnamed_taxa) %>%
    unlist() %>% unique()

  
  #### Generate CCA ####
  samples.dist <- cca(formula = allData[, unnamed.taxa.to.include] ~ sulfide_uM + Mn_uM + year, allData)
  sumPlot <- summary(samples.dist)
  
  
  #### Set up for CCA plotting ####
  plotData <- as.data.frame(sumPlot$sites)
  plotData$redoxState <- meta$redoxState
  plotData$Year <- as.character(meta$year)
  # Metadata scores
  sulf.scores <- as.data.frame(t(sumPlot$biplot["sulfide_uM",]))
  mn.scores <- as.data.frame(t(sumPlot$biplot["Mn_uM",]))
  year.scores <- as.data.frame(t(sumPlot$biplot["year",]))
  
  
  #### Generate plot ####
  brayPlot <- ggplot(plotData, aes(x = CCA1, y = CCA2, col = redoxState, shape = Year)) +
    #place groups in aes!
    geom_point(size = 1, aes(colour = redoxState, shape = Year)) +
    #geom_text(aes(label = data.scores$sulfide),hjust=0, vjust=0, size = 4)+ #this adds labels as the Mn_uM to each point
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
          axis.text.x = element_text(colour = "black", face = "bold", size = 12),
          legend.text = element_text(size = 12, face ="bold", colour ="black"),
          legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
          legend.title = element_text(size = 14, colour = "black", face = "bold"),
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
          legend.key=element_blank()) +
    geom_segment(data = sulf.scores, aes(x = 0, xend = CCA1, y = 0, yend = CCA2), #geom_segment() plots arrows. put proper envfit object inside
                 arrow = arrow(length = unit(0.25, "cm")), colour = "black", inherit.aes = FALSE) +
    geom_text(data = sulf.scores, #labels the environmental variable arrows
              aes(x=CCA1,y=CCA2,label="[sulfide]"), size = 5, inherit.aes = FALSE, hjust=-0.1) +
    geom_segment(data = mn.scores, aes(x = 0, xend = CCA1, y = 0, yend = CCA2), #geom_segment() plots arrows. put proper envfit object inside
                 arrow = arrow(length = unit(0.25, "cm")), colour = "black", inherit.aes = FALSE) +
    geom_text(data = mn.scores, #labels the environmental variable arrows
              aes(x=CCA1,y=CCA2,label="[Mn]"), size = 5, inherit.aes = FALSE, hjust=-0.1) +
    geom_segment(data = year.scores, aes(x = 0, xend = CCA1, y = 0, yend = CCA2), #geom_segment() plots arrows. put proper envfit object inside
                 arrow = arrow(length = unit(0.25, "cm")), colour = "black", inherit.aes = FALSE) +
    geom_text(data = year.scores, #labels the environmental variable arrows
              aes(x=CCA1,y=CCA2,label="year"), size = 5, inherit.aes = FALSE, hjust=-0.1) +
    stat_ellipse(aes(lty = Year)) + #stat_ellipse is adding the elipses
    # #colors here, make sure correct amount
    scale_colour_manual(values = color.vector) +
    ggtitle(title.to.use)
  
  #plotting brayplot and envfit object
  brayPlot
  
}



#### Make plots ####
ordination.by.lineage.0percent <- plot.ordination(group.by.taxa = TRUE,
                                                  percentWanted = 0.00,
                                                  taxaWanted = "lineage")
ordination.by.phyla.0percent <- plot.ordination(group.by.taxa = TRUE,
                                                percentWanted = 0.00,
                                                taxaWanted = "phylum")

ordination.by.lineage.1percent <- plot.ordination(group.by.taxa = TRUE,
                                                  percentWanted = 0.01,
                                                  taxaWanted = "lineage")
ordination.by.phyla.1percent <- plot.ordination(group.by.taxa = TRUE,
                                                percentWanted = 0.01,
                                                taxaWanted = "phylum")

(ordination.by.lineage.0percent + ordination.by.lineage.1percent) / (ordination.by.phyla.0percent + ordination.by.phyla.1percent)


#Plots for Figure
ordination.by.phyla.1percent <- plot.ordination(group.by.taxa = TRUE,
                                                percentWanted = 0.01,
                                                taxaWanted = "phylum")

ordination.by.order.1percent <- plot.ordination(group.by.taxa = TRUE,
                                                percentWanted = 0.01,
                                                taxaWanted = "order")

ordination.by.lineage.1percent <- plot.ordination(group.by.taxa = TRUE,
                                                  percentWanted = 0.01,
                                                  taxaWanted = "lineage")

(ordination.by.phyla.1percent + ordination.by.order.1percent + ordination.by.lineage.1percent)

