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
library(tidyverse)

rm(list = ls())
setwd("~/Box/Mendota_16S/")

#### Read in data ####
#OTU table
OTU = read.table("dataEdited/16S_filtering/5M_cleaner.count", header=TRUE, sep="\t", stringsAsFactors = FALSE)
colnames(OTU) <-  gsub("X5M", "5M", colnames(OTU))

#metadata
meta = read.csv("dataEdited/16S_filtering/16s_metadata.csv", sep=",",header=TRUE, stringsAsFactors = FALSE)


#### Clean up OTU data ####

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
OTUs.filtered <- OTUs.filtered %>% select(meta[,"sampleID"])

# Remove the entire total column
OTUs.filtered <- OTUs.filtered %>% select(-total)



#### Normalize the OTU table ####
OTUs.transposed <- t(OTUs.filtered)
rel.OTUs.filtered <- OTUs.transposed /rowSums(OTUs.transposed) 
rel.OTUs.filtered <- t(rel.OTUs.filtered)
rel.OTUs.filtered <- as.data.frame(rel.OTUs.filtered) 
rel.OTUs.filtered <- t(rel.OTUs.filtered)


#### Calculate Shannon diversity ####
d <- diversity(rel.OTUs.filtered, index = "shannon")
shannonData <- scores(d)
shannonData <- as.data.frame(shannonData)
sampIDs <- rownames(shannonData)
rownames(shannonData) <- NULL
shannonData$sampleID <- sampIDs

#### Prepare data for plotting ####
allData <- right_join(shannonData, meta) %>%
  mutate(status = paste(redoxState, year, sep = "-"),
         status = fct_relevel(status,
                              c("oxy-2016", "oxy-2017", "sub-2016",
                                "sub-2017", "eux-2016", "eux-2017")),
         year = as.character(year)) %>%
  select(sampleID, year, redoxState, status, Dim1)


#### Prep vector for year and status ####
color.vector <- c("#009E73", "#E69F00", "#8A2BE2")
names(color.vector) <- c("oxy", "sub", "eux")
point.vector <- c(17, 16)
names(point.vector) <- c(2016, 2017)


#### Generate figure ####
stripPlot <- allData %>%
  ggplot(aes(x = status,
             y = Dim1,
             group = status,
             col = redoxState,
             shape = year)) +
  geom_boxplot() +
  geom_jitter() +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12),
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "none", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank(), plot.title = element_text(face = "bold", size = 16, colour = "black", hjust = 0.5)) +
  labs(x = "",y = "Shannon Diversity Index (H)", colour = "Redox State and Year") + 
  # ggtitle(paste("Shannon Diversity", sep = "")) + 
  #colors here, make sure correct amount
  scale_colour_manual(values = color.vector) +
  scale_shape_manual(values = point.vector) #+


#### Save out plot ####
pdf("results/alpha_diversity/shannon_diversity.pdf",
    height = 5,
    width = 8)
stripPlot
dev.off()


#### Statistical analysis of differences in diversity ####
shannon.aov <- aov(Dim1 ~ year * redoxState,
                   data = allData)
# Check normality
par(mfrow = c(1,2))
hist(shannon.aov$residuals,
     breaks = 20)
plot(density(shannon.aov$residuals),
     main="Density plot of residuals",
     ylab="Density",
     xlab="Residuals")
shapiro.test(shannon.aov$residuals)
# QQ-normal plot
qqnorm(shannon.aov$residuals)
qqline(shannon.aov$residuals)
# Residuals are fairly close to normally distributed, but not quite.
# Little bit of left skew, according to Q-Q plot and density plot

# Let's go ahead and look at the output
summary(shannon.aov)
# Major impact of redox status on Shannon diversity.
# Smaller but still significant effect of year.
# No interactive effect, suggests interannual effects
# are independent of redox status and vice versa.
