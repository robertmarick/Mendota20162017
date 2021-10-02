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
library(ggpubr)

rm(list = ls())
setwd("~/Box/Mendota_16S/dataEdited/16S_filtering") 

#metadata
meta = read.csv("metaData_withNOx.csv", sep=",",header=TRUE, stringsAsFactors = FALSE)

waterChemistry <- function(desiredMonth = 8, desiredDay = 9, desiredYear = 2016, legendGen = FALSE) {
  desired <- filter(meta, month == desiredMonth & year == desiredYear & day == desiredDay)
  
  #For Sulfide line plot 
  desired <- group_by(desired, depth) %>%
    mutate(avgSulfide_uM = mean(sulfide_uM))
  
  #For DO line plot 
  desired <- group_by(desired, depth) %>%
    mutate(avgDO = mean(DO_mg_L))
  
  #For Mn line plot 
  desired <- group_by(desired, depth) %>%
    mutate(avguMMn = mean(Mn_uM))
  
  # #For NOx line plot 
  desired <- group_by(desired, depth) %>%
    mutate(avgNOx = mean(NOX_ppb))
  
  oxyData <- filter(desired, redoxState == "oxy")
  euxData <- filter(desired, redoxState == "eux")
  subData <- filter(desired, redoxState == "sub")
  oxyStart <- 0
  oxyEnd <- max(oxyData$depth) 
  euxStart <- min(euxData$depth)
  euxEnd <- 25
  subStart <- min(subData$depth)
  subEnd <- max(subData$depth) 
  redoxColor <- c("Oxy"="#009E73","Sub"="#E69F00","Eux"="#8A2BE2")
  cols <- c("Sulfide" = "black", "Dissolved Oxygen" = "red", "Manganese" = "purple", "NOx" = "cyan")
  lineT <- c("Sulfide" = 1, "Dissolved Oxygen" = 2, "Manganese" = 3, "NOx" = 4)
  
  if (legendGen == TRUE) {
    #create legend 
    figlegend <- ggplot(desired) +
      geom_line(aes(y = desired$avgSulfide_uM, x = desired$depth, color = "Sulfide")) +
      geom_point(aes(y = desired$avgSulfide_uM, x = desired$depth, shape = "Sulfide", color = "Sulfide"), size = 1.5) +
      geom_rect(mapping = aes(xmin = oxyStart, xmax = (oxyEnd+subStart)/2, ymin = -Inf, ymax = Inf, fill = "Oxy"), alpha = 0.25) +
      geom_rect(mapping = aes(xmin = (oxyEnd+subStart)/2, xmax = (subEnd+euxStart)/2, ymin = -Inf, ymax = Inf, fill = "Sub"), alpha =0.25) +
      geom_rect(mapping = aes(xmin = (subEnd+euxStart)/2, xmax = euxEnd, ymin = -Inf, ymax = Inf, fill = "Eux"), alpha = 0.25) +
      ##the 10 multipled is for plot scaling
      geom_line(aes(y = desired$avgDO*20, x = desired$depth, color = "Dissolved Oxygen")) +
      geom_point(aes(y = desired$avgDO*20, x = desired$depth, shape = "Dissolved Oxygen", color = "Dissolved Oxygen"), size = 1.5) +
      geom_line(aes(y = desired$avguMMn*20, x = desired$depth, color = "Manganese")) +
      geom_point(aes(y = desired$avguMMn*20, x = desired$depth, shape = "Manganese", color = "Manganese"), size = 1.5) +
      geom_line(aes(y = desired$avgNOx, x = desired$depth, color = "NOx")) +
      geom_point(aes(y = desired$avgNOx, x = desired$depth, shape = "NOx", color = "NOx"), size = 1.5) +
      #remove geom_rects if not using redoxState
      theme(axis.text.y = element_text(colour = "black", face = "bold", size = 10),
            axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
            legend.text = element_text(size = 6, face ="bold", colour ="black"), 
            legend.position = "right", axis.title.y = element_text(face = "bold", size = 8), 
            axis.title.x = element_text(colour = "black", face = "bold", size = 8), 
            axis.line = element_line(colour = "black", size = 0.5),
            legend.title = element_blank(),
            panel.background = element_blank(),
            legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
            plot.title = element_text(colour = "black", size = 14, face = "bold", hjust = 0.5)) +
      coord_flip(xlim = c(24,0), ylim = c(0,200)) +
      scale_x_continuous(trans = "reverse") +
      scale_color_manual(name = "Water Chemistry", values = cols) +
      scale_shape_manual(name = "Water Chemistry", values = lineT) +
      #scale y axis to same number multiplied by desired$avgDO
      scale_y_continuous(sec.axis = sec_axis(~./20,name = "Dissolved Oxygen mg/L, µM Manganese, ppb NOx")) +
      ggtitle(paste(desiredMonth, "/", desiredDay, "/", desiredYear, sep = "")) +
      scale_fill_manual(name = "Redox State", values = redoxColor) +
      labs(y="µM Sulfide", x = "Depth (m)")
    return(figlegend)
  }
  if (min(subData$depth) == Inf) {
    results <- ggplot(desired) +
      geom_line(aes(y = desired$avgSulfide_uM, x = desired$depth, color = "Sulfide")) +
      geom_point(aes(y = desired$avgSulfide_uM, x = desired$depth, shape = "Sulfide", color = "Sulfide"), size = 1.5) +
      geom_rect(mapping = aes(xmin = oxyStart, xmax = (oxyEnd+euxStart)/2, ymin = -Inf, ymax = Inf, fill = "Oxy"), alpha = 0.01) +
      geom_rect(mapping = aes(xmin = (oxyEnd+euxStart)/2, xmax = euxEnd, ymin = -Inf, ymax = Inf, fill = "Eux"), alpha = 0.01) +
      ##the 10 multipled is for plot scaling
      geom_line(aes(y = desired$avgDO*20, x = desired$depth, color = "Dissolved Oxygen")) +
      geom_point(aes(y = desired$avgDO*20, x = desired$depth, shape = "Dissolved Oxygen", color = "Dissolved Oxygen"), size = 1.5) +
      geom_line(aes(y = desired$avguMMn*20, x = desired$depth, color = "Manganese")) +
      geom_point(aes(y = desired$avguMMn*20, x = desired$depth, shape = "Manganese", color = "Manganese"), size = 1.5) +
      geom_line(aes(y = desired$avgNOx, x = desired$depth, color = "NOx")) +
      geom_point(aes(y = desired$avgNOx, x = desired$depth, shape = "NOx", color = "NOx"), size = 1.5) +
      #remove geom_rects if not using redoxState
      theme(axis.text.y = element_text(colour = "black", face = "bold", size = 10),
            axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
            legend.text = element_text(size = 6, face ="bold", colour ="black"), 
            legend.position = "none", axis.title.y = element_text(face = "bold", size = 8), 
            axis.title.x = element_text(colour = "black", face = "bold", size = 8), 
            axis.line = element_line(colour = "black", size = 0.5),
            legend.title = element_blank(),
            panel.background = element_blank(),
            legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
            plot.title = element_text(colour = "black", size = 14, face = "bold", hjust = 0.5)) +
      coord_flip(xlim = c(24,0), ylim = c(0,200)) +
      scale_x_continuous(trans = "reverse") +
      scale_color_manual(name = "Water Chemistry", values = cols) +
      scale_shape_manual(name = "Water Chemistry", values = lineT) +
      #scale y axis to same number multiplied by desired$avgDO
      scale_y_continuous(sec.axis = sec_axis(~./20,name = "Dissolved Oxygen mg/L, µM Manganese, ppb NOx")) +
      ggtitle(paste(desiredMonth, "/", desiredDay, "/", desiredYear, sep = "")) +
      scale_fill_manual(name = "Redox State", values = redoxColor) +
      labs(y="µM Sulfide", x = "Depth (m)")
  }
  else {
    results <- ggplot(desired) +
      geom_line(aes(y = desired$avgSulfide_uM, x = desired$depth, color = "Sulfide")) +
      geom_point(aes(y = desired$avgSulfide_uM, x = desired$depth, shape = "Sulfide", color = "Sulfide"), size = 1.5) +
      geom_rect(mapping = aes(xmin = oxyStart, xmax = (oxyEnd+subStart)/2, ymin = -Inf, ymax = Inf, fill = "Oxy"), alpha = 0.01) +
      geom_rect(mapping = aes(xmin = (oxyEnd+subStart)/2, xmax = (subEnd+euxStart)/2, ymin = -Inf, ymax = Inf, fill = "Sub"), alpha = 0.01) +
      geom_rect(mapping = aes(xmin = (subEnd+euxStart)/2, xmax = euxEnd, ymin = -Inf, ymax = Inf, fill = "Eux"), alpha = 0.01) +
      ##the 10 multipled is for plot scaling
      geom_line(aes(y = desired$avgDO*20, x = desired$depth, color = "Dissolved Oxygen")) +
      geom_point(aes(y = desired$avgDO*20, x = desired$depth, shape = "Dissolved Oxygen", color = "Dissolved Oxygen"), size = 1.5) +
      geom_line(aes(y = desired$avguMMn*20, x = desired$depth, color = "Manganese")) +
      geom_point(aes(y = desired$avguMMn*20, x = desired$depth, shape = "Manganese", color = "Manganese"), size = 1.5) +
      geom_line(aes(y = desired$avgNOx, x = desired$depth, color = "NOx")) +
      geom_point(aes(y = desired$avgNOx, x = desired$depth, shape = "NOx", color = "NOx"), size = 1.5) +
      #remove geom_rects if not using redoxState
      theme(axis.text.y = element_text(colour = "black", face = "bold", size = 10),
            axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
            legend.text = element_text(size = 6, face ="bold", colour ="black"), 
            legend.position = "none", axis.title.y = element_text(face = "bold", size = 8), 
            axis.title.x = element_text(colour = "black", face = "bold", size = 8), 
            axis.line = element_line(colour = "black", size = 0.5),
            legend.title = element_blank(),
            panel.background = element_blank(),
            legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
            plot.title = element_text(colour = "black", size = 14, face = "bold", hjust = 0.5)) +
      coord_flip(xlim = c(24,0), ylim = c(0,200)) +
      scale_x_continuous(trans = "reverse") +
      scale_color_manual(name = "Water Chemistry", values = cols) +
      scale_shape_manual(name = "Water Chemistry", values = lineT) +
      #scale y axis to same number multiplied by desired$avgDO
      scale_y_continuous(sec.axis = sec_axis(~./20,name = "Dissolved Oxygen mg/L, µM Manganese, ppb NOx")) +
      ggtitle(paste(desiredMonth, "/", desiredDay, "/", desiredYear, sep = "")) +
      scale_fill_manual(name = "Redox State", values = redoxColor) +
      labs(y="µM Sulfide", x = "Depth (m)")
  }
  return(results)
}

#generate all graphs
g816 <- waterChemistry(desiredMonth = 8, desiredDay = 9, desiredYear = 2016, legendGen = FALSE)
g916 <- waterChemistry(desiredMonth = 9, desiredDay = 22, desiredYear = 2016, legendGen = FALSE)
g1016 <- waterChemistry(desiredMonth = 10, desiredDay = 4, desiredYear = 2016, legendGen = FALSE)
g817 <- waterChemistry(desiredMonth = 8, desiredDay = 11, desiredYear = 2017, legendGen = FALSE)
g917 <- waterChemistry(desiredMonth = 9, desiredDay = 8, desiredYear = 2017, legendGen = FALSE)
g1017 <- waterChemistry(desiredMonth = 10, desiredDay = 4, desiredYear = 2017, legendGen = FALSE)
figlegend <- waterChemistry(desiredMonth = 8, desiredDay = 9, desiredYear = 2016, legendGen = TRUE)

noLegend <- grid.arrange(g816, g916, g1016, g817, g917, g1017, ncol=3, nrow = 2)
grid.arrange(as_ggplot(get_legend(figlegend)), noLegend, ncol = 2, nrow = 1, widths = c(1,8))
