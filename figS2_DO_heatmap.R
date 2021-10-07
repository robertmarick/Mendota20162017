#### code/waterChem/DO_heatmap.R ####
# Benjamin D. Peterson


#### Set it up ####

rm(list = ls())
setwd("~/Box/Mendota_16S/")
library(akima)
# library(colorspace)
# library(ggplot2)
library(fields)
library(ggimage)
library(lubridate)
library(tidyverse)
library(viridis)


#### Read in data ####

sonde.data <- read.csv("dataEdited/waterChem/sonde_2016_2017.csv",
                       stringsAsFactors = FALSE)


#### Dates ####

date.vector <- c("May", "Jun", "Jul", "Aug",
                 "Sep", "Oct", "Nov")
date.naming.vector <- yday(mdy(paste(date.vector,
                                     ", 2017",
                                     sep = "")))
names(date.naming.vector) <- date.vector


#### Make common limits for both years ####
unique.yday <- yday(sonde.data$sample.date) %>%
  unique() %>%
  sort()
# Start June 1st
#start.date <- unique.yday[1]
start.date <- yday("2016-05-10")
end.date <- yday("2016-11-05")

#### Define function ####

DO.heatmap <- function(year) {
  sonde.data.year <- sonde.data %>%
    filter(Year == year) %>%
    mutate(sampleDate = yday(as.Date(sample.date))) %>%
    filter(!is.na(DO.mg.L))
  
  heatmap.data <- sonde.data.year %>%
    select(sampleDate, depth.m, DO.mg.L) %>%
    spread(key = depth.m,
           value = DO.mg.L) %>%
    gather(key = depth,
           value = DO.mg.L,
           -1) %>%
    arrange(sampleDate)
  
  heatmap.data.interp <- interp(x = heatmap.data$sampleDate,
                                y = -as.numeric(heatmap.data$depth),
                                z = heatmap.data$DO.mg.L)
  
  par(mar = c(3, 3.5, 2, 1),
      mgp=c(1.5,0.4,0),
      tck=-0.01)
  image.plot(heatmap.data.interp,
             axes = F,
             col = viridis(20),
             xlim = c(start.date,
                      end.date),
             zlim = c(0,15))

  axis(side = 2,
       at = seq(-20, 0,
                by = 5),
       labels = seq(20, 0,
                    by = -5))
  
  axis(side = 1,
       at = date.naming.vector,
       labels = names(date.naming.vector),
       las = 1)
  
  title(main = paste("DO profile for ",
                     year,
                     sep = ""),
        ylab = "Depth (m)")
  
  mtext("mg/L",
        side=4,
        line=3.3,
        cex=1.2)
  
  
}


#### Read out plots ####
pdf("results/waterChem/DO_heatmap_2016.pdf",
    width = 7.5,
    height = 5)
par(mfrow = c(1, 1))
DO.heatmap(2016)
dev.off()


pdf("results/waterChem/DO_heatmap_2017.pdf",
    width = 7.5,
    height = 5)
par(mfrow = c(1, 1))
DO.heatmap(2017)
dev.off()
