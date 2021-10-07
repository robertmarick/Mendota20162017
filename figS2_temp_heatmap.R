

#### Get set up ####

rm(list = ls())
setwd("~/Box/Mendota_16S/")
library(tidyverse)
library(gtools)
library(LakeMetabolizer)
library(rLakeAnalyzer)
library(stringr)
library(RcppRoll)
library(lubridate)
library(viridis)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")


#### Set variables of interest for testing ####

data.file.of.interest <- "dataRaw/2016_2017_sensor_mendota_lake_watertemp_hourly.csv"
year.of.interest <- 2017
starting.date <- "2017-05-01"
ending.date <- "2017-11-05"


#### Define function to plot heatmap of water temp with thermocline laid over top ####

temp.profile.thermocline <- function(data.file.of.interest,
                                     year.of.interest,
                                     starting.date,
                                     ending.date,
                                     temp.scale = c(4, 30)) {
  
  if (!exists("data.file.of.interest")) {
    print("Gotta supply a data file.")
    break()
  }
  
  if (!exists("year.of.interest")) {
    print("Gotta set year of interest.")
    break()
  }
  
  # Read in the data and clean the timestamps
  LTERdf <- read.csv(file = data.file.of.interest,
                     stringsAsFactors = FALSE) %>%
    filter(year4 == year.of.interest) %>%
    mutate(hour = str_pad(hour, 4, pad = "0")) %>%
    mutate(date_time = ymd_hm(paste(sampledate, hour))) %>%
    select(date_time, depth, wtemp)
  
  if (exists("starting.date")) {
    print(paste("Cutting off dates before ",
                starting.date,
                ".",
                sep = ""
                ))
    LTERdf <- LTERdf %>%
      filter(as.Date(date_time) >= as.Date(starting.date))
  }
  if (exists("ending.date")) {
    print(paste("Cutting off dates after ",
                ending.date,
                ".",
                sep = ""
    ))
    LTERdf <- LTERdf %>%
      filter(as.Date(date_time) <= ending.date)
  }
  
  
  
  # Set up depth and date vectors
  depths <- unique(LTERdf$depth)
  dates <- unique(LTERdf$date_time)
  depths.to.plot <- seq(20, 0, by = -5)
  dates.to.plot <- seq(ceiling_date(min(dates), "weeks"),
                       floor_date(max(dates), "weeks"),
                       by='months')
  
  
  LTERdf <- LTERdf %>%
    spread(key = depth,
           value = wtemp,
           fill = NA)

  
  # Find dates with NAs and exclude them
  # LTERdf.missing.vector <- LTERdf %>%
  #   select(-date_time) %>%
  #   rowSums() %>%
  #   is.na()
  # if (length(which(LTERdf.missing.vector)) > 0) {
  #   print(paste("NAs present on ",
  #               LTERdf[which(LTERdf.missing.vector), "date_time"],
  #               sep = ""))
  #   LTERdf <- LTERdf[!LTERdf.missing.vector, ]
  # }
  # 
  
  # Set up matrix
  LTERmt <- LTERdf %>%
    select(-date_time) %>%
    as.matrix()
  # LTERmt_top<-LTERmt[,1:15]
  
  
  
  #### Calculate thermocline depth ####
  names(LTERdf)[-1] <- paste("var_", names(LTERdf)[-1], sep = "")
  thermocline.depth.ts <- LTERdf %>%
    rename(datetime = date_time) %>%
    ts.thermo.depth() %>%
    group_by(as.Date(datetime)) %>%
    summarize(thermo.depth = median(thermo.depth,
                                    na.rm = TRUE)) %>%
    rename(datetime = `as.Date(datetime)`)
  
  

  #### Generate plot ####

  
  par(mar = c(3, 3.5, 2, 1),
      mgp=c(1.5,0.4,0),
      tck=-0.01)
  filled.contour(x = dates,
                 y = depths,
                 z = LTERmt,
                 ylim = c(max(depths), 0),
                 zlim = temp.scale,
                 nlevels = 40,
                 color.palette = colorRampPalette(c("violet", "blue", "cyan", "green3", "yellow", "orange", "red"),
                                                  bias = 1, space = "rgb"),
                 plot.title = { title(main = paste("Temperature profile and thermocline depth for ",
                                                   year.of.interest,
                                                   sep = "")) },
                 plot.axes = { axis(1,
                                    at = dates.to.plot,
                                    labels = month(dates.to.plot, label = TRUE));
                   axis(2,
                        at = depths.to.plot);
                   lines(as.POSIXct(thermocline.depth.ts$datetime),
                         thermocline.depth.ts$thermo.depth,
                         lwd = 1,
                         col = 'black') }
                 )
  # lines(x = as.POSIXct(thermocline.depth.ts$datetime),
  #       y = thermocline.depth.ts$thermo.depth)
  
  title(ylab = "Depth (m)")

}


#### Plot out heat maps ####

pdf("results/buoy_data/2016_heatmap_thermocline.pdf",
    width = 7.5,
    height = 5)
temp.profile.thermocline(data.file.of.interest = "dataRaw/2016_2017_sensor_mendota_lake_watertemp_hourly.csv",
                         year.of.interest = 2016,
                         starting.date = "2016-05-10",
                         ending.date = "2016-11-05")
dev.off()


pdf("results/buoy_data/2017_heatmap_thermocline.pdf",
    width = 7.5,
    height = 5)
temp.profile.thermocline(data.file.of.interest = "dataRaw/2016_2017_sensor_mendota_lake_watertemp_hourly.csv",
                         year.of.interest = 2017,
                         starting.date = "2017-05-10",
                         ending.date = "2017-11-05")
dev.off()
