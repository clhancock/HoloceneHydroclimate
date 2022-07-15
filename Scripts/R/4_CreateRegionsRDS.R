#Purpose--------------------------------------------------------------------------------
#goal: compile regional data
#in:   
#out:  rdata file list of regions including key figure and reconstruction data

#Load Packages--------------------------------------------------------------------------------

library(compositeR)
library(doParallel)
library(dplyr)
library(foreach)
library(geoChronR)
library(lipdR)
library(magrittr)
library(purrr)
library(tidyverse)


#Set up directories and names--------------------------------------------------------------------------------

dir  <- getwd()# '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate' #
var  <- 'HC'
save <- FALSE
saveDir <- file.path(dir,'Data','RegionComposites',var)

region
polygon shapefile
xy figure adjustment
temp reconstrution
hc reconstruction
LiPD_proxies 
plot timeseries (50%)
