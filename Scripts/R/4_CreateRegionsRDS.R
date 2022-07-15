#Purpose--------------------------------------------------------------------------------
#goal: compile regional data
#in:   
#out:  rdata file list of regions including key figure and reconstruction data

#Load Packages--------------------------------------------------------------------------------

library(dplyr)
library(geoChronR)
library(lipdR)
library(maptools)
library(proj4)
library(sf)
library(sp)
library(tidyverse)

#Set up directories and names--------------------------------------------------------------------------------

dir  <- getwd()# '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate' #
var  <- 'HC'
save <- FALSE
saveDir <- file.path(dir,'Data','RegionComposites',var)


PROJ <- '+proj=robin   +ellps=WGS84 +datum=WGS84 +no_defs +lon_0=0 +x_0=0 +y_0=0 +units=m'
PROJorig <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
load(url('https://github.com/SantanderMetGroup/ATLAS/blob/main/reference-regions/IPCC-WGI-reference-regions-v4_R.rda?raw=true'))
refregions <-  spTransform(IPCC_WGI_reference_regions_v4, CRSobj = PROJ)

lipdData <- readRDS(file.path(dir,'Data','Proxy','LiPD','lipdData.rds'))

out <- vector(mode='list')


reg <- 'EAS'
out[[reg]] <- vector(mode='list')

out[[reg]]$polygon <- refregions[refregions@data$Acronym == reg, ]
out[[reg]]$name <- as.character(out[[reg]]$polygon$Name)
out[[reg]]$type <- as.character(out[[reg]]$polygon$Type)
out[[reg]]$latitude  <- out[[reg]]$polygon@polygons[[1]]@labpt[[1]]
out[[reg]]$longitude <- out[[reg]]$polygon@polygons[[1]]@labpt[[2]]
out[[reg]]$xadjust   <- NA
out[[reg]]$yadjust   <- NA
for (var in c('T','HC')){
  out[[reg]][[var]] <- vector(mode='list')
  out[[reg]][[var]]$LiPD <- filterTs(lipdData[[var]],paste('geo_ipccRegion ==',reg))
  out[[reg]][[var]]$nproxy <- length(out[[reg]][[var]]$LiPD)
  out[[reg]][[var]]$compositeEnsemble <- read.csv(file.path(dir,'Data','RegionComposites',var,paste(reg,'.csv',sep='')))
  pltTimeAvail <- plotTimeAvailabilityTs(out[[reg]][[var]]$LiPD,age.range=c(0,12000),group.var ='CategorySpecific',step=200)$dat %>% 
    group_by(yvec) %>% summarise(count=sum(value),nPct=round(100*sum(value)/out[[reg]][[var]]$nproxy,1))
  pltTimeAvail50range <- c(max(0,which(pltTimeAvail$yvec<6000 & pltTimeAvail$nPct < 50)+1),
          min(nrow(pltTimeAvail),which(pltTimeAvail$yvec>6000 & pltTimeAvail$nPct < 50)-1))
}
