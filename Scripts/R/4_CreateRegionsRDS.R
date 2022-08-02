#Purpose--------------------------------------------------------------------------------
#goal: Create image settings for all figures and compile regional data
#in:   lipd RDS; ipcc region URL, 
#regionData:  regional RDS

#Load Packages--------------------------------------------------------------------------------

library(dplyr)
library(geoChronR)
library(lipdR)
library(maptools)
library(proj4)
library(sf)
library(sp)
library(tidyverse)

#Set up directories and load data--------------------------------------------------------------------------------

dir  <- getwd()# '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate' #

PROJ <- '+proj=robin   +ellps=WGS84 +datum=WGS84 +no_defs +lon_0=0 +x_0=0 +y_0=0 +units=m'
PROJorig <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
load(url('https://github.com/SantanderMetGroup/ATLAS/blob/main/reference-regions/IPCC-WGI-reference-regions-v4_R.rda?raw=true'))
refregions <-  spTransform(IPCC_WGI_reference_regions_v4, CRSobj = PROJ)

lipdData  <- readRDS(file.path(dir,'Data','Proxy','LiPD','lipdData.rds'))
#modelData <- readRDS(file.path(dir,'Data','Model','Region','lipdData.rds')) 

#Set up directories and names--------------------------------------------------------------------------------

regionData <- vector(mode='list')
regionSummary <- data.frame(region=names(regionData),
                            name=rep(NA,length(regionData)),
                            type=rep(NA,length(regionData)),
                            latitude=rep(NA,length(regionData)),
                            longitude=rep(NA,length(regionData)),
                            nproxy_HC=rep(NA,length(regionData)),
                            nproxy_T=rep(NA,length(regionData)),
                            composite_HC=rep(NA,length(regionData)),
                            composite_T=rep(NA,length(regionData))
                            )

for (reg in sort(as.vector(refregions$Acronym))){
  regionData[[reg]]$polygon <- refregions[refregions@data$Acronym == reg, ]
  regionData[[reg]]$name <- as.character(regionData[[reg]]$polygon$Name)
  regionData[[reg]]$type <- as.character(regionData[[reg]]$polygon$Type)
  regionData[[reg]]$latitude  <- regionData[[reg]]$polygon@polygons[[1]]@labpt[[1]]
  regionData[[reg]]$longitude <- regionData[[reg]]$polygon@polygons[[1]]@labpt[[2]]
  regionData[[reg]]$xadjust   <- NA
  regionData[[reg]]$yadjust   <- NA
  for (var in c('HC')){
    regionData[[reg]][[var]] <- vector(mode='list')
    regionData[[reg]][[var]]$LiPD   <- filterTs(lipdData[[var]],paste('geo_ipccRegion ==',reg))
    regionData[[reg]][[var]]$nproxy <- length(regionData[[reg]][[var]]$LiPD)
    if (regionData[[reg]][[var]]$nproxy > 0){
      regionData[[reg]][[var]]$SummaryDF <- read.csv(file.path(dir,'Data','Proxy',paste('proxyMetadata_',var,'.csv',sep=''))) %>% 
        filter(ipccReg==reg)
      pltTimeAvail <- plotTimeAvailabilityTs(regionData[[reg]][[var]]$LiPD,age.range=c(0,12000),group.var ='CategorySpecific',step=100)$dat %>% 
        group_by(yvec) %>% summarise(count=sum(value),nPct=round(100*sum(value)/regionData[[reg]][[var]]$nproxy,1))
      regionData[[reg]][[var]]$pltTimeAvail50range <- seq(max(0,which(pltTimeAvail$yvec<6000 & pltTimeAvail$nPct < 50)+1),
                                min(nrow(pltTimeAvail),which(pltTimeAvail$yvec>6000 & pltTimeAvail$nPct < 50)-1))
    } else{
      regionData[[reg]][[var]]$SummaryDF <- NA
    }
    tryCatch(
      expr = {
        regionData[[reg]][[var]]$compositeEnsemble <- read.csv(file.path(dir,'Data','RegionComposites',var,paste(reg,'.csv',sep='')))},
      error = function(e){
        regionData[[reg]][[var]]$compositeEnsemble <- NA}
    )
    #Model Data?
  }
}

#Nudge values for global figure--------------------------------------------------------------------------------

#Americas
regionData[['GIC']]$xadjust  <-  0
regionData[['NWN']]$xadjust  <- -0.002
regionData[['NWN']]$yadjust  <- -0.02
regionData[['NEN']]$yadjust  <- -0.028
regionData[['NEN']]$xadjust  <- -0.002
regionData[['WNA']]$xadjust  <- -0.053
regionData[['CNA']]$xadjust  <-  0.007
regionData[['ENA']]$xadjust  <-  0.060
regionData[['NAS']]$xadjust  <- -0.04
regionData[['SCA']]$xadjust  <- -0.015
regionData[['CAR']]$xadjust  <-  0.03
regionData[['CAR']]$yadjust  <-  0.045
regionData[['NWS']]$xadjust  <- -0.01
regionData[['NSA']]$xadjust  <-  0.025
regionData[['NSA']]$yadjust  <-  0.048
regionData[['NES']]$yadjust  <-  0.022
regionData[['SAM']]$yadjust  <- -0.03
#Eruope/Africa
regionData[['NEU']]$xadjust  <-  0.015  
regionData[['NEU']]$yadjust  <-  0.01  
regionData[['WCE']]$xadjust  <- -0.005  
regionData[['MED']]$xadjust  <- -0.02   
regionData[['CAF']]$xadjust  <- -0.054   
regionData[['NEAF']]$yadjust <-  0.002   
regionData[['SEAF']]$yadjust <- -0.002   
regionData[['WSAF']]$xadjust <- -0.032  
regionData[['ESAF']]$xadjust <-  0.03   
#Asia/Australasia  
regionData[['WSB']]$xadjust <-  -0.04
regionData[['ESB']]$yadjust <-  0.014  
regionData[['ESB']]$yadjust <-  0.014  
regionData[['ESB']]$xadjust <- -0.01   
regionData[['RFE']]$xadjust <- 0.04
regionData[['RFE']]$yadjust <- -0.02   
regionData[['WCA']]$xadjust <- -0.021  
regionData[['ECA']]$xadjust <- -0.0035  
regionData[['ECA']]$yadjust <-  0.015   
regionData[['TIB']]$xadjust <-  0.003   
regionData[['TIB']]$yadjust <- -0.007   
regionData[['EAS']]$xadjust <-  0.021   
regionData[['SAS']]$yadjust <- -0.02  
regionData[['SAS']]$xadjust <-  0.0  
regionData[['SEA']]$xadjust <- -0.017  
regionData[['SAU']]$xadjust <- -0.044
regionData[['NZ']]$xadjust  <- -0.017
#Ocean  
regionData[['EPO']]$xadjust <- -0.02   
regionData[['ARO']]$xadjust <- -0.18  
regionData[['ARO']]$yadjust <- -0.03  
regionData[['NAO']]$yadjust <-  0.106
regionData[['NAO']]$xadjust <-  0.03  


#Save--------------------------------------------------------------------------------

saveRDS(regionData,file.path(dir,'Data','FigureSettings','regionData.rds'))








