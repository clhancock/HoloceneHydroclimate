#---
#goal: create regional composits based on ipcc regions. works with either temp or hc
#input: LiPD files, csv file from python script to filter TSids with
#output: individual csv for each region where columns repersent iterations. Also summary with each column as the regional median timeseries
#---
#Changed minN from 8 to 3 in standardize function to accomodate lake deposits
#   this is the number of measurments neaded to be located withing the search range
library(geoChronR)
library(lipdR)
#library(scales)
#library(cowplot)
library(purrr)
library(dplyr)
library(magrittr)
#library(ggplot2)
library(compositeR)
library(foreach)
library(doParallel)
library(tidyverse)
#library(abind)
#library(ncdf4)

#Set up directories and names
githubDir <- getwd()
climVar <- 'HC'
#

###Load Data
lipdData <- readRDS(file.path(githubDir,'Data','LiPD','lipdData.rds'))
lipdTSO <- lipdData[[climVar]]
regionNames <- sort(unique(as.character(pullTsVariable(lipdTSO,'geo_ipccRegion'))))
climVar <- 'HC'

save=FALSE
#set.seed(#) Set same sets of records which will make completely reproducable
#
#Set variables for composite code
nens          <- 5000  #make low to run quickly, set high to get large ensemble range (variation from standardization search range and order )
binsize       <- 100 #years
ageMin        <- -100 #age BP
ageMax        <- 12400 #age BP
searchDuration<- 3500 #yrs
minNrecords   <- 6 #num of records
samplePct <- 0.75
#
binvec   <- seq(ageMin-binsize/2, to = ageMax+binsize/2, by = binsize)
binYears <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))

#Set up data to add once
compositeEns      <- vector(mode="list")
medianCompositeTS <- data_frame(time=binYears)
#Do the compositing (by region)
for (region in regionNames) {
  #Filter the TS by cluster name and make sure have enough values
  lipdRegion <- filterTs(lipdTSO,paste('geo_ipccRegion ==',region))
  #Skip if number of records is too few
  if(length(lipdRegion)<minNrecords|sum(pullTsVariable(lipdRegion,'archiveType')!="LakeDeposits")<=2)next
  #
  timeAvail <- plotTimeAvailabilityTs(lipdRegion,age.range=c(0,12000),
                                      group.var ='Category',step=binsize)$dat %>%
    group_by(yvec) %>% 
    summarise(count=sum(value),countPct=sum(value)/length(lipdRegion))
  idx <- which(timeAvail$countPct<=0.5)
  idx_MH <- which(timeAvail$yvec==6000)
  if (length(idx > 0)){
    idx <- c(max(which(timeAvail$countPct[1:idx_MH]<=0.5),
                 which(timeAvail$yvec==0)),
             min(which(timeAvail$countPct[idx_MH:length(timeAvail$countPct)]<=0.5)+idx_MH,
                 which(timeAvail$yvec==10000)))
  } else{idx <- c(which(timeAvail$yvec==0),which(timeAvail$yvec==10000))}
  searchAgeMin  <- timeAvail$yvec[idx[1]] #age BP
  searchAgeMax  <- timeAvail$yvec[idx[2]] #age BP
  #setup and run ensemble ########## This is the main part of the code to edit ##########
  compEns <- matrix(NA,nrow = length(binYears),ncol=nens)
  for (i in 1:nens){
    weights <- pullTsVariable(lipdRegion,'ageRange') / pullTsVariable(lipdRegion,'ageResPlus')
    weights[which(is.na(weights))] <- mean(weights,na.rm=TRUE)
    lipdRegionSample <- sample(x    = pullTsVariable(lipdRegion,'paleoData_TSid'),
                               size = length(lipdRegion)*samplePct,
                               prob = weights / sum(weights))
    lipdRegionSample <- lipdRegion[which(pullTsVariable(lipdRegion,'paleoData_TSid') %in% lipdRegionSample)]
    tc <- compositeR::compositeEnsembles(lipdRegionSample,
                                         binvec,
                                         stanFun = standardizeMeanIteratively,
                                         minN = 3,
                                         ageVar  = "age",
                                         alignInterpDirection = TRUE,
                                         spread      = TRUE,
                                         duration    = searchDuration,
                                         searchRange = c(searchAgeMin,searchAgeMax),
                                         normalizeVariance = TRUE,
                                         scope = "climate",
                                         binFun = simpleBinTs) #sampleEnsembleThenBinTs
    compEns[,i] <- tc$composite
  }
  #
  # Return reconstruction and additional data for plotting
  compositeEns[[region]]      <- compEns
  medianCompositeTS[[region]] <- apply(compEns,1,median,na.rm=TRUE)
}

if(save){
  write.csv(medianCompositeTS, row.names = FALSE,
            file = file.path(githubDir,'Data','RegionComposites',climVar,'MedianTSbyRegion.csv'))
  for (region in names(compositeEns)){
    write.csv(compositeEns[[region]], row.names = FALSE,
              file = file.path(githubDir,'Data','RegionComposites',climVar,paste(region,'.csv',sep='')))
  }
  print("csv files saved")
}




