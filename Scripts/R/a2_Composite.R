#Purpose-----
#goal: create regional composits based on ipcc regions using either temp or hc
#in:   LiPD files
#out:  individual csv for each region where columns repersent iterations. 
#      summary table with each column as the regional median timeseries
#---
library(compositeR)
library(doParallel)
library(dplyr)
library(foreach)
library(geoChronR)
library(lipdR)
library(magrittr)
library(purrr)
library(tidyverse)


#Set up directories and names
dataDir <- getwd()
climVar <- 'T'
#

###Load Data
lipdData <- readRDS(file.path(dataDir,'Data','LiPD','lipdData.rds'))
lipdTSO  <- lipdData[[climVar]]
if(climVar == 'T'){
  lipdTSO <- filterTs(lipdTSO,'paleoData_units == degC')
  lipdTSO <- filterTs(lipdTSO,'paleoData_datum == abs')
}
regNames <- sort(unique(as.character(pullTsVariable(lipdTSO,'geo_ipccRegion'))))

save=TRUE
set.seed(5) #make reproducible#

#Set variables for composite code
nens          <- 1000    #lower = faster
binsize       <- 100   #years (median resolution = 107yrs)
ageMin        <- -100  #age BP
ageMax        <- 12400 #age BP
searchDur     <- 3500  #yrs (for 3 lake deposit data points)
nThresh       <- 6     #minimum # of records, else skip #arbitrary
samplePct     <- 0.75  #arbitrary

#Set bin vectors
binvec   <- seq(ageMin-binsize/2, to = ageMax+binsize/2, by = binsize)
binYears <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))

#Set up data to add once regional composite is calculated
compositeEns      <- vector(mode="list")
medianCompositeTS <- data_frame(time=binYears)

#Loop to composite (by region)
for (reg in regNames) {
  #Filter the TS by cluster name and make sure have enough values
  lipdReg  <- filterTs(lipdTSO,paste('geo_ipccRegion ==',reg))
  regCount <- length(lipdReg)
  #Skip if number of records is too few
  if(regCount < nThresh) next
  for (i in 1:regCount){
    if (lipdReg[[i]]$climateInterpretation1_interpDirection == 'negative'){
      lipdReg[[i]]$paleoData_values <- lipdReg[[i]]$paleoData_values*-1
    }
  }
  #
  #Calculate data density throughout the Holocene 
  timeN <- plotTimeAvailabilityTs(lipdReg,
                                  age.range=c(0,12000),
                                  group.var ='Category',
                                  step=binsize)$dat %>%
    group_by(yvec) %>% 
    summarise(count=sum(value),countPct=sum(value)/regCount)
  #Determine search range based on 0-10ka with >50% proxy coverage
  MH <- which(timeN$yvec==6000)
  if (length(which(timeN$countPct<=0.5) > 0)){
    idx <- c(max(which(timeN$countPct[1:MH]<=0.5),
                 which(timeN$yvec==0)),
             min(which(timeN$countPct[MH:length(timeN$countPct)]<=0.5)+MH,
                 which(timeN$yvec==12000)))
  } else{
    idx <- c(which(timeN$yvec==0),which(timeN$yvec==12000))
  }
  searchMin  <- timeN$yvec[idx[1]] #age BP
  searchMax  <- timeN$yvec[idx[2]] #age BP
  #
  #setup and run ensemble ##### This is the main part of the code to edit ##### 
  compEns <- matrix(NA,nrow = length(binYears),ncol=nens)
  for (i in 1:nens){
    #Sample from full region
    wght <- pullTsVariable(lipdReg,'ageRange')/pullTsVariable(lipdReg,'ageResPlus')
    wght[which(is.na(wght))] <- mean(wght,na.rm=TRUE)
    lipdRegSample <- sample(x    = pullTsVariable(lipdReg,'paleoData_TSid'),
                            #prob = wght / sum(wght),
                            size = length(lipdReg) * samplePct)
    #sample <- which(pullTsVariable(lipdReg,'paleoData_TSid') %in% pullTsVariable(lipdReg,'paleoData_TSid'))
    sample <- which(pullTsVariable(lipdReg,'paleoData_TSid') %in% lipdRegSample)
    lipdRegSample <- lipdReg#[sample]
    #Composite
    if (climVar == 'HC'){
      std <- TRUE
    } else{std <- FALSE}
    tc <- compositeR::compositeEnsembles(fTS      = lipdRegSample,
                                         binvec   = binvec,
                                         stanFun  = standardizeMeanIteratively,
                                         binFun   = simpleBinTs,
                                         ageVar   = "age",
                                         alignInterpDirection = FALSE,
                                         spread   = TRUE,
                                         duration = searchDur,
                                         searchRange = c(searchMin,searchMax),
                                         normalizeVariance    = std,
                                         scope    = "climate",
                                         minN     = 3) 
    compEns[,i] <- tc$composite
  }
  #
  # Return reconstruction and additional data for plotting
  compositeEns[[reg]]      <- compEns[which(binYears==0):which(binYears==12000),]
  medianCompositeTS[[reg]] <- apply(compEns,1,median,na.rm=TRUE)
}



if(save){
  write.csv(medianCompositeTS, row.names = FALSE,
            file = file.path(dataDir,'Data','RegionComposites',climVar,'MedianTSbyRegion.csv'))
  for (region in names(compositeEns)){
    write.csv(compositeEns[[region]], row.names = FALSE,
              file = file.path(dataDir,'Data','RegionComposites',climVar,paste(region,'.csv',sep='')))
  }
  print("csv files saved")
}

