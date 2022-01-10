#---
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
githubDir <- getwd()
climVar <- 'HC'
#

###Load Data
lipdData <- readRDS(file.path(githubDir,'Data','LiPD','lipdData.rds'))
lipdTSO  <- lipdData[[climVar]]
regNames <- sort(unique(as.character(pullTsVariable(lipdTSO,'geo_ipccRegion'))))
climVar  <- 'HC'

save=TRUE
set.seed(5) #make reproducible
#
#Set variables for composite code
nens          <- 500    #lower = faster
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
  lipdReg <- filterTs(lipdTSO,paste('geo_ipccRegion ==',reg))
  regCount   <- length(lipdReg)
  #Skip if number of records is too few
  if(regCount < nThresh) next
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
                 which(timeN$yvec==10000)))
  } else{
    idx <- c(which(timeN$yvec==0),which(timeN$yvec==10000))
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
    lipdRegSample  <- pullTsVariable(lipdReg,'paleoData_TSid')
   # lipdRegSample <- sample(x    = pullTsVariable(lipdReg,'paleoData_TSid'),
     #                       size = length(lipdReg) * samplePct,
     #                       prob = wght / sum(wght))
    sample <- which(pullTsVariable(lipdReg,'paleoData_TSid') %in% lipdRegSample)
    lipdRegSample <- lipdReg[sample]
    #Composite
    tc <- compositeR::compositeEnsembles(fTS      = lipdRegSample,
                                         ageVar   = "age",
                                         scope    = "climate",
                                         spread   = TRUE,
                                         binvec   = binvec,
                                         binFun   = simpleBinTs,
                                         stanFun  = standardizeMeanIteratively,
                                         duration = searchDur,
                                         minN     = 3,
                                         searchRange = c(searchMin,searchMax),
                                         alignInterpDirection = TRUE,
                                         normalizeVariance    = TRUE) 
    compEns[,i] <- tc$composite
  }
  #
  # Return reconstruction and additional data for plotting
  compositeEns[[reg]]      <- compEns
  medianCompositeTS[[reg]] <- apply(compEns,1,median,na.rm=TRUE)
}



#Save
if(save){
  write.csv(medianCompositeTS, row.names = FALSE,
            file = file.path(githubDir,'Data','RegionComposites',climVar,'MedianTSbyRegion.csv'))
  for (region in names(compositeEns)){
    write.csv(compositeEns[[region]], row.names = FALSE,
              file = file.path(githubDir,'Data','RegionComposites',climVar,paste(region,'.csv',sep='')))
  }
  print("csv files saved")
}



####



reg <- 'GIC'
#Filter the TS by cluster name and make sure have enough values
lipdReg <- filterTs(lipdTSO,paste('geo_ipccRegion ==',reg))
regCount   <- length(lipdReg)
#Skip if number of records is too few
#if(regCount < nThresh) next
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
               which(timeN$yvec==10000)))
} else{
  idx <- c(which(timeN$yvec==0),which(timeN$yvec==10000))
}
searchMin  <- timeN$yvec[idx[1]] #age BP
searchMax  <- timeN$yvec[idx[2]] #age BP

compositeEns      <- vector(mode="list")
medianCompositeTS <- data_frame(time=binYears)
nens <- 100
for (n in c(3:10,regCount)){
  #
  #setup and run ensemble ##### This is the main part of the code to edit ##### 
  compEns <- matrix(NA,nrow = length(binYears),ncol=nens)
  #Sample from full region
  #wght <- pullTsVariable(lipdReg,'ageRange')/pullTsVariable(lipdReg,'ageResPlus')
  #wght[which(is.na(wght))] <- mean(wght,na.rm=TRUE)
  lipdRegSample <- sample(x    = pullTsVariable(lipdReg,'paleoData_TSid'),
                          size = n)
                   #       prob = wght / sum(wght))
  sample <- which(pullTsVariable(lipdReg,'paleoData_TSid') %in% lipdRegSample)
  lipdRegSample <- lipdReg[sample]
  for (i in 1:nens){
    #Composite
    tc <- compositeR::compositeEnsembles(fTS      = lipdRegSample,
                                         ageVar   = "age",
                                         scope    = "climate",
                                         spread   = TRUE,
                                         binvec   = binvec,
                                         binFun   = simpleBinTs,
                                         stanFun  = standardizeMeanIteratively,
                                         duration = searchDur,
                                         minN     = 3,
                                         searchRange = c(searchMin,searchMax),
                                         alignInterpDirection = TRUE,
                                         normalizeVariance    = TRUE) 
    compEns[,i] <- tc$composite
  }
  compositeEns[[n]] <- compEns
}
  #

for (n in c(3:10)){ 
  corout <- corEns(time.1   = binYears,
                   values.1 = compositeEns[[n]],
                   time.2   = binYears,
                   values.2 = compositeEns[[length(lipdReg)]],
                   bin.step = 100,
                   max.ens  = 100,
                   isopersistent  = TRUE,
                   isospectral    = TRUE,
                   gaussianize    = TRUE)
  print(n)
  print(mean(corout[["cor.stats"]][["r"]]))
}
  