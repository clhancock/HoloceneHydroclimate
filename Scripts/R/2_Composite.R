#Purpose--------------------------------------------------------------------------------
#goal: create regional composite based on IPCC regions using either temp or HC
#in:   LiPD files
#out:  individual csv for each region where columns represent iterations. 
#      summary table with each column as the regional median time series

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
save <- TRUE
saveDir <- file.path(dir,'Data','RegionComposites',var)


#Load Data without winter+ or summer+ seasonality--------------------------------------------------------------------------------

lipdData <- readRDS(file.path(dir,'Data','Proxy','LiPD','lipdData.rds'))[[var]]
lipdTSO  <- lipdData[-which(pullTsVariable(lipdData,"climateInterpretation1_seasonalityGeneral") %in% c('Summer+','Winter+'))]

if(var == 'T'){
  lipdTSO <- filterTs(lipdTSO,'paleoData_units == degC')
  lipdTSO <- filterTs(lipdTSO,'paleoData_datum == abs')
  std <- FALSE      #Use calibrated data so no need to normalize variance
} else{std <- TRUE} #Normalize variance because data recorded with different units


#Set variables for composite code--------------------------------------------------------------------------------

nens          <- 200     #Ensemble numbers (lower = faster)
binsize       <- 100     #years (median resolution = 107yrs)
ageMin        <- 0       #age BP
ageMax        <- 12400   #age BP
searchDur     <- 4000    #yrs (for 3 lake deposit data points)
nThresh       <- 6       #minimum no. of records, else skip 

#Set bin vectors
binvec   <- seq(ageMin-binsize/2, to = ageMax+binsize/2, by = binsize)
binYears <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))

#ID regions to reconstruct based on number of records (nThresh)
regNames <- data.frame(name=pullTsVariable(lipdTSO,'geo_ipccRegion')) %>% 
  group_by(name) %>% 
  summarise(n = n()) %>% 
  filter(n >= nThresh)
regNames <- c(as.character(regNames$name),'EAN','SSA') #Add 2 SH regions with fewer records to gain global coverage


#Calculate reconstructions--------------------------------------------------------------------------------

#Set up data to add once regional composite is calculated
compositeEnsemble <- vector(mode='list')
medianCompositeTS <- data_frame(time=binYears[1:which(binYears==12000)])

#Loop to composite (by region)
for (reg in c(regNames)) {
  lipdReg  <- filterTs(lipdTSO,paste('geo_ipccRegion ==',reg))
  for (i in 1:length(lipdReg)){
       if (lipdReg[[i]]$climateInterpretation1_interpDirection == 'negative'){
         lipdReg[[i]]$paleoData_values <- lipdReg[[i]]$paleoData_values*-1
       }
  }
  set.seed(5) #Reproducibility
  ensOut <- foreach(i = 1:nens) %dopar% {
    tc <- compositeEnsembles(fTS                  = lipdReg,
                             binvec               = binvec,
                             stanFun              = standardizeMeanIteratively,
                             binFun               = simpleBinTs,
                             ageVar               = "age",
                             alignInterpDirection = FALSE,
                             spread               = TRUE,
                             duration             = searchDur,
                             searchRange          = c(1000,9000),
                             normalizeVariance    = std,
                             minN                 = 3)
    return(list(composite = tc$composite,count = tc$count))
  }
  regionComposite           <- as.matrix(purrr::map_dfc(ensOut,magrittr::extract,"composite"))
  rownames(regionComposite) <- binYears
  regionComposite           <- regionComposite[1:which(binYears==12000),]
  compositeEnsemble[[reg]]  <- regionComposite
  medianCompositeTS[[reg]]  <- apply(regionComposite,1,median,na.rm=TRUE)
  #plot region to confirm that everything looks good
  plotTimeseriesEnsRibbons(X = binYears[1:which(binYears==12000)],Y = compositeEnsemble[[reg]])+
    scale_x_continuous(name = "age (yr BP)",         oob = scales::squish)+
    scale_y_continuous(name = "Standardized Anomaly",oob = scales::squish)+
    theme_bw()+
    ggtitle(paste(reg,"Composite Ensemble"))
}

if(save){
  write.csv(medianCompositeTS, row.names=FALSE, file=file.path(saveDir,'MedianTS_byRegion.csv'))
  for (reg in names(compositeEnsemble)){
    write.csv(compositeEnsemble[[reg]], row.names=FALSE,file = file.path(saveDir,paste(reg,'.csv',sep='')))
  }
  print("csv files saved")
}







#Save--------------------------------------------------------------------------------

#Code for sampleing and weighting but not really worth it. Widens the spread but minimal impact on median/trends. And code is clunky

# #Loop to composite (by region)-----
#   #Filter the TS by cluster name and make sure have enough values
#   regCount <- length(lipdReg)
#   #Skip if number of records is too few
#   if(regCount < nThresh) next
#   for (i in 1:regCount){
#     if (lipdReg[[i]]$climateInterpretation1_interpDirection == 'negative'){
#       lipdReg[[i]]$paleoData_values <- lipdReg[[i]]$paleoData_values*-1
#     }
#   }
#   #
#   #Calculate data density throughout the Holocene 
#   timeN <- plotTimeAvailabilityTs(lipdReg,
#                                   age.range=c(0,12000),
#                                   group.var ='Category',
#                                   step=binsize)$dat %>%
#     group_by(yvec) %>% 
#     summarise(count=sum(value),countPct=sum(value)/regCount)
#   #Determine search range based on 0-10ka with >50% proxy coverage
#   MH <- which(timeN$yvec==6000)
#   if (length(which(timeN$countPct<=0.5) > 0)){
#     idx <- c(max(which(timeN$countPct[1:MH]<=0.5),
#                  which(timeN$yvec==0)),
#              min(which(timeN$countPct[MH:length(timeN$countPct)]<=0.5)+MH,
#                  which(timeN$yvec==12000)))
#   } else{
#     idx <- c(which(timeN$yvec==0),which(timeN$yvec==12000))
#   }
#   searchMin  <- timeN$yvec[idx[1]] #age BP
#   searchMax  <- timeN$yvec[idx[2]] #age BP
#   #
#   #setup and run ensemble ##### This is the main part of the code to edit ##### 
#   compEns <- matrix(NA,nrow = length(binYears),ncol=nens)
#   for (i in 1:nens){
#     #Sample from full region
#     #wght <- pullTsVariable(lipdReg,'ageRange')/pullTsVariable(lipdReg,'ageResPlus')
#     #wght[which(is.na(wght))] <- mean(wght,na.rm=TRUE)
#     #set.seed(5)
#     #lipdRegSample <- sample(x    = pullTsVariable(lipdReg,'paleoData_TSid'),
#                             #prob = wght / sum(wght),
#                             #size = length(lipdReg) * samplePct)
#     #sample <- which(pullTsVariable(lipdReg,'paleoData_TSid') %in% pullTsVariable(lipdReg,'paleoData_TSid'))
#     #sample <- which(pullTsVariable(lipdReg,'paleoData_TSid') %in% lipdRegSample)
#     #lipdRegSample <- lipdReg#[sample]
#     #Composite
#     if (climVar == 'T'){std <- FALSE} else{std <- TRUE}
#     set.seed(5)
#     tc <- compositeR::compositeEnsembles(fTS                  = lipdRegSample,
#                                          binvec               = binvec,
#                                          stanFun              = standardizeMeanIteratively,
#                                          binFun               = simpleBinTs,
#                                          ageVar               = "age",
#                                          alignInterpDirection = FALSE,
#                                          spread               = TRUE,
#                                          duration             = searchDur,
#                                          searchRange          = c(searchMin,searchMax),
#                                          normalizeVariance    = std,
#                                          scope                = "climate",
#                                          minN                 = 3) 
#     compEns[,i] <- tc$composite
#   }
#   mean(apply(compEns,1,median,na.rm=TRUE))
#   #
#   # Return reconstruction and additional data for plotting
#   compositeEns[[reg]]      <- compEns[which(binYears==0):which(binYears==12000),]
#   medianCompositeTS[[reg]] <- apply(compEns,1,median,na.rm=TRUE)
# }




