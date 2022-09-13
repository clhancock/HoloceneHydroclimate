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
lipdData <- lipdData[which(pullTsVariable(lipdData,'geo_ipccRegion') %in% c('WNA','ENA','CNA','NWN','NEN'))]  
lipdData <- lipdData[which(between(pullTsVariable(lipdData,'geo_latitude'),30,50))]



lipdTSO  <- vector(mode='list')
lipdTSO$Annual <- lipdData[which(pullTsVariable(lipdData,"climateInterpretation1_seasonalityGeneral") %in% c('Annual'))]
lipdTSO$Summer <- lipdData[which(pullTsVariable(lipdData,"climateInterpretation1_seasonalityGeneral") %in% c('Summer','Summer+'))]
lipdTSO$Winter <- lipdData[which(pullTsVariable(lipdData,"climateInterpretation1_seasonalityGeneral") %in% c('Winter','Winter+'))]
for (i in 1:length(lipdTSO)){
  print(length(lipdTSO[[i]]))
}
  

if(var == 'T'){
  lipdTSO <- filterTs(lipdTSO,'paleoData_units == degC')
  lipdTSO <- filterTs(lipdTSO,'paleoData_datum == abs')
  std <- FALSE      #Use calibrated data so no need to normalize variance
} else{std <- TRUE} #Normalize variance because data recorded with different units


#Set variables for composite code--------------------------------------------------------------------------------

nens          <- 500     #Ensemble numbers (lower = faster)
binsize       <- 100     #years (median resolution = 107yrs)
ageMin        <- 0       #age BP
ageMax        <- 12400   #age BP
searchDur     <- 3500    #yrs (for 3 lake deposit data points)
nThresh       <- 6       #minimum no. of records, else skip 

#Set bin vectors
binvec   <- seq(ageMin-binsize/2, to = ageMax+binsize/2, by = binsize)
binYears <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))


names <- c('Annual','Summer','Winter')

#Calculate reconstructions--------------------------------------------------------------------------------

#Set up data to add once regional composite is calculated
compositeEnsemble <- vector(mode='list')
medianCompositeTS <- data_frame(time=binYears[1:which(binYears==12000)])

#Loop to composite (by region)
for (reg in c(names)) {
  print(reg)
  lipdReg  <- lipdTSO[[reg]]
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
                             searchRange          = c(1000,10000),
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


