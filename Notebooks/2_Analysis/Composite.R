#Purpose--------------------------------------------------------------------------------
#goal: create regional composite based on IPCC regions using either temp or HC
#in:   LiPD files
#out:  individual csv for each region where columns represent iterations. 
#      summary table with each column as the regional median time series

#Load Packages--------------------------------------------------------------------------------


#library(compositeR)
devtools::install("/Users/chrishancock/Library/CloudStorage/OneDrive-NorthernArizonaUniversity/Research/Manuscript/HoloceneHydroclimate/compositeR")
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

wd = '/Users/chrishancock/Library/CloudStorage/OneDrive-NorthernArizonaUniversity/Research/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate'
var  <- 'HC'
save <- FALSE

#Load Data without winter+ or summer+ seasonality--------------------------------------------------------------------------------

lipdData <- readRDS(file.path(wd,'Data','Proxy','lipdData.rds'))[[var]]
lipdTSO  <- lipdData[-which(pullTsVariable(lipdData,"climateInterpretation1_seasonalityGeneral") %in% c('winter+','summer','Summer+','Winter+'))]

if(var == 'T'){
  lipdTSO <- filterTs(lipdTSO,'paleoData_units == degC')
  lipdTSO <- filterTs(lipdTSO,'paleoData_datum == abs')
  std <- FALSE       #Use calibrated data for T so no need to normalize variance #Change to true so that a more 1:1: comparison 4/24/23
} else{std <- TRUE} #Normalize HC variance because data recorded with different units


#Set variables for composite code--------------------------------------------------------------------------------

nens          <- 500     #Ensemble numbers (lower = faster)
binsize       <- 100     #years (median resolution = 109yrs)
ageMin        <- 0       #age BP 
ageMax        <- 12400   #age BP
searchDur     <- 3500    #yrs (for 3 lake deposit data points)
nThresh       <- 6       #minimum no. of records, else skip 

#Set bin vectors
binvec   <- seq(ageMin-binsize/2, to = ageMax+binsize/2, by = binsize)
binYears <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))

#ID regions to reconstruct based on number of records (nThresh)
regNames <- data.frame(name=pullTsVariable(lipdTSO,'geo_ipccRegion')) %>% 
  group_by(name) %>% 
  dplyr::summarise(n = n()) %>% 
  filter(n >= nThresh)

if (var != 'T'){
  regNames <- c(as.character(regNames$name),'EAN','SSA') #Add 2 SH regions with fewer records to gain global coverage
} else{regNames <- c(as.character(regNames$name))}

#Calculate reconstructions--------------------------------------------------------------------------------

#Set up data to add once regional composite is calculated
compositeEnsemble <- vector(mode='list')
medianCompositeTS <- data_frame(time=binYears[1:which(binYears==12000)])

#Loop to composite (by region)
for (reg in c('ENA')) { #regNames
  lipdReg  <- filterTs(lipdTSO,paste('geo_ipccRegion ==',reg))
  regN     <- length(lipdReg)
  # If sufficient data, use a sample for each iteration
  if (regN < 8){        pct <- 1.00
  } else if (regN < 30){pct <- 0.75
  }else{  pct <- 0.666}
  #
  #
  for (i in 1:length(lipdReg)){
    if (lipdReg[[i]]$climateInterpretation1_interpDirection == 'negative'){
      lipdReg[[i]]$paleoData_values <- lipdReg[[i]]$paleoData_values*-1
    }
  }
  #
  #
  set.seed(5) #Reproducibility
  ensOut <- foreach(i = 1:nens) %dopar% {
    tc <- compositeEnsembles(fTS                  = lipdReg[sample(seq(1,length(lipdReg)),length(lipdReg)*pct)],
                             binvec               = binvec,
                             stanFun              = standardizeMeanIteratively,
                             binFun               = simpleBinTs,
                             ageVar               = 'age',
                             alignInterpDirection = FALSE,
                             spread               = TRUE,
                             duration             = searchDur,
                             searchRange          = c(1000,10000),
                             normalizeVariance    = std,
                             minN                 = 3)
    return(list(composite = tc$composite,count = tc$count))
  }
  # Reformat Data
  regionComposite           <- as.matrix(purrr::map_dfc(ensOut,magrittr::extract,"composite"))
  rownames(regionComposite) <- binYears
  #Only Holocene 
  regionComposite           <- regionComposite[1:which(binYears==12000),]
  # Rescale for full timeseries
  if(std){regionComposite           <- (regionComposite-mean(regionComposite,na.rm=TRUE))/(sd(regionComposite,na.rm=TRUE))}
  # Save to matrix
  compositeEnsemble[[reg]]  <- regionComposite
  medianCompositeTS[[reg]]  <- apply(regionComposite,1,median,na.rm=TRUE)
}

#plot an example region to confirm that everything looks good
# plotTimeseriesEnsRibbons(X = binYears[1:which(binYears==12000)],Y = compositeEnsemble[['WNA']] )+
#   scale_x_reverse(name = "age (yr BP)",         oob = scales::squish)+
#   scale_y_continuous(name = "Standardized Anomaly",oob = scales::squish)+
#   theme_bw()#+ggtitle(paste("S Asia (using align interpretation = TRUE"))


#Save--------------------------------------------------------------------------------

if(save){
  write.csv(medianCompositeTS, row.names=FALSE, 
            file=file.path('Data','RegionComposites',var,
                           'MedianTS_byRegion.csv'))
  for (reg in names(compositeEnsemble)){
    write.csv(compositeEnsemble[[reg]], row.names=FALSE,
              file = file.path('Data','RegionComposites',var,
                               paste(reg,'.csv',sep='')))
  }
  print("csv files saved")
}




