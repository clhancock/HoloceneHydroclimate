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

wd  <- '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/2021_HoloceneHydroclimate/2021_HoloceneHydroclimate' #
var  <- 'HC'
save <- FALSE


#Load Data without winter+ or summer+ seasonality--------------------------------------------------------------------------------

lipdData <- readRDS(file.path(wd,'Data','Proxy','lipdData.rds'))[[var]]
lipdTSO  <- lipdData[-which(pullTsVariable(lipdData,"climateInterpretation1_seasonalityGeneral") %in% c('winter+','summer','Summer+','Winter+'))]

if(var == 'T'){
  lipdTSO <- filterTs(lipdTSO,'paleoData_units == degC')
  lipdTSO <- filterTs(lipdTSO,'paleoData_datum == abs')
  std <- FALSE      #Use calibrated data for T so no need to normalize variance
} else{std <- TRUE} #Normalize HC variance because data recorded with different units


#Set variables for composite code--------------------------------------------------------------------------------

nens          <- 5     #Ensemble numbers (lower = faster)
binsize       <- 100     #years (median resolution = 107yrs)
ageMin        <- -100       #age BP
ageMax        <- 12400   #age BP
searchDur     <- 3500    #yrs (for 3 lake deposit data points)
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

reg<-'ECA'
lipdReg  <- filterTs(lipdTSO,paste('geo_ipccRegion ==',reg))
ts<-lipdReg[[5]]
age<-ts$age
value<-ts$paleoData_values
newAge = NA
spreadBy = abs(mean(diff(binvec)))/10
maxGap = NA
maxPct = 0.75
minAge = -69 


#Loop to composite (by region)
for (reg in c(regNames)) {
  lipdReg  <- filterTs(lipdTSO,paste('geo_ipccRegion ==',reg))
  regN     <- length(lipdReg)
  set.seed(5) #Reproducibility
  ensOut <- foreach(i = 1:nens) %dopar% {
    tc <- compositeEnsembles(fTS                  = lipdReg,
                             binvec               = binvec,
                             stanFun              = standardizeMeanIteratively,
                             binFun               = simpleBinTs,
                             ageVar               = 'age',
                             alignInterpDirection = TRUE,
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
}

#plot region to confirm that everything looks good
plotTimeseriesEnsRibbons(X = binYears[1:which(binYears==12000)],Y = compositeEnsemble[['EAS']])+
  scale_x_reverse(name = "age (yr BP)",         oob = scales::squish)+
  scale_y_continuous(name = "Standardized Anomaly",oob = scales::squish)+
  theme_bw()+
  ggtitle(paste("E. Asia","Composite Ensemble"))


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




