library(compositeR)
?library(doParallel)
library(dplyr)
library(foreach)
library(geoChronR)
library(lipdR)
library(magrittr)
library(purrr)
library(tidyverse)


#Set up directories and names--------------------------------------------------------------------------------

wd  <- '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/2021_HoloceneHydroclimate/2021_HoloceneHydroclimate' #
var <- 'HC'

#Load Data without winter+ or summer+ seasonality--------------------------------------------------------------------------------

lipdData <- readRDS(file.path(wd,'Data','Proxy','lipdData.rds'))[[var]]
lipdTSO  <- lipdData[-which(pullTsVariable(lipdData,"climateInterpretation1_seasonalityGeneral") %in% 
                              c('winter+','summer','Summer+','Winter+'))]

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
nThresh       <- 6   
#Set bin vectors
binvec   <- seq(ageMin-binsize/2, to = ageMax+binsize/2, by = binsize)+1
binYears <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))

#Show issue with compositeR--------------------------------------------------------------------------------
reg<-'ECA'
lipdReg  <- filterTs(lipdTSO,paste('geo_ipccRegion ==',reg))

tc <- compositeEnsembles(fTS                  = lipdReg,
                         binvec               = binvec,
                         stanFun              = standardizeMeanIteratively,
                         binFun               = simpleBinTs,
                         ageVar               = 'age',
                         alignInterpDirection = TRUE,
                         spread               = TRUE,
                         duration             = 3500,
                         searchRange          = c(1000,10000),
                         normalizeVariance    = TRUE,
                         minN                 = 3)

#binMatR <- as.matrix(purrr::map_dfc(fTS, binFun, binvec, 
                                #    ageVar = ageVar, spread = spread, gaussianizeInput = gaussianizeInput, 
                                #    alignInterpDirection = alignInterpDirection, scope = scope))
#if (spread) {
 # sp <- spreadPaleoData(age = ts[[ageVar]], value = ts$paleoData_values, 
  #                      spreadBy = spreadBy, maxGap = as.numeric(quantile(abs(diff(ts[[ageVar]])), 
   #                                                                       probs = 0.75, na.rm = TRUE)))
  #age <- sp$spreadAge
  #vals <- sp$spreadVal
#}

#Example timeseries with issue--------------------------------------------------------------------------------

ts    <-lipdReg[[3]]
age   <-ts$age
value <-ts$paleoData_values

#Function defaults
newAge   = NA
spreadBy = 10
maxGap   = 6
maxPct   = 0.75
minAge   = -69 

#Function
if (length(age) == 0) {
  return(list(spreadAge = age, spreadVal = value))
}
good <- which(is.finite(age))
if (length(good) < 2) {
  return(list(spreadAge = matrix(NA, ncol = length(age)), 
              spreadVal = matrix(NA, ncol = length(age))))
}
age <- age[good]
value <- value[good]
hasNas <- FALSE
if (all(is.na(newAge))) {
  newAge <- seq(ceiling(min(age)), floor(max(age)), by = spreadBy)
} else {
  spreadBy <- median(newAge, na.rm = TRUE)
}
newVals <- pracma::interp1(as.vector(age), as.vector(value), 
                           xi = newAge, method = "nearest")
nv1 <- which(newVals == value[1])
if (any(diff(nv1) != 1)) {
  nv1 <- nv1[1:(min(which(diff(nv1) != 1)))]
}
f <- min(newAge) - min(c(0, length(nv1) - 1)) * spreadBy
t <- min(newAge) - spreadBy
if (f < t) {
  begAge <- seq(from = f, to = t, by = spreadBy)
  begVal <- rep(newVals[1], times = length(begAge))
} else {
  begAge <- c()
  begVal <- c()
}
end <- length(newAge)
nv2 <- which(newVals == value[length(value)])
if (any(diff(nv2) != 1)) {
  nv2 <- nv2[max(which(diff(nv2) != 1)):length(nv2)]
}
f <- max(newAge) + spreadBy
t <- max(newAge) + max(c(0, length(nv2) - 1)) * spreadBy
if (f < t) {
  endAge <- seq(from = f, to = t, by = spreadBy)
  endVal <- rep(newVals[length(newVals)], times = length(endAge))
} else {
  endAge <- c()
  endVal <- c()
}
newAgeOut <- c(begAge, newAge, endAge)
newValsOut <- c(begVal, newVals, endVal)

d2n <- sign(diff(map_dbl(newAgeOut, function(x) min(abs(x - age)))))
age
newAgeOut
d2n
locmaxi <- which(diff(sign(diff(d2n))) == -2) + 1
locmaxi
locmax <- pracma::interp1(x = as.vector(c(newAgeOut[1], 
                                          newAgeOut[locmaxi], 
                                          newAgeOut[length(newAgeOut)])), 
                          as.vector(c(d2n[locmaxi[1]], d2n[locmaxi], d2n[locmaxi[length(locmaxi)]])), 
                          xi = newAgeOut, method = "nearest")


sort(age)[1:20]
newAgeOut[1:20]
map_dbl(newAgeOut[1:20], function(x) min(abs(x - age)))
diff(map_dbl(newAgeOut[1:20], function(x) min(abs(x - age))))
sign(diff(map_dbl(newAgeOut[1:20], function(x) min(abs(x - age)))))
diff(sign(diff(map_dbl(newAgeOut[1:20], function(x) min(abs(x - age))))))

lmpct <- d2n/locmax
if (!is.na(maxGap)) {
  newValsOut[d2n > maxGap] <- NA
}
if (!is.na(maxPct)) {
  newValsOut[lmpct > maxPct] <- NA
}
ty <- which(newAgeOut < minAge)
if (length(ty) > 0) {
  newAgeOut <- newAgeOut[-ty]
  newValsOut <- newValsOut[-ty]
}

