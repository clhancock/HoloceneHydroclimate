#Purpose-----
#goal: Create r.data files for temp12k and hc12k to refrence from. 
#Also assign regional designations and additional metadata. 
#input: LiPD files; refRegions
#output: rdata for Temp12k and HC12k ts
#author: chris hancock
#Load Packages-----
library(lipdR)
library(maptools)
library(proj4)
library(sf)
library(sp)
library(tidyverse)
#Set up directories and names-----
githubDir <- getwd()
#Which versions of datesets to use
tempVers <- '1_0_2'
hcVers   <- '0_4_0'
#Cut off for lake deposite numbers
lakeDeposNo <- 9
#Load ipcc region spatial data
PROJ <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
PROJorig <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
refregions <- readShapePoly(file.path(githubDir,'Data','IPCC_refRegions','IPCC-WGI-reference-regions-v4.shp'),
                            proj4string=CRS(PROJorig))
refregions <-  spTransform(refregions, CRSobj = PROJ)
#Load LiPD Data-----
#Load Lipd Files
D_temp <- readLipd(paste("http://lipdverse.org/Temp12k/",tempVers,"/Temp12k",tempVers,".zip",sep=''))
D_hc   <- readLipd(paste("http://lipdverse.org/HoloceneHydroclimate/",hcVers,"/HoloceneHydroclimate",hcVers,".zip",sep=''))
D_new  <- readLipd(file.path(githubDir,'Data','LiPD','new'))
#Assign tsids for data compilations based on within correct dataset and version -----

#Combine LiPD files and extract data from 2 sources without duplicates
TS_temp <- extractTs(D_temp)
TS_hc   <- c(extractTs(D_hc),extractTs(D_new))
#TS_hc   <- extractTs(D_hc)
TS_temp <- splitInterpretationByScope(TS_temp)
TS_hc   <- splitInterpretationByScope(TS_hc)
TS_all  <- c(TS_temp,TS_hc)

tsidListHC   <- c('WEBb3fd19e6', 'WEBd3eaa693', 'WEBb6841ae1','WEB4fba7605', 'WEBa1c23512', 'WEB8ffcaee7')
tsidListTemp <- c()
for (ts in TS_all){
  compNo <- 1
  while (!is.null(ts[[paste('inCompilationBeta',compNo,'_compilationName',sep='')]])){
    compName   <- ts[[paste('inCompilationBeta',compNo,'_compilationName',sep='')]]
    compVers   <- ts[[paste('inCompilationBeta',compNo,'_compilationVersion',sep='')]]
    if (compName == 'Temp12k'){
      if (any(compVers == tempVers)){
        tsidListTemp <- c(tsidListTemp,ts$paleoData_TSid)
      }
    } else if (compName == 'HoloceneHydroclimate'){
      if (any(compVers == hcVers)){
        if (sum(!is.na(ts$paleoData_values[which(ts$age < 12000)])) >= lakeDeposNo){
          tsidListHC <- c(tsidListHC,ts$paleoData_TSid)
        }
      }
    }
    compNo <- compNo+1
  }
}

tsidListHC <- tsidListHC[-which(tsidListHC == 'WEBaf733834')]
lipdData <- list(Temp = TS_temp[which(pullTsVariable(TS_temp,"paleoData_TSid") %in% tsidListTemp)], 
                 HC   = TS_hc[  which(pullTsVariable(TS_hc  ,"paleoData_TSid") %in% tsidListHC)])

print(paste("Temp:",length(lipdData$Temp)))
print(paste("HC:",length(lipdData$HC)))



#Add metadata-----
for (climVar in names(lipdData)){
  lipd <- lipdData[[climVar]]
  for (ts in 1:length(lipd)){
    lipd[[ts]]$age <- as.numeric(lipd[[ts]]$age)
    lipd[[ts]]$paleoData_values <- as.numeric(lipd[[ts]]$paleoData_values)
    #Make sure name of climate interp field is standardized 
    if (is.null(lipd[[ts]]$climateInterpretation1_direction) == FALSE){
      lipd[[ts]]$climateInterpretation1_interpDirection <- lipd[[ts]]$climateInterpretation1_direction
    }
    #
    #Add Ipcc regions to metdata
    if (climVar == 'HC' & lipd[[ts]]$geo_latitude==-57.5627 & lipd[[ts]]$geo_longitude==-18.0918){
      lipd[[ts]]$geo_latitude <--18.0918
      lipd[[ts]]$geo_longitude<--57.5627
    }
    if (lipd[[ts]]$paleoData_TSid == 'WEBeab5d1e0') {lipd[[ts]]$age <- lipd[[ts]]$age*1000}
    pointData <- data.frame(longitude=c(lipd[[ts]]$geo_longitude),latitude=c(lipd[[ts]]$geo_latitude))
    pointData <- SpatialPointsDataFrame(coords=pointData,data = pointData, 
                                        proj4string = CRS(PROJorig))
    pointData <-  spTransform(pointData, CRSobj = PROJ)
    lipd[[ts]]$geo_ipccRegion <- as.character(over(pointData, refregions)$Acronym)
    #
    #Add age information
    tso <- data.frame(age = as.numeric(lipd[[ts]]$age), values = as.numeric(lipd[[ts]]$paleoData_values)) %>%
      filter(between(age,0,12000)) %>%
      filter(!is.na(values)) %>%
      arrange(age)
    lipd[[ts]]$ageMin     <- min(tso$age)
    lipd[[ts]]$ageMax     <- max(tso$age)
    lipd[[ts]]$ageRange   <- diff(range(tso$age))
    lipd[[ts]]$ageRes     <- median(diff(tso$age))
    lipd[[ts]]$ageResPlus <- median(diff(tso$age[which(diff(tso$values) != 0)]))
    #
    if (climVar == 'HC'){
      if (is.null(lipd[[ts]]$climateInterpretation1_seasonalityGeneral)){
        lipd[[ts]]$climateInterpretation1_seasonalityGeneral <- 'Annual'
      } else if (grepl('ummer',lipd[[ts]]$climateInterpretation1_seasonalityGeneral)){
        lipd[[ts]]$climateInterpretation1_seasonalityGeneral <- 'Summer'
      } else if (grepl('inter',lipd[[ts]]$climateInterpretation1_seasonalityGeneral)){
        lipd[[ts]]$climateInterpretation1_seasonalityGeneral <- 'Winter'
      } else {
        lipd[[ts]]$climateInterpretation1_seasonalityGeneral <- 'Annual'
      }
      if (lipd[[ts]]$climateInterpretation1_variable == 'M'){
        lipd[[ts]]$climateInterpretation1_variable <- 'P-E'
      }
    }
    #
    #Add category information
    if (climVar == 'HC'){
      archive  <- lipd[[ts]]$archiveType
      proxy    <- lipd[[ts]]$paleoData_proxy
      unit     <- lipd[[ts]]$paleoData_units
      if (is.null(proxy) | is.null(archive)){
        lipd[[ts]]$Category         <- 'Other'
        lipd[[ts]]$CategorySpecific <- 'Other (not calibrated)'
      } else if (archive == 'Speleothem'){
        lipd[[ts]]$Category           <- 'Speleothem'
        if (proxy == 'd18O' | proxy ==  'd13C'){
          lipd[[ts]]$CategorySpecific <- paste(archive,' (','\u3B4',substring(proxy, 2),')',sep='')
        } else{
          lipd[[ts]]$CategorySpecific <- 'Speleothem (other)'
        }
      } else if (archive == 'LakeDeposits'){
        lipd[[ts]]$Category           <- 'Lake Deposits'
        lipd[[ts]]$CategorySpecific   <- 'Lake Deposits'
      } else if (archive == 'GlacierIce'){
        lipd[[ts]]$Category           <- 'Glacier Ice'
        lipd[[ts]]$CategorySpecific   <- 'Glacier Ice'
      } else if (archive == 'LakeSediment' & proxy == 'd18O'){
        lipd[[ts]]$Category           <- paste('Lake Sediment (','18O)',sep="\u3B4")
        lipd[[ts]]$CategorySpecific   <- paste('Lake Sediment (','18O)',sep="\u3B4")
      } else if (proxy == 'dDwax'){
        lipd[[ts]]$Category           <- paste('Leaf Wax (','D)',sep="\u3B4")
        lipd[[ts]]$CategorySpecific   <- paste('Leaf Wax (','D)',sep="\u3B4")
      } else if (proxy == 'pollen'){
        lipd[[ts]]$Category           <- 'Pollen'
        if (is.null(unit)){
          lipd[[ts]]$CategorySpecific <- 'Pollen (not calibrated)'
        } else if (unit == 'mm' | unit == 'mm/a'){
          lipd[[ts]]$CategorySpecific <- 'Pollen (calibrated)'
        } else {
          lipd[[ts]]$CategorySpecific <- 'Pollen (not calibrated)'
        }
      } else {
        lipd[[ts]]$Category           <- 'Other'
        if (is.null(unit)){
          lipd[[ts]]$CategorySpecific <- 'Other (not calibrated)'
        } else if (unit == 'mm' | unit == 'mm/a'){
          lipd[[ts]]$CategorySpecific <- 'Other (calibrated)'
        } else {
          lipd[[ts]]$CategorySpecific <- 'Other (not calibrated)'
        }
      }
    } else{lipd[[ts]]$Category <- lipd[[ts]]$paleoData_proxyGeneral}
    #
    if (climVar == 'HC'){
      if (is.null(lipd[[ts]]$createdBy)){lipd[[ts]]$createdBy <- ''}
      if (is.null(lipd[[ts]]$originalDataUrl)){lipd[[ts]]$originalDataUrl <- ''}
      if (is.null(lipd[[ts]]$calibration_method)){lipd[[ts]]$calibration_method <- ''}
      if (is.null(lipd[[ts]]$pub1_doi)){lipd[[ts]]$pub1_doi <- ''}
      if (is.null(lipd[[ts]]$pub2_doi)){lipd[[ts]]$pub2_doi <- ''}
      if (lipd[[ts]]$createdBy == 'http://github.com/nickmckay/oxfordLakeStatus2Lipd'){
        lipd[[ts]]$Source = 'Oxford Lake Levels Database'
      }
      #  dataSource.append('Oxford Lake Status')
      else if (lipd[[ts]]$createdBy =='sisal2lipd'){
        lipd[[ts]]$Source = 'SISAL'
      }
       # dataSource.append('SISAL (Comas-Bru et al., 2020)')
      else if (lipd[[ts]]$originalDataUrl == 'wNAm'){
        lipd[[ts]]$Source = 'wNA'
      }
       # dataSource.append('wNA')
      else if (lipd[[ts]]$originalDataUrl == 'geochange.er.usgs.gov/midden/'){
        lipd[[ts]]$Source = 'wNA'
      }
      else if (lipd[[ts]]$calibration_method == 'JM18_MAT'){
        lipd[[ts]]$Source = 'Marsicek et al. (2018)'
      }
      else if (grepl('gov/paleo/study/15444',lipd[[ts]]$originalDataUrl) | '10.5194/cp-10-1605-2014' == lipd[[ts]]$pub2_doi){
        lipd[[ts]]$Source = 'Arctic Holocene'
      }
      else if (substr('lipd[[ts]]$dataSetName',1,2) == 'LS'){
        lipd[[ts]]$Source = 'iso2k'
      }
      else if (lipd[[ts]]$originalDataUrl == 'https://essd.copernicus.org/articles/12/2261/2020/'){
        lipd[[ts]]$Source = 'iso2k'
      }
      else if (grepl('10.25921/4RY2-G808',lipd[[ts]]$originalDataUrl) | grepl('/paleo/study/27330',lipd[[ts]]$originalDataUrl)){
        lipd[[ts]]$Source = 'Temp12k'
      } else{lipd[[ts]]$Source = 'Other'}
    } else{lipd[[ts]]$Source <- 'Temp12k'}
    #
  }
  #require >0 change points for lakedeposits
  lipd                <- lipd[which(!is.na(pullTsVariable(lipd,"ageResPlus")))]
  lipdData[[climVar]] <- lipd[which(pullTsVariable(lipd,"climateInterpretation1_seasonalityGeneral") %in% c('summer+','winter+') == FALSE)]
}

print(paste("Temp:",length(lipdData$Temp))) #810
print(paste("HC:",length(lipdData$HC))) #664

#Save-----
saveRDS(lipdData, file.path(githubDir,'Data','LiPD','lipdData.rds'))

