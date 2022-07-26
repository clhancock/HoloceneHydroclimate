#Purpose--------------------------------------------------------------------------------

#goal:   Create r.data files for temp12k and hc12k to reference from. 
 #       Also assign regional designations and additional metadata. 
#input:  url for LiPD files
#output: rdata file for Temp12k and HC12k timeseries lists
 #          
#author: chris hancock

#'WEB5aca062f',#paleodata to fill in data which was NA for Eilandvlei.Wuendsch.2018.lpd
# 'WEBeab5d1e0',#New ages for LagunaLaGaiba.Fornace.2016.lpd
# tsidListHC <- tsidListHC[-which(tsidListHC ==  'GH129d82be')] #New ages
# lipdData$HC <- lipdData$HC[-(which(pullTsVariable(lipdData$HC,"paleoData_TSid")=='WEBeab5d1e0')[1])]

#Load Packages--------------------------------------------------------------------------------

library(lipdR)
library(maptools)
library(proj4)
library(sf)
library(sp)
library(tidyverse)


#Set up directories and names--------------------------------------------------------------------------------

dir <- getwd()# '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate' #

tempVers <- '1_0_2' #Which versions of datesets to use
hcVers   <- '0_4_0' #Which versions of datesets to use


#Load Data--------------------------------------------------------------------------------

#Load ipcc region spatial data
PROJ <- '+proj=robin   +ellps=WGS84 +datum=WGS84 +no_defs +lon_0=0 +x_0=0 +y_0=0 +units=m'
PROJorig <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
load(url('https://github.com/SantanderMetGroup/ATLAS/blob/main/reference-regions/IPCC-WGI-reference-regions-v4_R.rda?raw=true'))
refregions <-  spTransform(IPCC_WGI_reference_regions_v4, CRSobj = PROJ)

#Load Lipd Files
D__t <- readLipd(paste('http://lipdverse.org/Temp12k/',tempVers,'/Temp12k',tempVers,'.zip',sep=''))
D_hc <- readLipd(paste('http://lipdverse.org/HoloceneHydroclimate/',hcVers,'/HoloceneHydroclimate',hcVers,'.zip',sep=''))


#Assign tsids for data compilations based on within correct dataset and version--------------------------------------------------------------------------------

TS__t <- splitInterpretationByScope(extractTs(D__t))
TS_hc <- splitInterpretationByScope(extractTs(D_hc))

getTSfromLiPDs <- function(input,compilationName,versNo){
  out <- vector(mode='list')
  for (tso in input){
    tsNo <- 1
    if (sum(!is.na(tso$paleoData_values[which(tso$age <= 12400)])) < 10){next}
    while (!is.null(tso[[paste('inCompilationBeta',tsNo,'_compilationName',sep='')]])){
      tsName <- tso[[paste('inCompilationBeta',tsNo,'_compilationName',sep='')]]
      tsVers <- tso[[paste('inCompilationBeta',tsNo,'_compilationVersion',sep='')]]
      if (tsName == compilationName & any(tsVers == versNo)){
        out[[length(out)+1]] <- tso
      }
      tsNo <- tsNo+1
    }
  }
  return(out)
}

TS__t <- getTSfromLiPDs(TS__t,'Temp12k',tempVers)
TS_hc <- getTSfromLiPDs(TS_hc,'HoloceneHydroclimate',hcVers)

print(paste("Temp:",length(TS__t)))
print(paste("HC:",  length(TS_hc)))


#Adjust metadata (mostly for HC)--------------------------------------------------------------------------------

lipdData <- list(T = TS__t, HC = TS_hc)

for (var in names(lipdData)){
  for (ts in 1:length(lipdData[[var]])){
    #
    tso                  <- lipdData[[var]][[ts]]
    tso$age              <- as.numeric(tso$age)
    tso$paleoData_values <- as.numeric(tso$paleoData_values)
    #
    #Make sure name of climate interp field is standardized 
    #
    if (is.null(tso$climateInterpretation1_direction) == FALSE){
      tso$climateInterpretation1_interpDirection <- tso$climateInterpretation1_direction
    }
    #
    #Add Ipcc regions to metadata
    #
    pointData <- data.frame(longitude=c(tso$geo_longitude),latitude=c(tso$geo_latitude))
    pointData <- spTransform(SpatialPointsDataFrame(coords=pointData,data = pointData, proj4string = CRS(PROJorig)),CRSobj = PROJ)
    tso$geo_ipccRegion <- as.character(over(pointData, refregions)$Acronym)
    #
    #Add Holocene age information
    #
    df <- data.frame(age = as.numeric(tso$age), 
                     values = as.numeric(tso$paleoData_values)) %>%
      filter(between(age,-100,12400)) %>%
      filter(!is.na(values)) %>%
      arrange(age)
    tso$ageMin     <- max(min(df$age),0)
    tso$ageMax     <- min(max(df$age),12000)
    tso$ageRange   <- min(diff(range(df$age)),12000)
    tso$ageRes     <- median(diff(df$age))
    tso$ageResPlus <- median(diff(df$age[which(diff(df$values) != 0)]))
    tso$paleoData_HoloceneValues <- df
    #
    #Revise seasonality & interpretation variable (M->P-E)
    #
    if (var == 'HC'){
      szn <- tso$climateInterpretation1_seasonalityGeneral
      if (is.null(szn)){                                szn <- 'Annual'
      } else if        (grepl('ummer',szn)){ 
        if (substr(szn, nchar(szn), nchar(szn)) == '+'){szn <- 'Summer+'
        } else{                                         szn <- 'Summer'}
      } else if (grepl('inter',szn)){ 
        if (substr(szn, nchar(szn), nchar(szn)) == '+'){szn <- 'Winter+'} 
        else{                                           szn <- 'Winter'}
      } else {                                          szn <- 'Annual'} 
      tso$climateInterpretation1_seasonalityGeneral <- szn
      #
      if (is.null(tso$climateInterpretation1_variable)){ 
        tso$climateInterpretation1_variable <- 'P-E'
      } else if (tso$climateInterpretation1_variable == 'M'){
        tso$climateInterpretation1_variable <- 'P-E'
      }
    }
    #
    #Add category information
    #
    if (var == 'HC'){
      archive  <- tso$archiveType
      proxy    <- tso$paleoData_proxy
      unit     <- tso$paleoData_units
      if (is.null(proxy) | is.null(archive)){
        tso$Category         <- 'Other'
        tso$CategorySpecific <- 'Other (not calibrated)'
      } else if (archive == 'Speleothem'){
        tso$Category           <- 'Speleothem'
        if (proxy == 'd18O' | proxy ==  'd13C'){
          tso$CategorySpecific <- paste(archive,' (','δ',substring(proxy, 2),')',sep='')
        } else{
          tso$CategorySpecific <- 'Speleothem (other)'
        }
      } else if (archive == 'LakeDeposits'){
        tso$Category           <- 'Shoreline'
        tso$CategorySpecific   <- 'Shoreline (Lake Level)'
      } else if (archive == 'GlacierIce'){
        tso$Category           <- 'Glacier Ice'
        tso$CategorySpecific   <- 'Glacier Ice'
      } else if (archive == 'LakeSediment' & proxy == 'd18O'){
        tso$Category           <- 'Lake Sediment (δ18O)'
        tso$CategorySpecific   <- 'Lake Sediment (δ18O)'
      } else if (proxy == 'dDwax'){
        tso$Category           <- 'Leaf Wax (δD)'
        tso$CategorySpecific   <- 'Leaf Wax (δD)'
      } else if (proxy == 'pollen'){
        tso$Category           <- 'Pollen'
        if (is.null(unit)){
          tso$CategorySpecific <- 'Pollen (not calibrated)'
        } else if (grepl('mm/',unit)){ 
          tso$CategorySpecific <- 'Pollen (calibrated)'
        } else {
          tso$CategorySpecific <- 'Pollen (not calibrated)'
        }
      } else {
        tso$Category           <- 'Other'
        if (is.null(unit)){
          tso$CategorySpecific <- 'Other (not calibrated)'
        } else if (grepl('mm/',unit)){ 
          tso$CategorySpecific <- 'Other (calibrated)'
        } else {
          tso$CategorySpecific <- 'Other (not calibrated)'
        }
      }
    } else{
      tso$Category         <- tso$paleoData_proxyGeneral
      tso$CategorySpecific <- tso$paleoData_proxyGeneral
    }
    #
    #Add source
    #
    if (var == 'HC'){
      if (is.null(tso$createdBy)){tso$createdBy <- ''}
      if (is.null(tso$originalDataUrl)){tso$originalDataUrl <- ''}
      if (is.null(tso$calibration_method)){tso$calibration_method <- ''}
      if (is.null(tso$pub1_doi)){tso$pub1_doi <- ''}
      if (is.null(tso$pub2_doi)){tso$pub2_doi <- ''}
      #
      if (grepl(tso$createdBy,'oxfordLakeStatus2Lipd')){
        tso$Source = 'Oxford Lake Level Database'
      }
      else if (grepl(tso$createdBy,'paleoDiver2lipd')){
        tso$Source = 'Liefert and Shuman, 2020'
      }
      else if (tso$createdBy =='sisal2lipd'){
        tso$Source = 'SISALv2'
      }
      else if (tso$originalDataUrl == 'wNAm'){
        tso$Source = 'wNA'
      }
      else if (tso$originalDataUrl == 'geochange.er.usgs.gov/midden/'){
        tso$Source = 'wNA'
      }
      else if (grepl(tso$createdBy,'LegacyClimate2LiPD')){
        tso$Source = 'Legacy Climate v1.0'
      }
      else if (tso$calibration_method == 'JM18_MAT'){
        tso$Source = 'Marsicek et al. (2018)'
      }
      else if (grepl('gov/paleo/study/15444',tso$originalDataUrl) | '10.5194/cp-10-1605-2014' == tso$pub2_doi){
        tso$Source = 'Arctic Holocene'
      }
      else if (substr(tso$dataSetName,1,2) == 'LS'){
        tso$Source = 'iso2k'
      }
      else if (tso$originalDataUrl == 'https://essd.copernicus.org/articles/12/2261/2020/'){
        tso$Source = 'iso2k'
      }
      else if (grepl('10.25921/4RY2-G808',tso$originalDataUrl) | grepl('/paleo/study/27330',tso$originalDataUrl)){
        tso$Source = 'Temp12k'
      } else{tso$Source = 'Other'}
    } else{tso$Source <- 'Temp12k'}
    #
    lipdData[[var]][[ts]] <- tso
  }
}

print(paste("Temp:",length(lipdData$T)))  #810
print(paste("HC:",  length(lipdData$HC))) #663


#Save--------------------------------------------------------------------------------

saveRDS(lipdData, file.path(dir,'Data','Proxy','LiPD','lipdData.rds'))


#Create and Save Summary Table for Data--------------------------------------------------------------------------------
for (var in names(lipdData)){
  lipdTSO <- lipdData[[var]]
  proxyDf <- tibble(dataset       = pullTsVariable(lipdTSO,'dataSetName'),
                    tsid          = pullTsVariable(lipdTSO,'paleoData_TSid'),
                    longitude     = pullTsVariable(lipdTSO,'geo_longitude'),
                    latitude      = pullTsVariable(lipdTSO,'geo_latitude'),
                    ipccReg       = pullTsVariable(lipdTSO,'geo_ipccRegion'),
                    archive       = pullTsVariable(lipdTSO,'archiveType'),
                    proxy         = pullTsVariable(lipdTSO,'paleoData_proxy'),
                    Category      = pullTsVariable(lipdTSO,'Category'),
                    CategorySpec  = pullTsVariable(lipdTSO,'CategorySpecific'),
                    minAge        = pullTsVariable(lipdTSO,'ageMin'),
                    maxAge        = pullTsVariable(lipdTSO,'ageMax'),
                    ageRange      = pullTsVariable(lipdTSO,'ageRange'),
                    ageRes        = pullTsVariable(lipdTSO,'ageRes'),
                    ageResPlus    = pullTsVariable(lipdTSO,'ageResPlus'),
                    season        = pullTsVariable(lipdTSO,'climateInterpretation1_seasonalityGeneral'),
                    climInterp    = pullTsVariable(lipdTSO,'climateInterpretation1_variable'),
                    source        = pullTsVariable(lipdTSO,'Source'),
                    direction     = pullTsVariable(lipdTSO,'climateInterpretation1_interpDirection'),
                    ka_0.5        = rep(NA, length(lipdTSO)),
                    ka_4          = rep(NA, length(lipdTSO)),
                    ka_6          = rep(NA, length(lipdTSO)),
                    ka_8          = rep(NA, length(lipdTSO)),
                    ka_10         = rep(NA, length(lipdTSO)),
                    maxValAge     = rep(NA, length(lipdTSO)),
                    minValAge     = rep(NA, length(lipdTSO)),
  ) #
  for (tso in lipdTSO){
    i = which(proxyDf$tsid == tso$paleoData_TSid)
    vals = tso$paleoData_values[which(between(tso$age,0,12000))]
    ages =              tso$age[which(between(tso$age,0,12000))]
    for (ka in c(0.5,seq(4,10,2))){
      bounds <- which(between(ages, 1000*ka-500, 1000*ka+500))
      proxyDf[i,paste('ka',ka,sep='_')] <- round(mean(vals[bounds],na.rm=TRUE),3)
    }
    if (!is.na(proxyDf[i,'direction']) & proxyDf[i,'direction']=='negative'){ vals  <- vals*-1}
    proxyDf[i,'maxValAge'] <- mean(ages[which(vals==max(vals,na.rm=TRUE))])
    proxyDf[i,'minValAge'] <- mean(ages[which(vals==min(vals,na.rm=TRUE))])
  }
  write.csv(proxyDf,file=file.path(dir,'Data','Proxy',paste('proxyMetadata_',var,'.csv',sep='')))
}


#Other --------------------------------------------------------------------------------
#tablelist <- vector(mode='list')
#for (archive in unique(proxyDf$archive)){
#  summary <- proxyDf[which(proxyDf$archive == archive),]%>% group_by(proxy) %>% summarise(n = n())
#  tablelist[[paste(archive,'---------------------------')]] = as.data.frame(summary)
#}



#Save summary tables
#capture.output(tablelist, file = file.path(dir,'Data','Proxy','ProxyArchiveCount_HoloceneHC.txt'))
#write.csv(proxyDf,file=file.path(dir,'Data','Proxy',paste('proxyMetaData_',climVar,'.csv',sep='')))

# library(janitor)
# z <- read.table("/Users/chrishancock/Downloads/Herzschuh-etal_2021_climate-recons_Asia (3).tab",
#                 header = F,sep="\t")  %>% 
#   row_to_names(row_number = 1) 
# 
# View(z)
# View(sort(unique(z$Site)))
# 
# View(z[which(grepl('Hidden',z$Site,ignore.case=TRUE)),])





