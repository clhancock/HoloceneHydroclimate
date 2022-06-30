#Purpose-----
#goal: Create r.data files for temp12k and hc12k to reference from. 
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
dir <- getwd()# '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate' #
#Which versions of datesets to use
tempVers <- '1_0_2'
hcVers   <- '0_4_0'
#Cut off for lake deposite numbers
LakeDepositsNo <- 9
#Load ipcc region spatial data
PROJ <- "+proj=robin   +ellps=WGS84 +datum=WGS84 +no_defs +lon_0=0 +x_0=0 +y_0=0 +units=m"
load(url("https://github.com/SantanderMetGroup/ATLAS/blob/main/reference-regions/IPCC-WGI-reference-regions-v4_R.rda?raw=true"))
refregions <-  spTransform(IPCC_WGI_reference_regions_v4, CRSobj = PROJ)

#Load LiPD Data-----
#Load Lipd Files
D_t    <- readLipd(paste("http://lipdverse.org/Temp12k/",tempVers,"/Temp12k",tempVers,".zip",sep=''))
D_hc   <- readLipd(paste("http://lipdverse.org/HoloceneHydroclimate/",hcVers,"/HoloceneHydroclimate",hcVers,".zip",sep=''))
D_new  <- readLipd(file.path(dir,'Data','Proxy','LiPD','new'))

#Assign tsids for data compilations based on within correct dataset and version -----
#Combine LiPD files and extract data from 2 sources without duplicates
TS_t   <- splitInterpretationByScope(extractTs(D_t))
TS_hc  <- splitInterpretationByScope(c(extractTs(D_hc),extractTs(D_new)))
TS_all <- c(TS_t,TS_hc)

#Sort data into list of TSIDs for appropriate variables. May not be needed after hc12k fully settled
tsidListHC <- c('WEBb3fd19e6', 'WEBd3eaa693', 'WEBb6841ae1','WEB4fba7605', 'WEBa1c23512', 'WEB8ffcaee7',
                'WEB5aca062f',#paleodata to fill in data which was NA for Eilandvlei.Wuendsch.2018.lpd
                'WEBeab5d1e0',#New ages for LagunaLaGaiba.Fornace.2016.lpd
                'WEB-ef183-b6de-44d7-8b46-61a47')#New ages for Alley.GISP2.2000
tsidListT  <- c()
for (tso in TS_all){
  tsNo <- 1
  while (!is.null(tso[[paste('inCompilationBeta',tsNo,'_compilationName',sep='')]])){
    tsName <- tso[[paste('inCompilationBeta',tsNo,'_compilationName',sep='')]]
    tsVers <- tso[[paste('inCompilationBeta',tsNo,'_compilationVersion',sep='')]]
    if (tsName == 'Temp12k'){
      if (any(tsVers == tempVers)){tsidListT <- c(tsidListT,tso$paleoData_TSid)}
    } else if (tsName == 'HoloceneHydroclimate'){
      if (any(tsVers == hcVers)){
        if (sum(!is.na(tso$paleoData_values[which(tso$age <= 12000)])) >= LakeDepositsNo){
          tsidListHC <- c(tsidListHC,tso$paleoData_TSid)
        } else{print(paste('Exclude TSid',tso$paleoData_TSid,'(',tso$archiveType,')'))}
      }
    }
    tsNo <- tsNo+1
  }
}

tsidListHC <- tsidListHC[-which(tsidListHC == 'WEBaf733834')] #Record length too short
tsidListHC <- tsidListHC[-which(tsidListHC == 'LPD6d6b5db7')] #Record length too short
tsidListHC <- tsidListHC[-which(tsidListHC == 'LPD2999f647')] #Non-linear
tsidListHC <- tsidListHC[-which(tsidListHC ==  'GH129d82be')] #New ages

lipdData <- list(T  = TS_t[ which(pullTsVariable(TS_t, "paleoData_TSid") %in% tsidListT)], 
                 HC = TS_hc[which(pullTsVariable(TS_hc,"paleoData_TSid") %in% tsidListHC)])

print(paste("Temp:",length(lipdData$T)))
print(paste("HC:",  length(lipdData$HC)))

#Fix metadata-----
lipdData$HC <- lipdData$HC[-(which(pullTsVariable(lipdData$HC,"paleoData_TSid")=='WEBeab5d1e0')[1])]
lipdData$HC[[which(pullTsVariable(lipdData$HC,"paleoData_TSid")=='WEBed535db2')]]$geo_latitude  <- -18.0918
lipdData$HC[[which(pullTsVariable(lipdData$HC,"paleoData_TSid")=='WEBed535db2')]]$geo_longitude <- -57.5627


#Add metadata (mostly for HC)-----
for (climVar in names(lipdData)){
  lipd <- lipdData[[climVar]]
  for (ts in 1:length(lipd)){
    lipd[[ts]]$age              <- as.numeric(lipd[[ts]]$age)
    lipd[[ts]]$paleoData_values <- as.numeric(lipd[[ts]]$paleoData_values)
    #Make sure name of climate interp field is standardized 
    if (is.null(lipd[[ts]]$climateInterpretation1_direction) == FALSE){lipd[[ts]]$climateInterpretation1_interpDirection <- lipd[[ts]]$climateInterpretation1_direction}
    #
    #Add Ipcc regions to metadata
    #
    pointData <- data.frame(longitude=c(lipd[[ts]]$geo_longitude),latitude=c(lipd[[ts]]$geo_latitude))
    pointData <- spTransform(SpatialPointsDataFrame(coords=pointData,data = pointData, proj4string = CRS(PROJorig)),CRSobj = PROJ)
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
      szn <- lipd[[ts]]$climateInterpretation1_seasonalityGeneral
      if (is.null(szn)){                                                        szn <- 'Annual'
      } else if        (grepl('ummer',szn)){ 
        if (substr(szn, nchar(szn), nchar(szn)) == '+'){szn  <- 'Summer+'} else{szn  <- 'Summer'}
      } else if (grepl('inter',szn)){ 
        if (substr(szn, nchar(szn), nchar(szn)) == '+'){szn  <- 'Winter+'} else{szn  <- 'Winter'}
      } else {                                                                  szn <- 'Annual'} 
      lipd[[ts]]$climateInterpretation1_seasonalityGeneral <- szn
      if (lipd[[ts]]$climateInterpretation1_variable == 'M'){lipd[[ts]]$climateInterpretation1_variable <- 'P-E'}
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
        lipd[[ts]]$SourceURL          <- 'https://researchdata.reading.ac.uk/256/'
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
    } else{
      lipd[[ts]]$Category         <- lipd[[ts]]$paleoData_proxyGeneral
      lipd[[ts]]$CategorySpecific <- lipd[[ts]]$paleoData_proxyGeneral
    }
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
  #Save TSid of removed files
  removedData <- pullTsVariable(lipd,"dataSetName")[which(!is.na(pullTsVariable(lipd,"ageResPlus")))]
  #require >0 change points for lakedeposits
  lipd                <- lipd[which(!is.na(pullTsVariable(lipd,"ageResPlus")))]
  lipdData[[climVar]] <- lipd#[which(pullTsVariable(lipd,"climateInterpretation1_seasonalityGeneral") %in% c('Summer+','Winter+') == FALSE)]
}

print(paste("Temp:",length(lipdData$T)))  #810
print(paste("HC:",  length(lipdData$HC))) #663

#Save-----
#saveRDS(lipdData, file.path(dir,'Data','Proxy','LiPD','lipdData.rds'))

#Save Summary Table for Data-----
climVar <- 'HC'
lipdTSO <- lipdData[[climVar]]
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
                  recordRange   = pullTsVariable(lipdTSO,'ageRange'),
                  recordRes     = pullTsVariable(lipdTSO,'ageRes'),
                  recordResPlus = pullTsVariable(lipdTSO,'ageResPlus'),
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
)

for (tso in lipdTSO){
  i = which(proxyDf$tsid == tso$paleoData_TSid)
  vals = tso$paleoData_values[which(between(tso$age,0,12000))]
  ages =              tso$age[which(between(tso$age,0,12000))]
  for (ka in c(0.5,seq(4,10,2))){
    proxyDf[i,paste('ka',ka,sep='_')] <- round(mean(vals[which(between(ages,1000*ka-500,1000*ka+500))],na.rm=TRUE),3)}
  if (!is.na(proxyDf[i,'direction']) & proxyDf[i,'direction']=='negative'){vals=vals*-1}
  proxyDf[i,'maxValAge'] <- mean(ages[which(vals==max(vals,na.rm=TRUE))])
  proxyDf[i,'minValAge'] <- mean(ages[which(vals==min(vals,na.rm=TRUE))])
}

tablelist <- vector(mode='list')
for (archive in unique(proxyDf$archive)){
  summary <- proxyDf[which(proxyDf$archive == archive),]%>% group_by(proxy) %>% summarise(n = n())
  tablelist[[paste(archive,'---------------------------')]] = as.data.frame(summary)
}

#Save summary tables
#capture.output(tablelist, file = file.path(dir,'Data','Proxy','ProxyArchiveCount_HoloceneHC.txt'))
#write.csv(proxyDf,file=file.path(dir,'Data','Proxy',paste('proxyMetaData_',climVar,'.csv',sep='')))

library(janitor)
z <- read.table("/Users/chrishancock/Downloads/Herzschuh-etal_2021_climate-recons_Asia (3).tab",
                header = F,sep="\t")  %>% 
  row_to_names(row_number = 1) 

View(z)
View(sort(unique(z$Site)))

View(z[which(grepl('Hidden',z$Site,ignore.case=TRUE)),])





