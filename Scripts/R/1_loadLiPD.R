#Purpose--------------------------------------------------------------------------------

#goal:   Create r.data files for temp12k and hc12k to reference from. 
 #       Also assign regional designations and additional metadata. 
#input:  url for LiPD files
#output: rdata file for Temp12k and HC12k timeseries lists
#        csv summarizing data in an easy to share table    
#        pdf of timeseries dashboards          
#
#author: chris hancock

#'WEB5aca062f',#paleodata to fill in data which was NA for Eilandvlei.Wuendsch.2018.lpd
#'WEBeab5d1e0',#New ages for LagunaLaGaiba.Fornace.2016.lpd
#tsidListHC <- tsidListHC[-which(tsidListHC ==  'GH129d82be')] #New ages
#lipdData$HC <- lipdData$HC[-(which(pullTsVariable(lipdData$HC,"paleoData_TSid")=='WEBeab5d1e0')[1])]

#Load Packages--------------------------------------------------------------------------------

library(cowplot)
library(ggstar)
library(ggplot2)
library(gsheet)
library(lipdR)
library(maptools)
library(proj4)
library(RCurl)
library(sf)
library(sp)
library(tidyverse)


#Set up directories and names--------------------------------------------------------------------------------

dir <-  '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/2021_HoloceneHydroclimate/2021_HoloceneHydroclimate' #
#getwd()
tempVers <- '1_0_2' #Which versions of temp12k datesets to use
hcVers   <- '0_7_0' #Which versions of Holocene hydroclimate datesets to use


#Load Data--------------------------------------------------------------------------------

#Load IPCC Region Spatial Data
PROJ <- '+proj=robin   +ellps=WGS84 +datum=WGS84 +no_defs +lon_0=0 +x_0=0 +y_0=0 +units=m'
PROJorig <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
load(url('https://github.com/SantanderMetGroup/ATLAS/blob/main/reference-regions/IPCC-WGI-reference-regions-v4_R.rda?raw=true'), verbose = TRUE)
refregions <-  spTransform(IPCC_WGI_reference_regions_v4, CRSobj = PROJ)

#Load Temperature and Hydroclimate Lipd Files
#D__t <- readLipd(paste0('https://lipdverse.org/Temp12k/',tempVers,'/Temp12k',tempVers,'.zip'))
#D_hc <- readLipd(paste0('https://lipdverse.org/HoloceneHydroclimate/',hcVers,'/HoloceneHydroclimate',hcVers,'.zip'))
#D_hc <- readLipd(file.path(dir,"Data","Proxy","LiPD",paste0("HoloceneHydroclimate",hcVers)))

#Assign tsids for data compilations based on within correct dataset and version--------------------------------------------------------------------------------
TS__t <- splitInterpretationByScope(extractTs(D__t))
TS_hc <- splitInterpretationByScope(extractTs(D_hc))

TS_hc[[which(pullTsVariable(TS_hc, "paleoData_TSID") == "S2LRbRVYOu2hWu")]][["age"]] <- TS_hc[[which(pullTsVariable(TS_hc, "paleoData_TSID") == "S2LRbRVYOu2hWu")]][["age"]]/1000

getTSfromLiPDs <- function(input,compilationName,versNo){
  out <- vector(mode='list')
  for (tso in input){
    tsNo <- 1
    while (!is.null(tso[[paste('inCompilationBeta',tsNo,'_compilationName',sep='')]])){
      tsName <- tso[[paste('inCompilationBeta',tsNo,'_compilationName',sep='')]]
      tsVers <- tso[[paste('inCompilationBeta',tsNo,'_compilationVersion',sep='')]]
      if (tsName == compilationName & any(tsVers == versNo)){
        if (sum(!is.na(tso$paleoData_values[which(tso$age <= 12400)])) < 10){print(tso$paleoData_TSid)}
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

#qcsheet <- gsheet2tbl('docs.google.com/spreadsheets/d/1rhYoL0B5OfE5A-rNwuQZfnmjI3Vj3NCX07r-Mif3Ncs') %>%
#  filter(inThisCompilation==TRUE)
#TS_hc<-TS_hc[which(pullTsVariable(TS_hc,"paleoData_TSid")%in%qcsheet$TSid)]
#print(paste("HC:",  length(TS_hc)))


TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB6e660f3d")]]$originalDataUrl <-"https://doi.org/10.1016/j.geomorph.2021.107896"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB8b32a2e0")]]$originalDataUrl <-"https://www.ncei.noaa.gov/access/paleo-search/study/5451"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB-a7bc6-3c3f-4541-b2f0-c02b6")]]$originalDataUrl <-"http://ncdc.noaa.gov/paleo/study/16677"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB8d973712")]]$originalDataUrl <- 'https://www.ncei.noaa.gov/access/paleo-search/study/32033'
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEBa5fcc0ed")]]$originalDataUrl <- 'https://doi.org/10.17632/zkn8rs76hy.1'
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB22a46599")]]$originalDataUrl <-'https://doi.org/10.17632/zkn8rs76hy.1'
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB8c80a27c")]]$originalDataUrl <- 'https://www.ncei.noaa.gov/access/paleo-search/study/13378'
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="LPD292174f8xxx")]]$originalDataUrl <- 'https://www.ncei.noaa.gov/access/paleo-search/study/18355'
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="NAm2kHydro096")]]$originalDataUrl <- 'https://www.ncei.noaa.gov/access/paleo-search/study/8640'
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="NAm2kHydro208")]]$originalDataUrl <- 'https://www.ncei.noaa.gov/access/paleo-search/study/23072'
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="NAm2kHydro215")]]$originalDataUrl <- 'https://www.ncei.noaa.gov/access/paleo-search/study/23072'
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEBe297c3e8")]]$originalDataUrl <- 'https://www.ncei.noaa.gov/access/paleo-search/study/33654'
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEBbf48ea2b")]]$originalDataUrl <- 'https://doi.org/10.6084/m9.figshare.12480344'
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB-49aac-c148-4e76-809a-362e4")]]$originalDataUrl <- 'https://doi.org/10.1594/PANGAEA.832385'
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB-84c26-cf61-4730-a2fb-5101f")]]$originalDataUrl <-'https://doi.org/10.1594/PANGAEA.921255'
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB21852084")]]$originalDataUrl <- 'https://doi.org/10.1594/PANGAEA.111890'
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="GH28d0af42")]]$originalDataUrl <- 'https://www.ncei.noaa.gov/access/paleo-search/study/11931'
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB-49d07-0b65-4966-88d1-10bad")]]$originalDataUrl <- 'https://www.ncei.noaa.gov/access/paleo-search/study/29692'
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB-af994-241c-4196-a6fa-3589e")]]$originalDataUrl <- 'https://www.ncei.noaa.gov/access/paleo-search/study/35393'
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB29098c59")]]$originalDataUrl <- "https://www.iceandclimate.nbi.ku.dk/data/"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB-c74a1-de1f-477c-becc-15b90")]]$originalDataUrl <- 'https://figshare.com/s/b4b5431fd9577afd95ef'
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB-8358f-5546-4fa4-9716-b98f3")]]$originalDataUrl <- 'https://doi.org/10.5880/GFZ.4.3.2021.005'
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="GH31325f74")]]$originalDataUrl <- 'https://www.ncei.noaa.gov/access/paleo-search/study/13097'
  
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB-a7bc6-3c3f-4541-b2f0-c02b6")]]$pub1_doi <- "10.1016/j.palaeo.2014.04.014"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB-7cf3f-ffa5-4f6d-8f4c-250b0")]]$pub1_doi <- "10.1016/j.quascirev.2021.107178"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEBa5fcc0ed")]]$pub1_doi <- "10.1016/j.quascirev.2021.106825"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB22a46599")]]$pub1_doi <- "10.1016/j.quascirev.2021.106825"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB15b6e09f")]]$pub1_doi <- "10.1016/j.palaeo.2017.09.032"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB-49aac-c148-4e76-809a-362e4")]]$pub1_doi <- "10.1016/j.quascirev.2014.04.006"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB-84c26-cf61-4730-a2fb-5101f")]]$pub1_doi <- "10.1016/j.quascirev.2014.04.006"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB21852084")]]$pub1_doi <- "10.1126/science.1080325"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="LPD17ad91b2")]]$pub1_doi <- "10.1016/j.epsl.2015.12.014"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="LPD277e0621")]]$pub1_doi <- "10.1016/j.epsl.2015.12.014"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB279cd5b2")]]$pub1_doi <- "10.1016/j.epsl.2015.12.014"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEBb00c1138")]]$pub1_doi <- "10.1016/j.epsl.2021.117148"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB19b9d003")]]$pub1_doi <- "10.1016/j.epsl.2021.117148"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB7baf4af4")]]$pub1_doi <- "10.1016/j.quascirev.2018.05.030"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB-49d07-0b65-4966-88d1-10bad")]]$pub1_doi <- "10.5194/cp-16-1097-2020"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB-af994-241c-4196-a6fa-3589e")]]$pub1_doi <- "10.1029/2021GL096611"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB-f9320-c765-469f-b4ee-6276b")]]$pub1_doi <- "10.1029/2020GL089183"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB-747c6-735e-4b08-a565-465fb")]]$pub1_doi <- "10.1002/jqs.2759"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB-85b98-f377-4cb2-a58b-dc9c3")]]$pub1_doi <- "10.1016/j.epsl.2005.02.025"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="LPD1d1e6750")]]$pub1_doi <- "10.1038/nature13196"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="LPD1756fdf4")]]$pub1_doi <- "10.1038/nature13196"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB-c74a1-de1f-477c-becc-15b90")]]$pub1_doi <- "10.1038/s41598-019-38626-3"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB-8358f-5546-4fa4-9716-b98f3")]]$pub1_doi <- "10.1038/s43247-022-00368-y"
TS_hc[[which(pullTsVariable(TS_hc,'paleoData_TSid')=="WEB-9119d-899e-4084-82c1-1fac9")]]$pub1_doi <- "10.1016/j.epsl.2012.03.016"
  
  




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
    tso$ageRange   <- tso$ageMax-tso$ageMin
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
    #Add Chron MetaData 
    #
    if (var == 'HC'){
      Dselect <- D_hc[[tso$dataSetName]]
      n_chron <- length(Dselect[["chronData"]][[1]][["measurementTable"]])
      n_paleo <- max(length(Dselect[["paleoData"]]),length(Dselect[["paleoData"]][[1]][["measurementTable"]]))
      if (n_chron == 0){                                         i_chron <- NA
      } else if (n_chron*n_paleo == 1){                          i_chron <- 1
      } else if (grepl("Composite",tso$paleoData_variableName)){
        offset<-c()
        for (j in 1:n_chron){
          chronAges <- Dselect[["chronData"]][[1]][["measurementTable"]][[j]][["age"]]$values
          offset <- c(offset,max(abs(tso$ageMin-min(chronAges,na.rm=TRUE)),
                                 abs(tso$ageMax-max(chronAges,na.rm=TRUE))))
        }
        i_chron <- which(offset==min(offset))
      } else if (n_chron != n_paleo & (n_chron+1 != n_paleo | tso$archiveType!='Speleothem')){      
        offset<-c()
        for (j in 1:n_chron){
          chronAges <- Dselect[["chronData"]][[1]][["measurementTable"]][[j]][["age"]]$values
          offset <- c(offset,max(abs(tso$ageMin-min(chronAges,na.rm=TRUE)),
                                 abs(tso$ageMax-max(chronAges,na.rm=TRUE))))
        }
        i_chron <- which(offset==min(offset))
      } else{
        for (i in 1:n_paleo){
          if (!is.null(Dselect[["paleoData"]][[1]][["measurementTable"]][[i]][[tso$paleoData_variableName]]$TSid)){
            if (Dselect[["paleoData"]][[1]][["measurementTable"]][[i]][[tso$paleoData_variableName]]$TSid == tso$paleoData_TSid){
              i_chron <- i
            }
          }
        }
      }
      if (tso$paleoData_TSid %in% c('lcRpbhIeSqLqRxKnBNF','lcRseG41EyZPFn7vBuL','LPD1498ea89','LPD03892251')){i_chron<-NA}
      if (!is.na(i_chron)){
        tso$chronData_table <- Dselect[["chronData"]][[1]][["measurementTable"]][[i_chron]]
        names <- names(tso$chronData_table)[grepl(c('age'),names(tso$chronData_table),ignore.case=TRUE)]
        for (name in c('type','error','min','max','hi','radio','14','lo','comment','use',"up","old","young","std","reservoir","σ","low","Rejected","±","comment","err","uncertainty","unc")){
          names <- names[which(grepl(name,names,ignore.case=TRUE)==FALSE)]
        }
        ageColName <- NA
        for (name in c("Th230/Th232","Median","calib.14C","14C.raw","14C Age","14C age BP","Cal yr chosen","age14C","14C age (yr BP)","C14 age dated",
                       tail(names,1),"corrected 230Th Age",
                       "Median Year BP",'age','Age','calAge','CalAge','CalibratedAge','Calibrated Age')){
          if (name %in% names(tso$chronData_table)){
            ageColName <- name
          } 
        }
        tso$chronData_ageName <- ageColName
        tso$chronData_ages <- as.numeric(tso$chronData_table[[tso$chronData_ageName]]$values)
        if (grepl("Composite",tso$paleoData_variableName)){
          tso$chronData_ages <- c()
          for (j in 1:n_chron){
            tso$chronData_ages <- c(tso$chronData_ages,as.numeric(Dselect[["chronData"]][[1]][["measurementTable"]][[j]][["age"]]$values))
          }
        }
        if ((length(which(tso$chronData_ages<13000))<=1) | is.na(sum(tso$chronData_ages)) | (sum(!is.na(tso$chronData_ages))<length(tso$chronData_ages)*0.5)){
          tso$chronData_table             <- NA
          tso$chronData_ageName           <- NA
          tso$chronData_ages              <- NA
          tso$chronData_ages_12k          <- NA
          tso$chronData_agesN_12k         <- NA
          tso$chronData_agesMaxGap_12k    <- NA
          tso$chronData_agesMedianGap_12k <- NA
        }else{
          tso$chronData_ages <- tso$chronData_ages[which(!is.na(tso$chronData_ages))]
          if(mean(tso$chronData_ages,na.rm=TRUE)<0){
            tso$chronData_ages<-tso$chronData_ages*-1 #
          }
          if(grepl("ka",tso$chronData_table[[tso$chronData_ageName]]$units)){
            tso$chronData_ages<-tso$chronData_ages*1000 #d 
          } else if(max(tso$chronData_ages,na.rm=TRUE)<300){
            tso$chronData_ages<-tso$chronData_ages*1000 #d 
          }
          tso$chronData_ages_12k          <- tso$chronData_ages[which(tso$chronData_ages<13000)]
          tso$chronData_agesN_12k         <- length(which(diff(tso$chronData_ages_12k)>0))+1
          tso$chronData_agesMaxGap_12k    <- max(diff(sort(c(tso$ageMin,tso$chronData_ages_12k,tso$ageMax))))
          tso$chronData_agesMedianGap_12k <- median(abs(diff(sort(c(tso$ageMin,tso$chronData_ages_12k,tso$ageMax)))))
        }
      } else{
        tso$chronData_table             <- NA
        tso$chronData_ageName           <- NA
        tso$chronData_ages              <- NA
        tso$chronData_ages_12k          <- NA
        tso$chronData_agesN_12k         <- NA
        tso$chronData_agesMaxGap_12k    <- NA
        tso$chronData_agesMedianGap_12k <- NA
      }
    }
    #
    #Add category information
    #
    if (var == 'HC'){
      archive  <- tso$archiveType
      proxy    <- tso$paleoData_proxy
      unit     <- tso$paleoData_units
      if(is.null(tso$climateInterpretation1_variableDetail)){tso$climateInterpretation1_variableDetail<-'Blank'} 
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
      } else if (tolower(tso$climateInterpretation1_variableDetail) == 'lakelevel@surface'){
        tso$Category           <- 'Shoreline'
        tso$CategorySpecific   <- 'Shoreline (Lake Level)'
      } else if (archive == 'GlacierIce'){
        tso$Category           <- 'Glacier Ice'
        tso$CategorySpecific   <- 'Glacier Ice (Accumulation)'
      } else if (archive == 'LakeSediment' & proxy == 'd18O'){
        tso$Category           <- 'Lake Sediment (δ18O)'
        tso$CategorySpecific   <- 'Lake Sediment (δ18O)'
      } else if (proxy == 'dDwax'){
        tso$Category           <- 'Leaf Wax'
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
      if (grepl('oxfordLakeStatus2Lipd',tso$createdBy)){
        tso$Source = 'Oxford Lake Level Database'
      }
      else if (grepl('paleoDiver2lipd',tso$createdBy)){
        tso$Source = 'Liefert and Shuman, 2020'
      }
      else if (tso$createdBy =='sisal2lipd'){
        tso$Source = 'SISALv2'
      }
      else if (grepl('LegacyClimate2LiPD',tso$createdBy)){
        tso$Source = 'Legacy Climate v1.0'
      }
      else if (tso$originalDataUrl == 'geochange.er.usgs.gov/midden/'){
        tso$Source = 'wNA'
      }
      else if (grepl('/study/30535',tso$originalDataUrl)){
        tso$Source = 'wNA'
      }
      else if (grepl('10.25921/bnxb-1n90',tso$originalDataUrl)){
        tso$Source = 'Mid-Latitude Holocene'
      }
      else if (grepl('/study/15444',tso$originalDataUrl)){
        tso$Source = 'Arctic Holocene'
      }
      else if (grepl('10.5194/essd-12-2261-2020',tso$originalDataUrl) | (substr(tso$dataSetName,1,2) == 'LS')){
        tso$Source = 'Iso2k'
      }
      else if (grepl('10.25921/4RY2-G808',tso$originalDataUrl) | grepl('/study/27330',tso$originalDataUrl)){
        tso$Source = 'Temp12k'
      } else{tso$Source = 'Other'}
    } else{tso$Source <- 'Temp12k'}
    #
    lipdData[[var]][[ts]] <- tso
  }
}

print(paste("Temp:",length(lipdData$T)))  #1319
print(paste("HC:",  length(lipdData$HC))) #817


#Save--------------------------------------------------------------------------------

saveRDS(lipdData, file.path(dir,'Data','Proxy','LiPD','lipdData.rds'))


#Create and Save Summary Table for Data--------------------------------------------------------------------------------
for (var in 'HC'){
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
                    ageCtrlN      = pullTsVariable(lipdTSO,'chronData_agesN_12k'),
                    ageCtrlMax    = pullTsVariable(lipdTSO,'chronData_agesMaxGap_12k'),
                    ageCtrlMedian = pullTsVariable(lipdTSO,'chronData_agesMedianGap_12k'),
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


proxyDf %>% 
  group_by(climInterp) %>% 
  summarize(count = n(),
            range = median(ageRange)/1000,
            resolution = median(ageRes))

#Other --------------------------------------------------------------------------------

var <- 'HC'

countries  <- spTransform(rworldmap::getMap("less islands"), CRSobj = PROJ)

proxyDf <- read.csv(file=file.path(dir,'Data','Proxy',paste('proxyMetadata_',var,'.csv',sep='')))
proxyDf <- as.data.frame(spTransform(SpatialPointsDataFrame(proxyDf[,c("longitude", "latitude")], proxyDf, proj4string=CRS(PROJorig)), CRSobj = PROJ))

lipdTSO <- lipdData[[var]]

plotSettings <- vector(mode='list')
plotSettings$names <- sort(unique(proxyDf$CategorySpec))
#
plotSettings$color <- as.character(plotSettings$names)
plotSettings$color[which(plotSettings$names=="Glacier Ice (Accumulation)")]             <- "powder blue"
plotSettings$color[which(plotSettings$names=="Shoreline (Lake Level)")]  <- "corn flower blue"
plotSettings$color[which(plotSettings$names=="Lake Sediment (δ18O)")]    <- "dark blue"
plotSettings$color[which(plotSettings$names=="Leaf Wax (δD)")]           <- "dark orchid" #δ
plotSettings$color[which(plotSettings$names=="Other (calibrated)")]      <- "grey40"
plotSettings$color[which(plotSettings$names=="Other (not calibrated)")]  <- "grey"
plotSettings$color[which(plotSettings$names=="Pollen (calibrated)")]     <- "forest green"
plotSettings$color[which(plotSettings$names=="Pollen (not calibrated)")] <- "yellowgreen"
plotSettings$color[which(plotSettings$names=="Speleothem (other)")]      <- "darkorange"
plotSettings$color[which(plotSettings$names=="Speleothem (δ13C)")]       <- "lightcoral"
plotSettings$color[which(plotSettings$names=="Speleothem (δ18O)")]       <- "firebrick"
#
plotSettings$shape <- as.character(plotSettings$names) 
plotSettings$shape[which(plotSettings$names=="Glacier Ice (Accumulation)")]             <- 12
plotSettings$shape[which(plotSettings$names=="Shoreline (Lake Level)")]  <- 21
plotSettings$shape[which(plotSettings$names=="Lake Sediment (δ18O)")]    <- 15
plotSettings$shape[which(plotSettings$names=="Leaf Wax (δD)")]           <- 5
plotSettings$shape[which(plotSettings$names=="Other (calibrated)")]      <- 6
plotSettings$shape[which(plotSettings$names=="Other (not calibrated)")]  <- 13
plotSettings$shape[which(plotSettings$names=="Pollen (calibrated)")]     <- 14
plotSettings$shape[which(plotSettings$names=="Pollen (not calibrated)")] <- 1
plotSettings$shape[which(plotSettings$names=="Speleothem (other)")]      <- 17
plotSettings$shape[which(plotSettings$names=="Speleothem (δ13C)")]       <- 23
plotSettings$shape[which(plotSettings$names=="Speleothem (δ18O)")]       <- 11 

plotSettings$shape <- as.numeric(plotSettings$shape)

for (reg in as.character(refregions@data[["Acronym"]])){
  n <- which(as.character(refregions@data[["Acronym"]])==reg)
  regTSO <- lipdTSO[which(pullTsVariable(lipdTSO,'geo_ipccRegion')==reg)]
  if (length(regTSO) == 0){next}
  refrenceRegShp <- subset(refregions, Acronym ==reg)
  regionDf <- proxyDf[which(proxyDf$ipccReg==reg),]
  idx <- which(plotSettings$names %in% regionDf$CategorySpec)
  regMap <-  ggplot() +
    geom_map(data=countries, map=fortify(countries),
             aes(x=long, y=lat, group=group, map_id=id), 
             fill = "grey80",color="grey90",size=0.4) +
    geom_map(data=refrenceRegShp, map=fortify(refrenceRegShp), 
             aes(x=long, y=lat, group=group, map_id=id),
             fill=NA, alpha=0.75, size=0.5, color='black') +
    geom_point(data=proxyDf,aes(x=longitude.1 ,y=latitude.1),size=1,color='Black') + 
    geom_star(data=regionDf,aes(x=longitude.1 ,y=latitude.1,fill=CategorySpec,starshape=CategorySpec),
              size=3,color='Black',starstroke=0.5) + 
    coord_fixed(xlim= range(fortify(refrenceRegShp)$lon), ylim=range(fortify(refrenceRegShp)$lat)) + 
    scale_fill_manual(     values=plotSettings$color[idx]) +
    scale_starshape_manual(values=plotSettings$shape[idx]) +
    theme_void() + 
    theme(legend.position = 'none') 
  #
  print(paste(n,reg,sep=". "))
  tsn <- 0
  for (ts in arrange(regionDf, desc(latitude), longitude)$tsid){
    tsn<-tsn+1
    siteDf  <- regionDf[which(regionDf$tsid==ts),]
    tso     <- regTSO[[which(pullTsVariable(regTSO,'paleoData_TSid')==ts)]]
    col <- plotSettings$color[which(plotSettings$names == tso$CategorySpecific)]
    shp <- plotSettings$shape[which(plotSettings$names == tso$CategorySpecific)]
    df <- data.frame(
      age    = as.numeric(tso$age), 
      values = as.numeric(tso$paleoData_values)) %>%
      filter(between(age,-100,12400)) %>%
      arrange(age) # tso$paleoData_HoloceneValues
    ages<-sort(c(df$age[1],df$age[-1]-diff(df$age)/2,df$age[-nrow(df)]+diff(df$age)/2,df$age[nrow(df)]))
    values <- rep(df$values, each=2)
    plt <- ggplot() +
      geom_hline(yintercept=mean(df$values,na.rm=TRUE)) +
      geom_star(aes(x=df$age,y=df$values),fill=col,color='black',starshape=shp,size=1,starstroke=0.5,alpha=0.7) + 
      geom_path(aes(x=ages,y=values),color=col,alpha=0.7) + 
      scale_x_reverse(name = "Age (yr BP)", limits=c(12000,0),expand=c(0,0),n.breaks=7,oob=scales::squish) +
      theme_bw() +
      ggtitle(paste0(reg,": ",tso$dataSetName," (",tso$paleoData_TSid,")")) 
    if (tso$climateInterpretation1_interpDirection == 'negative'){
      plt <- plt + 
        geom_point(aes(x=tso$chronData_ages_12k,y=tso$chronData_ages_12k*0+max(df$values,na.rm=TRUE)), color="black",size=2,shape=17)+
        scale_y_reverse(name=paste(tso$paleoData_variableName,' (',tso$paleoData_units,')',sep=''))
    } else{
      plt <- plt + 
        geom_point(aes(x=tso$chronData_ages_12k,y=tso$chronData_ages_12k*0+min(df$values,na.rm=TRUE)), color="black",size=2,shape=17)+
        scale_y_continuous(name=paste(tso$paleoData_variableName,' (',tso$paleoData_units,')',sep=''))
    }
    #
    map <- regMap + geom_star(data=siteDf,
                              aes(x=longitude.1, y=latitude.1), starshape=shp, fill=col,
                              size=3,color='gold',alpha=1,starstroke=2)
    h <- 0
    txt <- ggplot(siteDf) 
    for (name in c('geo_latitude','geo_longitude','geo_elevation',
                   'archiveType','Category','CategorySpecific',
                   'paleoData_proxyGeneral','paleoData_proxy','paleoData_proxyDetail','paleoData_variableName',
                   'climateInterpretation1_seasonalityGeneral',
                   'climateInterpretation1_variable','paleoData_units',
                   "chronData_agesN_12k" , "chronData_agesMaxGap_12k",
                   'Source','pub1_title','pub1_doi','pub2_doi','originalDataUrl')){
      h <- h-1
      #if (name %in% names(tso) == FALSE){tso[[name]]<- NA}
      txt <- txt + 
        annotate("text", x = 0, y = h, label = paste(name,': ',sep=''),hjust = 0,size=2.5) +
        annotate("text", x = 0.25, y = h, label = unique(tso[[name]]),hjust = 0,size=2.5)    
    }
    txt <- txt + scale_x_continuous(limits=c(0,0.6)) + theme_void() 
    bkg <- ggplot()+
      theme_void() +
      theme(plot.background = element_rect(fill = "White",color='White'))
    summary <- ggdraw(bkg) + 
      draw_plot(plt, x = 0,   y = 0.5, width = 1,   height = 0.5) +
      draw_plot(map, x = 0,   y = 0,   width = 0.5, height = 0.5) +
      draw_plot(txt, x = 0.5, y = 0,   width = 0.5, height = 0.5) 
    ggsave(file.path(dir,"Figures","Proxy","Dashboard",
                     paste(n,'_',reg,'_',tsn,'_',tso$paleoData_TSid,'.pdf',sep='')),device='pdf',
           plot=summary,width=10,height=6,units='in')
  }
}












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



