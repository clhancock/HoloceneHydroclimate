library(lipdR)
library(maptools)
library(ncdf4)
library(rworldmap)
githubDir <- getwd()
lipdData  <- readRDS(file=file.path(githubDir,'Data','LiPD','lipdData.rds'))

PROJ       <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
PROJorig   <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
refregions <- readShapePoly(file.path(githubDir,'Data','IPCC_refRegions','IPCC-WGI-reference-regions-v4.shp'),
                            proj4string=CRS(PROJorig))
refregions <-  spTransform(refregions, CRSobj = PROJ)

climVar <- 'HC'
lipdTSO <- lipdData[[climVar]]

modelData <- vector(mode='list')
for (model in c('hadcm','trace')){
  modelData[[model]] <- readRDS(file=file.path(githubDir,'Data','Model_transient','rData','raw',paste(model,'.rds',sep='')))
}  


#LAND MASK
traceMask <- nc_open(file.path(githubDir,'Data','Model_transient','netcdf',
                               'trace','trace.01-36.22000BP.cam2.LANDFRAC.22000BP_decavg_400BCE.nc'))
traceMask0 <- ncvar_get(traceMask,'LANDFRAC')[,,2204]
traceMask0[which(traceMask0>=0.5)] <- 1
traceMask0[which(traceMask0< 0.5)] <- NA

hadcmMask <- nc_open(file.path(githubDir,'Data','Model_transient','netcdf',
                               'hadcm','deglh.vn1_0.ht_mm_srf.monthly.ANN.100yr.nc'))
hadcmMaskMask0 <- ncvar_get(hadcmMask,'ht_mm_srf')[,,250]
hadcmMaskMask0[which(hadcmMaskMask0> 0)] <- 1
hadcmMaskMask0[which(hadcmMaskMask0<=0)] <- NA


binYears <- seq(0,12000,200)
seasons <- c('ANN','DJF','JJA')
vars    <- c('precip','pe','temp')
names   <- apply(expand.grid(seasons, vars ), 1, paste, collapse="_")

modelComposites <- vector(mode='list')
for (model in names(modelData)){
  modelComposites[[model]] <- vector(mode='list')
  md <- modelData[[model]][['ANN']]
  regionNamesMatrix <- matrix(NA,ncol=length(md$lat),nrow=length(md$lon))
  if (model == 'trace'){
    mask <- traceMask0
  } else{mask <- hadcmMaskMask0}
  for (i in 1:length(md$lon)){
    lon <- md$lon[i]
    pointData <- data.frame(longitude= rep(md$lon[i]*0.99999,length(md$lat)),
                            latitude = md$lat)
    pointData <- SpatialPointsDataFrame(coords=pointData,data = pointData, 
                                        proj4string = CRS(PROJorig))
    pointData <-  spTransform(pointData, CRSobj = PROJ)
    regionNamesMatrix[i,] <- as.character(over(pointData, refregions)$Acronym)
  }
  for (region in (unique(pullTsVariable(lipdTSO,'geo_ipccRegion')))){
    print(region)
    modelComposites[[model]][[region]] <- vector(mode='list')
    regionTS    <- lipdTSO[which(pullTsVariable(lipdTSO,'geo_ipccRegion') == region)]
    #
    regionMask <- regionNamesMatrix
    regionMask[which(regionMask!=region)] <- NA
    regionMask[which(regionMask==region)] <- 1
    regionMask <- apply(regionMask,c(1,2),as.numeric)
    #
    regionProxyArray <- array(dim=c(length(binYears),c(length(names),length(regionTS))),
                         dimnames=list('age'=binYears,
                                       'variable'=names,
                                       'TSid'=pullTsVariable(regionTS,'paleoData_TSid')))
    regionModelArray <- array(dim=c(length(binYears),c(length(names),1)),
                              dimnames=list('age'=binYears,
                                            'variable'=names,
                                            '3'=1))
    #
    regionModelSDs <- c()
    for (season in seasons){
      md <- modelData[[model]][[season]]
      #
      for (var in vars){
        w <- which(colnames(regionProxyArray)==paste(paste(season,var,sep='_')))
        #proxynetwork
        for (i in 1:length(regionTS)){
          ilat <- min(which(abs(regionTS[[i]]$geo_latitude -md$lat)==min(abs(regionTS[[i]]$geo_latitude -md$lat))))
          ilon <- min(which(abs(regionTS[[i]]$geo_longitude-md$lon)==min(abs(regionTS[[i]]$geo_longitude-md$lon))))
          if (is.null(regionTS[[i]]$ageMin)) {regionTS[[i]]$ageMin<-0}
          if (is.null(regionTS[[i]]$ageMax)) {regionTS[[i]]$ageMin<-12000}
          tMin <- min(which(abs(regionTS[[i]]$ageMin-md$time)==min(abs(regionTS[[i]]$ageMin-md$time))))
          tMax <- min(which(abs(regionTS[[i]]$ageMax-md$time)==min(abs(regionTS[[i]]$ageMax-md$time))))
          siteData <- md[[var]][ilon,ilat,]
          siteData <- (siteData-mean(siteData))/sd(siteData)
          siteData[-c(tMin:tMax)] <- NA
          regionProxyArray[,w,i] <- siteData
        }
        #fulRegion
        fullRegionArray <- array(dim=dim(md[[var]]))
        for (i in 1:dim(md[[var]])[3]){
          fullRegionArray[,,i] <- md[[var]][,,i]*regionMask
        }
        regionModelArray[,w,1] <- apply(fullRegionArray,3,median,na.rm=TRUE)
        regionModelSDs <- c(regionModelSDs,mean(apply(fullRegionArray,3,sd,na.rm=TRUE),na.rm=TRUE))
      }
    }
    proxyMatrix <- matrix(nrow=length(binYears),ncol=length(regionTS))
    modelMatrix <- matrix(nrow=length(binYears),ncol=length(regionTS))
    for (i in 1:length(regionTS)){
      ts <- regionTS[[i]]
      season <- "NA"
      season <- ts$climateInterpretation1_seasonalityGeneral
      if (is.null(season)){
        season <- 'ANN'
      }else if (tolower(season) == 'summeronly'){
        if (ts$geo_latitude > 0){
               season <- 'JJA'
        } else{season <- 'DJF'}
      } else if (tolower(season)  == 'winteronly'){
        if (ts$geo_latitude > 0){
               season <- 'DJF'
        } else{season <- 'JJA'}
      } else{season <- 'ANN'}
      #
      var <- ts$climateInterpretation1_variable
      if (climVar == 'T'){
        var <- 'temp'
      } else if (var == 'P'){
        var <- 'precip'
      } else{var <- 'pe'}
      proxyMatrix[,i] <- regionProxyArray[,paste(season,var,sep='_'),ts$paleoData_TSid]
    }
    colnames(proxyMatrix) <- pullTsVariable(regionTS,'paleoData_TSid')
    modelComposites[[model]][[region]][['proxySiteArraySpread']] <- regionProxyArray
    modelComposites[[model]][[region]][['proxySiteMedian']] <- apply(regionProxyArray,c(1,2),median,na.rm=TRUE)
    modelComposites[[model]][[region]][['proxySiteByInterp']] <- proxyMatrix
    modelComposites[[model]][[region]][['fullRegionMedian']] <- regionModelArray[,,1]
    modelComposites[[model]][[region]][['fullRegionStDev']] <- regionModelSDs
  }
}

saveRDS(modelComposites,file.path(githubDir,'Data','Model_transient','rData','regionalComposites',
                                  paste(climVar,'modelComposite.rds',sep='_')))

