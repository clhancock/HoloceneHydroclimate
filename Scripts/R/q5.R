PROJ <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
PROJorig <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
refregions <- readShapePoly(file.path(githubDir,'Data','IPCC_refRegions','IPCC-WGI-reference-regions-v4.shp'),
                            proj4string=CRS(PROJorig))
refregions <-  spTransform(refregions, CRSobj = PROJ)


pointData <- data.frame(longitude=c(lipd[[ts]]$geo_longitude),latitude=c(lipd[[ts]]$geo_latitude))
pointData <- SpatialPointsDataFrame(coords=pointData,data = pointData, 
                                    proj4string = CRS(PROJorig))
pointData <-  spTransform(pointData, CRSobj = PROJ)
lipd[[ts]]$geo_ipccRegion <- as.character(over(pointData, refregions)$Acronym)
#


for (model in c('trace')){
  modelData <- readRDS(file=file.path(githubDir,'Data','Model_transient','rData','raw',paste(model,'rds',sep='.')))\
  for (season in c('ANN')){
    pointData <- data.frame(longitude= modelData[[season]]$lon,
                            latitude = modelData[[season]]$lat)
    pointData <- SpatialPointsDataFrame(coords=pointData,data = pointData, 
                                        proj4string = CRS(PROJorig))
    pointData <-  spTransform(pointData, CRSobj = PROJ)
    for (variable in c('precip')){
      data <- modelData[[season]][[variable]]
      rast <- raster(data[,,1],crs=PROJorig)
    }
  }
}
#model_season_variable_full
#model_season_variable_proxy

z <- 
plot(refregions)