library(geoChronR)

githubDir <- getwd()
climVar <- 'HC'
compositeData <- read.csv(file.path(githubDir,'Data','RegionComposites',climVar,'MedianTSbyRegion.csv'))
#proxyData     <- readRDS(file=file.path(githubDir,'Data','LiPD','lipdData.rds'))[[climVar]]

#Calculate data to plot and connect it to region data. 
#More complicated than needed because part of larger code to make series of figures for each region 
pcaData <-vector(mode='list')
for (region in names(compositeData)[-1]){ #First column is the time variable
  pcaData[[region]]           <- vector(mode="list")
  pcaData[[region]][['Name']] <- region
  #Load region composite matrix csv
  pcaData[[region]]$compMatix <- read.csv(file.path(githubDir,'Data','RegionComposites',climVar,
                                                     paste(region,'.csv',sep='')))
}

#Set up PCA input
pcaInList <- vector(mode = "list")
ageRes <- 200
ageMax <- 12000 
for (region in names(compositeData)[-1]){
  pcaInList[[region]]   <- binEns(seq(0,12000,200), #Time sequence of data
                                  as.matrix(pcaData[[region]]$compMatix), #Output of composite code
                                  seq(0-(ageRes/2),12000+(ageRes/2),ageRes), #bin only over this range change to 8000 to make work
                                  bin.fun = mean, 
                                  max.ens = NA)
}
print(length(pcaInList))


pcaOut <- geoChronR::pcaEns(pcaInList)
pcaLoadingsMean <- apply(pcaOut$loadings,c(1,2),mean)
pcaLoadingsMean <- as.data.frame(pcaLoadingsMean,row.names = names(pcaInList))
colnames(pcaLoadingsMean) <- paste('PC',as.character(c(1:ncol(pcaLoadingsMean))),sep='')


#Save Output
saveRDS(pcaOut,file.path(githubDir,'Data','RegionComposites','PCA','pcaHC12k.rds'))
saveRDS(pcaLoadingsMean,file.path(githubDir,'Data','RegionComposites','PCA','pcaLoadingsMean.rds'))
