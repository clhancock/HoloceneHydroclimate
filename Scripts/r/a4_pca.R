library(geoChronR)

#githubDir <- getwd()
climVar <- 'HC'
compositeData <- read.csv(file.path(githubDir,'Data','RegionComposites',
                                    climVar,'MedianTSbyRegion.csv'))
lipdTSO <- readRDS(file=file.path(githubDir,'Data','LiPD','lipdData.rds'))[[climVar]]

#Calculate data to plot and connect it to region data. 
#More complicated than needed because part of larger code to make series of figures for each region 
pcaData <-vector(mode='list')
for (reg in names(compositeData)[-1]){ #First column is the time variable
  pcaData[[reg]]           <- vector(mode="list")
  pcaData[[reg]][['Name']] <- reg
  #Load region composite matrix csv
  pcaData[[reg]]$compMatix <- read.csv(file.path(githubDir,'Data',
                                                 'RegionComposites',climVar,
                                                 paste(reg,'csv',sep='.')))
  #
  #
  lipdReg    <- filterTs(lipdTSO,paste('geo_ipccRegion ==',reg))
  regCount   <- length(lipdReg)
  #Skip if number of records is too few
  #
  #Calculate data density throughout the Holocene 
  timeN <- plotTimeAvailabilityTs(lipdReg,
                                  age.range=c(0,12000),
                                  group.var ='Category',
                                  step=100)$dat %>%
    group_by(yvec) %>% 
    summarise(count=sum(value),countPct=sum(value)/regCount)
  #Determine search range based on 0-10ka with >50% proxy coverage
  MH <- which(timeN$yvec==6000)
  if (length(which(timeN$countPct<=0.5) > 0)){
    idx <- c(max(which(timeN$countPct[1:MH]<=0.5),
                 which(timeN$yvec==0)),
             min(which(timeN$countPct[MH:length(timeN$countPct)]<=0.5)+MH,
                 which(timeN$yvec==12000)))
  } else{
    idx <- c(which(timeN$yvec==0),which(timeN$yvec==12000))
  }
  searchMin  <- which(compositeData[,1]==timeN$yvec[idx[1]]) #age BP
  searchMax  <- which(compositeData[,1]==timeN$yvec[idx[2]]) #age BP
  print(c(searchMin,searchMax))
  pcaData[[reg]]$compMatix[1:searchMin,] <- NA
  pcaData[[reg]]$compMatix[searchMax:length(compositeData[,1]),] <- NA
  
  #
  
}

#Set up PCA input
pcaInList <- vector(mode = "list")
for (reg in names(compositeData)[-1]){
  pcaInList[[reg]]   <- binEns(compositeData[,1], #Time sequence of data
                               as.matrix(pcaData[[reg]]$compMatix), #Output of composite code
                               seq(0,12000,100), #bin only over this range change to 8000 to make work
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
