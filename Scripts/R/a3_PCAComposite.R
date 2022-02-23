library(geoChronR)
library(cowplot)
librar
#githubDir <- getwd()
climVar <- 'HC'
compositeData <- read.csv(file.path(dataDir,'Data','RegionComposites',
                                    climVar,'MedianTSbyRegion.csv'))
lipdTSO <- readRDS(file=file.path(dataDir,'Data','LiPD','lipdData.rds'))[[climVar]]

#Calculate data to plot and connect it to region data. 
#More complicated than needed because part of larger code to make series of figures for each region 
pcaData <-vector(mode='list')
for (reg in names(compositeData)[-1]){ #First column is the time variable
  pcaData[[reg]]           <- vector(mode="list")
  pcaData[[reg]][['Name']] <- reg
  #Load region composite matrix csv
  pcaData[[reg]]$compMatix <- read.csv(file.path(dataDir,'Data',
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
saveRDS(pcaOut,file.path(dataDir,'Data','RegionComposites','PCA','pcaHC12k.rds'))
saveRDS(pcaLoadingsMean,file.path(dataDir,'Data','RegionComposites','PCA','pcaLoadingsMean.rds'))


pca <- readRDS(file.path(dataDir,'Data','RegionComposites','PCA','pcaHC12k.rds'))
regNames   <- names(pcaInList)
PROJ       <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
PROJorig   <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
refregions <- readShapePoly(file.path(dataDir,'Data','IPCC_refRegions','IPCC-WGI-reference-regions-v4.shp'),proj4string=CRS(PROJorig))
refregions <-  spTransform(refregions, CRSobj = PROJ)
countries  <- getMap("less islands")
countries  <- spTransform(countries,  CRSobj = PROJ)
refrenceSubset <- subset(refregions, Acronym %in% regNames)
basemap <- ggplot() +
  #Set Border around plot - probably not the best way to do this
  borders(aggregate(refregions, FUN=length), fill=NA, colour='black', size=3) +
  geom_map(data=refregions, map=fortify(refregions),
           aes(x=long, y=lat, group=group, map_id=id), 
           fill="white", colour="white", size=1)+
  #Add Country data (basemap)
  geom_map(data=countries, map=fortify(countries),
           aes(x=long, y=lat, group=group, map_id=id), 
           fill = "grey80",color="grey90",size=0.2) +
  borders(database = refrenceSubset, fill=NA, colour='grey20') + 
  coord_fixed(1) + 
  theme_void() 

pcPlt <- ggdraw() 
for (pc in c(1,2)){
  dataTable <- fortify(refrenceSubset)
  dataTable$Count     <- NA
  for (reg in regNames){
    i = which(refregions@data[["Acronym"]]==reg)
    loading <- apply(pca[["loadings"]][,pc,],1,median,na.rm=TRUE)[which(regNames==reg)]
    dataTable$Count[which(dataTable$id==as.character(i-1))] <- loading
  }
  pcMapplot <- basemap + 
    #Add refrence regions boundaries
    geom_map(data=refrenceSubset, map=dataTable, alpha=0.75, size=0.5, color='black' , 
             aes(x=long, y=lat, group=group, map_id=id,fill=dataTable$Count)) +
    scale_fill_binned(type = "viridis",name =paste('Loading (PC',pc,')',sep='')) + 
    theme(text = element_text(family='sans',size=12),
          plot.background = element_rect(fill = 'white',color='White'),
          plot.margin = unit(c(0.25, 0.25, 0.2, 0.2), "in"),
          legend.position = c(0.5,0.12),
          legend.direction='horizontal',
          legend.box.background=element_rect(fill = 'white',color='Black'),
          legend.box.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.key.width=unit(20, 'points'),
          legend.key.height=unit(4, 'points'))
  pcTSplot <- ggplot() + 
    geom_line(aes(x=pca[["age"]]/1000,y=apply(pca[["PCs"]][,pc,],1,median,na.rm=TRUE)),
              color='Black',size=1)+
    labs(x='Age (yr BP)',y='Dry <    > Wet',subtitle=paste('PC',pc,sep='')) +
    theme_bw() +
    scale_x_reverse(name = "Age (ka BP)", 
                    limits=c(12,0), expand=c(0,0), n.breaks=7) +
    theme(text = element_text(family='sans',size=8),
          plot.background = element_rect(fill = 'white',color='White'),
          plot.margin = unit(c(0.25, 0.25, 0.2, 0.2), "in"))
  pcPlt <- pcPlt + 
    draw_plot(pcMapplot, x = 0,   y = ((2-pc)/2), width = 0.6, height = 0.5) + 
    draw_plot(pcTSplot,  x = 0.6, y = ((2-pc)/2), width = 0.4, height = 0.5)
}
pcPlt 
  
ggsave(plot=pcPlt, width = 6.5, height = 4.3, dpi = 600,
       filename = paste(file.path(dataDir,'Figures',paste('PCA_',climVar,'.png',sep=''))))
