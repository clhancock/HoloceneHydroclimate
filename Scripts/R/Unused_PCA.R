library(geoChronR)
library(cowplot)
#githubDir <- getwd()
var <- 'HC'
compositeData <- read.csv(file.path(dir,'Data','RegionComposites',
                                    var,'MedianTS_byRegion.csv'))
lipdTSO <- readRDS(file=file.path(dir,'Data','Proxy','LiPD','lipdData.rds'))[[var]]

#Calculate data to plot and connect it to region data. 
#More complicated than needed because part of larger code to make series of figures for each region 
pcaData <-vector(mode='list')
for (reg in names(compositeData)[-1]){ #First column is the time variable
  pcaData[[reg]]           <- vector(mode="list")
  pcaData[[reg]][['Name']] <- reg
  #Load region composite matrix csv
  pcaData[[reg]]$compMatix <- read.csv(file.path(dir,'Data',
                                                 'RegionComposites',var,
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
                               seq(0,10000,100), #bin only over this range change to 8000 to make work
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

plotScreeEns(pcaOut)
pcs <- 2
pcPlt <- ggdraw() 
for (pc in 1:pcs){
  dataTable <- fortify(regionsSelect$composite)
  dataTable$Count     <- NA
  for (reg in regNames){
    i = which(refregions@data[["Acronym"]]==reg)
    loading <- apply(pca[["loadings"]][,pc,],1,median,na.rm=TRUE)[which(regNames==reg)]
    dataTable$Count[which(dataTable$group==paste(reg,".1",sep=""))] <- loading
  }
  pcMapplot <- basemap + 
    #Add refrence regions boundaries
    geom_map(data=regionsSelect$composite , map=dataTable, alpha=0.75, size=0.5, color='black' , 
             aes(x=long, y=lat, group=group, map_id=id,fill=dataTable$Count)) +
    scale_fill_binned(type = "viridis",name =paste('Loading (PC',pc,')',sep='')) + 
    theme(text = element_text(family='sans',size=8),
          plot.background = element_rect(fill = 'white',color='White'),
          plot.margin = unit(c(0.25, 0.25, 0.2, 0.2), "in"),
          legend.position = c(0.5,0),
          legend.direction='horizontal',
          legend.box.background=element_rect(fill = 'white',color='Black'),
          legend.box.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.key.width=unit(30, 'points'),
          legend.key.height=unit(4, 'points'))
  pcTSplot <-   plotTimeseriesEnsRibbons(X=pca[["age"]]/1000,Y=pca[["PCs"]][,pc,],
                                         probs = c(0.005, 0.1, 0.5, 0.8, 0.995),
                                         line.width = 0.4)+
    labs(x='Age (yr BP)',y='',subtitle=paste('PC',pc,sep='')) +
    theme_bw() +
    scale_x_reverse(name = "Age (ka BP)", 
                    limits=c(10,0), expand=c(0,0), n.breaks=7) +
    theme(text = element_text(family='sans',size=8),
          plot.background = element_rect(fill = 'white',color='White'),
          plot.margin = unit(c(0.25, 0.25, 0.2, 0.2), "in"),
          panel.grid = element_blank())
  pcPlt <- pcPlt + 
    draw_plot(pcTSplot,  x = 0.55, y = ((pcs-pc)/pcs), width = 0.45, height = 1/pcs)+
    draw_plot(pcMapplot, x = 0,   y = ((pcs-pc)/pcs), width = 0.6, height = 1/pcs) 
}
pcPlt 
ggsave(plot=pcPlt, width = 6.5, height = 4.3, dpi = 600,
       filename = paste(file.path(dir,'Figures',paste('PCA_',var,'.png',sep=''))))
