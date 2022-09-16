#Purpose--------------------------------------------------------------------------------
#goal: create regional composite based on IPCC regions using either temp or HC
#in:   LiPD files
#out:  individual csv for each region where columns represent iterations. 
#      summary table with each column as the regional median time series

#Load Packages--------------------------------------------------------------------------------

library(compositeR)
library(doParallel)
library(dplyr)
library(foreach)
library(geoChronR)
library(lipdR)
library(magrittr)
library(purrr)
library(tidyverse)
library(cowplot)
library(dplyr) #and dplyr for data.frame manipulation
library(egg)
library(geoChronR) #for plotting mostly
library(ggplot2)
library(ggstar)
library(lipdR) #to read and interact with LiPD data
library(maptools)
library(proj4)
library(rworldmap)
library(sp)

#Set up directories and names--------------------------------------------------------------------------------

dir  <- getwd()# '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate' #
var  <- 'HC'
save <- FALSE
saveDir <- file.path(dir,'Data','RegionComposites',var)


#Load Data without winter+ or summer+ seasonality--------------------------------------------------------------------------------

lipdData <- readRDS(file.path(dir,'Data','Proxy','LiPD','lipdData.rds'))[[var]]
lipdData <- lipdData[which(pullTsVariable(lipdData,'geo_ipccRegion') %in% c('NCA','WNA','ENA','CNA','NWN','NEN'))]  
lipdData <- lipdData[which(between(pullTsVariable(lipdData,'geo_latitude'),30,50))]



lipdTSO  <- vector(mode='list')
lipdTSO$'Annual (All)' <- lipdData[which(pullTsVariable(lipdData,"climateInterpretation1_seasonalityGeneral") %in% c('Annual'))]
lipdTSO$'Annual (Pollen)'    <- lipdTSO$'Annual (All)'[which(pullTsVariable(lipdTSO$'Annual (All)','Category') == 'Pollen')]
lipdTSO$'Annual (Shoreline)' <- lipdTSO$'Annual (All)'[which(pullTsVariable(lipdTSO$'Annual (All)','Category') == 'Shoreline')]
lipdTSO$'Annual (Other)'     <- lipdTSO$'Annual (All)'[which(pullTsVariable(lipdTSO$'Annual (All)','Category') %in% c('Pollen','Shoreline')==FALSE)]
lipdTSO$Summer <- lipdData[which(pullTsVariable(lipdData,"climateInterpretation1_seasonalityGeneral") %in% c('Summer','Summer+'))]
lipdTSO$Winter <- lipdData[which(pullTsVariable(lipdData,"climateInterpretation1_seasonalityGeneral") %in% c('Winter','Winter+'))]
for (i in 1:length(lipdTSO)){
  print(length(lipdTSO[[i]]))
}
  

if(var == 'T'){
  lipdTSO <- filterTs(lipdTSO,'paleoData_units == degC')
  lipdTSO <- filterTs(lipdTSO,'paleoData_datum == abs')
  std <- FALSE      #Use calibrated data so no need to normalize variance
} else{std <- TRUE} #Normalize variance because data recorded with different units


#Set variables for composite code--------------------------------------------------------------------------------

nens          <- 200     #Ensemble numbers (lower = faster)
binsize       <- 100     #years (median resolution = 107yrs)
ageMin        <- 0       #age BP
ageMax        <- 12400   #age BP
searchDur     <- 3500    #yrs (for 3 lake deposit data points)
nThresh       <- 4       #minimum no. of records, else skip 

#Set bin vectors
binvec   <- seq(ageMin-binsize/2, to = ageMax+binsize/2, by = binsize)
binYears <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))


names <- names(lipdTSO)

#Calculate reconstructions--------------------------------------------------------------------------------

#Set up data to add once regional composite is calculated
compositeEnsemble <- vector(mode='list')
medianCompositeTS <- data_frame(time=binYears[1:which(binYears==12000)])

#Loop to composite (by region)
for (reg in names) {
  print(reg)
  lipdReg  <- lipdTSO[[reg]]
  for (i in 1:length(lipdReg)){
    if (lipdReg[[i]]$climateInterpretation1_interpDirection == 'negative'){
      lipdReg[[i]]$paleoData_values <- lipdReg[[i]]$paleoData_values*-1
    }
  }
  set.seed(5) #Reproducibility
  ensOut <- foreach(i = 1:nens) %dopar% {
    tc <- compositeEnsembles(fTS                  = lipdReg,
                             binvec               = binvec,
                             stanFun              = standardizeMeanIteratively,
                             binFun               = simpleBinTs,
                             ageVar               = "age",
                             alignInterpDirection = FALSE,
                             spread               = TRUE,
                             duration             = searchDur,
                             searchRange          = c(1000,8000),
                             normalizeVariance    = std,
                             minN                 = 3)
    return(list(composite = tc$composite,count = tc$count))
  }
  regionComposite           <- as.matrix(purrr::map_dfc(ensOut,magrittr::extract,"composite"))
  rownames(regionComposite) <- binYears
  regionComposite           <- regionComposite[1:which(binYears==12000),]
  compositeEnsemble[[reg]]  <- regionComposite
  medianCompositeTS[[reg]]  <- apply(regionComposite,1,median,na.rm=TRUE)
  #plot region to confirm that everything looks good
  plotTimeseriesEnsRibbons(X = binYears[1:which(binYears==12000)],Y = compositeEnsemble[[reg]])+
    scale_x_continuous(name = "age (yr BP)",         oob = scales::squish)+
    scale_y_continuous(name = "Standardized Anomaly",oob = scales::squish)+
    theme_bw()+
    ggtitle(paste(reg,"Composite Ensemble"))
}

if(save){
  write.csv(medianCompositeTS, row.names=FALSE, file=file.path(saveDir,'MedianTS_byRegion.csv'))
  for (reg in names(compositeEnsemble)){
    write.csv(compositeEnsemble[[reg]], row.names=FALSE,file = file.path(saveDir,paste(reg,'.csv',sep='')))
  }
  print("csv files saved")
}

binvec<-binvec[1:which(binYears==12000)]
colors <- vector(mode='list')
colors$Summer <- RColorBrewer::brewer.pal(9,'Reds')[c(2,5,9)]
colors$Winter <- RColorBrewer::brewer.pal(9,'Blues')[c(2,5,9)]
colors$'Annual (All)'      <- RColorBrewer::brewer.pal(9,'Greys')[c(2,5,9)]
colors$'Annual (Pollen)'  <- RColorBrewer::brewer.pal(9,'Greens')[c(2,5,9)]
colors$'Annual (Shoreline)'<- RColorBrewer::brewer.pal(9,'BuPu')[c(2,5,9)]
colors$'Annual (Other)'     <- RColorBrewer::brewer.pal(9,'RdPu')[c(2,5,9)]

#RColorBrewer::display.brewer.pal(9,'BuPu')

pltVector <- vector(mode='list')
for (name in names){
  sznTso   <- lipdTSO[[name]]
  #Standardize Ensemble Composite Values
  composite <-compositeEnsemble[[name]]
  composite <- as.matrix(composite - as.numeric(apply(composite[which(between(binvec,0,1000)),],2,mean,na.rm=TRUE)))
  #Create matrix for portion of timeseries with >50% data coverage
  plt <- ggdraw(ggplot()+theme_void()+theme(plot.background= element_rect(colour='White',fill='White')))
  #Plot timeseries for region
  compBands <- vector(mode = 'list')
  plotlim <- c(-4,2)
  labs <- max(100,length(sznTso))
  labs <- c(0,round(labs/2),labs) 
  compBands <-  plotTimeseriesEnsRibbons(ggplot()+geom_hline(yintercept=0,size=0.05,color='black'),
                                         X=binvec, Y=composite, alp=0.9,line.width=0.1,
                                         color.low = colors[[name]][[1]],
                                         color.high= colors[[name]][[2]],
                                         color.line= colors[[name]][[3]]) + 
    geom_hline(yintercept=0,size=0.2,color='black') +
    scale_x_reverse(limits=c(12000,0), expand=c(0,0), n.breaks=7)+ 
    scale_y_continuous(name=name,limits=c(-1000,1000),breaks=seq(plotlim[1],plotlim[2],2),labels=seq(plotlim[1],plotlim[2],2),position="left", expand=c(0,0),
                       sec.axis = sec_axis( trans=~(.+plotlim[1]*-1)*labs[3]/diff(plotlim), name="Count",breaks=labs,labels=labs))+
    coord_cartesian(xlim=c(12000,0), ylim=c(plotlim[1],plotlim[2])) +
    #ggtitle(name) + 
    theme_bw() +
    theme(panel.background=element_rect(colour='Black',fill=NA),
          panel.border    =element_rect(colour='Black',fill=NA),
          plot.background =element_rect(colour='White',fill=NA),
          axis.ticks      =element_line(color = 'black', size=0.4), 
          axis.ticks.length.x = unit(-3,"pt"),
          axis.ticks.length.y = unit(2,"pt"),
          axis.text.x = element_blank(),
          axis.text.y =element_text(family='sans',size=8,color='Black'),
          plot.title = element_text(hjust = 0,vjust=-0.3,family='sans',size=8),
          plot.caption.position = "plot",
          axis.title.x      = element_blank(),
          panel.grid      = element_blank(),
          plot.margin     = unit(c(0.05, 0.05, 0.05, 0.05), "in"),
          text = element_text(family='sans',size=8),
          legend.position='none') +
          theme(plot.background= element_rect(colour='White',fill=NA),
          panel.border    =element_rect(colour='Black',fill=NA))
  
  #Load Data for Region
  sznDf <- proxyDf[which(proxyDf$tsid %in% pullTsVariable(sznTso,"paleoData_TSid")),]
  sznDf <- as.data.frame(spTransform(SpatialPointsDataFrame(sznDf[,c("longitude", "latitude")], sznDf, proj4string=CRS(PROJorig)), CRSobj = PROJ))
  idx      <- which(plotSettings$names %in% sort(unique(sznDf$CategorySpec)))
  #RegionMap
  RegShp   <- IPCC_WGI_reference_regions_v4[IPCC_WGI_reference_regions_v4@data$Acronym %in% c('ENA','CNA','WNA','NWN','NEN'), ]#regionData[[reg]]$polygon
  regMap   <-  map +
    geom_star(data=sznDf,size=1.5,color='Black',alpha=1,starstroke=0.4,
              aes(x=longitude,y=latitude,starshape=CategorySpec,fill=CategorySpec)) +####lat lons to change when prj
    coord_fixed(xlim=c(-125,-55),ylim=c(30,50))+
    #geom_hline(yintercept=30)+
    #geom_hline(yintercept=50)+
    #geom_rect(aes(xmin=-180,xmax=0,ymin=0,ymax=30),fill='grey',alpha=0.7) +
    #geom_rect(aes(xmin=-180,xmax=0,ymin=50,ymax=90),fill='grey',alpha=0.7) +
    scale_fill_manual(values=plotSettings$color[idx],name= 'Proxy Category') +
    scale_starshape_manual(values=plotSettings$shape[idx],name= 'Proxy Category') +
    theme_void() + 
    theme(panel.border    = element_rect(colour='Black',fill=NA,size=0.75),
          plot.background = element_rect(colour='White',fill='White'),
          panel.background = element_rect(colour='Black',fill='White'),
          plot.margin     = unit(c(0.05, 0.05, 0.05, 0.05), "in"),
          legend.position = 'none') 
  #Plot Time Availability for region
  #if(length(sznTso)<10){labs <- c('0',paste('0',length(sznTso),sep=''))
  #} else{labs <- c(-3,)}
  pltTime <- plotTimeAvailabilityTs(sznTso,age.range=c(0,12000),group.var ='CategorySpecific',step=100)
  pltTime <- ggplot(pltTime$data,aes(yvec,value))+
    geom_area(aes(fill=group),color=NA,alpha=0.3,size=0.2) +
    scale_fill_manual(values=rep('grey',length(plotSettings$color[idx])),name= 'Proxy Category') +
    scale_x_reverse(limits=c(12000,0), expand=c(0,0), n.breaks=7)+ 
    scale_y_continuous(limits=c(0,labs[3],expand=c(0,0)),
                       labels=c(-3,0,3),breaks=c(0,length(lipdData)/2,length(lipdData)),expand=c(0,0),name=name,
                       sec.axis = sec_axis( trans=~(.+plotlim[1]*-1)*labs[3]/diff(plotlim), name="Count",breaks=labs,labels=labs))+
    coord_cartesian(xlim=c(12000,0), ylim=c(0,labs[3])) +
    theme_bw()+ 
    theme(panel.background=element_rect(colour='Black',fill=NA),
          panel.border= element_rect(colour='Black',fill=NA),
          panel.grid  = element_blank(),
          axis.title.x  = element_blank(),
          axis.ticks  = element_line(color = 'Black',size=0.4), 
          axis.text.x = element_blank(),
          axis.text.y =element_text(family='sans',size=8,color='White'),
          axis.ticks.length.y=unit(2,"pt"),
          axis.ticks.length.x=unit(-3,"pt"),
          plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "in"),
          text = element_text(family='sans',size=8),
          legend.position='none')
  
  pltVector[[name]] <- ggdraw(ggplot() + theme(plot.background= element_rect(colour='White',fill='White'),
                                            panel.background = element_rect(colour='White',fill='White')))+
    draw_plot(pltTime,  x = 0,   y = 0, width = 0.7, height = 1)+
    draw_plot(compBands,  x = 0,   y = 0., width = 0.7, height = 1)+
    draw_plot(regMap,  x = 0.7,   y = 0, width = 0.3, height = 1)
}

z <- ggdraw(ggplot() + theme(plot.background= element_rect(colour='White',fill='White'),
                        panel.background = element_rect(colour='White',fill='White')))+
  draw_plot(pltVector[['Summer']],            x = 0,   y = 5/6, width = 1, height = 1/6)+
  draw_plot(pltVector[['Winter']],            x = 0,   y = 4/6, width = 1, height = 1/6)+
  draw_plot(pltVector[['Annual (All)']],      x = 0,   y = 3/6, width = 1, height = 1/6)+
  draw_plot(pltVector[['Annual (Pollen)']],   x = 0,   y = 2/6, width = 1, height = 1/6)+
  draw_plot(pltVector[['Annual (Shoreline)']],x = 0,   y = 1/6, width = 1, height = 1/6)+
  draw_plot(pltVector[['Annual (Other)']],    x = 0,   y = 0/6, width = 1, height = 1/6)

ggsave(plot=z, width = 6.5, height = 7, dpi = 600,
       filename = file.path(dir,'Figures','RegionComposites','NA_compositeSummaryStack.png'))

