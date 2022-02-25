#Load Packages----
library(lipdR) #to read and interact with LiPD data
library(geoChronR) #for plotting mostly
library(dplyr) #and dplyr for data.frame manipulation
library(ggplot2)
library(cowplot)
library(ggstar)
library(maptools)
library(proj4)
library(sp)
dataDir <- '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate'
var     <- 'HC'
project = FALSE
#Load Data----
#Projections
PROJ     <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
PROJorig <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
#Refrence Region Shapefiles
refReg   <- readShapePoly(file.path(dataDir,'Data','IPCC_refRegions','IPCC-WGI-reference-regions-v4.shp'),
                          proj4string=CRS(PROJorig))
if (project){refReg   <-  spTransform(refReg, CRSobj = PROJ)}
#Countries for basemap
countries  <- getMap("less islands")
if (project){countries  <- spTransform(countries,  CRSobj = PROJ)}
#Load Proxy Information
regPlts  <- read.csv(file.path(dataDir,'Data','RegionComposites',var,'MedianTSbyRegion.csv'))
regnames <- as.character(refReg@data[["Acronym"]])[which(refReg@data[["Acronym"]]%in%names(regPlts))]
regPlts  <- vector(mode='list')
lipdTso  <- readRDS(file.path(dataDir,'Data','LiPD','lipdData.rds'))[[var]]

#Figure Settings----
#Climate variable Settings for HC vs T
settings    <- vector(mode='list')
settings$T  <- vector(mode='list')
settings$HC <- vector(mode='list')
settings$T$Pal  <- 'RdBu'
settings$T$Dir  <- 01
settings$T$Col  <- c("#FFEBEE","#FFCDD2","#EF9A9A") #reds
settings$T$Lab  <- "Cool <- (Standardized Anomoly) -> Warm"
settings$T$NAv  <- 'light grey'
settings$T$Seq  <- seq(-2.2,2.2,0.4)
settings$HC$Pal <- 'BrBG'
settings$HC$Dir <- 01
settings$HC$Col <- c("#E8EAF6","#303F9F","#1A237E") #indigo
settings$HC$Lab <- "Dry <- (Standardized Anomoly) -> Wet"
settings$HC$Lim <-  c(-2.5,2.5)
settings$HC$NAv <- 'light grey'
settings$HC$Seq <- seq(-2.2,2.2,0.4)
#Proxy Value Settings
CatColor <- c("powder blue","corn flower blue","dark blue","dark orchid",
              "grey40","grey",
              "forest green","yellowgreen",
              "firebrick","lightcoral","darkorange")
CatNames <- c("Glacier Ice","Lake Deposits",'Lake Sediment (δ18O)','Leaf Wax (δD)',
              'Other (calibrated)','Other (not calibrated)',
              'Pollen (calibrated)','Pollen (not calibrated)',
              'Speleothem (δ18O)','Speleothem (δ13C)','Speleothem (other)')
CatShape <- c(12,21,15,5, 6,13, 14,1, 23,11,17)

#Create Plots-----
sample = 1
if (var == 'HC'){
  if (sample == 1){regNames <- c('NWN','WNA','NEN','CNA','GIC','ENA','NEU','ESB','WCE','ECA','MED','WCA','SAU','NZ')
  } else{          regNames <- c('NCA','SCA','SAH','EAS','NEAF','TIB','SEAF','SAS','ESAF','SEA','WSAF','NES','NWS','SAM')}
} else{
  if (sample == 1){regNames <- regnames[1:(length(regnames)/2)]
  } else{          regNames <- regnames[((length(regnames)/2)+1):length(regnames)]}
}

for (reg in regNames){
  #Load Data for Region
  regTso <- filterTs(lipdTso,paste('geo_ipccRegion ==',reg))
  name   <- as.character(refReg@data[["Name"]])[which(refReg@data[["Acronym"]]==reg)]
  df <- tibble(dataset       = pullTsVariable(regTso,'dataSetName'),
               tsid          = pullTsVariable(regTso,'paleoData_TSid'),
               longitude     = pullTsVariable(regTso,'geo_longitude'),
               latitude      = pullTsVariable(regTso,'geo_latitude'),
               ipccReg       = pullTsVariable(regTso,'geo_ipccRegion'),
               Category      = pullTsVariable(regTso,'Category'),
               CategorySpec  = pullTsVariable(regTso,'CategorySpecific'))
  df <- SpatialPointsDataFrame(df[,c("longitude", "latitude")], df, proj4string=CRS(PROJorig))
  if (project){df <- spTransform(df, CRSobj = PROJ)}
  df <- as.data.frame(df)
  idx <- which(CatNames %in% unique(df$CategorySpec))
  regEnsNA <- read.csv(file.path(dataDir,'Data','RegionComposites',var,paste(reg,'.csv',sep='')))
  #Standardize mean at 0
  regEnsNA <- as.matrix(regEnsNA - as.numeric(apply(regEnsNA,2,mean,na.rm=TRUE)))
  regEnsNA <- regEnsNA / as.numeric(apply(regEnsNA,2,sd,na.rm=TRUE))
  regEns   <- matrix(NA,nrow(regEnsNA),ncol(regEnsNA))
  #RegionMap
  RegShp   <- subset(refReg, Acronym ==reg)
  regMap <-  ggplot() +
    geom_map(data=countries, map=fortify(countries),fill = "grey80",color="grey30",size=0.4,
             aes(x=long, y=lat, group=group, map_id=id)) +
    geom_map(data=RegShp, map=fortify(RegShp),fill=NA, alpha=0.75, size=0.5, color='black',
             aes(x=long, y=lat, group=group, map_id=id)) +
    geom_star(data=df,size=1.7,color='Black',alpha=0.8,starstroke=0.4,
              aes(x=longitude.1,y=latitude.1,starshape=CategorySpec,fill=CategorySpec)) +
    coord_fixed(xlim=range(fortify(RegShp)$lon)+80*c(-0.1,0.1),#diff(range(fortify(RegShp)$lon))
                ylim=range(fortify(RegShp)$lat)+80*c(-0.2,0.2))+
    scale_fill_manual(values=CatColor[idx],name= 'Proxy Category') +
    scale_starshape_manual(values=CatShape[idx],name= 'Proxy Category') +
    theme_void() + 
    theme(panel.border    = element_rect(colour='Black',fill=NA,size=0.75),
          plot.background = element_rect(colour='White',fill='White'),
          panel.background = element_rect(colour='Black',fill='aliceblue'),
          plot.margin     = unit(c(0.05, 0, 0.05,0.05), "in"),
          legend.position = 'none') 
  #Plot Time Availability for region
  pltT <- plotTimeAvailabilityTs(regTso,age.range=c(0,12000),group.var ='CategorySpecific',step=100)
  pltTime <- ggplot(pltT$data,aes(yvec,value))+
    geom_area(alpha=1,aes(fill=group),color='Black',size=0.2) +
    scale_fill_manual(values=CatColor[idx],name= 'Proxy Category') +
    scale_x_reverse(limits=c(12000,0),expand=c(0,0),n.breaks=7)+
    scale_y_continuous(limits=c(0,max(length(regTso)*1.25,10),expand=c(0,0)),
                       labels=c(0,length(regTso)),breaks=c(0,length(regTso)),expand=c(0,0))+
    theme_bw()+ 
    theme(panel.background=element_rect(colour='Black',fill=NA),
          panel.border    =element_rect(colour='Black',color=,fill=NA),
          panel.grid.major=element_line(colour='light Grey'),
          axis.title  = element_blank(),
          axis.ticks  = element_line(color = 'black',size=0.4), 
          axis.text.x = element_blank(),
          axis.ticks.length.y=unit(-3,"pt"),
          axis.ticks.length.x=unit(3,"pt"),
          plot.margin = unit(c(0, 0.05, 0.05, 0.05), "in"),
          text = element_text(family='sans',size=8),
          legend.position='none')
  #
  timeN <- pltT$dat %>% group_by(yvec)  %>% summarise(count=sum(value),nPct=sum(value)/length(regTso))
  if (length(which(timeN$nPct<=0.5) > 0)){
    mh <- which(timeN$yvec==6000)
    idx <- c(max(which(timeN$nPct[1:mh]<=0.5),1),
             min(which(timeN$nPct[mh:length(timeN$yvec)]<=0.5)+mh,length(timeN$yvec)))
    regEns[idx[1]:idx[2],] <- regEnsNA[idx[1]:idx[2],]                  
  } else{regEns<-regEnsNA}
  plotlimit_set <- 5
  compBands <- vector(mode = 'list')
  compBands$na <-  plotTimeseriesEnsRibbons(X=timeN$yvec, Y=regEnsNA, alp=0.8,line.width=0.2,
                                           color.low='grey90',
                                           color.high='grey50',
                                           color.line='grey20')
  compBands$ts <- plotTimeseriesEnsRibbons(X=timeN$yvec, Y=regEns, alp=0.8,line.width=0.2,
                                           color.low=settings[[var]]$Col[1],
                                           color.high=settings[[var]]$Col[2],
                                           color.line=settings[[var]]$Col[3])
  for (plt in names(compBands)){
    compBands[[plt]] <- compBands[[plt]] + 
      geom_hline(yintercept=0,size=0.2,color='black') +
      scale_x_reverse(limits=c(12000,0), expand=c(0,0), n.breaks=7)+ 
      scale_y_continuous(limits=c(plotlimit_set*c(-1000,1000)),breaks=seq(-4,4,2),labels=seq(-4,4,2))+
      coord_cartesian(xlim=c(12000,0), ylim=c(plotlimit_set*c(-1,1))) +
      ggtitle(paste(name,' (',reg,')',sep='')) + 
      theme_bw() +
      theme(panel.background=element_rect(colour='Black',fill=NA),
            panel.border    =element_rect(colour='Black',fill=NA),
            plot.background =element_rect(colour='White',fill=NA),
            axis.ticks      = element_line(color = 'black', size=0.4), 
            axis.ticks.length = unit(-3,"pt"),
            axis.text.x = element_blank(),
            plot.title = element_text(hjust = 0.5,vjust=-0.3,family='sans',size=8),
            axis.title      = element_blank(),
            panel.grid      = element_blank(),
            plot.margin     = unit(c(0.05, 0.05, 0, 0.05), "in"),
            text = element_text(family='sans',size=8),
            legend.position='none')
  }
  compBands$na <- compBands$na +
    theme(plot.background= element_rect(colour='White',fill='White'),
        panel.border    =element_rect(colour='Black',fill=NA),
        panel.grid.major=element_line(colour='light Grey',size=0.2))
  regPlts[[reg]] <- ggdraw(ggplot() + theme(plot.background= element_rect(colour='White',fill='White'),
                                            panel.background = element_rect(colour='White',fill='White')))+
    draw_plot(ggarrange(compBands$na, pltTime, nrow = 2,heights=c(0.7,0.3)), x = 0.25, y = 0, width = 0.7, height = 1) +
    draw_plot(ggarrange(compBands$ts, pltTime, nrow = 2,heights=c(0.7,0.3)), x = 0.25, y = 0, width = 0.7, height = 1) + 
    draw_plot(regMap,  x = 0,   y = 0, width = 0.25, height = 0.8)
}


#Save Summary Plot-----
plt <- ggplot() + theme(plot.background = element_rect(colour='White',fill='White'),
                        panel.background = element_rect(colour='White',fill='White'))
plt <- ggdraw(plt) 
v <- 1
dv <- 2*v/length(regNames)
for (i in 1:length(regNames)){
  h <- (i%%2-1)/-2
  if (i%%2 == 1){v <- v-dv}
  plt <- plt + draw_plot(regPlts[[regNames[[i]]]], x = h, y = v, width = 0.5, height = dv)
}
ggsave(plot=plt, width = 6.5, height = 8.5, dpi = 200,
       filename = file.path(dataDir,'Figures',paste('compositeSummaryGrid_',var,'_',sample,'.png',sep='')))

#-----