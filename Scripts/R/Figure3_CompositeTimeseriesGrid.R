#Load Packages----
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
dir <-  getwd()#'/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate'
var <- 'HC'

basemapMercator <- ggplot() +
  #Set Border around plot - probably not the best way to do this
  borders(aggregate(IPCC_WGI_reference_regions_v4, FUN=length), fill='white', colour='white', size=2) +
  geom_map(data=IPCC_WGI_reference_regions_v4, map=fortify(IPCC_WGI_reference_regions_v4),
           aes(x=long, y=lat, group=group, map_id=id), fill="#F4F7FF", colour='NA', size=1)+
  #Add Country data (basemap)
  geom_map(data=rworldmap::getMap("less islands"), map=fortify(rworldmap::getMap("less islands")),
           aes(x=long, y=lat, group=group, map_id=id), fill = "grey80",color="grey90",size=0.2) +
  coord_fixed(1) + 
  theme_void()  

basemapMercator

#Load Data----
regionData <- readRDS(file.path(dir,'Data','FigureSettings','regionData.rds'))
binvec     <- read.csv(file.path(dir,'Data','RegionComposites',var,'MedianTS_byRegion.csv'))$time
plotSettings <- plotSettings #From Figure_Settings.Rmd
map          <- basemapMercator      #From Figure_Settings.Rmd
Csettings    <- Csettings    #From Figure_Settings.Rmd
alph<-1

sample = 1
if (1==1){
  if (sample == 1){regNames <- c('NEN','GIC','NWN','WNA','CNA','ENA','NEU','WSB','WCE','ESB','MED','RFE','SAU','NZ','SSA','EAN')
                   position <- c(letters[(length(regNames)-1):26],'aa','ab','ac','ad')
  } else{          regNames <- c('WCA','ECA','TIB','EAS','NCA','SAS','SCA','SEA','NWS','SAH','SAM','NEAF','NES','SEAF','WSAF','ESAF')
                   position <- c(letters[1:length(regNames)])
  }
} else{
  if (sample == 1){regNames <- regnames[1:(length(regnames)/2)]
  } else{          regNames <- regnames[((length(regnames)/2)+1):length(regnames)]}
}

figHeight <- (8/9)*length(regNames)/2
  

regPlts <- vector(mode='list')

for (reg in regNames){
  #Load Data for Region
  regTso   <- regionData[[reg]][[var]][["LiPD"]]
  #Standardize Ensemble Composite Values
  regEnsNA <- read.csv(file.path(dir,'Data','RegionComposites',var,paste(reg,'.csv',sep='')))
  regEnsNA <- as.matrix(regEnsNA - as.numeric(apply(regEnsNA,2,mean,na.rm=TRUE))) / as.numeric(apply(regEnsNA,2,sd,na.rm=TRUE))
  #Create matrix for portion of timeseries with >50% data coverage
  regEns <- matrix(NA,nrow(regEnsNA),ncol(regEnsNA))
  regEns[regionData[[reg]][[var]][["pltTimeAvail50range"]],] <- regEnsNA[regionData[[reg]][[var]][["pltTimeAvail50range"]],]
  regPlt <- ggdraw(ggplot()+theme_void()+theme(plot.background= element_rect(colour='White',fill='White')))
  #Plot timeseries for region
  compBands <- vector(mode = 'list')
  compBands$na <-  plotTimeseriesEnsRibbons(ggplot()+geom_hline(yintercept=0,size=0.05,color='black'),
                                            X=binvec, Y=regEnsNA, alp=alph,line.width=0.1,
                                            color.low='grey95',
                                            color.high='grey80',
                                            color.line='grey40')
  compBands$ts <- plotTimeseriesEnsRibbons(X=binvec, Y=regEns, alp=alph-0.2,line.width=0.1,
                                           color.low=Csettings[1],
                                           color.high=Csettings[2],
                                           color.line=Csettings[3])
  for (plt in names(compBands)){
    compBands[[plt]] <- compBands[[plt]] + 
      geom_hline(yintercept=0,size=0.2,color='black') +
      scale_x_reverse(limits=c(12000,0), expand=c(0,0), n.breaks=7)+ 
      scale_y_continuous(limits=c(-1000,1000),breaks=seq(-4,4,2),labels=c(-4,'',0,'',4),position="right", expand=c(0,0))+
      coord_cartesian(xlim=c(12000,0), ylim=c(-4.5,4.5)) +
      ggtitle(paste('(',position[which(regNames==reg)],') ',regionData[[reg]]$name,' (',reg,')',sep='')) + 
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
            axis.title      = element_blank(),
            panel.grid      = element_blank(),
            plot.margin     = unit(c(0.05, 0.05, 0, 0.05), "in"),
            text = element_text(family='sans',size=8),
            legend.position='none')
  }
  compBands$na <- compBands$na +
    theme(plot.background= element_rect(colour='White',fill='White'),
          panel.border    =element_rect(colour='Black',fill=NA))
  #Load Data for Region
  regionDf <- regionData[[reg]]$HC$SummaryDF
  regionDf <- as.data.frame(spTransform(SpatialPointsDataFrame(regionDf[,c("longitude", "latitude")], regionDf, proj4string=CRS(PROJorig)), CRSobj = PROJ))
  idx      <- which(plotSettings$names %in% sort(unique(regionDf$CategorySpec)))
  #RegionMap
  RegShp   <- IPCC_WGI_reference_regions_v4[IPCC_WGI_reference_regions_v4@data$Acronym == reg, ]#regionData[[reg]]$polygon
  latrange <- range(fortify(RegShp)$lat)
  lonrange <- range(fortify(RegShp)$long)
  range <- max(diff(latrange),diff(lonrange))/2
  range <- range*1.2
  if(diff(latrange)>diff(lonrange)){
    latrange <- mean(latrange) + c(range*-1,range)*1.2
    lonrange <- mean(lonrange) + c(range*-1,range)*1.8
  } else{
    latrange <- mean(latrange) + c(range*-1,range)*0.8
    lonrange <- mean(lonrange) + c(range*-1,range)*1.2
  }
  #if (range<20){range <- 20}
  if (lonrange[1]< -180){
    lonrange[2] <- lonrange[2] + (abs(lonrange[1])-180)
    lonrange[1] <- -180
  }
  if (lonrange[2]> 180){
    lonrange[1] <- lonrange[1] - (abs(lonrange[2])-180)
    lonrange[2] <- 180
  }
  if (latrange[1]< -90){
    latrange[2] <- latrange[2] + (abs(latrange[1])-90)
    latrange[1] <- -90
  }
  if (latrange[2]> 90){
    latrange[1] <- latrange[1] - (abs(latrange[2])-90)
    latrange[2] <- 90
  }
  if (reg == 'EAN'){
    lonrange <- c(-120,150)
    latrange <- c(-90,90)
  }
  regMap   <-  map +
    geom_map(data=RegShp, map=fortify(RegShp),fill=NA, alpha=0.75, size=0.35, color='black',
             aes(x=long, y=lat, group=group, map_id=id)) +
    geom_star(data=regionDf,size=1.5,color='Black',alpha=1,starstroke=0.4,
              aes(x=longitude,y=latitude,starshape=CategorySpec,fill=CategorySpec)) +####lat lons to change when prj
    coord_fixed(xlim=lonrange,ylim=latrange)+
    scale_fill_manual(values=plotSettings$color[idx],name= 'Proxy Category') +
    scale_starshape_manual(values=plotSettings$shape[idx],name= 'Proxy Category') +
    theme_void() + 
    theme(panel.border    = element_rect(colour='Black',fill=NA,size=0.75),
          plot.background = element_rect(colour='White',fill='White'),
          panel.background = element_rect(colour='Black',fill='White'),
          plot.margin     = unit(c(0.05, 0, 0.05,0.05), "in"),
          legend.position = 'none') 
  #Plot Time Availability for region
  if(length(regTso)<10){labs <- c('0',paste('0',length(regTso),sep=''))
  } else{labs <- c(0,length(regTso))}
  pltTime <- plotTimeAvailabilityTs(regTso,age.range=c(0,12000),group.var ='CategorySpecific',step=100)
  pltTime <- ggplot(pltTime$data,aes(yvec,value))+
    geom_area(aes(fill=group),color='Black',size=0.2) +
    scale_fill_manual(values=plotSettings$color[idx],name= 'Proxy Category') +
    scale_x_reverse(limits=c(12000,0),expand=c(0,0),n.breaks=7,labels=seq(0,12,2))+
    scale_y_continuous(limits=c(0,max(length(regTso)*1.2),expand=c(0,0)),
                       labels=labs,breaks=c(0,length(regTso)),expand=c(0,0))+
    theme_bw()+ 
    theme(panel.background=element_rect(colour='Black',fill=NA),
          panel.border= element_rect(colour='Black',color=,fill=NA),
          panel.grid  = element_blank(),
          axis.title  = element_blank(),
          axis.ticks  = element_line(color = 'Black',size=0.4), 
          axis.text.x = element_blank(),
          axis.text.y =element_text(family='sans',size=8,color='Black'),
          axis.ticks.length.y=unit(2,"pt"),
          axis.ticks.length.x=unit(3,"pt"),
          plot.margin = unit(c(0, 0.05, 0.05, 0.05), "in"),
          text = element_text(family='sans',size=8),
          legend.position='none')
  #
  regPlts[[reg]] <- ggdraw(ggplot() + theme(plot.background= element_rect(colour='White',fill='White'),
                                            panel.background = element_rect(colour='White',fill='White')))+
    draw_plot(ggarrange(compBands$na, pltTime, nrow = 2,heights=c(0.7,0.3)), x = 0, y = 0, width = 0.7, height = 1) +
    draw_plot(ggarrange(compBands$ts, pltTime, nrow = 2,heights=c(0.7,0.3)), x = 0, y = 0, width = 0.7, height = 1) + 
    draw_plot(regMap,  x = 0.7,   y = 0, width = 0.28, height = 0.9)
}

scale <- ggplot() + geom_point(aes(x=0,y=0),size=0,color='white') +
  scale_x_reverse('Age (ka BP)',limits=c(12,0),expand=c(0,0),n.breaks=7)+
  scale_y_continuous(limits=c(-4,10), breaks=seq(0,10,10),sec.axis = sec_axis(~ .,breaks=seq(-4,4,4))) +
  theme_bw()+ 
  theme(panel.background=element_rect(colour='White',fill='White'),
        panel.border    =element_rect(colour='White',fill='White'),
        panel.grid.major=element_line(colour='White'),
        axis.title.y  =  element_blank(),
        axis.text.y =element_text(family='sans',size=8,color='White'),
        axis.line.x = element_line(color = 'black'),
        axis.ticks.x  = element_line(color = 'Black',size=0.4), 
        axis.ticks.y  = element_line(color = 'White',size=0.4), 
        axis.text = element_text(family='sans',size=8),
        axis.title.x = element_text(family='sans',size=8),
        axis.ticks.length.y=unit(2,"pt"),
        axis.ticks.length.x=unit(3,"pt"),
        plot.margin = unit(c(0, 0.05, 0.05, 0.05), "in"),
        text = element_text(family='sans',size=8),
        legend.position='none')
scale <- ggdraw(ggplot() + theme(plot.background= element_rect(colour='White',fill='White'),
                          panel.background = element_rect(colour='White',fill='White')))+
  draw_plot(scale, x = 0, y = 0, width = 0.7, height = 1) 


#Save Summary Plot-----
plt <- ggplot() + theme(plot.background = element_rect(colour='White',fill='White'),
                        panel.background = element_rect(colour='White',fill='White'))
plt <- ggdraw(plt) 
v <- 1
dv <- 2*v/(length(regNames)+1)
for (i in 1:length(regNames)){
  print(i)
  h <- (i%%2-1)/-2
  if (i%%2 == 1){v <- v-dv}
  plt <- plt + draw_plot(regPlts[[regNames[[i]]]], x = h, y = v, width = 0.5, height = dv)
}
for (i in 1:2){
  h <- (i%%2-1)/-2
  if (i%%2 == 1){v <- v-dv/2}
  plt <- plt + draw_plot(scale, x = h, y = v, width = 0.5, height = dv/2)
}
ggsave(plot=plt, width = 6.5, height = figHeight, dpi = 600,
       filename = file.path(dir,'Figures',paste('compositeSummaryGrid_',var,'_',sample,'.png',sep='')))

#-----