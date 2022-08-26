library(lipdR) #to read and interact with LiPD data
library(geoChronR) #for plotting mostly
library(dplyr) #and dplyr for data.frame manipulation
library(ggplot2)
library(cowplot)
library(ggstar)
library(maptools)
library(proj4)
library(rworldmap)
library(sp)

var      <- 'HC'
modelVar <- "p-e_ann"

if (var == "HC"){
  regions <- c('WCA','ECA','TIB','EAS','NCA','SAS','SCA','SEA','NWS','SAH','SAM','NEAF','NES','SEAF','WSAF','ESAF',
                'NEN','GIC','NWN','WNA','CNA','ENA','NEU','WSB','WCE','ESB','MED','RFE','SAU','NZ','SSA','EAN')
  postion <- as.character(c(paste("(",letters,")",sep=""),(paste("(a",letters[1:(length(regions)-26)],")",sep=""))))
  postion <- c(letters,"aa","ab","ac","ad","ae","af")
}


if (var == 'HC'){geo <-'all'
} else{geo<-'all'}

#Load Data----
Data <- vector(mode='list')
Data$proxy <- read.csv(file.path(dir,'Data','RegionComposites',var,'MedianTS_byRegion.csv'))
if (!is.na(modelVar)){
  for (model in c('trace','hadcm','cmip6')){
    Data[[model]]<- read.csv(file.path(dir,'Data','Model','RegionalTS',
                                       paste('regional_',modelVar,'_',model,'_',geo,'.csv',sep='')))
  }
}
regNames <-names(Data$proxy)[-1]
binvec   <- Data$proxy$time

regionData <- readRDS(file.path(dir,'Data','FigureSettings','regionData.rds'))
xSize <- 0.1
ySize <- 0.07

#Figure Settings----
#Climate variable Settings for HC vs T

Csettings <- Csettings
Chadcm <- '#DDAA33'
Ctrace <- '#BB5566'

alph <- 0.8
map <- ggdraw(basemap) #+ borders(database = regionsSelect$composite, fill=NA, colour='grey40',size=0.1)) 

#Plot---- 

for (reg in regNames){ 
  #Load Data for Region
  regTso   <- regionData[[reg]][[var]][["LiPD"]]
  regEnsNA <- read.csv(file.path(dir,'Data','RegionComposites',var,paste(reg,'.csv',sep='')))
  scaleVal <- 1
  if (var == 'HC'){scaleVal <- 30}
  if (!is.na(modelVar)){
    hadcmVals <- (Data[['hadcm']][[reg]]-mean(Data[['hadcm']][1:5,reg],na.rm=TRUE))[1:121]*scaleVal
    traceVals <- (Data[['trace']][[reg]]-mean(Data[['trace']][1:5,reg],na.rm=TRUE))[1:121]*scaleVal
    cmip6Vals <- Data[['cmip6']][[reg]]*scaleVal
  }
  #Standardize mean at 0
  regEnsNA <- as.matrix(regEnsNA - as.numeric(apply(regEnsNA[which(between(binvec,0,500)),],2,mean,na.rm=TRUE)))
  if (var == 'HC'){
    regEnsNA <- regEnsNA / as.numeric(apply(regEnsNA,2,sd,na.rm=TRUE))
    if (!is.na(modelVar)){regEnsNA <- regEnsNA * mean(sd(traceVals,na.rm=TRUE),sd(hadcmVals,na.rm=TRUE))}
  }
  regEns   <- matrix(NA,nrow(regEnsNA),ncol(regEnsNA))
  regEns[regionData[[reg]][[var]][["pltTimeAvail50range"]],] <- regEnsNA[regionData[[reg]][[var]][["pltTimeAvail50range"]],]
  plotlimit_set <- max(abs(regEns),na.rm=TRUE)*c(-1,1)
  regPlt <- ggdraw(ggplot()+theme_void()+theme(plot.background= element_rect(colour='White',fill='White')))
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
  if (!is.na(modelVar)){
    plotlimit_set <- range(c(traceVals,hadcmVals,cmip6Vals,
                             apply(regEns,1,mean,na.rm=TRUE)+apply(regEns,1,sd,na.rm=TRUE),apply(regEns,1,mean,na.rm=TRUE)-apply(regEns,1,sd,na.rm=TRUE)
    ),na.rm=TRUE) 
    plotlimit_set <- plotlimit_set + diff(range(plotlimit_set))*c(-0.1,0.1)
    #plotlimit_set <- c(max(abs(plotlimit_set))*c(-1.02,1.01))
    compBands$ts <- compBands$ts + geom_hline(yintercept=0,size=0.05,color='black') +
      geom_line(aes(x=binvec[which(between(binvec,0,12000))],y=hadcmVals),color=Chadcm,size=0.3,alpha=alph)+
      geom_line(aes(x=binvec[which(between(binvec,0,12000))],y=traceVals),color=Ctrace,size=0.3,alpha=alph)+
      geom_boxplot(aes(x=6000,y=cmip6Vals),width=1000,size=0.1,alpha=alph,
                   outlier.size=0.5,outlier.stroke = 0.15,outlier.alpha=1,outlier.colour='Black')
  } else{
    plotlimit_set <- quantile(regEns,c(0.001,0.999),na.rm=TRUE)
    plotlimit_set<- plotlimit_set+diff(range(plotlimit_set))*0.1*c(-1,1)
  }
  for (plt in names(compBands)){
    compBands[[plt]] <- compBands[[plt]] + 
      scale_x_reverse(limits=c(12100,-100), expand=c(0,0), n.breaks=7,sec.axis = dup_axis())+ 
      scale_y_continuous(limits=c(-1000,1000),breaks=seq(-100,100,2),sec.axis = dup_axis())+
      coord_cartesian(xlim=c(12000,0), ylim=c(plotlimit_set),expand	=FALSE) +
      theme_void() +
      theme(panel.border    = element_rect(color='Black',fill=NA,size=0.5),
            axis.ticks      = element_line(color='Black',size=0.1), 
            axis.ticks.length = unit(-1,"pt"),
            plot.margin       = unit(c(0, 0, 0, 0), "in"), legend.position='none')
    if (plt == "ts" & var == "HC"){
      if (abs(plotlimit_set)[1] < abs(plotlimit_set)[2]){
        yval=plotlimit_set[2]-diff(range(plotlimit_set))/5
      } else{
        yval=plotlimit_set[1]+diff(range(plotlimit_set))/5
      }
      compBands$ts <- compBands$ts +  
        geom_label(aes(x = 1000, y = yval,
                       label = paste("(",postion[which(regions==reg)],")",sep="")), 
                   fill = "white",family='sans',size = 1.5,label.padding = unit(0.05, "lines"),
                   alpha = 0.75,label.size=NA)
    }
    regPlt <- regPlt + draw_plot(compBands[[plt]], x = 0, y = 0, width = 1, height = 1)
  }
  plotLat <- regionData[[reg]][["latitude"]]
  plotLon <- regionData[[reg]][["longitude"]]
  if (is.na(regionData[[reg]]$xadjust)){regionData[[reg]]$xadjust <- 0}
  if (is.na(regionData[[reg]]$yadjust)){regionData[[reg]]$yadjust <- 0}
  map <- map + draw_plot(regPlt,width = xSize, height = ySize,
                         x = regionData[[reg]]$xadjust+(0.5-xSize/2)+0.5*plotLat/(9050504*2), 
                         y = regionData[[reg]]$yadjust+(0.5-ySize/2)+0.92*plotLon/(8625155*2), 
                         )
}




#plotTimeseriesEnsRibbons(X=timeN$yvec, Y=regEns, alp=0.7,line.width=0.1,color.low='#80DBF1',color.high='#253DA1',color.line='#000137')

scale <- ggplot() +
  scale_x_reverse('Age (ka BP)',limits=c(12,0),expand=c(0,0),breaks=seq(0,12,3),labels=c('0','','6','','12'))+
  #geom_segment(aes(x=11.5,xend=6,y=1,yend=1),size=2,color='Black') +
  #geom_segment(aes(x=11.5,xend=6,y=2,yend=2),size=2,color='Black') +
  #geom_segment(aes(x=11.5,xend=6,y=3,yend=3),size=2,color='Black') +
  geom_segment(aes(x=11.7,xend=8.5,y=1,yend=1),size=1,color=Chadcm,alpha=alph) +
  geom_segment(aes(x=11.7,xend=8.5,y=2,yend=2),size=1,color=Ctrace,alpha=alph) +
  geom_segment(aes(x=11.7,xend=8.5,y=3,yend=3),size=1,color=Csettings[2],alpha=alph) +
  annotate("text",label=paste("HadCM (",toupper(substr(modelVar,5,7)),")",sep=""), x = 4, y = 1,family='sans',color='Black',size = 1.7)+
  annotate("text",label=paste("TraCE (",toupper(substr(modelVar,5,7)),")",sep=""), x = 4, y = 2,family='sans',color='Black',size = 1.7) + 
  annotate("text",label="Proxy", x = 4, y = 3,family='sans',color='Black',size = 1.7) + 
  scale_y_continuous(limits=c(0,3.7),expand=c(0,0))+
  theme_void()+ 
  theme(panel.background=element_rect(colour='White',fill='White'),
        plot.background    =element_rect(colour='white',fill='White'),
        axis.line.x = element_line(color = 'black',size=0.3),
        axis.ticks.x  = element_line(color = 'Black',size=0.3), 
        axis.text.x = element_text(family='sans',size=5),
        axis.title.x = element_text(family='sans',size=5),
        axis.ticks.length.x=unit(3,"pt"),
        #plot.margin = unit(c(0.1, 0.1, 0.5,0.1), "in"),
        text = element_text(family='sans',size=8))


if (var == 'HC'){
  map2 <- map +
    theme(plot.margin = unit(c(0,rep(0,4)), "in"))+
    annotate("text",label=paste("(a) Hydroclimate (",toupper(substr(modelVar,1,3)),")",sep=""), x = 0.11, y = 0.93,family='sans',color='Black',size = 2)
} else{
  map2 <- map +annotate("text",label="(b) Temperature", x = 0.11, y = 0.93,family='sans',color='Black',size = 2.2)
  
}
map3 <- map2 + draw_plot(scale,
                        x = 0.14, 
                        y = 0.22, 
                        width = xSize, height = ySize*2.5)

ggsave(plot=map3, width = 6.5, height = 3.38, dpi = 600,
       filename = paste(file.path(dir,'Figures','RegionComposites','global_'),modelVar,'_compBandPlt.png',sep=''))
