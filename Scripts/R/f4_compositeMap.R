library(lipdR) #to read and interact with LiPD data
library(geoChronR) #for plotting mostly
library(dplyr) #and dplyr for data.frame manipulation
library(ggplot2)
library(cowplot)
library(ggstar)
library(maptools)
library(proj4)
library(sp)

dataDir <- getwd()# '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate'
var      <- 'HC'
modelvar <- 'p-e_ANN'
project = TRUE

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
lipdTso  <- readRDS(file.path(dataDir,'Data','LiPD','lipdData.rds'))[[var]]

Data <- vector(mode='list')
Data$proxy <- read.csv(file.path(dataDir,'Data','RegionComposites',var,'MedianTSbyRegion.csv'))
for (model in c('trace','hadcm','cmip6')){
  Data[[model]]<- read.csv(file.path(dataDir,'Data','Model','RegionTS',
                                     paste('regional_',modelvar,'_',model,'.csv',sep='')))
}
regNames <-names(Data$proxy)[-1]
binvec <- Data$proxy$time
#ageMin <- min(binvec)
#ageMax <- max(binvec)

#Figure Settings----
#Climate variable Settings for HC vs T


xSize <- 0.085
ySize <- 0.06
nudgeVals <- data.frame(region=regNames,x=rep(0,length(regNames)),y=rep(0,length(regNames)))
nudgeVals[which(regNames=='ARO'),'x'] <- 0.025
nudgeVals[which(regNames=='ARO'),'y'] <- -0.01
nudgeVals[which(regNames=='CAF'),'x'] <- -0.04
nudgeVals[which(regNames=='CAR'),'x'] <- 0.007
nudgeVals[which(regNames=='CAR'),'y'] <- 0.04
nudgeVals[which(regNames=='NWN'),'x'] <- 0.0115
nudgeVals[which(regNames=='NWN'),'y'] <- 0.007
nudgeVals[which(regNames=='NEN'),'y'] <- -0.022
nudgeVals[which(regNames=='NEN'),'x'] <- 0.007
nudgeVals[which(regNames=='GIC'),'x'] <- 0.003
nudgeVals[which(regNames=='WNA'),'x'] <- -0.042
nudgeVals[which(regNames=='CNA'),'x'] <- 0.0035
nudgeVals[which(regNames=='ENA'),'x'] <- 0.042
nudgeVals[which(regNames=='SAM'),'y'] <- -0.02
nudgeVals[which(regNames=='NES'),'y'] <- 0.02
nudgeVals[which(regNames=='NSA'),'x'] <- 0.025
nudgeVals[which(regNames=='NSA'),'y'] <- 0.04
nudgeVals[which(regNames=='NAO'),'y'] <- -0.03
nudgeVals[which(regNames=='NAO'),'x'] <- 0.004
nudgeVals[which(regNames=='NEU'),'x'] <- 0.004
nudgeVals[which(regNames=='NEU'),'y'] <- 0.01
nudgeVals[which(regNames=='WCE'),'x'] <- -0.005
nudgeVals[which(regNames=='WSAF'),'x'] <- -0.022
nudgeVals[which(regNames=='ESAF'),'x'] <- 0.022
nudgeVals[which(regNames=='WCA'),'x'] <- -0.01
nudgeVals[which(regNames=='ESB'),'y'] <- 0.004
nudgeVals[which(regNames=='ESB'),'x'] <- -0.01
nudgeVals[which(regNames=='ECA'),'y'] <- 0.013
nudgeVals[which(regNames=='RFE'),'x'] <- -0.04
nudgeVals[which(regNames=='EPO'),'x'] <- -0.04
nudgeVals[which(regNames=='SAS'),'y'] <- -0.002
nudgeVals[which(regNames=='SAS'),'x'] <- -0.005
nudgeVals[which(regNames=='EAS'),'x'] <- 0.003
nudgeVals[which(regNames=='SEA'),'x'] <- -0.017
nudgeVals[which(regNames=='SAU'),'x'] <- -0.036
nudgeVals[which(regNames=='NZ'),'x'] <- -0.025
nudgeVals[which(regNames=='NZ'),'y'] <- 0.02

#Basemap----

if (var=='T'){ Csettings  <- c("#fddbc7","#d6604d","#b2182b") #blues
}else{
  Csettings <- c("#ebfaeb","#145214","#0a290a") #greens
  Csettings <- c("#f6e8c3","#bf812d","#8c510a") #yellows
}

#Plot----

map<- ggdraw(ggplot() +
               #Set Border around plot - probably not the best way to do this
               borders(aggregate(refReg, FUN=length), fill=NA, colour='black', size=1) +
               geom_map(data=refReg, map=fortify(refReg),aes(x=long, y=lat, group=group, map_id=id), 
                        fill="white", colour="white", size=0.5)+
               #Add Country data (basemap)
               geom_map(data=countries, map=fortify(countries), aes(x=long, y=lat, group=group, map_id=id), 
                        fill = "grey80",color="grey90",size=0.2) +
               borders(database = subset(refReg, Acronym %in% regNames), fill=NA, colour='grey40',size=0.25) + 
               coord_fixed(1) + 
               theme_void())

for (reg in regNames){ 
  #Load Data for Region
  regTso <- filterTs(lipdTso,paste('geo_ipccRegion ==',reg))
  regEnsNA <- read.csv(file.path(dataDir,'Data','RegionComposites',var,paste(reg,'.csv',sep='')))
  scaleVal <- 1
  if (var == 'HC'){scaleVal <- 30}
  hadcmVals <- (Data[['hadcm']][[reg]]-mean(Data[['hadcm']][1:5,reg],na.rm=TRUE))[1:121]*scaleVal
  traceVals <- (Data[['trace']][[reg]]-mean(Data[['trace']][1:5,reg],na.rm=TRUE))[1:121]*scaleVal
  nudgeVals[which(nudgeVals$region==reg),'x']
  #Standardize mean at 0
  regEnsNA <- as.matrix(regEnsNA - as.numeric(apply(regEnsNA[which(between(binvec,0,500)),],2,mean,na.rm=TRUE)))
  if (var == 'HC'){
    regEnsNA <- regEnsNA / as.numeric(apply(regEnsNA,2,sd,na.rm=TRUE))
    if (!is.na(modelvar)){regEnsNA <- regEnsNA * mean(sd(traceVals,na.rm=TRUE),sd(hadcmVals,na.rm=TRUE))}
  }
  regEns   <- matrix(NA,nrow(regEnsNA),ncol(regEnsNA))
  #RegionMap
  timeN <-  plotTimeAvailabilityTs(regTso,age.range=c(0,12000),step=100)$dat %>% 
    group_by(yvec) %>% 
    summarise(count=sum(value),nPct=sum(value)/length(regTso))
  if (length(which(timeN$nPct<=0.5) > 0)){
    mh <- which(timeN$yvec==6000)
    idx <- c(max(which(timeN$nPct[1:mh]<=0.5),1),
             min(which(timeN$nPct[mh:length(timeN$yvec)]<=0.5)+mh,length(timeN$yvec)))
    regEns[idx[1]:idx[2],] <- regEnsNA[idx[1]:idx[2],]                  
  } else{regEns<-regEnsNA}
  plotlimit_set <- 5
  regPlt <- ggdraw(ggplot()+theme_void()+theme(plot.background= element_rect(colour='White',fill='White')))
  compBands <- vector(mode = 'list')
  compBands$na <-  plotTimeseriesEnsRibbons(ggplot()+geom_hline(yintercept=0,size=0.2,color='black'),
                                            X=timeN$yvec, Y=regEnsNA, alp=0.8,line.width=0.1,
                                            color.low='grey90',
                                            color.high='grey50',
                                            color.line='grey20')
  compBands$ts <- plotTimeseriesEnsRibbons(X=timeN$yvec, Y=regEns, alp=0.4,line.width=0.1,
                                           color.low=Csettings[1],
                                           color.high=Csettings[2],
                                           color.line=Csettings[3])
  if (!is.na(modelvar)){
    plotlimit_set <- max(abs(c(traceVals,hadcmVals,Data[['cmip6']][[reg]],apply(regEns,1,mean))),na.rm=TRUE)
    compBands$ts <- compBands$ts + geom_hline(yintercept=0,size=0.2,color='black') +
      geom_line(aes(x=binvec[which(between(binvec,0,12000))],y=hadcmVals),color='#4527A0',size=0.3,alpha=0.6)+
      geom_line(aes(x=binvec[which(between(binvec,0,12000))],y=traceVals),color='#2E7D32',size=0.3,alpha=0.6)+
      geom_boxplot(aes(x=6000,y=Data[['cmip6']][[reg]]*scaleVal),width=1000,size=0.1,alpha=0.8,
                   outlier.size=0.5,outlier.stroke = 0.15,outlier.alpha=1,outlier.colour='Black')
  }
  for (plt in names(compBands)){
    compBands[[plt]] <- compBands[[plt]] + 
      scale_x_reverse(limits=c(12100,-100), expand=c(0,0), n.breaks=7,sec.axis = dup_axis())+ 
      scale_y_continuous(limits=c(plotlimit_set*c(-1000,1000)),breaks=seq(-100,100,2),sec.axis = dup_axis())+
      coord_cartesian(xlim=c(12000,0), ylim=c(plotlimit_set*c(-1,1))) +
      theme_void() +
      theme(panel.border    = element_rect(color='Black',fill=NA,size=0.5),
            axis.ticks      = element_line(color='Black',size=0.1), 
            axis.ticks.length = unit(-1,"pt"),
            plot.margin       = unit(c(0, 0, 0, 0), "in"), legend.position='none')
    regPlt <- regPlt + draw_plot(compBands[[plt]], x = 0, y = 0, width = 1, height = 1)
  }
  plotLat <- refReg@polygons[[which(refReg@data[["Acronym"]]==reg)]]@Polygons[[1]]@labpt[2]
  plotLon <- refReg@polygons[[which(refReg@data[["Acronym"]]==reg)]]@Polygons[[1]]@labpt[1]
  map <- map + draw_plot(regPlt,
                         x = nudgeVals[which(nudgeVals$region==reg),'x']+(0.5-xSize/2)+0.5*plotLon/(9050504*2), 
                         y = nudgeVals[which(nudgeVals$region==reg),'y']+(0.5-ySize/2)+0.92*plotLat/(8625155*2), 
                         width = xSize, height = ySize)
}



#plotTimeseriesEnsRibbons(X=timeN$yvec, Y=regEns, alp=0.6,line.width=0.1,color.low='#f6e8c3',color.high='#dfc27d',color.line='#8c510a')

scale <- ggplot() +
  scale_x_reverse('Age (ka BP)',limits=c(12,0),expand=c(0,0),breaks=seq(0,12,6))+
  geom_segment(aes(x=11.7,xend=6.3,y=1,yend=1),size=2,color='#4527A0',alpha=0.6,label='HadCM') +
  geom_segment(aes(x=11.7,xend=6.3,y=2,yend=2),size=2,color='#2E7D32',alpha=0.6,label='Trace') +
  annotate("text",label="HadCM", x = 2.7, y = 1,size=2,family='sans',color='#4527A0') + 
  annotate("text",label="TraCE", x = 2.7, y = 2,size=2,family='sans',color='#2E7D32') + 
  scale_y_continuous(limits=c(0,2.5),expand=c(0,0))+
  theme_void()+ 
  theme(panel.background=element_rect(colour='White',fill='White'),
        plot.background    =element_rect(colour='white',fill='White'),
        axis.line.x = element_line(color = 'black'),
        axis.ticks.x  = element_line(color = 'Black',size=0.4), 
        axis.text.x = element_text(family='sans',size=8),
        axis.title.x = element_text(family='sans',size=8),
        axis.ticks.length.x=unit(3,"pt"),
        #plot.margin = unit(c(0.1, 0.1, 0.5,0.1), "in"),
        text = element_text(family='sans',size=8))

map2 <- map + draw_plot(scale,
                        x = 0.15, 
                        y = 0.22, 
                        width = xSize, height = ySize*2.5)

ggsave(plot=map2, width = 6.5, height = 6.5*0.5072, dpi = 400,
       filename = paste(file.path(dataDir,'Figures','global_'),modelvar,'_compBandPlt.png',sep=''))
