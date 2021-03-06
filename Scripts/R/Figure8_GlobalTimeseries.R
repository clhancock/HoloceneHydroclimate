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

dir <-  '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/'
climVar  <- 'HC'
modelVar <- NA
project = TRUE

#Load Data----
#Projections
PROJ     <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
PROJorig <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
#Refrence Region Shapefiles
refReg   <- readShapePoly(file.path(dir,'Data','IPCC_refRegions','IPCC-WGI-reference-regions-v4.shp'),
                          proj4string=CRS(PROJorig))
if (project){refReg   <-  spTransform(refReg, CRSobj = PROJ)}
#Countries for basemap
countries  <- getMap("less islands")
if (project){countries  <- spTransform(countries,  CRSobj = PROJ)}
#Load Proxy Information
lipdTso  <- readRDS(file.path(dir,'Data','Proxy','LiPD','lipdData.rds'))[[climVar]]

Data <- vector(mode='list')
Data$proxy <- read.csv(file.path(dir,'Data','RegionComposites',climVar,'_MedianTS_byRegion.csv'))
for (model in c('trace','hadcm','cmip6')){
  Data[[model]]<- read.csv(file.path(dir,'Data','Model','RegionTS',
                                     paste('regional_',modelVar,'_',model,'.csv',sep='')))
}
regNames <-names(Data$proxy)[-1]
binvec <- Data$proxy$time
#ageMin <- min(binvec)
#ageMax <- max(binvec)

#Figure Settings----
#Climate variable Settings for HC vs T


xSize <- 0.1
ySize <- 0.07

#Basemap----

if (climVar=='T'){ Csettings  <- c("#fddbc7","#d6604d","#b2182b") #reds
}else{
  Csettings <- c("#ebfaeb","#145214","#0a290a") #greens
  Csettings <- c("#f6e8c3","#bf812d","#8c510a") #yellows
  Csettings <- c("#E1E6EA","#8599AB",'#434D55') #Blues
  
}
Csettings <- c("#80DBF1","#253DA1",'#000137') #Blues
Chadcm <- '#AB8599'
Ctrace <- '#99AB85'
alph<-0.8
map<- ggdraw(ggplot() +
               #Set Border around plot - probably not the best way to do this
               borders(aggregate(refReg, FUN=length), fill=NA, colour='black', size=1) +
               geom_map(data=refReg, map=fortify(refReg),aes(x=long, y=lat, group=group, map_id=id), 
                        fill="white", colour="white", size=0.3)+
               #Add Country data (basemap)
               geom_map(data=countries, map=fortify(countries), aes(x=long, y=lat, group=group, map_id=id), 
                        fill = "grey80",color="grey90",size=0.1) +
               borders(database = subset(refReg, Acronym %in% regNames), fill=NA, colour='grey40',size=0.1) + 
               coord_fixed(1) + 
               theme_void())


#Plot---- 
               

for (reg in regNames){ 
  #Load Data for Region
  regTso <- filterTs(lipdTso,paste('geo_ipccRegion ==',reg))
  regEnsNA <- read.csv(file.path(dir,'Data','RegionComposites',climVar,paste(reg,'.csv',sep='')))
  scaleVal <- 1
  if (climVar == 'HC'){scaleVal <- 30}
  hadcmVals <- (Data[['hadcm']][[reg]]-mean(Data[['hadcm']][1:5,reg],na.rm=TRUE))[1:121]*scaleVal
  traceVals <- (Data[['trace']][[reg]]-mean(Data[['trace']][1:5,reg],na.rm=TRUE))[1:121]*scaleVal
  nudgeVals[which(nudgeVals$region==reg),'x']
  #Standardize mean at 0
  regEnsNA <- as.matrix(regEnsNA - as.numeric(apply(regEnsNA[which(between(binvec,0,500)),],2,mean,na.rm=TRUE)))
  if (climVar == 'HC'){
    regEnsNA <- regEnsNA / as.numeric(apply(regEnsNA,2,sd,na.rm=TRUE))
    if (!is.na(modelVar)){regEnsNA <- regEnsNA * mean(sd(traceVals,na.rm=TRUE),sd(hadcmVals,na.rm=TRUE))}
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
  compBands$na <-  plotTimeseriesEnsRibbons(ggplot()+geom_hline(yintercept=0,size=0.05,color='black'),
                                            X=timeN$yvec, Y=regEnsNA, alp=alph,line.width=0.1,
                                            color.low='grey90',
                                            color.high='grey50',
                                            color.line='grey20')
  compBands$ts <- plotTimeseriesEnsRibbons(X=timeN$yvec, Y=regEns, alp=alph-0.2,line.width=0.1,
                                           color.low=Csettings[1],
                                           color.high=Csettings[2],
                                           color.line=Csettings[3])
  if (!is.na(modelVar)){
    plotlimit_set <- max(abs(c(traceVals,hadcmVals,Data[['cmip6']][[reg]],apply(regEns,1,mean))),na.rm=TRUE)*1.2
    compBands$ts <- compBands$ts + geom_hline(yintercept=0,size=0.05,color='black') +
      geom_line(aes(x=binvec[which(between(binvec,0,12000))],y=hadcmVals),color=Chadcm,size=0.3,alpha=alph)+
      geom_line(aes(x=binvec[which(between(binvec,0,12000))],y=traceVals),color=Ctrace,size=0.3,alpha=alph)+
      geom_boxplot(aes(x=6000,y=Data[['cmip6']][[reg]]*scaleVal),width=1000,size=0.1,alpha=alph,
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



plotTimeseriesEnsRibbons(X=timeN$yvec, Y=regEns, alp=0.7,line.width=0.1,color.low='#80DBF1',color.high='#253DA1',color.line='#000137')

scale <- ggplot() +
  scale_x_reverse('Age (ka BP)',limits=c(12,0),expand=c(0,0),breaks=seq(0,12,3),labels=c('0','','6','','12'))+
  #geom_segment(aes(x=11.5,xend=6,y=1,yend=1),size=2,color='Black') +
  #geom_segment(aes(x=11.5,xend=6,y=2,yend=2),size=2,color='Black') +
  #geom_segment(aes(x=11.5,xend=6,y=3,yend=3),size=2,color='Black') +
  #geom_segment(aes(x=11.4,xend=6.1,y=1,yend=1),size=1.5,color=Chadcm,alpha=alph) +
  #geom_segment(aes(x=11.4,xend=6.1,y=2,yend=2),size=1.5,color=Ctrace,alpha=alph) +
  #geom_segment(aes(x=11.4,xend=6.1,y=3,yend=3),size=1.5,color=Csettings[2],alpha=alph) +
  #annotate("text",label="HadCM", x = 3, y = 1,family='sans',color='Black',size = 1.7)+
  #annotate("text",label="TraCE", x = 3, y = 2,family='sans',color='Black',size = 1.7) + 
  #annotate("text",label="Proxy", x = 3, y = 3,family='sans',color='Black',size = 1.7) + 
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


if (climVar == 'HC'){
  map2 <- map +
    theme(plot.margin = unit(c(0,rep(0,4)), "in"))
    #annotate("text",label="(a) Hydroclimate", x = 0.11, y = 0.93,family='sans',color='Black',size = 2)
} else{
  map2 <- map +annotate("text",label="(b) Temperature", x = 0.11, y = 0.93,family='sans',color='Black',size = 2.2)
  
}
map3 <- map2 + draw_plot(scale,
                        x = 0.14, 
                        y = 0.22, 
                        width = xSize, height = ySize*2.5)

ggsave(plot=map3, width = 6.5, height = 3.4, dpi = 600,
       filename = paste(file.path(dir,'Figures','global_'),modelVar,'_compBandPlt.png',sep=''))
