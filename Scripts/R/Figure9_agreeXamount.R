dir <- getwd()#'/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate' #
var <- 'pre'

#Load Data
proxyDataAgree <- read.csv(file.path(dir,'Data','Proxy','proxyMetaData_HC.csv'))
proxyRegionTS <- read.csv(file.path(dir,'Data','RegionComposites','HC','MedianTS_byRegion.csv'))
hadcmRegionTS_ann <- read.csv(file.path(dir,'Data','Model','RegionalTS',paste('regional_',var,'_ANN_hadcm_land.csv',sep='')))
traceRegionTS_ann <- read.csv(file.path(dir,'Data','Model','RegionalTS',paste('regional_',var,'_ANN_trace_land.csv',sep='')))
hadcmRegionTS_jja <- read.csv(file.path(dir,'Data','Model','RegionalTS',paste('regional_',var,'_JJA_hadcm_land.csv',sep='')))
traceRegionTS_jja <- read.csv(file.path(dir,'Data','Model','RegionalTS',paste('regional_',var,'_JJA_trace_land.csv',sep='')))
hadcmRegionTS_djf <- read.csv(file.path(dir,'Data','Model','RegionalTS',paste('regional_',var,'_DJF_hadcm_land.csv',sep='')))
traceRegionTS_djf <- read.csv(file.path(dir,'Data','Model','RegionalTS',paste('regional_',var,'_DJF_trace_land.csv',sep='')))

cmipRegionTS_ann <- read.csv(file.path(dir,'Data','Model','RegionalTS',paste('regional_',var,'_ANN_cmip6_land.csv',sep='')))

#modelRegionTS <- 

#Create Empty DF
regNames <-  names(proxyRegionTS[-1])
NAs <- rep(NA,length(regNames))
data <- data.frame(regions=regNames,
                   regType=NAs,
                   proxyAgreePct=NAs,
                   #proxyhadcmCor=NAs,
                   #proxytraceCor=NAs,
                   #hadcmVariance=NAs,
                   #traceVariance=NAs,
                   hadcmRange=NAs,
                   traceRange=NAs,
                   averageRange=NAs,
                   hadcmSummerPct=NAs,
                   traceSummerPct=NAs,
                   averageSummerPct=NAs)

#Populate DF
for (reg in regNames){
  row <- which(regNames==reg)
  #
  proxySub  <- proxyDataAgree[which(proxyDataAgree$ipccReg == reg),]
  proxyVals <- proxySub$ka_6-proxySub$ka_0.5
  for (i in (which(proxySub$direction=='negative'))){proxyVals[i] <- proxyVals[i]*-1}
  proxyVals <- proxyVals[which(!is.na(proxyVals))]
  proxyVals <- proxyVals[which(proxyVals != 0)]
  data[row,'regions'] <- reg
  data[row,'proxyAgreePct'] <- abs(round(100*length(which(proxyVals>0))/length(proxyVals),1)-50)+50
  if (reg %in% c('EAS','TIB','SAS','NEAF','SAH','ESAF','SAM')){data[row,'regType'] <- 'Tropics (Monsoon)'
  } else if (reg %in% c('NCA','SCA','NES','SEAF','SEA','WSAF','NWS')){      data[row,'regType'] <- 'Tropics (Non-Monsoon)'
  } else if (reg %in% c('SAU','NZ','SSA')){                   data[row,'regType'] <- 'Mid-Latitude (Southern)'
  } else if (reg %in% c('NEN','NWN','ENA','CNA','WNA','NEU','WCE','MED','EEU','WCA','ECA','WSB','ESB','RFE')){data[row,'regType'] <- 'Mid-Latitude (Northern)'
  }else{data[row,'regType'] <-'Polar'}
  data[row,'cmip6Range'] <- mean(cmipRegionTS_ann[,reg])
  data[row,'hadcmRange'] <- -1*(mean(hadcmRegionTS_ann[which(between(hadcmRegionTS_ann$X,0,1000)),reg],na.rm=TRUE)-mean(hadcmRegionTS_ann[which(between(hadcmRegionTS_ann$X,5500,6500)),reg],na.rm=TRUE)) #diff(range(hadcmRegionTS_ann[,reg],na.rm=TRUE))
  data[row,'traceRange'] <- -1*(mean(traceRegionTS_ann[which(between(traceRegionTS_ann$X,0,1000)),reg],na.rm=TRUE)-mean(traceRegionTS_ann[which(between(traceRegionTS_ann$X,5500,6500)),reg],na.rm=TRUE))#diff(range(traceRegionTS_ann[,reg],na.rm=TRUE))
  data[row,'proxyAgreePct'] <- (round(100*length(which(proxyVals>0))/length(proxyVals),1)-50)*sign(mean(data[row,'hadcmRange'],data[row,'traceRange']))+50
  #data[row,'proxyAgreePct'] <- (round(100*length(which(proxyVals>0))/length(proxyVals),1)-50)*sign(median(data[row,'cmip6Range']))+50
  #Nick/Michael comments
  data[row,'hadcmRange'] <- abs(data[row,'hadcmRange'])
  data[row,'traceRange'] <- abs(data[row,'traceRange'])
  #data[row,'hadcmRange'] <- 100*data[row,'hadcmRange']/mean(mean(hadcmRegionTS_ann[which(between(hadcmRegionTS_ann$X,0,1000)),reg],na.rm=TRUE),mean(hadcmRegionTS_ann[which(between(hadcmRegionTS_ann$X,5500,6500)),reg],na.rm=TRUE)) #diff(range(hadcmRegionTS_ann[,reg],na.rm=TRUE))
  #data[row,'traceRange'] <- 100*data[row,'traceRange']/mean(mean(traceRegionTS_ann[which(between(traceRegionTS_ann$X,0,1000)),reg],na.rm=TRUE),mean(traceRegionTS_ann[which(between(traceRegionTS_ann$X,5500,6500)),reg],na.rm=TRUE)) #diff(range(traceRegionTS_ann[,reg],na.rm=TRUE))
  data[row,'averageRange'] <- mean(c(data[row,'traceRange'],data[row,'hadcmRange']))
  if (reg %in% c('NWS','SAM','NES','SSA','WSAF','ESAF','SAU','NZ','EAF')){
    data[row,'hadcmSummerPct']<-round(100*mean(traceRegionTS_djf[,reg],na.rm=TRUE)/(mean(traceRegionTS_ann[,reg],na.rm=TRUE)*4))
    data[row,'traceSummerPct']<-round(100*mean(hadcmRegionTS_djf[,reg],na.rm=TRUE)/(mean(hadcmRegionTS_ann[,reg],na.rm=TRUE)*4))
  } else{
    data[row,'hadcmSummerPct']<-round(100*mean(traceRegionTS_jja[,reg],na.rm=TRUE)/(mean(traceRegionTS_ann[,reg],na.rm=TRUE)*4))
    data[row,'traceSummerPct']<-round(100*mean(hadcmRegionTS_jja[,reg],na.rm=TRUE)/(mean(hadcmRegionTS_ann[,reg],na.rm=TRUE)*4))
  }
  data[row,'averageSummerPct'] <- mean(c(data[row,'hadcmSummerPct'],data[row,'traceSummerPct']))
}


pointData <- data.frame(longitude=c(proxyDf$longitude),latitude=c(proxyDf$latitude))
pointData <- spTransform(SpatialPointsDataFrame(coords=pointData,data = pointData, proj4string = CRS(PROJorig)),CRSobj = PROJ)
proxyDf$lonsPrj <- pointData@coords[,1]
proxyDf$latsPrj <- pointData@coords[,2]

dataTable <- fortify(regionsSelect$composite)
dataTable$regType     <- NA
dataTable$regType[which(dataTable$id %in% c('NEN','NWN','ENA','CNA','WNA','NEU','WCE','MED','EEU','WCA','ECA','WSB','ESB','RFE'))] <- '#1D6996'
dataTable$regType[which(dataTable$id %in% c('SAU','NZ','SSA'))] <- '#0F8554'
dataTable$regType[which(dataTable$id %in% c('GIC','EAN'))] <- '#6F4070'
dataTable$regType[which(dataTable$id %in% c('EAS','TIB','SAS','NEAF','SAH','ESAF','SAM'))] <- '#CC503E'
dataTable$regType[which(dataTable$id %in% c('NCA','SCA','NES','SEAF','SEA','WSAF','NWS'))] <- '#EDAD08'
  

proxyMapRegions <- basemap + 
  #Add refrence regions boundaries
  geom_map(data=regionsSelect$composite, map=dataTable, alpha=0.75, size=0.1, color='black' , 
           aes(x=long, y=lat, group=group, map_id=id,fill=as.factor(dataTable$regType))) +
  scale_fill_manual(values=c('#0F8554','#1D6996','#6F4070','#CC503E','#EDAD08'),)+
  #Add proxy sites
  geom_point(data=as.data.frame(proxyDf), aes(x=lonsPrj , y=latsPrj), shape=1,size=0.5,stroke=0.3) +
  #Format Legend
  theme(plot.background = element_rect(fill = NA,color=NA),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "in"),
        legend.position = 'none')
proxyMapRegions
library(plyr )
df <- na.omit(data)
find_hull <- function(df) df[chull(df$proxyAgreePct, df$averageRange), ]
hulls <- ddply(df, "regType", find_hull)

scatter <- ggplot(data=data,aes(x=averageRange,y=proxyAgreePct,shape=regType,fill=regType,group = regType))+
  geom_hline(yintercept = 50,size=0.2)+
  geom_polygon(data = hulls, alpha = 0.5) +
  geom_point(size=2.5,color='Black')+
  geom_point(size=2.5,color='Black')+
  annotate("text",label="50% Proxy \n Even Wet/Dry Split", x = 0.25, y = 50,family='sans',color='grey40',size=1.5)+
  annotate("text",label="Proxy-Proxy Agreement & \n Proxy-Model Agreement", x = 0.25, y = 90,family='sans',color='grey40',size=1.5)+
  annotate("text",label="Proxy-Proxy Agreement but \n Proxy-Model Disagreement", x = 0.25, y = 10,family='sans',color='grey40',size=1.5)+
  annotate("segment", x = 0.25, xend = 0.25, y = 60, yend = 80, colour = "grey40", size=0.5, arrow=arrow(type='closed',length = unit(0.05, "inches")))+
  annotate("segment", x = 0.25, xend = 0.25, y = 40, yend = 20, colour = "grey40", size=0.5, arrow=arrow(type='closed',length = unit(0.05, "inches")))+
  scale_fill_manual(values=c('#1D6996','#0F8554','#6F4070','#CC503E','#EDAD08'))+
  scale_shape_manual(values=c(24,25,23,21,22))+
  theme_bw()+
  #scale_x_continuous(limits=c(0,1000))+
  coord_cartesian(xlim=c(0,0.3), ylim=c(0,100),expand	=FALSE) +
  labs(y='% of Proxy Records with the Same Sign 6-0 ka Anomaly \n as the Simulated 6-0 ka Anomaly',x='Absolute Value of Simulated 6-0 ka Anomalies (mm/day)')+
  #labs(y='% of Proxy Records With the Same Sign 6-0 ka Anomaly',x='Absolute Value of Simulated 6-0 ka Anomalies \n _______________________________________ \n \n Mean of Values At 0ka and 6ka')+
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        text = element_text(family=figFont,size=6),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "in"),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        plot.background= element_blank(),
        legend.key.height=unit(0.15,"in"),
        legend.key.width=unit(0.15,"in"),
        legend.background	= element_rect(fill='white',color='black'))+    
  
  guides(fill = guide_legend(ncol = 2))
scatter

z <- cowplot::ggdraw(ggplot() + 
                       coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand=FALSE)+
                       annotate("text",label="(a)", x = 0.15, y = 0.97,family='sans',color='black',size=3)+
                       annotate("text",label="(b)", x = 0.15, y = 0.67,family='sans',color='black',size=3)+
                       theme_void()+
                       theme(plot.background= element_rect(colour='White',fill='White'),
                             panel.background = element_rect(colour='White',fill='White')))+
  draw_plot(scatter,  x = 0,   y = 0, width = 1, height = 0.66)+
  draw_plot(proxyMapRegions,  x = 0.05,   y = 0.62, width = 0.95, height = 0.4)

if (save) {
  ggsave(plot=z, width = 3.25, height = 4.5, dpi = 600,
         filename = file.path(dir,"Figures","Model",'ProxyAgreeByModelAnom3.png'))
}
save<-TRUE