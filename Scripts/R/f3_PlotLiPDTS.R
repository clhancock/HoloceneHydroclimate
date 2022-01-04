library(lipdR) #to read and interact with LiPD data
library(geoChronR) #for plotting mostly
library(magrittr) #we'll be using the magrittr pipe ( %>% ) for simplicity
library(dplyr) #and dplyr for data.frame manipulation
library(ggplot2)
library(cowplot)
library(geoChronR)
library(ggplot2)
library(ggstar)
library(lipdR)
library(maptools)
library(proj4)
library(RColorBrewer)
library(rworldmap)
library(scales)
library(sp)
library(tidyverse)
githubDir <- '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate' #getwd()
githubDir <- getwd()
#Settings 
plotSettings <- plotVals$specific
climVar <- 'HC'
lipdTSO <- readRDS(file.path(githubDir,'Data','LiPD','lipdData.rds'))[[climVar]]
##############################
proxyMetaData <- tibble(dataset       = pullTsVariable(lipdTSO,'dataSetName'),
                        tsid          = pullTsVariable(lipdTSO,'paleoData_TSid'),
                        longitude     = pullTsVariable(lipdTSO,'geo_longitude'),
                        latitude      = pullTsVariable(lipdTSO,'geo_latitude'),
                        ipccReg       = pullTsVariable(lipdTSO,'geo_ipccRegion'),
                        Category      = pullTsVariable(lipdTSO,'Category'),
                        CategorySpec  = pullTsVariable(lipdTSO,'CategorySpecific'),
                        recordRange   = pullTsVariable(lipdTSO,'ageRange'),
                        recordRes     = pullTsVariable(lipdTSO,'ageRes'),
                        recordResPlus = pullTsVariable(lipdTSO,'ageResPlus'),
                        season        = pullTsVariable(lipdTSO,'climateInterpretation1_seasonalityGeneral'),
                        climInterp    = pullTsVariable(lipdTSO,'climateInterpretation1_variable'),
                        source        = pullTsVariable(lipdTSO,'Source'))

#Catergory settings from f1 rmd file
plotVals <- vector(mode='list')
plotVals$general  <- vector(mode='list')
plotVals$specific <- vector(mode='list')

plotVals$general$names  <- sort(unique(proxyMetaData$Category))
plotVals$general$colors <- c("powder blue","corn flower blue",
                             "dark blue","dark orchid","grey",
                             "forest green","firebrick")
plotVals$general$shapes <- c(12,21,15,5,13,14,11)

plotVals$specific$names <- sort(unique(proxyMetaData$CategorySpec))
plotVals$specific$colors<- c("powder blue","corn flower blue",
                             "dark blue","dark orchid","grey40","grey",
                             "forest green","yellowgreen","lightcoral","firebrick","darkorange")
plotVals$specific$shapes<- c(12,21,15,5,6,13,14,1,23,11,17)

plotSettings <- plotVals$specific

PROJ       <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
PROJorig   <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
refregions <- readShapePoly(file.path(githubDir,'Data','IPCC_refRegions','IPCC-WGI-reference-regions-v4.shp'),
                            proj4string=CRS(PROJorig))
refregions <-  spTransform(refregions, CRSobj = PROJ)

countries  <- getMap("less islands")
countries  <- spTransform(countries,  CRSobj = PROJ)


for (regNo in as.numeric(refregions@data[["Acronym"]])){
  reg <- levels(refregions@data[["Acronym"]])[regNo]
  print(reg)
  count <- which(as.numeric(refregions@data[["Acronym"]])==regNo)
  #
  #Stack Plot
  #
  stack.df <- filterTs(lipdTSO,paste('geo_ipccRegion ==',reg))
  if (length(stack.df) == 0){next}
  stack.df <- tidyTs(stack.df,age.var = "age") %>% 
    filter(between(age,0,12000)) %>% #only years from Holocene
    group_by(paleoData_TSid) %>%     #group by column
    arrange(geo_latitude)            #North to south
  stackColors <- c()
  #Get colors for categories within region
  for (category in unique(stack.df$CategorySpecific)){
    stackColors <- c(stackColors,
                     plotSettings$colors[which(plotSettings$names == category)])
  }
  stackPlot <- plotTimeseriesStack(stack.df,
                                   time.var = "age",
                                   color.var =  "CategorySpecific",
                                   invert.var = 'climateInterpretation1_interpDirection',
                                   scale.factor= 0.1,
                                   color.ramp = stackColors) + 
    scale_x_reverse(name = "Age (yr BP)", limits=c(12000,0),expand=c(0,0),
                    n.breaks=7,sec.axis = sec_axis(~.,labels=NULL)) 
  #
  #Region map
  #
  proxyDataPrj <- SpatialPointsDataFrame(stack.df[,c("geo_longitude", "geo_latitude")], 
                                         stack.df, proj4string=CRS(PROJorig))
  proxyDataPrj <- spTransform(proxyDataPrj, CRSobj = PROJ)
  proxyDFprj   <- as.data.frame(proxyDataPrj)
  idx <- which(plotSettings$names %in% unique(proxyDFprj$CategorySpecific))
  refrenceRegShp <- subset(refregions, Acronym ==reg)
  regMap <-  ggplot() +
     geom_map(data=countries, map=fortify(countries),
             aes(x=long, y=lat, group=group, map_id=id), 
             fill = "grey80",color="grey90",size=0.4) +
    coord_fixed(1,xlim= range(fortify(refrenceRegShp)$lon),
                ylim=range(fortify(refrenceRegShp)$lat)) + 
    geom_map(data=refrenceRegShp, map=fortify(refrenceRegShp), 
             aes(x=long, y=lat, group=group, map_id=id),
             fill=NA, alpha=0.75, size=0.5, color='black') +
    geom_star(data=proxyDFprj,
              aes(x=geo_longitude.1 , y=geo_latitude.1,
                  starshape=CategorySpecific, fill=CategorySpecific),
              size=4,color='Black',alpha=1,starstroke=0.5) + 
    scale_fill_manual(values=plotSettings$colors[idx],name= 'Proxy Category') +
    scale_starshape_manual(values=plotSettings$shapes[idx],name= 'Proxy Category') +
    theme_void() + 
    theme(legend.position = 'none') 
  ggsave(file.path(githubDir,'Figures','Dashboard','_Summary',
                   paste(count,'_',reg,'__regMap.png',sep='')),
         plot=regMap,device='png',width=5,height=5,units='in')
  ggsave(file.path(githubDir,'Figures','Dashboard','_Summary',
                   paste(count,'_',reg,'__stackPlot.png',sep='')),
         plot=stackPlot,device='png',width=8.5,height=11,units='in')
  #
  #
  #
  for (i in 1:length(unique(proxyDFprj$paleoData_TSid))){
    df <- proxyDFprj %>% filter(paleoData_TSid==unique(proxyDFprj$paleoData_TSid)[i])
    ts <- lipdTSO[[which(pullTsVariable(lipdTSO,'paleoData_TSid')==unique(proxyDFprj$paleoData_TSid)[i])]]
    col <- plotSettings$colors[which(plotSettings$names == ts$CategorySpecific)]
    shp <- plotSettings$shapes[which(plotSettings$names == ts$CategorySpecific)]
    plt <- ggplot(df,aes(x=age,y=paleoData_values)) +
      geom_hline(yintercept=mean(df$paleoData_values,na.rm=TRUE)) +
      geom_star(fill=col,color='black',starshape=shp,size=1,starstroke=0.5,alpha=0.7) + 
      geom_line(color=col,alpha=0.7) + 
      scale_x_reverse(name = "Age (yr BP)", limits=c(12000,0),expand=c(0,0),
                      n.breaks=7,sec.axis = sec_axis(~.,labels=NULL)) +
      theme_bw() +     
      ggtitle(paste(unique(df[,'geo_ipccRegion']),': ',unique(df[,'dataSetName']),sep='')) + 
      theme(text = element_text(family='sans',size=8),
            plot.title = element_text(hjust = 0.5,family='sans',size=10,face='bold'))
    if (ts$climateInterpretation1_interpDirection == 'negative'){
      plt <- plt + scale_y_reverse(name=paste(ts$paleoData_variableName,
                                              ' (',ts$paleoData_units,')',sep=''))
    } else{
      plt <- plt + scale_y_continuous(name=paste(ts$paleoData_variableName,
                                                 ' (',ts$paleoData_units,')',sep=''))
    }
    map <- regMap + 
      geom_star(data=df,aes(x=geo_longitude.1 , y=geo_latitude.1,
                            starshape=CategorySpecific, fill=CategorySpecific),
                size=4,color='gold',alpha=1,starstroke=4) +
      theme(text = element_text(family='sans',size=8))
    h <- 0
    txt <- ggplot(df) 
    for (name in c('paleoData_TSid','dataSetName',
                   'geo_latitude','geo_longitude','geo_elevation',
                   'archiveType','Category','CategorySpecific',
                   'paleoData_proxyGeneral','paleoData_proxy','paleoData_proxyDetail','paleoData_variableName',
                   'climateInterpretation1_seasonalityGeneral',
                   'climateInterpretation1_variable','paleoData_units',
                   'Source','pub1_title','pub1_doi','pub2_title','pub2_doi')){
      h <- h-1
      if (name %in% names(df) == FALSE){df[,name] <- NA}
      txt <- txt + 
        annotate("text", x = 0, y = h, label = paste(name,': ',sep=''),hjust = 0,size=2.5) +
        annotate("text", x = 0.25, y = h, label = unique(df[,name]),hjust = 0,size=2.5)    
    }
    txt <- txt + scale_x_continuous(limits=c(0,0.6)) + theme_void() 
    bkg <- ggplot()+
      theme_void() +
      theme(plot.background = element_rect(fill = NA,color='Black'))
    summary <- ggdraw() + 
      draw_plot(plt, x = 0, y = 0.5, width = 1, height = 0.5) +
      draw_plot(map, x = 0, y = 0, width = 0.5, height = 0.5) +
      draw_plot(txt, x = 0.5, y = 0, width = 0.5, height = 0.5) +
      draw_plot(bkg, x = 0, y = 0, width = 1, height = 1) 
    ggsave(file.path(githubDir,'Figures','Dashboard',
                     paste(count,'_',
                           unique(df[,'geo_ipccRegion']),'_',
                           unique(df[,'geo_latitude']),'_',
                           unique(df[,'dataSetName']),'_',
                           unique(df[,'paleoData_TSid']),
                           '.png',sep='')),
           plot=summary,device='png',width=10,height=5,units='in')
  }
}
