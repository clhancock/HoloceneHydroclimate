#Script to create individual dashboard pdfs for each record in the Holocene Hydroclimate dataset

#Load Packages--------------------------------------------------------------------------------

library(cowplot)
library(ggstar)
library(ggplot2)
library(gsheet)
library(lipdR)
library(maptools)
library(proj4)
library(RCurl)
library(sf)
library(sp)
library(tidyverse)

print("Packages Loaded")

#Set Working Directory
wd = '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/2021_HoloceneHydroclimate/2021_HoloceneHydroclimate/'

#### Load Data

var     <- 'HC'
lipdTSO <- readRDS(file.path(wd,'Data','Proxy','lipdData.rds'))[[var]]
proxyDf <- read.csv(file=file.path(wd,'Data','Proxy',paste0('proxyMetaData_',var,'.csv')))
print("Proxy data loaded ")

#Load IPCC region data
load(url('https://github.com/SantanderMetGroup/ATLAS/blob/main/reference-regions/IPCC-WGI-reference-regions-v4_R.rda?raw=true'), verbose = TRUE)

#Set Projections
PROJ     <- '+proj=robin   +ellps=WGS84 +datum=WGS84 +no_defs +lon_0=0 +x_0=0 +y_0=0 +units=m'
PROJorig <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
#Transform projection
refregions <-  spTransform(IPCC_WGI_reference_regions_v4, CRSobj = PROJ)
#Countries for basemap
countries  <- rworldmap::getMap("less islands")
countries  <- sp::spTransform(countries,  CRSobj = PROJ)
#Transform proxy projections
proxyDf <- as.data.frame(spTransform(SpatialPointsDataFrame(proxyDf[,c("longitude", "latitude")], 
                                                            proxyDf, proj4string=CRS(PROJorig)), CRSobj = PROJ))


#### Figure Settings
save     <- FALSE
specific <- TRUE 

if (save){ print(paste0("save ",var," figs"))
} else{    print(paste0("plot ",var," figs"))}

figFont <- 'Times New Roman'
figText <- 10
figSize <- c(6.5,3)

#Colors and Shapes--------------------------------------------------------------------------------

plotSettings <- vector(mode='list')
Csettings <- c("#92C5DE","#4393C3",'#2166AC')
if (var == 'HC'){
  if(specific){
    plotSettings$names <- sort(unique(proxyDf$CategorySpec))
    #https://carto.com/carto-colors/
    plotSettings$color <- as.character(plotSettings$names)
    plotSettings$color[which(plotSettings$names=="Glacier Ice (Accumulation)")] <- "#5F4690" #"powder blue"
    plotSettings$color[which(plotSettings$names=="Shoreline (Lake Level)")]  <- "#38A6A5" #"corn flower blue"
    plotSettings$color[which(plotSettings$names=="Lake Sediment (δ18O)")]    <- "#1D6996" #"dark blue"
    plotSettings$color[which(plotSettings$names=="Leaf Wax (δD)")]           <- "#94346E" # "dark orchid" #δ
    plotSettings$color[which(plotSettings$names=="Other (calibrated)")]      <- "grey40" #"grey40"
    plotSettings$color[which(plotSettings$names=="Other (not calibrated)")]  <- "grey" #"grey"
    plotSettings$color[which(plotSettings$names=="Pollen (calibrated)")]     <- "#0F8554" #"forest green"
    plotSettings$color[which(plotSettings$names=="Pollen (not calibrated)")] <- "#73AF48" #"" #"yellowgreen"
    plotSettings$color[which(plotSettings$names=="Speleothem (other)")]      <- "#EDAD08" #"darkorange"
    plotSettings$color[which(plotSettings$names=="Speleothem (δ13C)")]       <- "#E17C05" #"lightcoral"
    plotSettings$color[which(plotSettings$names=="Speleothem (δ18O)")]       <- "#CC503E" #"firebrick"
    #
    plotSettings$shape <- as.character(plotSettings$names) 
    plotSettings$shape[which(plotSettings$names=="Glacier Ice (Accumulation)")] <- 12
    plotSettings$shape[which(plotSettings$names=="Shoreline (Lake Level)")]  <- 21
    plotSettings$shape[which(plotSettings$names=="Lake Sediment (δ18O)")]    <- 15
    plotSettings$shape[which(plotSettings$names=="Leaf Wax (δD)")]           <- 5
    plotSettings$shape[which(plotSettings$names=="Other (calibrated)")]      <- 6
    plotSettings$shape[which(plotSettings$names=="Other (not calibrated)")]  <- 13
    plotSettings$shape[which(plotSettings$names=="Pollen (calibrated)")]     <- 14
    plotSettings$shape[which(plotSettings$names=="Pollen (not calibrated)")] <- 1
    plotSettings$shape[which(plotSettings$names=="Speleothem (other)")]      <- 17
    plotSettings$shape[which(plotSettings$names=="Speleothem (δ13C)")]       <- 23
    plotSettings$shape[which(plotSettings$names=="Speleothem (δ18O)")]       <- 11  
  }
}
plotSettings$shape <- as.numeric(plotSettings$shape)



for (reg in as.character(refregions@data[["Acronym"]])){
  n <- which(as.character(refregions@data[["Acronym"]])==reg)
  regTSO <- lipdTSO[which(pullTsVariable(lipdTSO,'geo_ipccRegion')==reg)]
  if (length(regTSO) == 0){next}
  refrenceRegShp <- subset(refregions, Acronym ==reg)
  regionDf <- proxyDf[which(proxyDf$ipccReg==reg),]
  idx <- which(plotSettings$names %in% regionDf$CategorySpec)
  regMap <-  ggplot() +
    geom_map(data=countries, map=fortify(countries),
             aes(x=long, y=lat, group=group, map_id=id), 
             fill = "grey80",color="grey90",size=0.4) +
    geom_map(data=refrenceRegShp, map=fortify(refrenceRegShp), 
             aes(x=long, y=lat, group=group, map_id=id),
             fill=NA, alpha=0.75, size=0.5, color='black') +
    geom_point(data=proxyDf,aes(x=longitude.1 ,y=latitude.1),size=1,color='Black') + 
    geom_star(data=regionDf,aes(x=longitude.1 ,y=latitude.1,fill=CategorySpec,starshape=CategorySpec),
              size=3,color='Black',starstroke=0.5) + 
    coord_fixed(xlim= range(fortify(refrenceRegShp)$lon), ylim=range(fortify(refrenceRegShp)$lat)) + 
    scale_fill_manual(     values=plotSettings$color[idx]) +
    scale_starshape_manual(values=plotSettings$shape[idx]) +
    theme_void() + 
    theme(legend.position = 'none') 
  #
  print(paste(n,reg,sep=". "))
  tsn <- 0
  for (ts in arrange(regionDf, desc(latitude), longitude)$tsid){
    #ts <- arrange(regionDf, desc(latitude), longitude)$tsid#[3]
    tsn<-tsn+1
    siteDf  <- regionDf[which(regionDf$tsid==ts),]
    tso     <- regTSO[[which(pullTsVariable(regTSO,'paleoData_TSid')==ts)]]
    col <- plotSettings$color[which(plotSettings$names == tso$CategorySpecific)]
    shp <- plotSettings$shape[which(plotSettings$names == tso$CategorySpecific)]
    df <- data.frame(
      age    = as.numeric(tso$age), 
      values = as.numeric(tso$paleoData_values)) %>%
      filter(between(age,-100,12400)) %>%
      arrange(age) # tso$paleoData_HoloceneValues
    ages<-sort(c(df$age[1],df$age[-1]-diff(df$age)/2,df$age[-nrow(df)]+diff(df$age)/2,df$age[nrow(df)]))
    values <- rep(df$values, each=2)
    plt <- ggplot() +
      geom_hline(yintercept=mean(df$values,na.rm=TRUE)) +
      geom_star(aes(x=df$age,y=df$values),fill=col,color='black',starshape=shp,size=1,starstroke=0.5,alpha=0.7) + 
      geom_path(aes(x=ages,y=values),color=col,alpha=0.7) + 
      scale_x_reverse(name = "Age (yr BP)", limits=c(12000,0),expand=c(0,0),n.breaks=7,oob=scales::squish) +
      theme_bw() +
      ggtitle(paste0(reg,": ",tso$dataSetName," (",tso$paleoData_TSid,")")) 
    if (tso$climateInterpretation1_interpDirection == 'negative'){
      plt <- plt + 
        geom_point(aes(x=tso$chronData_ages_12k,y=tso$chronData_ages_12k*0+max(df$values,na.rm=TRUE)), color="black",size=2,shape=17)+
        scale_y_reverse(name=paste(tso$paleoData_variableName,' (',tso$paleoData_units,')',sep=''))
    } else{
      plt <- plt + 
        geom_point(aes(x=tso$chronData_ages_12k,y=tso$chronData_ages_12k*0+min(df$values,na.rm=TRUE)), color="black",size=2,shape=17)+
        scale_y_continuous(name=paste(tso$paleoData_variableName,' (',tso$paleoData_units,')',sep=''))
    }
    if (length(tso$paleoData_valuesMax)>0){
        plt <- plt + geom_ribbon(aes(x= as.numeric(tso$age), ymin=tso$paleoData_valuesMin, ymax=tso$paleoData_valuesMax), fill='grey10',color=NA,alpha=0.1)
    }
    #
    map <- regMap + geom_star(data=siteDf,
                              aes(x=longitude.1, y=latitude.1), starshape=shp, fill=col,
                              size=3,color='gold',alpha=1,starstroke=2)
    h <- 0
    txt <- ggplot(siteDf) 
    for (name in c('geo_latitude','geo_longitude','geo_elevation',
                   'archiveType','Category','CategorySpecific',
                   'paleoData_proxyGeneral','paleoData_proxy','paleoData_proxyDetail','paleoData_variableName',
                   'climateInterpretation1_seasonalityGeneral',
                   'climateInterpretation1_variable','paleoData_units',
                   "chronData_agesN_12k" , "chronData_agesMaxGap_12k",
                   'Source','pub1_title','pub1_doi','pub2_doi','originalDataUrl')){
      h <- h-1
      #if (name %in% names(tso) == FALSE){tso[[name]]<- NA}
      txt <- txt + 
        annotate("text", x = 0, y = h, label = paste(name,': ',sep=''),hjust = 0,size=2.5) +
        annotate("text", x = 0.25, y = h, label = unique(tso[[name]]),hjust = 0,size=2.5)    
    }
    txt <- txt + scale_x_continuous(limits=c(0,0.6)) + theme_void() 
    bkg <- ggplot()+
      theme_void() +
      theme(plot.background = element_rect(fill = "White",color='White'))
    summary <- ggdraw(bkg) + 
      draw_plot(plt, x = 0,   y = 0.5, width = 1,   height = 0.5) +
      draw_plot(map, x = 0,   y = 0,   width = 0.5, height = 0.5) +
      draw_plot(txt, x = 0.5, y = 0,   width = 0.5, height = 0.5) 
    ggsave(file.path(wd,"Figures","Proxy","Dashboard",
                     paste(n,'_',reg,'_',tsn,'_',tso$paleoData_TSid,'.pdf',sep='')),device='pdf',
           plot=summary,width=10,height=6,units='in')
  }
}









