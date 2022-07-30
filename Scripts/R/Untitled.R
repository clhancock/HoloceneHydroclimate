library(gsheet)
df <- gsheet2tbl('docs.google.com/spreadsheets/d/1rhYoL0B5OfE5A-rNwuQZfnmjI3Vj3NCX07r-Mif3Ncs') %>%
  filter(inThisCompilation==TRUE)


PROJ <- '+proj=robin   +ellps=WGS84 +datum=WGS84 +no_defs +lon_0=0 +x_0=0 +y_0=0 +units=m'
PROJorig <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
load(url('https://github.com/SantanderMetGroup/ATLAS/blob/main/reference-regions/IPCC-WGI-reference-regions-v4_R.rda?raw=true'))
refregions <-  spTransform(IPCC_WGI_reference_regions_v4, CRSobj = PROJ)

pointData <- data.frame(longitude=c(df$lon),latitude=c(df$lat))
pointData <- spTransform(SpatialPointsDataFrame(coords=pointData,data = pointData, proj4string = CRS(PROJorig)),CRSobj = PROJ)
df$lonsPrj <- pointData@coords[,1]
df$latsPrj <- pointData@coords[,2]
df$region <- as.character(over(pointData, refregions)$Acronym)

for (row in 1:nrow(df)){
      archive  <- df[row,'archiveType'] 
      proxy    <- df[row,'proxy'] 
      unit     <- df[row,'units'] 
      if(is.null(df[row,'climateVariableDetail'])){df[row,'climateVariableDetail'] <-'Blank'} 
      if (is.na(proxy) | is.na(archive)){
        Category         <- 'Other'
        CategorySpecific <- 'Other (not calibrated)'
      } else if (archive == 'Speleothem'){
        Category           <- 'Speleothem'
        if (proxy == 'd18O' | proxy ==  'd13C'){
          CategorySpecific <- paste(archive,' (','δ',substring(proxy, 2),')',sep='')
        } else{
          CategorySpecific <- 'Speleothem (other)'
        }
      } else if ( df[row,'climateVariableDetail'] == 'LakeLevel@surface'){
        Category           <- 'Shoreline'
        CategorySpecific   <- 'Shoreline (Lake Level)'
      } else if (archive == 'GlacierIce'){
        Category           <- 'Glacier Ice'
        CategorySpecific   <- 'Glacier Ice'
      } else if (archive == 'LakeSediment' & proxy == 'd18O'){
        Category           <- 'Lake Sediment (δ18O)'
        CategorySpecific   <- 'Lake Sediment (δ18O)'
      } else if (proxy == 'dDwax'){
        Category           <- 'Leaf Wax'
        CategorySpecific   <- 'Leaf Wax (δD)'
      } else if (proxy == 'pollen'){
        Category           <- 'Pollen'
        if (is.null(unit)){
          CategorySpecific <- 'Pollen (not calibrated)'
        } else if (grepl('mm/',unit)){ 
          CategorySpecific <- 'Pollen (calibrated)'
        } else {
          CategorySpecific <- 'Pollen (not calibrated)'
        }
      } else {
        Category           <- 'Other'
        if (is.null(unit)){
          CategorySpecific <- 'Other (not calibrated)'
        } else if (grepl('mm/',unit)){ 
          CategorySpecific <- 'Other (calibrated)'
        } else {
          CategorySpecific <- 'Other (not calibrated)'
        }
      }
      df[row,'Category']         <- Category
      df[row,'CategorySpecific'] <- CategorySpecific
}


proxyMapSites <- basemap +
  #geom_map(data=refrenceSubset, map=dataTable, alpha=0.75, size=0.3, color='black' , 
           #aes(x=long, y=lat, group=group, map_id=id),fill=NA,linetype = "dashed") + 
  geom_star(data= as.data.frame(df%>%filter(CategorySpecific  =="Pollen (calibrated)")),aes(x=lonsPrj , y=latsPrj,
                                                  starshape=CategorySpecific, fill=CategorySpecific),
            size=1.5,color='Black',alpha=1,starstroke=0.3) + 
  geom_star(data= as.data.frame(df%>%filter(CategorySpecific  !="Pollen (calibrated)")),aes(x=lonsPrj , y=latsPrj,
                                                                                            starshape=CategorySpecific, fill=CategorySpecific),
            size=1.5,color='Black',alpha=1,starstroke=0.3) + 
  scale_fill_manual(values=plotSettings$color,name= 'Proxy Category') +
  scale_starshape_manual(values=plotSettings$shape,name= 'Proxy Category') +
  theme(text = element_text(family=figFont,size=figText),
        plot.background = element_rect(fill = 'white',color='white'),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "in"),
        legend.background =  element_rect(fill = alpha('white', 0.85),size=0.5),#alpha=0.5),
        legend.margin = margin(c(0, 0.1, 0.05, 0.06), unit="in"),
        legend.key.height = unit(0.12, "in"),
        legend.key.width = unit(0.12, "in"),
        legend.position = c(0.15,0.28),
        legend.title = element_blank())#,legend.position = 'none')# legend.position = c(0.5,0.7))

proxyMapSites

save <- TRUE
if (save) {
  ggsave(plot=proxyMapSites, width = 6.5, height = (6/1.97+0.5), dpi = 600,
         filename = file.path(dir,"Figures",paste('PlotProxyTypeRegion_',var,'.png',sep='')))
}


#Plot Figure 2 count by regions
regionsSelect <- vector(mode='list')
regionsSelect$all       <- refregions
regionsSelect$select    <- subset(refregions, Acronym %in% df$region)

dataTable <- fortify(regionsSelect$select)
dataTable$Count     <- NA
dataTable$CountName <- NA
dataTable$Region    <- NA

#Set Breaks and color scheme with lowest group as grey
figBreaks <- c(5,10,15,20,30,50,80)
figlabels <- c(paste("<",figBreaks[1]+1),
               paste(as.character(figBreaks[-length(figBreaks)]+1),"-",
                     as.character(figBreaks[-1])),
               paste(">",figBreaks[length(figBreaks)]))
figPalette <- c("#CCCCCC",RColorBrewer::brewer.pal(n=7, name = "YlOrRd"))

for (reg in list(levels(regionsSelect$select@data[["Acronym"]]))[[1]]){
  i = which(refregions@data[["Acronym"]]==reg)
  count <- length(which(df$region==reg))
  if (count == 0){         countName <- NA
  } else if (count <= figBreaks[1]){ countName <- 0
  } else if (count <= figBreaks[2]){ countName <- 1
  } else if (count <= figBreaks[3]){ countName <- 2
  } else if (count <= figBreaks[4]){ countName <- 3
  } else if (count <= figBreaks[5]){ countName <- 4
  } else if (count <= figBreaks[6]){ countName <- 5
  } else if (count <= figBreaks[7]){ countName <- 6
  } else {countName <- 7}
  dataTable$Count[which(dataTable$id==reg)]     <- count
  dataTable$CountName[which(dataTable$id==reg)] <- countName
  dataTable$Region[which(dataTable$id==reg)]    <-reg
  print(reg)
  print(count)
}

#Plot
proxyMapRegions <- basemap + 
  #Add refrence regions boundaries
  geom_map(data=regionsSelect$select, map=dataTable, alpha=0.75, size=0.5, color='black' , 
           aes(x=long, y=lat, group=group, map_id=id,fill=as.factor(dataTable$CountName))) +
  scale_fill_manual(values=figPalette,labels=figlabels)+
  #Add proxy sites
  geom_point(data=as.data.frame(df), aes(x=lonsPrj , y=latsPrj), shape=1,size=0.8) +
  #Format Legend
  theme(plot.background = element_rect(fill = 'white',color='white'),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "in"),
        legend.position = c(0.5,-0.01),
        legend.direction = "horizontal",
        legend.title=element_blank(),
        legend.background=element_rect(colour='black',size=0.4),
        legend.margin=margin(t = 1, r = 1, b = 3, l = -4.5, unit = "pt"),
        legend.key.width  = (unit(35, 'pt')),
        legend.key.height = (unit(5, 'pt')),
        legend.text = element_text(family='sans',size=8))

proxyMapRegions

if (save) {
  ggsave(plot=proxyMapRegions, width = 6, height = (6/1.97+0.5), dpi = 600,
         filename = file.path(dir,"Figures",paste('PlotCountByRegion_',var,'.png',sep='')))
}

