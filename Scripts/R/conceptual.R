library(ncdf4) 
library(ggplot2)

dataDir <-getwd()
data <- vector(mode='list')

data("wrld_simpl")
PROJorig <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

#Load Insolation Data
data[['w/m2']] <- as.data.frame(read.table(file.path(dataDir,'Data','insolation_bein1.dat.txt')))
names(data[['w/m2']]) <-c('age','lat',c(1:12))
diff <- data[['w/m2']][which(data[['w/m2']]$age==-6),] - data[['w/m2']][which(data[['w/m2']]$age==0),]
data[['w/m2DF']] <- data.frame(lat=data[['w/m2']][which(data[['w/m2']]$age==0),'lat'],
                                ANN=apply(diff[,which(names(diff) %in% c(1:12))],1,sum),
                                JJA=apply(diff[,which(names(diff) %in% c(1:2,12))],1,sum),
                                DJF=apply(diff[,which(names(diff) %in% c(6:8))],1,sum))

proxyT <- read.csv(file.path(dataDir,'Data','temp12k_allmethods_percentiles.csv'))
data[['proxyT']] <- apply(proxyT[which(between(proxyT$ages,5500,6500)),],2,mean) -
                    apply(proxyT[which(between(proxyT$ages,   0,1000)),],2,mean)
data[['proxyT']] <- data.frame(lat=seq(-90,60,30)+15,
                               tas=rev(data[['proxyT']][which(grepl('median',names(data[['proxyT']])))[-1]]))

#Load Climate Data
cmip <- nc_open(file.path(dataDir,'Data','Model','cmip6_ANN.nc'))  
lats <- ncvar_get(cmip,paste('lat','regrid',sep='_'))
lons <- ncvar_get(cmip,paste('lon','regrid',sep='_'))
for (var in c('tas','p-e','pre')){
  data[[var]] <- vector(mode='list')
  vals <- ncvar_get(cmip,paste(var,'regrid',sep='_'))
  if (var == 'tas'){vals <- vals+273.15}
  valL <- (apply((vals[,,,1]-vals[,,,2]),c(1,2),mean))
  valL <- raster(apply(t(as.matrix(valL)),2,rev),crs=PROJorig)
  extent(valL) <- extent(c(0, 360, -90,90))
  valL <- rotate(valL)
  valL <- as.data.frame(rasterToPoints(mask(valL,wrld_simpl))) 
  valL <- valL %>% group_by(y) %>% summarize(v= mean(layer))
  diff <- apply((vals[,,,1] - vals[,,,2]),c(1,2),mean)
  df <- data.frame(lat=lats,ka0 =apply(vals[,,,2],2,mean), ka6 =apply(vals[,,,2],2,mean),
                   diff=(apply((vals[,,,1]-vals[,,,2]),2,mean)),
                   diffLand = rep(NA,length(lats)),
                   pctD=100*apply((vals[,,,1]-vals[,,,2])/((vals[,,,1]+vals[,,,2])/2),2,mean),
                   pctDLand = rep(NA,length(lats)))
  for (y in valL$y){df[which(df$lat == y),'diffLand'] = valL$v[which(valL$y==y)]}
  #df[-1] <- df[-1]-colMeans(df[-1],na.rm=TRUE)[col(df[-1])]
  data[[paste(var,'DF')]] <- df
}

base <- ggplot() +
  scale_x_continuous(limits=c(-90,90), expand=c(0,0),breaks=seq(-90,90,30)) +
  scale_y_continuous(limits=c(0,1),expand=c(0,0)) +
  theme_void() + 
  theme(legend.position = 'none',
        panel.grid.major.x = element_line(color='Grey',size=0.1),
        panel.border = element_rect(color='Black',fill=NA,size=1),
        plot.margin = unit(rep(0.05 ,4), "in")
  )


p<-vector(mode='list')
p[['insolation']] <- base + 
  geom_smooth(data=data[['w/m2DF']],aes(x=lat,y=ANN),color='Black',size=1.2,se=FALSE) + 
  geom_smooth(data=data[['w/m2DF']],aes(x=lat,y=JJA),color='Blue3',size=1.2,se=FALSE) + 
  geom_smooth(data=data[['w/m2DF']],aes(x=lat,y=DJF),color= 'Red3',size=1.2,se=FALSE) + 
  geom_hline(yintercept=0) + 
  scale_y_continuous(limits=c(-150,300),expand=c(0,0)) 


for (var in c('p-e','tas')){
  if (var == 'tas'){
    pal <- rev(RColorBrewer::brewer.pal(11, 'RdBu')[-c(1,5,7,11)])
    dir <- -1
  } else{
    pal <- RColorBrewer::brewer.pal(11, 'BrBG')[c(2,2,3,3,4,6,8,9,9,10,10)]
    dir <- 1
  } 
  p[[paste(var,'Diff')]] <- base + 
    geom_bar(data=data[[paste(var,'DF')]],aes(x=lat,y=1,fill=diff),stat="identity",width=2,alpha=1) +
    geom_bar(data=data[[paste(var,'DF')]],aes(x=lat,y=1,fill=diffLand),stat="identity",width=2,alpha=1) +
    scale_fill_gradientn(colors=(pal),
                         limits=c(-1,1)*max(abs(data[[paste(var,'DF')]]$diffLand),na.rm=TRUE),
                         na.value=NA)
  #p[[paste(var,'0ka')]] <- base + 
  #  geom_bar(data=data[[paste(var,'DF')]],aes(x=lat,y=1,fill=ka0),stat="identity",width=2,alpha=1) +
  #  scale_fill_distiller(palette=pal,direction=dir,na.value='grey90')
 # p[[paste(var,'6ka')]] <- base + 
  #  geom_bar(data=data[[paste(var,'DF')]],aes(x=lat,y=1,fill=ka6),stat="identity",width=2,alpha=1) +
  #  scale_fill_distiller(palette=pal,direction=dir,na.value='grey90')    
}



proxy <- TRUE
if (proxy){
  p[[paste('tas','Diff')]] <- p[[paste('tas','Diff')]] + 
    geom_bar(data=data[['proxyT']],aes(x=lat,y=0.5,fill=tas),stat="identity",width=30,alpha=1) + 
    geom_hline(yintercept = 0.5) 
}


valT <- apply(ncvar_get(cmip,paste('tas','regrid',sep='_'))[,,,2],c(1,2),mean)
valT <- raster(apply(t(as.matrix(valT)),2,rev),crs=PROJorig)
extent(valT) <- extent(c(0, 360, -90,90))
valT <- as.data.frame(rasterToPoints(rotate(valT))) 

pltT <- base + 
  scale_y_continuous(limits=c(-180,180),expand=c(0,0)) +
  scale_x_continuous(limits=c(-90,90), expand=c(0,0),breaks=seq(-90,90,30)) +
  geom_raster(data = valT, aes(x = y, y = x, fill = layer)) +
  scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(11, 'RdBu')[-c(1,5,7,11)]), 
                       limits=c(-1,1)*max(abs(valT[,3])),
                       na.value=NA)

valHC <- apply( ncvar_get(cmip,paste('p-e','regrid',sep='_'))[,,,2],c(1,2),mean)
valHC <- raster(apply(t(as.matrix(valHC)),2,rev),crs=PROJorig)
extent(valHC) <- extent(c(0, 360, -90,90))
valHC <- as.data.frame(rasterToPoints(mask(rotate(valHC),wrld_simpl)))

library(rworldmap)
data(coastsCoarse)


pltP <- base + scale_y_continuous(limits=c(-180,180),expand=c(0,0)) +
  geom_raster(data = valHC, aes(x = y, y = x, fill = layer)) +
  scale_x_continuous(limits=c(-90,90), expand=c(0,0),breaks=seq(-90,90,30)) +
  scale_fill_gradientn(colors=RColorBrewer::brewer.pal(11, 'BrBG')[c(2,2,3,3,4,6,8,9,9,10,10)],
                       limits=c(-1,1)*max(abs(valHC[,3])),
                       na.value=NA) 


pp <- ggplot() +
  geom_point(aes(x=0,y=0),color='white') + 
  scale_x_continuous(limits=c(-90,90), expand=c(0,0),breaks=seq(-90,90,30),labels=seq(-90,90,30),sec.axis = dup_axis()) +
  geom_hline(yintercept = 0.5) + scale_y_continuous(limits=c(0,1)) +
  theme_void() +
  theme(axis.line.x = element_line(color='black',size=1),
        axis.ticks.x = element_line(color='black',size=0.5),
        axis.ticks.length.x=unit(3,"pt"),
        axis.text.x = element_text(size=8, family='sans'),
        plot.margin = unit(rep(0.1 ,4), "in")
        )


library(cowplot)
map <- ggdraw(ggplot()+theme_void()+
                theme(plot.background=element_rect(colour='White',fill='White'),
                      plot.margin = unit(rep(0 ,4), "in"))) +
  draw_plot(pltT,x=0, y=0, width = 1, height = 1) +
  draw_plot(pltP,x=0, y=0, width = 1, height = 1)

plt <- ggarrange(pp,p[['p-e Diff']],p[['tas Diff']],p[['insolation']],
                 ncol=1,heights=c(6,1,1,2))
ggsave(plot=plt, width = 6, height = 4, dpi = 400,device='png',
       filename = paste(file.path(dataDir,'Figures','conceptual.png',sep='')))


