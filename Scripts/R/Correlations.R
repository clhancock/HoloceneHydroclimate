library(dplyr)
library(geoChronR)

dataDir<-getwd()
hcNames <- names(read.csv(file.path(dataDir,'Data','RegionComposites','HC','MedianTSbyRegion.csv')))
tNames  <- names(read.csv(file.path(dataDir,'Data','RegionComposites','T','MedianTSbyRegion.csv')))

time <- as.matrix(read.csv(file.path(dataDir,'Data','RegionComposites','HC','MedianTSbyRegion.csv'))[,'time'])
regNames <- hcNames[which(hcNames %in% tNames)][-1]

for (reg in regNames){
  HCseries <- read.csv(file.path(dataDir,'Data','RegionComposites','HC',paste(reg,'.csv',sep='')))
  T_series <- read.csv(file.path(dataDir,'Data','RegionComposites','T',paste(reg,'.csv',sep='')))
  HCseries <- as.matrix(HCseries[which(between(time,0,12000)),])
  T_series <- as.matrix(T_series[which(between(time,0,12000)),])
  time <- as.matrix(time[which(between(time,0,12000)),])
  z= corEns(time.1=time,values.1=HCseries,time.2=time,values.2=T_series,bin.step=100)#,seq(0,12000,100),100,1000)
  
  
}