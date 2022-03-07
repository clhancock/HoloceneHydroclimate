library(dplyr)
library(geoChronR)

dataDir<-getwd()
hcNames <- names(read.csv(file.path(dataDir,'Data','RegionComposites','HC','MedianTSbyRegion.csv')))
tNames  <- names(read.csv(file.path(dataDir,'Data','RegionComposites','T','MedianTSbyRegion.csv')))
time <- as.matrix(read.csv(file.path(dataDir,'Data','RegionComposites','HC','MedianTSbyRegion.csv'))[,'time'])
regNames <- hcNames[which(hcNames %in% tNames)][-1]

x = 250
Proxy_HC_T_Rvals = matrix(nrow= x*x,ncol=length(regNames))
for (reg in regNames){
  HCseries <- read.csv(file.path(dataDir,'Data','RegionComposites','HC',
                                 paste(reg,'.csv',sep='')))
  T_series <- read.csv(file.path(dataDir,'Data','RegionComposites','T',
                                 paste(reg,'.csv',sep='')))
  HCseries <- as.matrix(HCseries[which(between(time,0,12000)),])
  T_series <- as.matrix(T_series[which(between(time,0,12000)),])
  time <- as.matrix(time[which(between(time,0,12000)),])
  data <- corEns(time.1=time, values.1=HCseries[,1:x],
                 time.2=time, values.2=T_series[,1:x],
                 bin.step=100)#,seq(0,12000,100),100,1000)
  Proxy_HC_T_Rvals[,which(regNames==reg)] = as.numeric(data$cor.df$r)
  
}

z <- as.data.frame(Proxy_HC_T_Rvals)
names(z) <- regNames
write.csv(z,
          file=file.path(dataDir,'Data',
                         paste('HC_T_RegionalProxyEnsCorrelations.csv',sep='')))
