library(dplyr)
library(geoChronR)

wd<-getwd()
wd = '/Users/chrishancock/Library/CloudStorage/OneDrive-NorthernArizonaUniversity/Research/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate'

hcNames <- names(read.csv(file.path(wd,'Data','RegionComposites','HC','MedianTS_byRegion.csv')))
tNames  <- names(read.csv(file.path(wd,'Data','RegionComposites','T','MedianTS_byRegion.csv')))
time <- as.matrix(read.csv(file.path(wd,'Data','RegionComposites','HC','MedianTS_byRegion.csv'))[,'time'])
regNames <- hcNames[which(hcNames %in% tNames)][-1]

x = 500
Proxy_HC_T_Rvals = matrix(nrow= x*x,ncol=length(regNames))
for (reg in c('ENA')){
  print(reg)
  HCseries <- read.csv(file.path(wd,'Data','RegionComposites','HC',
                                 paste(reg,'.csv',sep='')))
  T_series <- read.csv(file.path(wd,'Data','RegionComposites','T',
                                 paste(reg,'.csv',sep='')))
  HCseries <- as.matrix(HCseries[which(between(time,0,12000)),])
  T_series <- as.matrix(T_series[which(between(time,0,12000)),])
  time <- as.matrix(time[which(between(time,0,12000)),])
  data <- corEns(time.1=time, values.1=HCseries[,1:x],
                 time.2=time, values.2=T_series[,1:x],
                 bin.step=100)#,seq(0,12000,100),100,1000)
  Proxy_HC_T_Rvals[,which(regNames==reg)] = as.numeric(data$cor.df$r)
  
}

correlationEnsemble <- as.data.frame(Proxy_HC_T_Rvals)
names(correlationEnsemble) <- regNames
write.csv(correlationEnsemble,
          file=file.path(wd,'Data',
                         paste('HC_T_RegionalProxyEnsCorrelations.csv',sep='')))

#read.csv()