dir <- '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate' #
var <- 'pre'

#Load Data
proxyDataAgree <- read.csv(file.path(dir,'Data','Proxy','proxyMetaData_HC.csv'))
#proxyRegionTS <- read.csv(file.path(dir,'Data','RegionComposites','HC','MedianTSbyRegion.csv'))
#hadcmRegionTS <- read.csv(file.path(dir,'Data','Model','RegionTS',paste('regional_',var,'_ANN_hadcm.csv',sep='')))
#traceRegionTS <- read.csv(file.path(dir,'Data','Model','RegionTS',paste('regional_',var,'_ANN_trace.csv',sep='')))
#modelRegionTS <- 

#Create Empty DF
regNames <-  names(proxyRegionTS[-1])
NAs <- rep(NA,length(regNames))
data <- data.frame(regions=regNames,
                   proxyAgreePct=NAs,
                   proxyhadcmCor=NAs,
                   proxytraceCor=NAs,
                   hadcmVariance=NAs,
                   traceVariance=NAs,
                   hadcmRange=NAs,
                   traceRange=NAs)

#Populate DF
for (reg in regNames){
  row <- which(regNames==reg)
  #
  proxySub  <- proxyDataAgree[which(proxyData$ipccReg == reg),]
  proxyVals <- proxySub$ka_6-proxySub$ka_0.5
  for (i in (which(proxySub$direction=='negative'))){proxyVals[i] <- proxyVals[i]*-1}
  proxyVals <- proxyVals[which(!is.na(proxyVals))]
  proxyVals <- proxyVals[which(proxyVals != 0)]
  data[row,'regions'] <- reg
  data[row,'proxyAgreePct'] <- abs(round(100*length(which(proxyVals>0))/length(proxyVals),1)-50)+50
  #proxy <- proxyRegionTS[2:122,reg]
  #data[row,'proxyhadcmCor'] <- cor(proxy[which(!is.na(proxy))],hadcmRegionTS[2:122,reg][which(!is.na(proxy))])
  #data[row,'proxytraceCor'] <- cor(proxy[which(!is.na(proxy))],traceRegionTS[2:122,reg][which(!is.na(proxy))])
  #data[row,'hadcmVariance'] <- var(hadcmRegionTS[2:122,reg])
  #data[row,'traceVariance'] <- var(traceRegionTS[2:122,reg])
  #data[row,'hadcmRange'] <- diff(range(hadcmRegionTS[2:122,reg],na.rm=TRUE))
  #data[row,'traceRange'] <- diff(range(traceRegionTS[2:122,reg],na.rm=TRUE))
  #data[row,'hadcmRange'] <- abs(mean(hadcmRegionTS[57:67,reg]) - mean(hadcmRegionTS[2:12,reg]))
  #data[row,'traceRange'] <- abs(mean(traceRegionTS[57:67,reg]) - mean(traceRegionTS[2:12,reg]))
}

plt1 <- ggplot(data=data) + 
  geom_point(aes(traceRange,proxyAgreePct,fill='TraCE'),color='Black',size=4,stroke=0.7,shape=21)+
  geom_point(aes(hadcmRange,proxyAgreePct,fill='HadCM'),color='Black',size=4,stroke=0.7,shape=21)+
  scale_fill_manual('Key',breaks=c('TraCE', 'HadCM'),values=c(Ctrace,Chadcm)) +
  theme_bw() +
  labs(y= "% of Proxy Records Indicating the Same Direction of Change", x = "Simulated Absolute Difference Between 6 and 0.5ka")+
  theme(text = element_text(family='sans',size=8),
        legend.position = c(0.8,0.2),
        legend.title = element_blank(),
        legend.box.background=element_rect(fill='White',color='Black'))
plt1



###### IGNORE BELOW



regNames <-  names(proxyRegionTS[-1])
NAs <- rep(NA,length(regNames))
data <- data.frame(regions=regNames,
                   proxyAgreePct=NAs,
                   proxyhadcmCor=NAs,
                   proxytraceCor=NAs,
                   hadcmVariance=NAs,
                   traceVariance=NAs,
                   hadcmRange=NAs,
                   traceRange=NAs)

for (reg in regNames){
 row <- which(regNames==reg)
 proxySub  <- proxyDataAgree[which(proxyData$ipccReg == reg),]
 proxyVals <- proxySub$ka_6-proxySub$ka_0.5
 for (i in (which(proxySub$direction=='negative'))){proxyVals[i] <- proxyVals[i]*-1}
 proxyVals <- proxyVals[which(!is.na(proxyVals))]
 proxyVals <- proxyVals[which(proxyVals != 0)]
 data[row,'regions'] <- reg
 data[row,'proxyAgreePct'] <- abs(round(100*length(which(proxyVals>0))/length(proxyVals),1)-50)+50
 proxy <- proxyRegionTS[2:122,reg]
 data[row,'proxyhadcmCor'] <- cor(proxy[which(!is.na(proxy))],hadcmRegionTS[2:122,reg][which(!is.na(proxy))])
 data[row,'proxytraceCor'] <- cor(proxy[which(!is.na(proxy))],traceRegionTS[2:122,reg][which(!is.na(proxy))])
 data[row,'hadcmVariance'] <- var(hadcmRegionTS[2:122,reg])
 data[row,'traceVariance'] <- var(traceRegionTS[2:122,reg])
 #data[row,'hadcmRange'] <- diff(range(hadcmRegionTS[2:122,reg],na.rm=TRUE))
 #data[row,'traceRange'] <- diff(range(traceRegionTS[2:122,reg],na.rm=TRUE))
 data[row,'hadcmRange'] <- abs(mean(hadcmRegionTS[57:67,reg]) - mean(hadcmRegionTS[2:12,reg]))
 data[row,'traceRange'] <- abs(mean(traceRegionTS[57:67,reg]) - mean(traceRegionTS[2:12,reg]))
}

Chadcm <- '#AB8599'
Ctrace <- '#99AB85'

plt1 <- ggplot(data=data) + 
 geom_point(aes(traceRange,proxyAgreePct,fill='TraCE'),color='Black',size=4,stroke=0.7,shape=21)+
 geom_point(aes(hadcmRange,proxyAgreePct,fill='HadCM'),color='Black',size=4,stroke=0.7,shape=21)+
 scale_fill_manual('Key',breaks=c('TraCE', 'HadCM'),values=c(Ctrace,Chadcm)) +
 theme_bw() +
 labs(y= "% of Proxy Records Indicating the Same Direction of Change", x = "Simulated Absolute Difference Between 6 and 0.5ka")+
 theme(text = element_text(family='sans',size=8),
       legend.position = c(0.8,0.2),
       legend.title = element_blank(),
       legend.box.background=element_rect(fill='White',color='Black'))
plt1
plt2 <- ggplot(data=data) + 
 geom_point(aes(traceVariance,proxytraceCor,fill='TraCE'),color='Black',size=4,stroke=0.7,shape=21)+
 geom_point(aes(hadcmVariance,proxyhadcmCor,fill='HadCM'),color='Black',size=4,stroke=0.7,shape=21)+
 scale_fill_manual('Key',breaks=c('TraCE', 'HadCM'),values=c(Ctrace,Chadcm)) +
 #scale_y_continuous(position = "right") +
 theme_bw() +
 labs(y= "Correlation Between Proxy and Model Time Series", x = "Simulated Variance") +
 theme(text = element_text(family='sans',size=8),
       legend.position = c(0.8,0.2),
       legend.title = element_blank(),
       legend.box.background=element_rect(fill='White',color='Black'))


plt <- ggdraw() + 
 draw_plot(plt1,x=0., y=0, width = 0.5, height =1) +
 draw_plot(plt2,x=0.5, y=0, width = 0.5, height =1) +
 draw_label( "(a)", x = 0.007, y = 0.97,fontfamily='sans',color='Black',size =8,fontface='bold',hjust = 0)+
 draw_label( "(b)", x = 0.507, y = 0.97,fontfamily='sans',color='Black',size =8,fontface='bold',hjust = 0)
  

ggsave(plot=plt, width = 6.5, height = 3.75, dpi = 600,device='png',
       filename = paste(file.path(dir,'Figures','magnitudePlot.png',sep='')))

