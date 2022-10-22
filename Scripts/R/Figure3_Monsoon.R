insolation <- as.data.frame(read.table(file.path(dir,'Data','insolation_bein1.dat.txt')))
insolation$V1 <- insolation$V1*-1000

insolation30 <- insolation[which(insolation$V2==30),]
insolation30 <- data.frame(age=insolation30$V1,
                           DJF=apply(insolation30[c(1,2,12)+2],1,sum),
                           JJA=apply(insolation30[c(6,7,8)+2],1,sum),
                           ANN=apply(insolation30[seq(1,12)+2],1,sum))[1:13,]
for (i in 1:ncol(insolation30)){
  insolation30[,i] <- insolation30[,i] - insolation30[1,i] 
}

insolationNeg30 <- insolation[which(insolation$V2==-30),]
insolationNeg30 <- data.frame(age=insolationNeg30$V1,
                              ANN=apply(insolationNeg30[seq(1,12)+2],1,sum),
                              DJF=apply(insolationNeg30[c(1,2,12)+2],1,sum),
                              JJA=apply(insolationNeg30[c(6,7,8)+2],1,sum))[1:13,]
for (i in 1:ncol(insolationNeg30)){
  insolationNeg30[,i] <- insolationNeg30[,i] - insolationNeg30[1,i] 
}

var <- 'HC'
RegionTS <- read.csv(file.path(dir,'Data','RegionComposites',var,'MedianTS_byRegion.csv'))
for (i in 2:ncol(RegionTS)){
  RegionTS[,i] <- (RegionTS[,i] - mean(RegionTS[1:10,i],na.rm=TRUE))#/sd(RegionTS[,i],na.rm=TRUE) 
}
NHmonsoon <- RegionTS[,c('EAS','TIB','SAH','NEAF','ECA','WCA')]
SHmonsoon <- RegionTS[,c('SAM','ESAF')]


scale <- 0.03
shift <- 0
#The patterns described in this paragraph (wet-dry-wet, opposite trends between hemisphere, and 6 ka shift) are difficult to see in the figure. I suggest adding a summary figure with all of the reconstructions plotted together along with the insolation forcing and ITCZ shift-DK
s<-0.75
fig<-ggplot(RegionTS)+
  geom_hline(yintercept =0,color='grey70')+
  geom_line(aes(x=time,y=EAS,color='EAS',linetype='North'),size=s) +
  geom_line(aes(x=time,y=TIB,color='TIB',linetype='North'),size=s) +
  geom_line(aes(x=time,y=SAH,color='SAH',linetype='North'),size=s) +
  geom_line(aes(x=time,y=NEAF,color='NEAF',linetype='North'),size=s) +
  geom_line(aes(x=time,y=SAM,color='SAM',linetype='South'),size=s) +
  geom_line(aes(x=time,y=ESAF,color='ESAF',linetype='South'),size=s) +
  geom_smooth(data=insolation30,aes(   x=age,y=(DJF)*scale+shift,fill="30\u00B0S DJF"),alpha=1,color="Grey50",linetype=4,size=s) +
  geom_smooth(data=insolationNeg30,aes(x=age,y=(JJA)*scale+shift,fill="30\u00B0N JJA"),alpha=1,color="Black",linetype=1,size=s) +
  scale_x_reverse(limits = c(30000,0),breaks = seq(0, 12000, 2000),labels=seq(0, 12, 2),name="Age (ka BP)") +
  scale_y_continuous(limits=c(-10,10),name="Proxy Composite Anomaly",sec.axis = sec_axis((~./scale),name=bquote('Insolation Anomaly (W m'^-2*")"))) + 
  coord_cartesian(xlim=c(12000,0), ylim=c(-4,4),expand	=FALSE) +
  scale_linetype_manual(values=c(1,4,1,1,1,1,4,1,1,4),name='Hemisphere') +
  scale_fill_manual(values=c("White","White",rep(NA,8)),name='Insolation') +
  scale_color_manual(values=c('#fcde9c','#faa476','#f0746e','#dc3977','#b9257a','#7c1d6f','black','grey'),name='Region')+#c(RColorBrewer::brewer.pal(12,"Paired")[seq(2,12,2)]),name='Region') +
  guides(fill  = guide_legend(override.aes = list(linetype = c(1,4),color=c('Black','Grey50'))),
         color = guide_legend(override.aes = list(linetype = c(1,4,1,1,4,1))),
         linetype = guide_legend(override.aes = list(linetype = c(1,4))))+
  theme_bw() + 
  theme(legend.key = element_rect(fill=NA,color = "white"),
        legend.key.height =  unit(0.15, 'in'),
        legend.spacing.y = unit(0.02, 'in'),
        panel.grid = element_blank(),
        text = element_text(family=figFont,size=8))
fig
#fcde9c,#faa476,#f0746e,#e34f6f,#dc3977,#b9257a,#7c1d6f
ggsave(plot=fig, width = 4.5, height = 3, dpi = 600,
       filename = file.path(dir,'Figures','RegionComposites',paste('Monsoon_',var,'.png',sep='')))
