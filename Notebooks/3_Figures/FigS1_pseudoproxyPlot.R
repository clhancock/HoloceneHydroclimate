dir = '/Users/chrishancock/Library/CloudStorage/OneDrive-NorthernArizonaUniversity/Research/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate'
var <- 'pre'
data <- read.csv(file.path(dir,'Data','Model','pseudoProxyCorr_pre.csv'))
data3 <- read.csv(file.path(dir,'Data','Model','pseudoProxyCorr_pre_byCount.csv'))
#data <- read.csv(file.path(dir,'Data','Model','pseudoProxyCorr_pre_NoSzn.csv'))
#data <- read.csv(file.path(dir,'Data','Model','pseudoProxyCorr_pre_OnlyGeo.csv'))


Chadcm <- '#DDAA33'
Ctrace <- '#BB5566'
    
lineX <- sort(unique(data$X1))
lineH<-c()
lineT<-c()
for (i in lineX){
    lineH <- append(lineH,median(data$X3[which(data$X1==i)]))
    lineT <- append(lineT,mean(data$X2[which(data$X1==i)]))  
}
data2=data.frame(lineX,lineT,lineH)
  
plt1 <- ggplot(data=data) + 
    #geom_smooth(aes(X1,apply(data[,4:5],1,mean)),color='black')+
    geom_line(data=data3,aes(X0,X1),color=Ctrace,size=1,linetype=1)+
    geom_line(data=data3,aes(X0,X2),color=Chadcm,size=1,linetype=1)+
    geom_line(data=data3,aes(X0,X3),color=Ctrace,size=1,linetype=6)+
    geom_line(data=data3,aes(X0,X4),color=Chadcm,size=1,linetype=6)+
    #geom_line(data=data2,aes(lineX,lineT),color=Ctrace)+  
    #geom_line(data=data2,aes(lineX,apply(data2[,2:3],1,mean)),color='black')+ 
    scale_fill_manual('Key',breaks=c('TraCE', 'HadCM'),values=c(Ctrace,Chadcm)) +
    geom_hline(yintercept=0,size=0.4)+
    geom_point(aes(X1,X2,fill='TraCE'),color='Black',size=4,stroke=0.7,shape=21)+
    geom_point(aes(X1,X3,fill='HadCM'),color='Black',size=4,stroke=0.7,shape=21)+
    theme_bw() +
    #geom_vline(xintercept=5.5) +
    annotate(geom = "rect", colour = NA, fill = "grey", xmin = 0, xmax = 5.5, ymin = -1, ymax = 1,alpha = 0.3) +
    scale_y_continuous("Pearson's Correlation Coefficient",limits=c(-1,1),expand=c(0,0),position = "left") +
    #ggtitle("Based on geographic distribution + season + p vs p-e")+
    scale_x_continuous("Number of Proxy Sites",limits=c(0,25),expand=c(0,0)) +
    theme(text = element_text(family='Times New Roman',size=10),
          legend.position = c(0.8,0.2),
          legend.title = element_blank(),
          legend.box.background=element_rect(fill='White',color='Black'))
plt1
  
  
ggsave(plot=plt1, width = 5.5, height = 4.5, dpi = 600,device='png',
       filename = paste(file.path(dir,'Figures','Model','pseudoproxyPlot.png',sep='')))
