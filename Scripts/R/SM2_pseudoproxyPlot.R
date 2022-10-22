dir <- '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/2021_HoloceneHydroclimate/2021_HoloceneHydroclimate' #
var <- 'pre'
data <- read.csv(file.path(dir,'Data','Model','pseudoProxyCorr_pre.csv'))

Chadcm <- '#DDAA33'
Ctrace <- '#BB5566'

plt1 <- ggplot(data=data) + 
 geom_point(aes(X1,X2,fill='TraCE'),color='Black',size=4,stroke=0.7,shape=21)+
 geom_point(aes(X1,X3,fill='HadCM'),color='Black',size=4,stroke=0.7,shape=21)+
 scale_fill_manual('Key',breaks=c('TraCE', 'HadCM'),values=c(Ctrace,Chadcm)) +
 theme_bw() +
 geom_vline(xintercept=6) +
 annotate(geom = "rect", colour = NA, fill = "grey", xmin = 0, xmax = 6, ymin = -1, ymax = 1,alpha = 0.3) +
 scale_y_continuous("Pearson's Correlation Coefficient",limits=c(-1,1),expand=c(0,0),position = "left") +
 scale_x_continuous("Number of Proxy Sites",limits=c(0,100),expand=c(0,0)) +
 theme(text = element_text(family='Times New Roman',size=10),
       legend.position = c(0.8,0.2),
       legend.title = element_blank(),
       legend.box.background=element_rect(fill='White',color='Black'))
plt1


ggsave(plot=plt1, width = 5.5, height = 4.5, dpi = 600,device='png',
       filename = paste(file.path(dir,'Figures','Model','pseudoproxyPlot.png',sep='')))

"“I don't really like this argument.  Greater model uncertainty doesn't mean that there's no systematic bias, if the bias is in the proxies.”-mpe"