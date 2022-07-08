dataDir <- '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate' #
var <- 'p-e'
data <- read.csv(file.path(dataDir,'Data','Model','pseudoProxyCorr_p-e.csv'))

Chadcm <- '#AB8599'
Ctrace <- '#99AB85'

plt1 <- ggplot(data=data) + 
 geom_point(aes(X1,X2,fill='TraCE'),color='Black',size=4,stroke=0.7,shape=21)+
 geom_point(aes(X1,X3,fill='HadCM'),color='Black',size=4,stroke=0.7,shape=21)+
 scale_fill_manual('Key',breaks=c('TraCE', 'HadCM'),values=c(Ctrace,Chadcm)) +
 theme_bw() +
 geom_vline(xintercept=6) +
 annotate(geom = "rect", colour = NA, fill = "grey", xmin = 0, xmax = 6, ymin = -1, ymax = 1,alpha = 0.3) +
 scale_y_continuous("Pearson's Correlation Coefficient",limits=c(-1,1),expand=c(0,0),position = "left") +
 scale_x_continuous("Number of Proxy Sites",limits=c(0,90),expand=c(0,0)) +
 theme(text = element_text(family='sans',size=8),
       legend.position = c(0.8,0.2),
       legend.title = element_blank(),
       legend.box.background=element_rect(fill='White',color='Black'))

plt2 <- ggplot(data=data) + 
 geom_point(aes(X1,X4,color='Max',fill='Max'),size=4,shape=21)+
 geom_point(aes(X1,X5,color='Average',fill='Average'),size=4,shape=1)+
 scale_color_manual('Key',breaks=c('Max', 'Average'),values=c('Black','Black')) +
 scale_fill_manual('Key',breaks=c('Max', 'Average'),values=c('Black','White')) +
 theme_bw() +
 geom_vline(xintercept=6) +
 annotate(geom = "rect", colour = NA, fill = "grey", xmin = 0, xmax = 6, ymin = -1, ymax = 1,alpha = 0.3) +
 scale_y_continuous("Pearson's Correlation Coefficient",limits=c(-1,1),expand=c(0,0),position = "right") +
 scale_x_continuous("Number of Proxy Sites",limits=c(0,90),expand=c(0,0)) +
 theme(text = element_text(family='sans',size=8),
       legend.position = c(0.8,0.2),
       legend.title = element_blank(),
       legend.box.background=element_rect(fill='White',color='Black'))



plt <- ggdraw() + 
 draw_plot(plt1,x=0., y=0, width = 0.48, height =1) +
 draw_plot(plt2,x=0.52, y=0, width = 0.48, height =1) +
 draw_label( "(a)", x = 0.007, y = 0.97,fontfamily='sans',color='Black',size =8,fontface='bold',hjust = 0)+
 draw_label( "(b)", x = 0.49, y = 0.97,fontfamily='sans',color='Black',size =8,fontface='bold',hjust = 0)


ggsave(plot=plt, width = 6.5, height = 3.75, dpi = 600,device='png',
       filename = paste(file.path(dataDir,'Figures','pseudoproxyPlot.png',sep='')))
