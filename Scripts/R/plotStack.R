library(dplyr)
library(geoChronR)
library(ggplot2)
library(lipdR)
library(magrittr)

dataDir <- '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate' #
#D <-  readLipd(file.path(dataDir,'Data','PlotStack'))
TSidList <- c('S2LR7FM5LvHegu','WEBc7a2aaa7','WEB0b46bdab','WEB55f96880')
              
TS <- extractTs(D)
TS <- TS[which(pullTsVariable(TS,"paleoData_TSid") %in% TSidList)]
TS[[which(pullTsVariable(TS,"paleoData_TSid") == 'WEB55f96880')]]$direction <- 'negative'

df <- tidyTs(TS,age.var = "age")
plot.df <-  df %>% 
 filter(between(age,0,12000)) %>% #only years from 1600 to 2000
 group_by(paleoData_TSid) %>% #group by column
 arrange(geo_longitude) #and sort by longitude
stackplot <- plotTimeseriesStack(plot.df,time.var ='age',
                                 scale.factor=0.3,
                                 scale.height = 0.38,
                                 lab.buff=0.04,
                                 lab.space = 2,
                                 lab.size	= 2.4,
                                 invert.var='direction',
                                 color.ramp = rep('Black',4)) + 
 scale_x_reverse(name="Age (ka BP)" ,
                 expand=c(0,0), limits= c(12000,0),
                 breaks=seq(0,12000,2000),labels=seq(0,12,2))+
 scale_fill_manual(values=c('grey','firebrick','dark orchid','grey'))+
 scale_color_manual(values=c('grey','firebrick','dark orchid','grey'))+
 theme_void() +
 theme(legend.position = 'none',
       text = element_text(family='sans',size=8),
       plot.background = element_rect(fill = 'white',color='white'),
       #panel.border = element_rect(fill = NA,color='grey20',size=0.6),
       plot.margin = unit(c(0.1,0,0.1,0), "in"),
       axis.ticks = element_line(color='Black'),
       axis.ticks.length.x =unit(5,"pt"),
       axis.text.x =element_text(family='sans',size=8,color='Black'),
       axis.title.x =element_text(family='sans',size=8,color='Black'),
       axis.line.x = element_line(color='Black',size=0.5),
       panel.grid.major.y=element_line(color='Black',size=0.5))


stackplot2 <- stackplot +
 scale_x_reverse(name="Age (ka BP)" ,
                 expand=c(0,0), limits= c(13100,-1100),
                 breaks=seq(0,12000,2000),labels=seq(0,12,2))+
 theme(plot.margin = unit(rep(0.1,4), "in"),
       panel.grid.major.y=element_blank(),
       plot.background = element_rect(fill = 'white',color='white'),
       panel.border = element_rect(fill = NA,color='white'),
       axis.ticks = element_line(color='White'),
       axis.text.x =element_text(family='sans',size=8,color='White'),
       axis.title.x =element_text(family='sans',size=8,color='White'),
       axis.line.x = element_line(color='White',size=0.5))

x <- 0.12

plt <- ggdraw() + 
 draw_plot(stackplot2,x=0, y=0, width = 1, height = 0.98) +
 draw_plot(stackplot,x=0.09, y=0, width = 0.82, height = 0.98) +
 draw_label( "(a) Lake Dobson (Tazmania)", x = x, y = 0.98,fontfamily='sans',color='Black',size =8,fontface='bold',hjust = 0)+
 draw_label("(b) Lake Voëlvlei (South Africa)" , x = x, y = 0.67,fontfamily='sans',color='Black',size = 8,fontface='bold',hjust = 0)+
 draw_label("(c) Botuver Cave (Brazil)", x = x, y = 0.53,fontfamily='sans',color='Black',size = 8,fontface='bold',hjust = 0)+
 draw_label( "(d) Lake Tamar (Chile)", x = x, y = 0.12,fontfamily='sans',color='Black',size = 8,fontface='bold',hjust = 0)


ggsave(plot=plt, width = 6.5, height = 7, dpi = 600,device='png',
     filename = paste(file.path(dataDir,'Figures','stackPlot.png',sep='')))


#Brazil   S2LR7FM5LvHegu   # positive = more extratropical vs amazon
#Chile    WEBc7a2aaa7      # Positive = more precipitation
#Africa
#Aus                       #Shift from Westerly to ENSO dominated cliamte through holocene
#NZ     #nd. A more negative δ18O will result from increased contribution to the annual rainfall budget from winter southwesterlies
