library(ggplot2)
library(oxcAAR)
quickSetupOxcal()
getCalMedian <- function(c14,unc){
 calout <- oxcAAR::oxcalCalibrate(c14,unc,)
 cp <- cumsum(calout$`1`$raw_probabilities$probabilities)
 cp <- cp/max(cp)
 calAge <- geoChronR::convertAD2BP(calout$`1`$raw_probabilities$dates[which(abs(cp-.5)==min(abs(cp-.5)))])
 return(calAge)
}


data <- read.csv("/Users/chrishancock/Desktop/book1.csv")
ages <- c()
for (n in 1:nrow(data)){
  ages<-c(ages,getCalMedian(as.numeric(data$age)[n],200))
}
data$newages <- ages

write.csv(data,"/Users/chrishancock/Desktop/book1.csv")

#setOxcalExecutablePath("/Users/chrishancock/Desktop/OxCal")

#Test
# z <- c(0,
#        0.5,
#        1,
#        1.5,
#        2,
#        3,
#        4,
#        5,
#        6,
#        7,
#        8,
#        9,
#        10,
#        11,
#        12)
# 
# i <- rep(200,length(z))
# 
# for (n in 1:length(z)){
#  print(getCalMedian(z[n]*1000,i[n]))
# }

# ages<-c()
# for (age in TS__t[[i]]$age){
#  ages<-c(ages,getCalMedian(age,200))
# }

