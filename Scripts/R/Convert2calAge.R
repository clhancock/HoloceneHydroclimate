getCalMedian <- function(c14,unc){
 calout <- oxcAAR::oxcalCalibrate(c14,unc,)
 cp <- cumsum(calout$`1`$raw_probabilities$probabilities)
 cp <- cp/max(cp)
 calAge <- geoChronR::convertAD2BP(calout$`1`$raw_probabilities$dates[which(abs(cp-.5)==min(abs(cp-.5)))])
 return(calAge)
}

#setOxcalExecutablePath("/Users/chrishancock/Desktop/OxCal")


z <- c(1040,
       1650,
       2410,
       3910,
       5340,
       6580,
       7010,
       7790,
       9100,
       9990,
       10290,
       11520,
       34980,
       34260,
       34080)

i <- c(180,
       150,
       130,
       130,
       150,
       140,
       210,
       180,
       190,
       100,
       140,
       110,
       690,
       580,
       650)

for (n in 1:length(z)){
 print(getCalMedian(z[n],i[n]))
}


data <- read.csv("/Users/chrishancock/Desktop/book1.csv")
ages <- c()
for (n in 1:nrow(data)){
 ages<-c(ages,getCalMedian(as.numeric(data$age)[n],200))
}
data$newages <- ages

write.csv(data,"/Users/chrishancock/Desktop/book1.csv")
