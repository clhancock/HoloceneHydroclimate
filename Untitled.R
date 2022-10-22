z<- c()
ids <- c()
ds<-c()
median<-c()
max<-c()
for (ts in 1:length(lipdData[[var]])){
  #
  id<-'WEB88fe721c'
  ts <- which(pullTsVariable(lipdData[[var]],"paleoData_TSid")==id)
  tso                  <- lipdData[[var]][[ts]]
  Dselect <- D_hc[[tso$dataSetName]]
  n_chron <- length(Dselect[["chronData"]][[1]][["measurementTable"]])
  n_paleo <- max(length(Dselect[["paleoData"]]),length(Dselect[["paleoData"]][[1]][["measurementTable"]]))
  if (n_chron == 0){                                         i_chron <- NA
  } else if (n_chron*n_paleo == 1){                          i_chron <- 1
  } else if (grepl("Composite",tso$paleoData_variableName)){
    offset<-c()
    for (j in 1:n_chron){
      paleoAges <- tso$age
      chronAges <- Dselect[["chronData"]][[1]][["measurementTable"]][[j]][["age"]]$values
      offset <- c(offset,max(abs(max(min(paleoAges,na.rm=TRUE),0)-min(chronAges,na.rm=TRUE)),
                             abs(min(max(paleoAges,na.rm=TRUE),12000)-max(chronAges,na.rm=TRUE))))
    }
    i_chron <- which(offset==min(offset))
  } else if (n_chron != n_paleo & (n_chron+1 != n_paleo | tso$archiveType!='Speleothem')){      
      offset<-c()
      for (j in 1:n_chron){
        paleoAges <- tso$age
        chronAges <- Dselect[["chronData"]][[1]][["measurementTable"]][[j]][["age"]]$values
        offset <- c(offset,max(abs(min(paleoAges,na.rm=TRUE)-min(chronAges,na.rm=TRUE)),
                               abs(max(paleoAges,na.rm=TRUE)-max(chronAges,na.rm=TRUE))))
      }
      i_chron <- which(offset==min(offset))
  } else{
    for (i in 1:n_paleo){
      if (!is.null(Dselect[["paleoData"]][[1]][["measurementTable"]][[i]][[tso$paleoData_variableName]]$TSid)){
        if (Dselect[["paleoData"]][[1]][["measurementTable"]][[i]][[tso$paleoData_variableName]]$TSid == tso$paleoData_TSid){
          i_chron <- i
        }
      }
    }
  }
  if (tso$paleoData_TSid %in% c('lcRpbhIeSqLqRxKnBNF','lcRseG41EyZPFn7vBuL','LPD1498ea89','LPD03892251')){i_chron<-NA}
  tso$ageMin     <- max(min(tso$age[which(tso$age<12400)],na.rm=TRUE),0)
  tso$ageMax     <- min(max(tso$age[which(tso$age<12400)],na.rm=TRUE),12000)
  if (!is.na(i_chron)){
    tso$chronData_table <- Dselect[["chronData"]][[1]][["measurementTable"]][[i_chron]]
    names <- names(tso$chronData_table)[grepl(c('age'),names(tso$chronData_table),ignore.case=TRUE)]
    for (name in c('type','error','min','max','hi','radio','14','lo','comment','use',"up","old","young","std","reservoir","σ","low","Rejected","±","comment","err","uncertainty","unc")){
      names <- names[which(grepl(name,names,ignore.case=TRUE)==FALSE)]
    }
    ageColName <- NA
    for (name in c("Th230/Th232","Median","calib.14C","14C.raw","14C Age","14C age BP","Cal yr chosen","age14C","14C age (yr BP)","C14 age dated",
                   tail(names,1),"corrected 230Th Age",
                   "Median Year BP",'age','Age','calAge','CalAge','CalibratedAge','Calibrated Age')){
      if (name %in% names(tso$chronData_table)){
        ageColName <- name
      } 
    }
    tso$chronData_ageName <- ageColName
    tso$chronData_ages <- as.numeric(tso$chronData_table[[tso$chronData_ageName]]$values)
    if (grepl("Composite",tso$paleoData_variableName)){
      tso$chronData_ages <- c()
      for (j in 1:n_chron){
        tso$chronData_ages <- c(tso$chronData_ages,as.numeric(Dselect[["chronData"]][[1]][["measurementTable"]][[j]][["age"]]$values))
      }
    }
    if ((length(which(tso$chronData_ages<13000))<=1) | is.na(sum(tso$chronData_ages)) | (sum(!is.na(tso$chronData_ages))<length(tso$chronData_ages)*0.5)){
      tso$chronData_table             <- NA
      tso$chronData_ageName           <- NA
      tso$chronData_ages              <- NA
      tso$chronData_ages_12k          <- NA
      tso$chronData_agesN_12k         <- NA
      tso$chronData_agesMaxGap_12k    <- NA
      tso$chronData_agesMedianGap_12k <- NA
    }else{
      tso$chronData_ages <- tso$chronData_ages[which(!is.na(tso$chronData_ages))]
      if(mean(tso$chronData_ages,na.rm=TRUE)<0){
        tso$chronData_ages<-tso$chronData_ages*-1 #
      }
      if(grepl("ka",tso$chronData_table[[tso$chronData_ageName]]$units)){
        tso$chronData_ages<-tso$chronData_ages*1000 #d 
      } else if(max(tso$chronData_ages,na.rm=TRUE)<300){
        tso$chronData_ages<-tso$chronData_ages*1000 #d 
      }
      tso$chronData_ages_12k          <- tso$chronData_ages[which(tso$chronData_ages<13000)]
      tso$chronData_agesN_12k         <- length(which(diff(tso$chronData_ages_12k)>0))+1
      tso$chronData_agesMaxGap_12k    <- max(diff(sort(c(tso$ageMin,tso$chronData_ages_12k,tso$ageMax))))
      tso$chronData_agesMedianGap_12k <- median(abs(diff(sort(c(tso$ageMin,tso$chronData_ages_12k,tso$ageMax)))))
    }
  } else{
    tso$chronData_table             <- NA
    tso$chronData_ageName           <- NA
    tso$chronData_ages              <- NA
    tso$chronData_ages_12k          <- NA
    tso$chronData_agesN_12k         <- NA
    tso$chronData_agesMaxGap_12k    <- NA
    tso$chronData_agesMedianGap_12k <- NA
  }
  if(tso$paleoData_TSid== "S2LR984dImqQzT"){print(tso$chronData_agesMaxGap_12k)} 
  z<- c(z,tso$chronData_agesN_12k)
  ids <- c(ids,tso$paleoData_TSid)
  ds<-c(ds,tso$dataSetName)
  median<-c(median,tso$chronData_agesMedianGap_12k)
  max<-c(max,tso$chronData_agesMaxGap_12k)
}

df <- data.frame(dataset=ds,tsid=ids,n=z,maxz=max,medianz=median)
View(df)
#View(df%>%filter(chron== "n chron != n paleo"))
#View(df%>%filter(chron== "Composite"))
