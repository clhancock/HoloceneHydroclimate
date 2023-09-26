#Load packages
library(rcrossref)
library(lipdR)
library(dplyr)
D_hc <- readLipd(paste0('https://lipdverse.org/HoloceneHydroclimate/0_8_0/HoloceneHydroclimate0_8_0.zip'))

#Load Data
style   <- "american-geophysical-union"
var     <- 'HC'
wd      <- '/Users/chrishancock/Library/CloudStorage/OneDrive-NorthernArizonaUniversity/Research/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate'
lipdTSO <- readRDS(file.path(wd,'Data','Proxy','lipdData.rds'))[[var]]

#Get dataframe of all DOIs
doi_df <- data.frame('Dataset' = pullTsVariable(lipdTSO,'dataSetName'),
                     'TSid'    = pullTsVariable(lipdTSO,'paleoData_TSid'))

for (i in c(1,2,3,4,5,6,7)){
  dois <- c()
  for (row in (1:length(lipdTSO))){
    doi <- lipdTSO[[row]][[paste0('pub',i,'_doi')]]
    if (length(doi) == 0){
      dois<-c(dois,NA)
    }else if (doi %in%  c("unpublished","","noPubOnRecord")){
      dois<-c(dois,NA)
    }else{
      dois<-c(dois,doi)
    }
  } 
  doi_df[paste0('DOI_',i)] <- dois
}


#Create empty table
doi_list <- data.frame('DOI'=sort(unique(unlist(doi_df[3:ncol(doi_df)]))))
doi_list$references <- doi_list$Author <- doi_list$Year <- NA

pb <- txtProgressBar(min=0,max = nrow(doi_list), style = 3)
#Get references from DOIs
for (i in 1:nrow(doi_list)){
  #Progress bar
  Sys.sleep(0.1) 
  setTxtProgressBar(pb, i)
  #
  doi <- doi_list$DOI[i]
  if (substr(doi,1,4)=="http"){doi<-sub(".*org/", "", doi)}
  #if (doi ){doi<-NA}
  if (substr(doi,1,3) %in% c("10.","htt") & grepl('10.',doi)){
    ref <- cr_cn(doi, format="bibentry")
    doi_list$references[i] <- cr_cn(doi, format="text", style=style)
    doi_list$Author[i] <- ref$author
    doi_list$Year[i]   <- as.numeric(ref$year)
  }else{
    dsn <- doi_df[which(doi_df==doi, arr.ind = TRUE)[1.1],1]
    if(length(D_hc[[dsn]]$pub)>0){
      for (p in 1:length(D_hc[[dsn]]$pub)){
        if (doi == D_hc[[dsn]]$pub[[p]]$doi){
          if(!is.null(D_hc[[dsn]]$pub[[p]]$author[[1]]$name)){
            doi_list$Author[i] <- D_hc[[dsn]]$pub[[p]]$author[[1]]$name
          } else{
            doi_list$Author[i] <- strsplit(doi, "[.,]")[[1]][1]
          }
          doi_list$Year[i] <- as.numeric(D_hc[[dsn]]$pub[[p]]$year)
        } else{
          doi_list$Author[i] <- strsplit(doi, "[.,]")[[1]][1]
          doi_list$Year[i] <- tail(strsplit(dsn, '[.]')[[1]],n=1)
        }
      } 
    }else{print(paste('NO REF INFO FOR DATASET',dsn))}
    doi_list$references[i] <- doi
  }
  # Add author and years for sorting
  if(!is.na(doi_list$references[i])){

  }
}

outDir <- file.path(wd,'Data','Proxy','References')
#List of references
doi_list <- doi_list %>% arrange(Author, Year) 
write.csv(doi_list[,c('references')], file.path(outDir,paste0("LiPD_referencesList_",var,".csv")))              

#Create New DF to output
doi_df_out <- doi_df
for (r in 1:nrow(doi_df_out)){
  for (c in 3:ncol(doi_df_out)){
    if(!is.na(doi_df_out[r,c])){
      doi_df_out[r,c] <- which(doi_list$DOI == doi_df[r,c])
    } else{
      if (c == 3){
        doi_df_out[r,c] <- 'noPubOnRecord'
      } else{
        doi_df_out[r,c] <- ''
      }
    }
  }
}
write.csv(doi_df_out, file.path(outDir,paste0("LiPD_referencesKey_",var,".csv")))                 # Apply writeLines function

#Write text file with count for each references
references <- doi_list$references
for (i in 1:length(references)){
  references[[i]] <- paste0(i,'. ',references[[i]])
}
writeLines(references, file.path(outDir,paste0("LiPD_referencesList_",var,".txt")))                 # Apply writeLines function

