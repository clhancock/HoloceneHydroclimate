library(maptools)
library(dplyr)
library(lipdR)
library(reactable)
library(reactablefmtr)
library(htmlTable)
library(htmltools)
library(htmlwidgets)
library(webshot)
dataDir <- '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate' #
climVar <- 'HC'

lipdData <- readRDS(file=file.path(dataDir,'Data','LiPD','lipdData.rds'))
lipdTSO <- lipdData[[climVar]]

tbl <- tibble(Region   = pullTsVariable(lipdTSO,'geo_ipccRegion'),
              Dataset     = pullTsVariable(lipdTSO,'dataSetName'),
              Archive  = pullTsVariable(lipdTSO,'archiveType'),
              Proxy    = pullTsVariable(lipdTSO,'paleoData_proxy'),
              #Lat      = round(pullTsVariable(lipdTSO,'geo_latitude'),2),
              #Lon      = round(pullTsVariable(lipdTSO,'geo_longitude'),2),
              Season   = pullTsVariable(lipdTSO,'climateInterpretation1_seasonalityGeneral'),
              Interp    = pullTsVariable(lipdTSO,'climateInterpretation1_variable'),
              Ages     = paste(round(pullTsVariable(lipdTSO,'ageMax')/1000,1),
                              round(pullTsVariable(lipdTSO,'ageMin')/1000,1),sep=' - '),
              Publication = pullTsVariable(lipdTSO,'pub1_doi'),
              Source   = pullTsVariable(lipdTSO,'originalDataUrl')) #Source

tbl$Proxy[which(tbl$Proxy == 'Ice Accumulation')] <- 'Acc.'

regList <- readShapePoly(file.path(dataDir,'Data','IPCC_refRegions','IPCC-WGI-reference-regions-v4.shp'))

size = 56

for (reg in regList$Acronym){
 print(reg)
 if (reg == regList$Acronym[1]){
  page  <- 1
  out <- vector(mode='list')
  freshdf <- data.frame(matrix(ncol = ncol(tbl)))[-1,]
  colnames(freshdf) <- colnames(tbl)
  df <- freshdf
 }
 regTbl <- tbl  %>% filter(tbl$Region==reg) %>% arrange(Dataset)
 if (nrow(df)+nrow(regTbl) <=  size){ #If it fits, 
  df <- rbind(df,regTbl)             #   it sits
 } else{
  if (nrow(df) ==  size){
   out[[page]] <- df
   page <- page + 1
   df <- freshdf 
  }
  #Add data to new pages as needed for number of records
  count <- size-nrow(df)
  divides <- c(0,count)
  while (count <= nrow(regTbl)){ #determine number of new pages
   count <- count + size
   divides <- c(divides,count)
  }
  divides <- divides[which(divides<count)]
  for (d in divides){ #fill data
   df <- rbind(df,regTbl[(d+1):min(d+size-nrow(df),nrow(regTbl)),])
   if (nrow(df) ==  size){
    out[[page]] <- df
    page <- page + 1
    df <- freshdf 
   }
  }
 }
}
out[[page]] <- df

n <- 0
for (i in 1:length(out)){
 print(nrow(out[[i]]))
 n <- n+nrow(out[[i]])
}


idx <- c(length(out[[1]]$Region))
for (i in 1:(length(out[[1]]$Region)-1)){
 if(out[[1]]$Region[i] != out[[1]]$Region[i+1]){
  idx <- c(idx,i)
 }
}
outTbl <- reactable(
 out[[1]],
 defaultPageSize=size,
 #compact=TRUE,
 wrap=FALSE,
 resizable=TRUE,
 borderless = TRUE,
 compact=TRUE,
 defaultColDef = colDef(width=50),
 theme = reactableTheme(
  borderWidth = 1.2,
  borderColor='Black',
  headerStyle = list(`border-bottom` = "double",`border-top` = "thin solid",
                     borderColor = "#555",borderWidth = 2,align = 'left')),
  style = list(fontFamily = 'sans',fontSize = 10),
 rowStyle = function(index) {
  if (index %in% idx) {list(`border-bottom` = "thin solid")}},
 columns = list(
  Region = colDef(
   align = 'left',
   style = JS("function(rowInfo, colInfo, Region) {
        const firstSorted = Region.sorted[0]
        // Merge cells if unsorted or sorting by school
        if (!firstSorted || firstSorted.id === 'Region') {
          const prevRow = Region.pageRows[rowInfo.viewIndex - 1]
          if (prevRow && rowInfo.row['Region'] === prevRow['Region']) {
            return { visibility: 'hidden' }
          }
        }
      }")),
  Dataset = colDef(width=150,align = 'left'),
  Archive  = colDef(width=80,align = 'left'),
  Proxy  = colDef(width=55,align = 'left'),
  #Lat  = colDef(width=40,align = 'left'),
  #Lon  = colDef(width=50,align = 'left'),
  Ages = colDef(width=55,align = 'left'),
  Interp = colDef(width=40,align = 'left'),
  Publication  = colDef(width=180,align = 'left'),
  Source  = colDef(width=230,align = 'left')
 ))  
outTbl

html_file <- "table.html"
saveWidget(widget = outTbl, file = html_file, selfcontained = TRUE)
webshot(url = html_file, file =  file.path(dataDir,'Figures','Table','table1.png'),
        delay = 0.1)
