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
              Season   = pullTsVariable(lipdTSO,'climateInterpretation1_seasonalityGeneral'),
              Interp    = pullTsVariable(lipdTSO,'climateInterpretation1_variable'),
              AgeRange      = paste(round(pullTsVariable(lipdTSO,'ageMax')/1000,1),
                               round(pullTsVariable(lipdTSO,'ageMin')/1000,1),sep=' - '),
              PublicationDOI = pullTsVariable(lipdTSO,'pub1_doi'),
              SourceURL  = pullTsVariable(lipdTSO,'originalDataUrl'),
              Lat      = round(pullTsVariable(lipdTSO,'geo_latitude'),2),
              Lon      = round(pullTsVariable(lipdTSO,'geo_longitude'),2)) #Source

tbl$Proxy[which(tbl$Proxy == 'Ice Accumulation')] <- 'Acc.'

regList <- readShapePoly(file.path(dataDir,'Data','IPCC_refRegions','IPCC-WGI-reference-regions-v4.shp'))

size = nrow(tbl)

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
 regTbl$Region <- rep(regList$Name[which(regList$Acronym==reg)],nrow(regTbl))
 df <- rbind(df,regTbl)             #   it sits
}
out[[1]] <- df

idx <- c(length(out[[1]]$Region))
idx2 <- c()
for (i in 1:(length(out[[1]]$Region)-1)){
 if(out[[1]]$Region[i] != out[[1]]$Region[i+1]){
  idx <- c(idx,i)
  idx2 <- c(idx2,i+1)
 }
}
outTbl <- reactable(
 out[[1]],
 groupBy = "Region",
 defaultPageSize=size,
 wrap=FALSE,
 resizable=TRUE,
 searchable=TRUE,
 filterable=TRUE,
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
  if (index %in% idx) {list(`border-bottom` = "thin solid")}
  else if (index %in% idx2) {list(`border-top` = "thin solid")}},
 columns = list(
  Region = colDef(width=160,align = 'left'),
  Dataset = colDef(width=150,align = 'left'),
  Archive  = colDef(width=80,align = 'left'),
  Proxy  = colDef(width=55,align = 'left'),
  Lat  = colDef(width=40,align = 'left'),
  Lon  = colDef(width=50,align = 'left'),
  AgeRange = colDef(width=60,align = 'left'),
  Interp = colDef(width=40,align = 'left'),
  PublicationDOI  = colDef(width=250,align = 'left'),
  SourceURL  = colDef(width=250,align = 'left')
 )) 
outTbl

html_file <- file.path(dataDir,'Figures','Table','appendix1.html')
saveWidget(widget = outTbl, file = html_file, selfcontained = TRUE)
#webshot(url = html_file, file =  file.path(dataDir,'Figures','Table','table1.png'),
  #      delay = 0.1)