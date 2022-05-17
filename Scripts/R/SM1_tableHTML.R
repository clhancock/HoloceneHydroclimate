library(maptools)
library(dplyr)
library(lipdR)
library(reactable)
library(reactablefmtr)
library(htmlTable)
library(htmltools)
library(htmlwidgets)
library(webshot)
dir     <- getwd() #'/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate' #
climVar <- 'HC'

#Load Data
lipdData <- readRDS(file=file.path(dir,'Data','LiPD','lipdData.rds'))
lipdTSO <- lipdData[[climVar]]
regList <- readShapePoly(file.path(dir,'Data','IPCC_refRegions','IPCC-WGI-reference-regions-v4.shp'))

#Create table with metadata from LiPD files
tbl <- tibble(Region         = pullTsVariable(lipdTSO,'geo_ipccRegion'),
              Dataset        = pullTsVariable(lipdTSO,'dataSetName'),
              #TSid           = pullTsVariable(lipdTSO,'paleoData_TSid'),
              Archive        = pullTsVariable(lipdTSO,'archiveType'),
              Category       = pullTsVariable(lipdTSO,'CategorySpecific'),
              Proxy          = pullTsVariable(lipdTSO,'paleoData_proxy'),
              Season         = pullTsVariable(lipdTSO,'climateInterpretation1_seasonalityGeneral'),
              Interp         = pullTsVariable(lipdTSO,'climateInterpretation1_variable'),
              Direction      = pullTsVariable(lipdTSO,'climateInterpretation1_direction'),
              AgeRange       = paste(round(pullTsVariable(lipdTSO,'ageMax')/1000,1),
                                     round(pullTsVariable(lipdTSO,'ageMin')/1000,1),sep=' - '),
              Resolution     = round(pullTsVariable(lipdTSO,'ageRes')),
              PublicationDOI = pullTsVariable(lipdTSO,'pub1_doi'),
              SourceURL      = pullTsVariable(lipdTSO,'originalDataUrl'),
              Lat            = round(pullTsVariable(lipdTSO,'geo_latitude'),2),
              Lon            = round(pullTsVariable(lipdTSO,'geo_longitude'),2)) #Source

#Modify url to include https
for (i in 1:nrow(tbl)){
 row <- tbl[i,]
 if(substr(row$PublicationDOI, 1, 3) == '10.'){row$PublicationDOI <- paste('https://doi.org/',row$PublicationDOI,sep='')}
 if(substr(row$SourceURL,      1, 3) == '10.'){row$SourceURL      <- paste('https://doi.org/',row$SourceURL,sep='')}
 tbl[i,] <- row
}

#Order data by region and dataset and replace region name from Acronym to full name
tbl2 <- tbl
tbl2[] <- NA
start <- 1
for (reg in regList$Acronym){
  regData <- tbl %>% filter(tbl$Region==reg) %>% arrange(Dataset)
  regData$Region <- as.character(regList$Name[which(regList$Acronym==reg)])
  if (nrow(regData) == 0){next}
  end <- start + nrow(regData) - 1 
  tbl2[start:end,] <- regData
  start <- end + 1
}


#Id Divides between regions
#idx <- c(lengthtbl2$Region))
#idx2 <- c()
#for (i in 1:(length(tbl2$Region)-1)){
 #if(tbl2$Region[i] != tbl2$Region[i+1]){
  #idx <- c(idx,i)
  #idx2 <- c(idx2,i+1)
 #}
#}

#Create Table
outTbl <- reactable(tbl2,
 groupBy = "Region",
 defaultPageSize = nrow(tbl2),
 height = 500,#), nrow(tbl2),
 wrap         = FALSE,
 resizable    = TRUE,
 searchable   = FALSE,
 filterable   = TRUE,
 sortable     = FALSE,
 showSortable = FALSE,
 borderless   = TRUE,
 compact      = TRUE,
 theme = reactableTheme(
  borderWidth = 1.2, borderColor='Black',
  headerStyle = list(`border-bottom` = "double",`border-top` = "thin solid", borderColor = "#555",borderWidth = 2,align = 'left')),
 style = list(fontFamily = 'sans',fontSize = 10),
 #rowStyle = function(index) {
  #if (index %in% idx) {list(      `border-bottom` = "thin solid")}
  #else if (index %in% idx2) {list(`border-top   ` = "thin solid")}},
 columns   = list(
  Region          = colDef(width=160,align = 'left'),
  Dataset         = colDef(width=150,align = 'left', cell = function(value, index){
    htmltools::tags$a(href = paste('http://lipdverse.org/HoloceneHydroclimate/current_version/',value,'.html',sep=''), target = "_blank", as.character(value))}),
  Archive         = colDef(width=80, align = 'left'),
  Category        = colDef(width=80, align = 'left'),
  Proxy           = colDef(width=80, align = 'left'),
  Interp          = colDef(width=45, align = 'center'),
  Season          = colDef(width=60, align = 'center'),
  Direction       = colDef(width=60, align = 'center'),
  AgeRange        = colDef(width=60, align = 'center', name='Age Range',filterable=FALSE),
  Resolution      = colDef(width=60, align = 'center',filterable=FALSE),
  PublicationDOI  = colDef(width=150,align = 'left', cell = function(value, index){
    htmltools::tags$a(href = value, target = "_blank", as.character(value))}),
  SourceURL       = colDef(width=150,align = 'left', cell = function(value, index){
    htmltools::tags$a(href = value, target = "_blank", as.character(value))}),
  Lat             = colDef(width=40, align = 'center'),
  Lon             = colDef(width=50, align = 'center')
  ))  %>%
  add_title("SM1. List of proxy records included in the Holocene Hydroclimate dataset.
             Data are grouped by geographical region which are ordered according to Iturbide et al., (2020). 
             Within each region, records are listed alphabetically according to their dataset name.
             Columns can be resized using the column boundaries in the header row. 
             The empty boxes below the heading allow you to search for specific records",
            font_size=12, font_weight='normal', font_style='italic')

outTbl


html_file <- file.path(dir,'Figures','Table','SM1.html')
saveWidget(widget = outTbl, file = html_file, selfcontained = TRUE)
