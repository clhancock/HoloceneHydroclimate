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

lipdData <- readRDS(file=file.path(dir,'Data','LiPD','lipdData.rds'))
lipdTSO <- lipdData[[climVar]]

tbl <- tibble(Region         = pullTsVariable(lipdTSO,'geo_ipccRegion'),
              Dataset        = pullTsVariable(lipdTSO,'dataSetName'),
              TSid           = pullTsVariable(lipdTSO,'paleoData_TSid'),
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
              Long           = round(pullTsVariable(lipdTSO,'geo_longitude'),2)) #Source

for (i in 1:nrow(tbl)){
 row <- tbl[i,]
 if (row$SourceURL == 'wNAm'){row$SourceURL <- 'https://lipdverse.org/wNAm/1_0_0/'}
 if (row$PublicationDOI == 'https://doi.org/10.1126/sciadv.1601503'){row$SourceURL <- 'https://www.ncdc.noaa.gov/paleo/study/21091'}
 if (row$TSid == 'LS14ZHJU01A'){row$SourceURL <- 'https://doi.org/10.1016/j.epsl.2014.07.013'}
 if (row$TSid == 'LPDc4757948'){row$SourceURL <- 'https://www.ncdc.noaa.gov/paleo/study/30932'}
 if (row$PublicationDOI == '10.1016/j.gca.2012.10.040'){row$SourceURL <- row$PublicationDOI}
 if (row$TSid == 'WEB29098c59'){row$SourceURL <- 'https://www.ncdc.noaa.gov/paleo/study/15615'}
 if (row$TSid == 'GHd57f7596'){row$SourceURL <- 'https://www.ncdc.noaa.gov/paleo/study/17855'}
 if (row$TSid == ''){row$SourceURL <- 'https://www.ncdc.noaa.gov/paleo/study/23084'}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if (row$TSid == ''){row$SourceURL <- ''}
 if(substr(row$PublicationDOI, 1, 3) == '10.'){row$PublicationDOI <- paste('https://doi.org/',row$PublicationDOI,sep='')}
 if(substr(row$SourceURL, 1, 3) == '10.'){row$SourceURL <- paste('https://doi.org/',row$SourceURL,sep='')}
 tbl[i,] <- row
 
}
tbl <- tbl[, !(colnames(tbl) %in% c('TSid'))]
regList <- readShapePoly(file.path(dir,'Data','IPCC_refRegions','IPCC-WGI-reference-regions-v4.shp'))

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
 defaultPageSize = size,
 wrap       = FALSE,
 resizable  = TRUE,
 searchable = FALSE,
 filterable = TRUE,
 sortable = FALSE,
 borderless = TRUE,
 compact    = TRUE,
 defaultColDef = colDef(width=50),
 theme = reactableTheme(
  borderWidth = 1.2,
  borderColor='Black',
  headerStyle = list(`border-bottom` = "double",`border-top` = "thin solid",
                      borderColor = "#555",borderWidth = 2,align = 'left')),
 style = list(fontFamily = 'sans',fontSize = 10),
 rowStyle = function(index) {
  if (index %in% idx) {list(      `border-bottom` = "thin solid")}
  else if (index %in% idx2) {list(`border-top   ` = "thin solid")}},
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
  AgeRange        = colDef(width=60, align = 'center'),
  Resolution      = colDef(width=60, align = 'center'),
  PublicationDOI  = colDef(width=150,align = 'left', cell = function(value, index){
    htmltools::tags$a(href = value, target = "_blank", as.character(value))}),
  SourceURL       = colDef(width=150,align = 'left', cell = function(value, index){
    htmltools::tags$a(href = value, target = "_blank", as.character(value))}),
  Lat             = colDef(width=40, align = 'center'),
  Long            = colDef(width=50, align = 'center')
  ))  %>%
  add_title("SM1. List of proxy records included in the Holocene Hydroclimate dataset.
             Data are grouped by geographical region which are ordered according to Iturbide et al., (2020). 
             Within each region, records are listed alphabetically according to their dataset name.
             Columns can be resized using the column boundaries in the header row. 
             The empty boxes below the heading allow you to search for specific rows based on metadata characteristics within the column.",
            font_size=12, font_weight='normal', font_style='italic')

outTbl


html_file <- file.path(dir,'Figures','Table','SM1.html')
saveWidget(widget = outTbl, file = html_file, selfcontained = TRUE)
